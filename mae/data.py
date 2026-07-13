"""HDF5 loading + mask reconstruction for the masked autoencoder.

`MAEDataset` reads features, masks, and continuous conditioning vars from an
HDF5 file, reconstructs all-zero masks for features that don't have a stored
mask, standardizes the continuous conditioning vars, and exposes the model
inputs as plain numpy arrays:

    ds = MAEDataset(
        data_path="catalog.h5",
        feature_names=[...],            # 51 names, defines column order of x
        masked_feature_names=[...],     # 48 names that have stored mask arrays
    )
    x_train, m_train, c_train = ds.train_arrays()
    x_val,   m_val,   c_val   = ds.val_arrays()
    model.train(x_train, m_train, c_train,
                val_x=x_val, val_mask=m_val, val_c=c_val)

Compared to `flow.FlowDataset`: the MAE encoder takes `(x, mask)` as inputs
and conditions on `c` separately via FiLM, so we expose the three arrays
independently rather than concatenating mask flags into `c`.

Imbalanced catalogs: `frac_lowsnr` / `frac_lowres` cap what share of the
loaded sample the low-S/N and poorly-resolved bulk may occupy (see
`MAEDataset` docstring) so it doesn't dominate training.
"""
from __future__ import annotations

import warnings
from typing import Optional, Sequence

import h5py
import numpy as np


DEFAULT_COND_VARS = ("mag_ellip", "fwhm", "sky_median", "sky_sigma", "sn_per_pixel")


class MAEDataset:
    """Load `(x, mask, c)` from HDF5 for the MAE.

    `x`     : (N, F)  continuous features, robust-scaled per column.
    `mask`  : (N, F)  binary; 1 where the value was imputed / bad.
    `c`     : (N, K)  continuous conditioning vars (standardized in place).

    Rebalancing (optional, composes with `n_subset`): most of the catalog is
    small, low-S/N, poorly-resolved sources; left alone they dominate
    training. Two caps bound their share of the *loaded* sample (train and
    val alike — the split happens downstream):

        frac_lowsnr : at most this fraction flagged by the lowsnr rule
        frac_lowres : at most this fraction flagged by the lowres rule

    A row is flagged when its variable falls on the *bad* side of the
    threshold, chosen by `lowsnr_side` / `lowres_side` ("below" or "above"):
    the lowsnr default is "below" (low S/N = small `sn_per_pixel`); the
    lowres default is "above" (bad seeing = large `fwhm`). Override the side
    together with the variable — e.g. a resolution-element count `nres`
    needs `lowres_side="below"` (poorly resolved = few elements).

    Thresholds are in raw catalog units (applied before standardization).
    `lowsnr_var` / `lowres_var` are read via `cond_path` and do not need to
    be in `cond_var_names`. A row flagged by both rules counts against both
    caps. Caps are "at most": hit exactly while flagged rows are plentiful;
    if unflagged rows run out before reaching `n_subset`, the sample shrinks
    (with a warning) rather than violating the caps. `balance_info` records
    totals, kept counts, and achieved fractions. Note `feat_median/iqr` and
    `cond_mu/sd` are then computed on the rebalanced sample — save and reuse
    them for any downstream scoring.
    """

    def __init__(
        self,
        data_path: str,
        feature_names: Sequence[str],
        masked_feature_names: Optional[Sequence[str]] = None,
        cond_var_names: Sequence[str] = DEFAULT_COND_VARS,
        *,
        feature_path: str = "training/{name}",
        mask_path: str = "masks/mask_{name}",
        cond_path: str = "condition/{name}",
        standardize_continuous_cond: bool = True,
        normalize_features: bool = True,
        val_frac: float = 0.1,
        n_subset: Optional[int] = None,
        frac_lowsnr: Optional[float] = None,
        lowsnr_thresh: Optional[float] = None,
        lowsnr_var: str = "sn_per_pixel",
        lowsnr_side: str = "below",
        frac_lowres: Optional[float] = None,
        lowres_thresh: Optional[float] = None,
        lowres_var: str = "fwhm",
        lowres_side: str = "above",
        seed: int = 0,
    ):
        for frac, thresh, side, label in (
            (frac_lowsnr, lowsnr_thresh, lowsnr_side, "lowsnr"),
            (frac_lowres, lowres_thresh, lowres_side, "lowres"),
        ):
            if side not in ("below", "above"):
                raise ValueError(f"{label}_side must be 'below' or 'above', got {side!r}")
            if frac is not None:
                if not (0.0 <= frac <= 1.0):
                    raise ValueError(f"frac_{label} must be in [0, 1], got {frac!r}")
                if thresh is None:
                    raise ValueError(f"frac_{label} requires {label}_thresh")
        self.data_path = data_path
        self.feature_names = list(feature_names)
        self.masked_feature_names = list(masked_feature_names or [])
        self.cond_var_names = list(cond_var_names)
        self.feature_path = feature_path
        self.mask_path = mask_path
        self.cond_path = cond_path
        self.standardize_continuous_cond = standardize_continuous_cond
        self.normalize_features = normalize_features
        self.val_frac = val_frac
        self.n_subset = n_subset
        self.frac_lowsnr = frac_lowsnr
        self.lowsnr_thresh = lowsnr_thresh
        self.lowsnr_var = lowsnr_var
        self.lowsnr_side = lowsnr_side
        self.frac_lowres = frac_lowres
        self.lowres_thresh = lowres_thresh
        self.lowres_var = lowres_var
        self.lowres_side = lowres_side
        self.seed = seed
        self.balance_info: Optional[dict] = None

        self._load()

    def _subset_indices(self, f: h5py.File) -> np.ndarray:
        """Sorted row indices to load: uniform `n_subset`, or the rebalanced
        sample honoring the `frac_lowsnr` / `frac_lowres` caps (see class
        docstring for the caps' semantics)."""
        n_total = f[self.feature_path.format(name=self.feature_names[0])].shape[0]
        rng = np.random.default_rng(self.seed)

        if self.frac_lowsnr is None and self.frac_lowres is None:
            return np.sort(rng.choice(n_total, size=self.n_subset, replace=False))

        def flagged(var, thresh, side):
            v = np.asarray(f[self.cond_path.format(name=var)][:])
            return v < thresh if side == "below" else v > thresh

        flag_snr = (flagged(self.lowsnr_var, self.lowsnr_thresh, self.lowsnr_side)
                    if self.frac_lowsnr is not None else np.zeros(n_total, dtype=bool))
        flag_res = (flagged(self.lowres_var, self.lowres_thresh, self.lowres_side)
                    if self.frac_lowres is not None else np.zeros(n_total, dtype=bool))
        frac_snr = self.frac_lowsnr or 0.0
        frac_res = self.frac_lowres or 0.0

        idx_snr   = np.flatnonzero(flag_snr & ~flag_res)   # flagged by one rule only
        idx_res   = np.flatnonzero(flag_res & ~flag_snr)
        idx_both  = np.flatnonzero(flag_snr & flag_res)    # counts against both caps
        idx_clean = np.flatnonzero(~flag_snr & ~flag_res)

        # Fixed point on the final size T (caps are fractions of T itself).
        # Singly-flagged rows are taken before doubly-flagged ones — the
        # latter consume both budgets, so spending them last maximizes how
        # much flagged data survives. T is non-increasing, so this terminates.
        T = self.n_subset if self.n_subset is not None else n_total
        while True:
            cap_snr, cap_res = int(frac_snr * T), int(frac_res * T)
            n_snr   = min(len(idx_snr), cap_snr)
            n_res   = min(len(idx_res), cap_res)
            n_both  = min(len(idx_both), cap_snr - n_snr, cap_res - n_res)
            n_clean = min(len(idx_clean), T - n_snr - n_res - n_both)
            T_new = n_snr + n_res + n_both + n_clean
            if T_new == T:
                break
            T = T_new

        if self.n_subset is not None and T < self.n_subset:
            warnings.warn(
                f"rebalanced subset has {T} rows, short of n_subset="
                f"{self.n_subset}: only {len(idx_clean)} rows pass both cuts "
                "and the frac_lowsnr/frac_lowres caps bound the rest."
            )

        parts = [rng.choice(pool, size=k, replace=False)
                 for pool, k in ((idx_snr, n_snr), (idx_res, n_res),
                                 (idx_both, n_both), (idx_clean, n_clean))
                 if k > 0]
        keep = (np.sort(np.concatenate(parts)) if parts
                else np.empty(0, dtype=np.int64))

        self.balance_info = {
            "n_total": int(n_total),
            "n_kept": int(T),
            "n_lowsnr_total": int(flag_snr.sum()),
            "n_lowsnr_kept": int(n_snr + n_both),
            "frac_lowsnr_kept": (n_snr + n_both) / max(T, 1),
            "n_lowres_total": int(flag_res.sum()),
            "n_lowres_kept": int(n_res + n_both),
            "frac_lowres_kept": (n_res + n_both) / max(T, 1),
        }
        return keep

    def _load(self):
        feats: list[np.ndarray] = []
        masks: list[np.ndarray] = []
        cond:  list[np.ndarray] = []

        with h5py.File(self.data_path, "r") as f:
            stored_masks = set(self.masked_feature_names)

            rebalance = self.frac_lowsnr is not None or self.frac_lowres is not None
            if self.n_subset is not None or rebalance:
                subset_idx = self._subset_indices(f)
                self.subset_idx = subset_idx
            else:
                subset_idx = slice(None)
                self.subset_idx = None

            for name in self.feature_names:
                feats.append(f[self.feature_path.format(name=name)][:][subset_idx])
            n = feats[0].shape[0]

            for name in self.feature_names:
                if name in stored_masks:
                    m = f[self.mask_path.format(name=name)][:][subset_idx]
                else:
                    m = np.zeros(n, dtype=np.uint8)
                masks.append(m.astype(np.float32))

            for name in self.cond_var_names:
                cond.append(f[self.cond_path.format(name=name)][:][subset_idx])

        self.x    = np.stack(feats, axis=1).astype(np.float32)
        self.mask = np.stack(masks, axis=1).astype(np.float32)
        self.c    = np.stack(cond,  axis=1).astype(np.float32)

        if not np.isfinite(self.x).all():
            raise ValueError(f"{(~np.isfinite(self.x)).sum()} non-finite values in x")
        if not np.isfinite(self.c).all():
            raise ValueError(f"{(~np.isfinite(self.c)).sum()} non-finite values in c")

        if self.normalize_features:
            f_dim = self.x.shape[1]
            self.feat_median = np.zeros(f_dim, dtype=np.float32)
            self.feat_iqr    = np.ones(f_dim,  dtype=np.float32)
            valid = self.mask < 0.5
            for j in range(f_dim):
                col = self.x[valid[:, j], j]
                if col.size == 0:
                    continue
                med = np.median(col)
                p5, p95 = np.percentile(col, [5, 95])
                self.feat_median[j] = med
                self.feat_iqr[j]    = max(float(p95 - p5), 1e-8)
            self.x = ((self.x - self.feat_median) / self.feat_iqr).astype(np.float32)
        else:
            f_dim = self.x.shape[1]
            self.feat_median = np.zeros(f_dim, dtype=np.float32)
            self.feat_iqr    = np.ones(f_dim,  dtype=np.float32)

        if self.standardize_continuous_cond:
            self.cond_mu = self.c.mean(axis=0)
            self.cond_sd = self.c.std(axis=0) + 1e-8
            self.c = (self.c - self.cond_mu) / self.cond_sd
        else:
            d = self.c.shape[1]
            self.cond_mu = np.zeros(d, dtype=np.float32)
            self.cond_sd = np.ones(d,  dtype=np.float32)

        rng = np.random.default_rng(self.seed)
        idx = rng.permutation(self.x.shape[0])
        n_val = int(self.x.shape[0] * self.val_frac)
        self.val_idx   = idx[:n_val]
        self.train_idx = idx[n_val:]

    # --- shapes ---
    @property
    def n_features(self) -> int:
        return self.x.shape[1]

    @property
    def n_cond(self) -> int:
        return self.c.shape[1]

    # --- model-ready slices ---
    def train_arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.x[self.train_idx], self.mask[self.train_idx], self.c[self.train_idx]

    def val_arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.x[self.val_idx], self.mask[self.val_idx], self.c[self.val_idx]
