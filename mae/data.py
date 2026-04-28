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
"""
from __future__ import annotations

from typing import Optional, Sequence

import h5py
import numpy as np


DEFAULT_COND_VARS = ("mag_ellip", "fwhm", "sky_median", "sky_sigma", "sn_per_pixel")


class MAEDataset:
    """Load `(x, mask, c)` from HDF5 for the MAE.

    `x`     : (N, F)  continuous features, expected pre-scaled to ~(-10, 10).
    `mask`  : (N, F)  binary; 1 where the value was imputed / bad.
    `c`     : (N, K)  continuous conditioning vars (standardized in place).
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
        val_frac: float = 0.1,
        n_subset: Optional[int] = None,
        seed: int = 0,
    ):
        self.data_path = data_path
        self.feature_names = list(feature_names)
        self.masked_feature_names = list(masked_feature_names or [])
        self.cond_var_names = list(cond_var_names)
        self.feature_path = feature_path
        self.mask_path = mask_path
        self.cond_path = cond_path
        self.standardize_continuous_cond = standardize_continuous_cond
        self.val_frac = val_frac
        self.n_subset = n_subset
        self.seed = seed

        self._load()

    def _load(self):
        feats: list[np.ndarray] = []
        masks: list[np.ndarray] = []
        cond:  list[np.ndarray] = []

        with h5py.File(self.data_path, "r") as f:
            stored_masks = set(self.masked_feature_names)

            if self.n_subset is not None:
                n_total = f[self.feature_path.format(name=self.feature_names[0])].shape[0]
                rng = np.random.default_rng(self.seed)
                subset_idx = np.sort(rng.choice(n_total, size=self.n_subset, replace=False))
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
