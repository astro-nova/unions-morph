"""Masked Autoencoder with FiLM conditioning.

Single-class sklearn-style API:

    model = MaskedAutoencoder(n_features=51, n_cond=5)
    model.train(x_tr, m_tr, c_tr, val_x=x_va, val_mask=m_va, val_c=c_va)
    z         = model.encode(x, m, c)              # (N, latent_dim)
    x_hat     = model.reconstruct(x, m, c)         # (N, n_features)
    score     = model.anomaly_score(x, m, c)       # (N,)
    per_feat  = model.per_feature_anomaly(x, m, c) # (N, n_features)
    model.save("ckpt.eqx")
    MaskedAutoencoder(n_features=51, n_cond=5).load("ckpt.eqx")

Architecture:
    encoder : (x ⊙ (1-m), m, c) → z         FiLM-MLP, depth `n_hidden_layers`, width `hidden_dim`
    decoder : (z, c)            → x_hat     FiLM-MLP, mirrored

Training:
    Each batch we draw an extra random mask `m_extra ~ Bernoulli(mask_ratio)`
    on top of the catalog mask `m_orig`. The encoder sees zeros at all hidden
    positions plus a binary `m_total` channel that flags them. The decoder
    reconstructs all `n_features`, but the loss is averaged only over
    positions that were *artificially* hidden and *originally* observed —
    standard MAE / BERT-style denoising objective. This forces the latent
    code to capture cross-feature structure rather than collapse to identity.

    Loss options (`loss_type`):
      - 'mse' : plain squared error on μ. Decoder out-dim = F.
      - 'nll' : Gaussian NLL with predicted log σ² per feature.
                Decoder out-dim = 2F (μ, log σ²); log σ² clipped to [-7, 7].
                Anomaly scoring then supports 'mahalanobis' = (x−μ)²/σ² and
                'nll' = 0.5(sq/σ² + log σ²) in addition to plain 'mse'.

The conditioning `c` (typically the 5 standardized observing-condition vars,
brightness first) modulates every hidden activation through FiLM:

    h ← (1 + γ(c)) ⊙ Linear(h) + β(c)

i.e. a per-channel affine transform whose scale and shift are linear in `c`.
At init, γ ≈ 0 / β ≈ 0 (small LeCun init on the FiLM head), so the network
starts as a plain MLP and learns brightness modulation gradually.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Optional

import numpy as np
import jax
import jax.numpy as jnp
import jax.random as jr
import equinox as eqx
import optax

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


# ---------------------------------------------------------------- modules
class FiLMBlock(eqx.Module):
    """Linear → FiLM(c) → GELU. FiLM head init is small so block ≈ Linear at start."""
    linear: eqx.nn.Linear
    film:   eqx.nn.Linear
    out_dim: int = eqx.field(static=True)

    def __init__(self, in_dim: int, out_dim: int, cond_dim: int, *, key):
        k1, k2 = jr.split(key, 2)
        self.linear = eqx.nn.Linear(in_dim, out_dim, key=k1)
        film = eqx.nn.Linear(cond_dim, 2 * out_dim, key=k2)
        # Shrink FiLM init so γ≈0, β≈0 → starts as identity modulation.
        film = eqx.tree_at(lambda l: l.weight, film, film.weight * 0.01)
        film = eqx.tree_at(lambda l: l.bias,   film, film.bias   * 0.0)
        self.film = film
        self.out_dim = out_dim

    def __call__(self, h, c):
        h = self.linear(h)
        gb = self.film(c)
        gamma, beta = gb[:self.out_dim], gb[self.out_dim:]
        h = (1.0 + gamma) * h + beta
        return jax.nn.gelu(h)


class _Encoder(eqx.Module):
    blocks: list
    head: eqx.nn.Linear

    def __init__(self, in_dim, hidden, depth, latent, cond_dim, *, key):
        keys = jr.split(key, depth + 1)
        blocks, d = [], in_dim
        for i in range(depth):
            blocks.append(FiLMBlock(d, hidden, cond_dim, key=keys[i]))
            d = hidden
        self.blocks = blocks
        self.head = eqx.nn.Linear(hidden, latent, key=keys[-1])

    def __call__(self, x_in, m, c):
        h = jnp.concatenate([x_in, m], axis=-1)
        for blk in self.blocks:
            h = blk(h, c)
        return self.head(h)


class _Decoder(eqx.Module):
    blocks: list
    head: eqx.nn.Linear

    def __init__(self, latent, hidden, depth, out_dim, cond_dim, *, key):
        keys = jr.split(key, depth + 1)
        blocks, d = [], latent
        for i in range(depth):
            blocks.append(FiLMBlock(d, hidden, cond_dim, key=keys[i]))
            d = hidden
        self.blocks = blocks
        self.head = eqx.nn.Linear(hidden, out_dim, key=keys[-1])

    def __call__(self, z, c):
        h = z
        for blk in self.blocks:
            h = blk(h, c)
        return self.head(h)


_LOG_VAR_CLIP = 7.0  # clip log σ² to [-7, 7] for NLL stability


class _MAEModule(eqx.Module):
    encoder: _Encoder
    decoder: _Decoder
    loss_type: str = eqx.field(static=True)
    n_features: int = eqx.field(static=True)

    def split_output(self, out):
        """Split decoder output into (mu, log_var). For mse, log_var is zeros."""
        if self.loss_type == "nll":
            mu = out[:self.n_features]
            log_var = jnp.clip(out[self.n_features:], -_LOG_VAR_CLIP, _LOG_VAR_CLIP)
        else:
            mu = out
            log_var = jnp.zeros_like(out)
        return mu, log_var


# ---------------------------------------------------------------- jitted kernels
def _forward_single(model, x_i, m_i, c_i):
    x_in = x_i * (1.0 - m_i)
    z = model.encoder(x_in, m_i, c_i)
    out = model.decoder(z, c_i)
    mu, log_var = model.split_output(out)
    return mu, log_var, z


@eqx.filter_jit
def _loss_fn(model, x, m_orig, c, key, mask_ratio):
    m_extra = (jr.uniform(key, x.shape) < mask_ratio).astype(x.dtype)
    m_total = jnp.maximum(m_orig, m_extra)

    def fwd(x_i, m_i, c_i):
        mu, log_var, _ = _forward_single(model, x_i, m_i, c_i)
        return mu, log_var

    mu, log_var = jax.vmap(fwd)(x, m_total, c)
    loss_mask = m_extra * (1.0 - m_orig)        # artificially-hidden & originally-observed
    sq = (mu - x) ** 2
    if model.loss_type == "nll":
        # 0.5 * ((x-μ)² / σ² + log σ²); constant 0.5*log(2π) dropped.
        per_elem = 0.5 * (sq * jnp.exp(-log_var) + log_var)
    else:
        per_elem = sq
    return (per_elem * loss_mask).sum() / jnp.clip(loss_mask.sum(), 1.0, None)


@eqx.filter_jit
def _step_fn(model, opt_state, x, m, c, key, mask_ratio, optimizer):
    loss, grads = eqx.filter_value_and_grad(_loss_fn)(model, x, m, c, key, mask_ratio)
    updates, opt_state = optimizer.update(grads, opt_state, model)
    model = eqx.apply_updates(model, updates)
    return model, opt_state, loss


@eqx.filter_jit
def _encode_fn(model, x, m, c):
    def per_row(x_i, m_i, c_i):
        x_in = x_i * (1.0 - m_i)
        return model.encoder(x_in, m_i, c_i)
    return jax.vmap(per_row)(x, m, c)


def _per_elem_score(sq, log_var, score_type):
    """Per-element score: 'mse' = sq; 'mahalanobis' = sq/σ²; 'nll' = 0.5(sq/σ² + log σ²)."""
    if score_type == "mahalanobis":
        return sq * jnp.exp(-log_var)
    if score_type == "nll":
        return 0.5 * (sq * jnp.exp(-log_var) + log_var)
    return sq


@eqx.filter_jit
def _reconstruct_fn(model, x, m, c):
    def per_row(x_i, m_i, c_i):
        mu, _, _ = _forward_single(model, x_i, m_i, c_i)
        return mu
    return jax.vmap(per_row)(x, m, c)


@eqx.filter_jit
def _predict_fn(model, x, m, c):
    """Per-row (mu, log_var). For loss_type='mse', log_var is zeros."""
    def per_row(x_i, m_i, c_i):
        mu, lv, _ = _forward_single(model, x_i, m_i, c_i)
        return mu, lv
    return jax.vmap(per_row)(x, m, c)


@eqx.filter_jit
def _self_recon_err_fn(model, x, m, c, score_type):
    """Per-source error on originally-observed positions (no extra masking)."""
    def per_row(x_i, m_i, c_i):
        mu, log_var, _ = _forward_single(model, x_i, m_i, c_i)
        sq = (mu - x_i) ** 2
        per_elem = _per_elem_score(sq, log_var, score_type)
        obs = 1.0 - m_i
        return (per_elem * obs).sum() / jnp.clip(obs.sum(), 1.0, None)
    return jax.vmap(per_row)(x, m, c)


@eqx.filter_jit
def _mc_recon_err_fn(model, x, m, c, key, mask_ratio, score_type):
    """Per-source error on positions we artificially hide (one MC draw)."""
    m_extra = (jr.uniform(key, x.shape) < mask_ratio).astype(x.dtype)
    m_total = jnp.maximum(m, m_extra)

    def per_row(x_i, m_total_i, m_extra_i, m_orig_i, c_i):
        mu, log_var, _ = _forward_single(model, x_i, m_total_i, c_i)
        loss_mask = m_extra_i * (1.0 - m_orig_i)
        sq = (mu - x_i) ** 2
        per_elem = _per_elem_score(sq, log_var, score_type)
        return (per_elem * loss_mask).sum() / jnp.clip(loss_mask.sum(), 1.0, None)

    return jax.vmap(per_row)(x, m_total, m_extra, m, c)


@eqx.filter_jit
def _per_feature_loo_fn(model, x, m, c):
    """For each (source, feature j), reconstruct x_j after masking only j.

    Returns (N, F) reconstructed values at the masked diagonal.
    """
    F = x.shape[1]
    eye = jnp.eye(F, dtype=x.dtype)

    def per_row(x_i, m_orig_i, c_i):
        def mask_j(j_onehot):
            m_total = jnp.maximum(m_orig_i, j_onehot)
            mu, _, _ = _forward_single(model, x_i, m_total, c_i)
            return mu
        x_hats = jax.vmap(mask_j)(eye)        # (F, F)
        return jnp.diagonal(x_hats)           # (F,)

    return jax.vmap(per_row)(x, m, c)


# ---------------------------------------------------------------- main class
class MaskedAutoencoder:
    """Masked autoencoder over `n_features` dims, FiLM-conditioned on `n_cond`."""

    def __init__(
        self,
        n_features: int,
        n_cond: int,
        *,
        # architecture
        latent_dim: int = 16,
        hidden_dim: int = 256,
        n_hidden_layers: int = 3,
        # training
        mask_ratio: float = 0.5,
        batch_size: int = 8192,
        learning_rate: float = 1e-3,
        n_epochs: int = 10,
        grad_clip: float = 1.0,
        loss_type: str = "mse",
        seed: int = 0,
    ):
        if loss_type not in ("mse", "nll"):
            raise ValueError(f"loss_type must be 'mse' or 'nll', got {loss_type!r}")
        self.n_features = n_features
        self.n_cond = n_cond
        self.latent_dim = latent_dim
        self.hidden_dim = hidden_dim
        self.n_hidden_layers = n_hidden_layers
        self.mask_ratio = mask_ratio
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.n_epochs = n_epochs
        self.grad_clip = grad_clip
        self.loss_type = loss_type
        self.seed = seed

        self.history: dict = {}
        self._build()

    def _build(self):
        key = jr.PRNGKey(self.seed)
        ke, kd = jr.split(key)
        encoder = _Encoder(
            in_dim=2 * self.n_features,        # [x_in | mask]
            hidden=self.hidden_dim,
            depth=self.n_hidden_layers,
            latent=self.latent_dim,
            cond_dim=self.n_cond,
            key=ke,
        )
        dec_out = 2 * self.n_features if self.loss_type == "nll" else self.n_features
        decoder = _Decoder(
            latent=self.latent_dim,
            hidden=self.hidden_dim,
            depth=self.n_hidden_layers,
            out_dim=dec_out,
            cond_dim=self.n_cond,
            key=kd,
        )
        self.model = _MAEModule(
            encoder=encoder,
            decoder=decoder,
            loss_type=self.loss_type,
            n_features=self.n_features,
        )

    # ----------------------------------------------------------------- batching
    def _iter_batches(self, x, m, c, *, shuffle, rng, batch_size=None):
        bs = batch_size or self.batch_size
        idx = np.arange(x.shape[0])
        if shuffle:
            idx = rng.permutation(idx)
        for i in range(0, len(idx), bs):
            b = idx[i:i + bs]
            yield jnp.asarray(x[b]), jnp.asarray(m[b]), jnp.asarray(c[b])

    # ----------------------------------------------------------------- training
    def train(
        self,
        x: np.ndarray,
        mask: np.ndarray,
        c: np.ndarray,
        *,
        val_x: Optional[np.ndarray] = None,
        val_mask: Optional[np.ndarray] = None,
        val_c: Optional[np.ndarray] = None,
        log_every: int = 0,
        progress: bool = True,
        checkpoint_dir: Optional[str | Path] = None,
    ):
        """Fit the MAE on `(x, mask, c)`. Optionally tracks validation loss.

        Populates `self.history`:
            step_loss          per-step training loss
            epoch_train_loss   mean training loss per epoch
            epoch_val_loss     mean validation loss per epoch (if val provided)
            epoch_time         seconds per epoch
            best_val_loss      lowest val loss observed
            best_epoch         epoch index of best val loss
        """
        if not (x.shape[0] == mask.shape[0] == c.shape[0]):
            raise ValueError(
                f"x, mask, c row-count mismatch: {x.shape[0]} / {mask.shape[0]} / {c.shape[0]}"
            )
        if x.shape[1] != self.n_features:
            raise ValueError(f"x has {x.shape[1]} features, expected {self.n_features}")
        if mask.shape[1] != self.n_features:
            raise ValueError(f"mask has {mask.shape[1]} columns, expected {self.n_features}")
        if c.shape[1] != self.n_cond:
            raise ValueError(f"c has {c.shape[1]} dims, expected {self.n_cond}")

        optimizer = optax.chain(
            optax.clip_by_global_norm(self.grad_clip),
            optax.adam(self.learning_rate),
        )
        opt_state = optimizer.init(eqx.filter(self.model, eqx.is_inexact_array))

        rng = np.random.default_rng(self.seed + 1)
        key = jr.PRNGKey(self.seed + 2)

        ckpt_dir = Path(checkpoint_dir) if checkpoint_dir is not None else None
        if ckpt_dir is not None:
            ckpt_dir.mkdir(parents=True, exist_ok=True)

        history = {
            "step_loss": [],
            "epoch_train_loss": [],
            "epoch_val_loss": [],
            "epoch_time": [],
            "best_val_loss": float("inf"),
            "best_epoch": -1,
        }
        has_val = val_x is not None and val_mask is not None and val_c is not None

        n_steps_per_epoch = (x.shape[0] + self.batch_size - 1) // self.batch_size
        use_tqdm = progress and tqdm is not None
        epoch_bar = tqdm(range(self.n_epochs), desc="epoch") if use_tqdm else range(self.n_epochs)

        for epoch in epoch_bar:
            t0 = time.time()

            train_losses = []
            step_iter = self._iter_batches(x, mask, c, shuffle=True, rng=rng)
            if use_tqdm:
                step_iter = tqdm(step_iter, total=n_steps_per_epoch,
                                 desc=f"epoch {epoch}", leave=False)

            for step_i, (xb, mb, cb) in enumerate(step_iter):
                key, sub = jr.split(key)
                self.model, opt_state, loss = _step_fn(
                    self.model, opt_state, xb, mb, cb, sub, self.mask_ratio, optimizer
                )
                lf = float(loss)
                train_losses.append(lf)
                history["step_loss"].append(lf)

                if use_tqdm:
                    recent = train_losses[-50:]
                    step_iter.set_postfix(loss=f"{lf:.4f}",
                                          mean=f"{np.mean(recent):.4f}")
                if log_every and step_i % log_every == 0:
                    print(f"  epoch {epoch} step {step_i:6d} loss={lf:.4f}")

            val_l = float("nan")
            if has_val:
                key, vsub = jr.split(key)
                val_losses = []
                for xb, mb, cb in self._iter_batches(val_x, val_mask, val_c,
                                                     shuffle=False, rng=rng):
                    vsub, kk = jr.split(vsub)
                    val_losses.append(float(_loss_fn(self.model, xb, mb, cb,
                                                     kk, self.mask_ratio)))
                val_l = float(np.mean(val_losses))
                history["epoch_val_loss"].append(val_l)
                if val_l < history["best_val_loss"]:
                    history["best_val_loss"] = val_l
                    history["best_epoch"] = epoch

            train_l = float(np.mean(train_losses))
            dt = time.time() - t0
            history["epoch_train_loss"].append(train_l)
            history["epoch_time"].append(dt)

            if use_tqdm:
                postfix = {"train": f"{train_l:.4f}", "dt": f"{dt:.1f}s"}
                if has_val:
                    postfix["val"] = f"{val_l:.4f}"
                epoch_bar.set_postfix(**postfix)
            else:
                msg = f"epoch {epoch}: train={train_l:.4f}"
                if has_val:
                    msg += f"  val={val_l:.4f}"
                msg += f"  ({dt:.1f}s)"
                print(msg)

            if ckpt_dir is not None:
                self.save(ckpt_dir / f"mae_epoch{epoch:03d}.eqx")

        self.history = history
        return self

    # --------------------------------------------------------------- inference
    def encode(
        self,
        x: np.ndarray,
        mask: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> np.ndarray:
        """Latent embedding `z`. No extra masking — uses the catalog mask only.

        Returns `(N, latent_dim)`.
        """
        out = np.empty((x.shape[0], self.latent_dim), dtype=np.float32)
        bs = batch_size or self.batch_size
        n_steps = (x.shape[0] + bs - 1) // bs
        rng = np.random.default_rng(0)
        it = self._iter_batches(x, mask, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="encode")
        pos = 0
        for xb, mb, cb in it:
            zb = _encode_fn(self.model, xb, mb, cb)
            n = zb.shape[0]
            out[pos:pos + n] = np.asarray(zb)
            pos += n
        return out

    def reconstruct(
        self,
        x: np.ndarray,
        mask: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> np.ndarray:
        """Decoded reconstruction `x_hat`. No extra masking. Returns `(N, n_features)`."""
        out = np.empty((x.shape[0], self.n_features), dtype=np.float32)
        bs = batch_size or self.batch_size
        n_steps = (x.shape[0] + bs - 1) // bs
        rng = np.random.default_rng(0)
        it = self._iter_batches(x, mask, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="reconstruct")
        pos = 0
        for xb, mb, cb in it:
            xh = _reconstruct_fn(self.model, xb, mb, cb)
            n = xh.shape[0]
            out[pos:pos + n] = np.asarray(xh)
            pos += n
        return out

    def anomaly_score(
        self,
        x: np.ndarray,
        mask: np.ndarray,
        c: np.ndarray,
        *,
        score_type: str = "mse",
        mc_seeds: int = 0,
        mask_ratio: Optional[float] = None,
        batch_size: Optional[int] = None,
        progress: bool = False,
        seed: int = 0,
    ) -> np.ndarray:
        """Per-source reconstruction-error score. Higher = more anomalous.

        `score_type`:
          - 'mse'         : raw squared error (μ − x)². Always available.
          - 'mahalanobis' : (μ − x)² / σ². Down-weights features the model
                            knows are noisy. Requires `loss_type='nll'`.
          - 'nll'         : 0.5 * ((μ − x)² / σ² + log σ²). The training
                            objective itself. Requires `loss_type='nll'`.

        Two MC modes:
          - `mc_seeds=0` (default): scored on originally observed positions,
            no extra masking. Cheapest; risk is that the bottleneck partly
            identity-bypasses high-information features.
          - `mc_seeds>0`: averages over `mc_seeds` random masking draws at
            `mask_ratio` (defaults to the training value). Scored only on
            artificially hidden, originally observed positions. Closer to
            the training objective; more robust to identity bypass.

        Returns `(N,)` numpy array.
        """
        if score_type not in ("mse", "mahalanobis", "nll"):
            raise ValueError(
                f"score_type must be 'mse', 'mahalanobis', or 'nll'; got {score_type!r}"
            )
        if score_type != "mse" and self.loss_type != "nll":
            raise ValueError(
                f"score_type={score_type!r} requires loss_type='nll' "
                f"(model has loss_type={self.loss_type!r}); σ² is not learned."
            )

        bs = batch_size or self.batch_size
        n_steps = (x.shape[0] + bs - 1) // bs
        rng = np.random.default_rng(0)
        out = np.empty(x.shape[0], dtype=np.float32)

        it = self._iter_batches(x, mask, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="anomaly")

        if mc_seeds <= 0:
            pos = 0
            for xb, mb, cb in it:
                eb = _self_recon_err_fn(self.model, xb, mb, cb, score_type)
                n = eb.shape[0]
                out[pos:pos + n] = np.asarray(eb)
                pos += n
        else:
            mr = self.mask_ratio if mask_ratio is None else mask_ratio
            key = jr.PRNGKey(seed)
            pos = 0
            for xb, mb, cb in it:
                acc = jnp.zeros(xb.shape[0], dtype=jnp.float32)
                for _ in range(mc_seeds):
                    key, sub = jr.split(key)
                    acc = acc + _mc_recon_err_fn(
                        self.model, xb, mb, cb, sub, mr, score_type
                    )
                eb = acc / mc_seeds
                n = eb.shape[0]
                out[pos:pos + n] = np.asarray(eb)
                pos += n
        return out

    def per_feature_anomaly(
        self,
        x: np.ndarray,
        mask: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = 256,
        progress: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Leave-one-out per-feature reconstruction.

        For each (source, feature j), masks only feature j and reconstructs
        it from the rest of the row. Returns:
            x_hat_loo : (N, F)  reconstructed value at the masked position
            sq_err    : (N, F)  (x - x_hat_loo)**2

        Cost is F forward passes per source, so use on a subset (this is the
        per-source diagnostic equivalent of the flow's `per_feature_anomaly_score`).
        """
        n, f = x.shape
        out = np.empty((n, f), dtype=np.float32)
        bs = batch_size or 256
        n_steps = (n + bs - 1) // bs
        rng = np.random.default_rng(0)
        it = self._iter_batches(x, mask, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="per_feature")
        pos = 0
        for xb, mb, cb in it:
            xh = _per_feature_loo_fn(self.model, xb, mb, cb)
            k = xh.shape[0]
            out[pos:pos + k] = np.asarray(xh)
            pos += k
        sq = (x - out) ** 2
        return out, sq.astype(np.float32)

    # --------------------------------------------------------------- I/O
    def save(self, path: str | Path):
        eqx.tree_serialise_leaves(str(path), self.model)

    def load(self, path: str | Path):
        """Load weights into the current architecture (in place). Returns self."""
        self.model = eqx.tree_deserialise_leaves(str(path), self.model)
        return self
