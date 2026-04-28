"""Conditional rational-quadratic spline coupling flow.

Single-class sklearn-style API:

    model = ConditionalFlow(n_features=51, n_cond=56)
    model.train(x_train, c_train, val_x=x_val, val_c=c_val)
    log_prob = model.evaluate(x, c)         # per-source log p(x|c)
    z        = model.transform(x, c)        # latent (Gaussianized)
    model.save("ckpt.eqx")
    model = ConditionalFlow(n_features=51, n_cond=56).load("ckpt.eqx")
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

from flowjax.flows import coupling_flow
from flowjax.bijections import RationalQuadraticSpline
from flowjax.distributions import Normal

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


# Module-level jitted kernels — argument signatures stay constant across
# instances, so JIT cache is reused.
@eqx.filter_jit
def _loss_fn(flow, x, c):
    return -flow.log_prob(x, condition=c).mean()


@eqx.filter_jit
def _step_fn(flow, opt_state, x, c, optimizer):
    loss, grads = eqx.filter_value_and_grad(_loss_fn)(flow, x, c)
    updates, opt_state = optimizer.update(grads, opt_state, flow)
    flow = eqx.apply_updates(flow, updates)
    return flow, opt_state, loss


@eqx.filter_jit
def _log_prob_fn(flow, x, c):
    return flow.log_prob(x, condition=c)


@eqx.filter_jit
def _latent_fn(flow, x, c):
    def single(x_i, c_i):
        z, _ = flow.bijection.inverse_and_log_det(x_i, condition=c_i)
        return z
    return jax.vmap(single)(x, c)


@eqx.filter_jit
def _split_logp_fn(flow, x, c):
    def single(x_i, c_i):
        z, log_det = flow.bijection.inverse_and_log_det(x_i, condition=c_i)
        log_pz = flow.base_dist.log_prob(z)
        return log_pz, log_det
    return jax.vmap(single)(x, c)


@eqx.filter_jit
def _score_fn(flow, x, c):
    def single_lp(x_i, c_i):
        return flow.log_prob(x_i, condition=c_i)
    return jax.vmap(jax.grad(single_lp))(x, c)


@eqx.filter_jit
def _jac_fn(flow, x, c):
    def single(x_i, c_i):
        z, _ = flow.bijection.inverse_and_log_det(x_i, condition=c_i)
        return z
    return jax.vmap(jax.jacfwd(single))(x, c)


@eqx.filter_jit
def _walk_fn(flow, x, c, u, alphas):
    def walk_one(x_i, c_i):
        z_i, _ = flow.bijection.inverse_and_log_det(x_i, condition=c_i)
        zs = z_i[None, :] + alphas[:, None] * u[None, :]
        xs = jax.vmap(lambda z: flow.bijection.transform(z, condition=c_i))(zs)
        return xs, zs
    return jax.vmap(walk_one)(x, c)


@eqx.filter_jit
def _local_axis_fn(flow, x, c, u):
    def single(x_i, c_i):
        z, _ = flow.bijection.inverse_and_log_det(x_i, condition=c_i)
        return z
    Js = jax.vmap(jax.jacfwd(single))(x, c)
    vs = jax.vmap(lambda J: jnp.linalg.solve(J, u))(Js)
    norms = jnp.linalg.norm(vs, axis=1)
    safe = jnp.where(norms > 0, norms, 1.0)
    return vs / safe[:, None], norms


class ConditionalFlow:
    """Conditional normalizing flow over `n_features` dims, conditioned on `n_cond`."""

    def __init__(
        self,
        n_features: int,
        n_cond: int,
        *,
        # architecture
        n_coupling_layers: int = 10,
        n_bins: int = 8,
        spline_bound: float = 10.0,
        nn_width: int = 256,
        nn_depth: int = 3,
        # training
        batch_size: int = 8192,
        learning_rate: float = 1e-3,
        n_epochs: int = 5,
        grad_clip: float = 1.0,
        seed: int = 0,
    ):
        self.n_features = n_features
        self.n_cond = n_cond
        self.n_coupling_layers = n_coupling_layers
        self.n_bins = n_bins
        self.spline_bound = spline_bound
        self.nn_width = nn_width
        self.nn_depth = nn_depth
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.n_epochs = n_epochs
        self.grad_clip = grad_clip
        self.seed = seed

        self.history: dict = {}
        self._build()

    def _build(self):
        key = jr.PRNGKey(self.seed)
        base = Normal(jnp.zeros(self.n_features))
        transformer = RationalQuadraticSpline(knots=self.n_bins, interval=self.spline_bound)
        self.flow = coupling_flow(
            key=key,
            base_dist=base,
            transformer=transformer,
            cond_dim=self.n_cond,
            flow_layers=self.n_coupling_layers,
            nn_width=self.nn_width,
            nn_depth=self.nn_depth,
        )

    # ----------------------------------------------------------------- batching
    def _iter_batches(self, x, c, *, shuffle, rng, batch_size=None):
        bs = batch_size or self.batch_size
        idx = np.arange(x.shape[0])
        if shuffle:
            idx = rng.permutation(idx)
        for i in range(0, len(idx), bs):
            b = idx[i:i + bs]
            yield jnp.asarray(x[b]), jnp.asarray(c[b])

    # ----------------------------------------------------------------- training
    def train(
        self,
        x: np.ndarray,
        c: np.ndarray,
        *,
        val_x: Optional[np.ndarray] = None,
        val_c: Optional[np.ndarray] = None,
        log_every: int = 0,
        progress: bool = True,
        checkpoint_dir: Optional[str | Path] = None,
    ):
        """Fit the flow on `(x, c)`. Optionally tracks validation loss.

        Populates `self.history`:
            step_loss          per-step training loss
            epoch_train_loss   mean training loss per epoch
            epoch_val_loss     mean validation loss per epoch (if val provided)
            epoch_time         seconds per epoch
            best_val_loss      lowest val loss observed
            best_epoch         epoch index of best val loss
        """
        if x.shape[0] != c.shape[0]:
            raise ValueError(f"x ({x.shape[0]}) and c ({c.shape[0]}) row count mismatch")
        if x.shape[1] != self.n_features:
            raise ValueError(f"x has {x.shape[1]} features, expected {self.n_features}")
        if c.shape[1] != self.n_cond:
            raise ValueError(f"c has {c.shape[1]} dims, expected {self.n_cond}")

        optimizer = optax.chain(
            optax.clip_by_global_norm(self.grad_clip),
            optax.adam(self.learning_rate),
        )
        opt_state = optimizer.init(eqx.filter(self.flow, eqx.is_inexact_array))

        rng = np.random.default_rng(self.seed + 1)

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
        has_val = val_x is not None and val_c is not None

        n_steps_per_epoch = (x.shape[0] + self.batch_size - 1) // self.batch_size
        use_tqdm = progress and tqdm is not None
        epoch_bar = tqdm(range(self.n_epochs), desc="epoch") if use_tqdm else range(self.n_epochs)

        for epoch in epoch_bar:
            t0 = time.time()

            train_losses = []
            step_iter = self._iter_batches(x, c, shuffle=True, rng=rng)
            if use_tqdm:
                step_iter = tqdm(step_iter, total=n_steps_per_epoch,
                                 desc=f"epoch {epoch}", leave=False)

            for step_i, (xb, cb) in enumerate(step_iter):
                self.flow, opt_state, loss = _step_fn(
                    self.flow, opt_state, xb, cb, optimizer
                )
                lf = float(loss)
                train_losses.append(lf)
                history["step_loss"].append(lf)

                if use_tqdm:
                    # running mean over the last 50 steps for a less twitchy display
                    recent = train_losses[-50:]
                    step_iter.set_postfix(loss=f"{lf:.4f}",
                                          mean=f"{np.mean(recent):.4f}")
                if log_every and step_i % log_every == 0:
                    print(f"  epoch {epoch} step {step_i:6d} loss={lf:.4f}")

            # validation
            val_l = float("nan")
            if has_val:
                val_losses = [
                    float(_loss_fn(self.flow, xb, cb))
                    for xb, cb in self._iter_batches(
                        val_x, val_c, shuffle=False, rng=rng
                    )
                ]
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
                self.save(ckpt_dir / f"flow_epoch{epoch:03d}.eqx")

        self.history = history
        return self

    # --------------------------------------------------------------- inference
    def evaluate(
        self,
        x: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> np.ndarray:
        """Per-source log p(x | c). Returns `(N,)` numpy array."""
        out = np.empty(x.shape[0], dtype=np.float32)
        pos = 0
        bs = batch_size or self.batch_size
        n_steps = (x.shape[0] + bs - 1) // bs
        rng = np.random.default_rng(0)  # unused; shuffle=False
        it = self._iter_batches(x, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="evaluate")
        for xb, cb in it:
            lp = _log_prob_fn(self.flow, xb, cb)
            n = lp.shape[0]
            out[pos:pos + n] = np.asarray(lp)
            pos += n
        return out

    def transform(
        self,
        x: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> np.ndarray:
        """Latent (Gaussianized) representation z. Returns `(N, n_features)`."""
        out = np.empty((x.shape[0], self.n_features), dtype=np.float32)
        pos = 0
        bs = batch_size or self.batch_size
        n_steps = (x.shape[0] + bs - 1) // bs
        rng = np.random.default_rng(0)
        it = self._iter_batches(x, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="transform")
        for xb, cb in it:
            zb = _latent_fn(self.flow, xb, cb)
            n = zb.shape[0]
            out[pos:pos + n] = np.asarray(zb)
            pos += n
        return out

    # --------------------------------------------------------------- jacobian diagnostics
    def split_log_prob(
        self,
        x: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Decompose log p(x|c) = log_pz(z) + log|det J_{x→z}|. Returns two `(N,)` arrays."""
        n = x.shape[0]
        log_pz = np.empty(n, dtype=np.float32)
        log_det = np.empty(n, dtype=np.float32)
        bs = batch_size or self.batch_size
        rng = np.random.default_rng(0)
        n_steps = (n + bs - 1) // bs
        it = self._iter_batches(x, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="split_log_prob")
        pos = 0
        for xb, cb in it:
            lpz, ld = _split_logp_fn(self.flow, xb, cb)
            k = lpz.shape[0]
            log_pz[pos:pos + k] = np.asarray(lpz)
            log_det[pos:pos + k] = np.asarray(ld)
            pos += k
        return log_pz, log_det

    def score(
        self,
        x: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> np.ndarray:
        """∇_x log p(x|c). Returns `(N, n_features)`. Norm is a useful outlier signal."""
        n, f = x.shape
        out = np.empty((n, f), dtype=np.float32)
        bs = batch_size or self.batch_size
        rng = np.random.default_rng(0)
        n_steps = (n + bs - 1) // bs
        it = self._iter_batches(x, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="score")
        pos = 0
        for xb, cb in it:
            s = _score_fn(self.flow, xb, cb)
            k = s.shape[0]
            out[pos:pos + k] = np.asarray(s)
            pos += k
        return out

    def jacobian(
        self,
        x: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> np.ndarray:
        """∂z/∂x. Returns `(N, n_features, n_features)`. Memory ~ N·F²·4 bytes."""
        n, f = x.shape
        bs = batch_size or 256
        out = np.empty((n, f, f), dtype=np.float32)
        rng = np.random.default_rng(0)
        n_steps = (n + bs - 1) // bs
        it = self._iter_batches(x, c, shuffle=False, rng=rng, batch_size=bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="jacobian")
        pos = 0
        for xb, cb in it:
            J = _jac_fn(self.flow, xb, cb)
            k = J.shape[0]
            out[pos:pos + k] = np.asarray(J)
            pos += k
        return out

    def walk_axis(
        self,
        x: np.ndarray,
        c: np.ndarray,
        u: np.ndarray,
        alphas: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Walk seed sources along unit vector `u` in latent space.

        For each seed `i`: z_i = f(x_i, c_i), then z_i + α u for each α, mapped
        back to data space (with the seed's own conditioning held fixed).

        Args:
            x:      (N, F) seed sources.
            c:      (N, n_cond) seed conditioning.
            u:      (F,) latent-space direction (typically unit-norm).
            alphas: (K,) step sizes.

        Returns:
            xs: (N, K, F) walked sources in data space.
            zs: (N, K, F) walked points in latent space.
        """
        u = jnp.asarray(u, dtype=jnp.float32)
        alphas = jnp.asarray(alphas, dtype=jnp.float32)
        xs, zs = _walk_fn(self.flow, jnp.asarray(x), jnp.asarray(c), u, alphas)
        return np.asarray(xs), np.asarray(zs)

    def local_axis_in_x(
        self,
        x: np.ndarray,
        c: np.ndarray,
        u: np.ndarray,
        *,
        batch_size: Optional[int] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Local feature-space direction corresponding to latent vector `u`.

        v(x_i) = J(x_i)^{-1} u, returned normalized to unit length, plus the
        pre-normalization length ‖J^{-1} u‖. Small length means a unit step in
        z costs little in x — the source sits on a steep part of the axis;
        large length means the axis is flat there.

        Returns:
            v:    (N, F) unit vectors.
            cost: (N,)   ‖J(x_i)^{-1} u‖ before normalization.
        """
        u = jnp.asarray(u, dtype=jnp.float32)
        n, f = x.shape
        bs = batch_size or 256
        v = np.empty((n, f), dtype=np.float32)
        cost = np.empty(n, dtype=np.float32)
        rng = np.random.default_rng(0)
        pos = 0
        for xb, cb in self._iter_batches(x, c, shuffle=False, rng=rng, batch_size=bs):
            vb, nb = _local_axis_fn(self.flow, xb, cb, u)
            k = vb.shape[0]
            v[pos:pos + k] = np.asarray(vb)
            cost[pos:pos + k] = np.asarray(nb)
            pos += k
        return v, cost

    # --------------------------------------------------------------- I/O
    def save(self, path: str | Path):
        eqx.tree_serialise_leaves(str(path), self.flow)

    def load(self, path: str | Path):
        """Load weights into the current architecture (in place). Returns self."""
        self.flow = eqx.tree_deserialise_leaves(str(path), self.flow)
        return self
