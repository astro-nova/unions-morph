"""Masked Autoencoder with FiLM conditioning.

Single-class sklearn-style API:

    model = MaskedAutoencoder(n_features=51, n_cond=5)
    model.train(x_tr, m_tr, c_tr, val_x=x_va, val_mask=m_va, val_c=c_va)
    z         = model.encode(x, m, c)              # (N, latent_dim)
    x_hat     = model.reconstruct(x, m, c)         # (N, n_features)
    x_hat     = model.decode(z, c)                 # (N, n_features) from z + c
    score     = model.anomaly_score(x, m, c)       # (N,)
    per_feat  = model.per_feature_anomaly(x, m, c) # (N, n_features)
    ks, err   = model.truncation_curve(x, m, c)    # recon error vs z[:k]
    model.save("ckpt.eqx")
    MaskedAutoencoder(n_features=51, n_cond=5).load("ckpt.eqx")

Architecture:
    conditioner : c                  → e      shared non-linear embedding of c
    encoder     : (x, m_orig, m_extra, e) → z  FiLM-MLP, depth `n_hidden_layers`
    decoder     : (z, e)             → x_hat  FiLM-MLP, mirrored

Two kinds of "masking", kept disentangled (input side vs. output side are
orthogonal — a position can be fed-but-not-scored, scored-but-corrupted,
both, or neither):

    Q = m_orig   intrinsic mask  — catalog quality flags. The value was never
                 measured, so it is *never* a loss target (no ground truth).
    R = m_extra  random pretext  — positions we artificially hide each step as
                 regularization, forcing the latent to learn cross-feature
                 structure. Drawn from the *observed* positions only, so
                 R ∩ Q = ∅ by construction.

Input side: masked entries get distinct learned tokens (`mask_token` for R,
`fill_token` for Q), and a single bit-encoded mask channel
`mask_code = m_orig + 2·m_extra ∈ {0=observed, 1=Q, 2=R}` is appended, so the
encoder input is `[x_filled, mask_code]` of width `2 * n_features`.

Output side: the decoder reconstructs all `n_features`. The objective here is
compression (50→latent_dim) with masking as regularization, so the loss runs
over *every observed* position O = ¬Q — the held-out R positions and the
still-visible positions weighted independently (R up-weighted by default).
Q is not scored by default, but `quality_loss_weight > 0` folds the catalog-
imputed values at Q positions into the objective (the same flag also weights
Q in the headline anomaly score). See `_loss_fn` / `_weighted_recon_loss`.

The conditioning `c` (typically the 5 standardized observing-condition vars,
brightness first) is first embedded by a small shared MLP `conditioner(c) = e`,
then every hidden activation is FiLM-modulated:

    h ← (1 + γ(e)) ⊙ Linear(h) + β(e)

i.e. a per-channel affine transform whose scale/shift come from a *non-linear*
function of c (the conditioner), so γ/β can bend with S/N and resolution
rather than depending on them linearly. At init the per-layer FiLM heads are
small (γ ≈ 0 / β ≈ 0), so the network starts as a plain MLP and learns the
conditioning gradually.

Optional adversarial disentanglement (`adv_lambda > 0`):
    A small MLP head predicts (a subset of) `c` from `z`. A gradient-reversal
    layer flips the sign of gradients flowing back into the encoder, so the
    encoder is pushed to make `z` *unpredictive* of the targeted condition
    vars while the head simultaneously learns to be a good predictor.

      total_loss = recon_loss + ‖adv_head(GRL(z, λ)) − c[:, target_dims]‖²

    The strength λ is ramped linearly from 0 to `adv_lambda` over
    `adv_warmup_steps` global steps, so the head becomes a competent predictor
    before its gradient meaningfully shapes `z`. `adv_target_dims` chooses which
    columns of `c` to push out of `z` (default: all of `c`). Set `adv_lambda=0`
    (the default) to disable entirely.

Optional latent-structure regularizers (both off by default, independent of
each other and of the adversary). Plain reconstruction pins down only the
latent *subspace*, not its basis — any invertible linear map of `z` can be
absorbed by the decoder's first layer — so these break that symmetry to make
individual dims meaningful:

    Nested dropout (`nested_dropout_p > 0`): each training row draws
    k ~ Geometric(p) clipped to `latent_dim`, and the decoder sees only
    `z[:k]` (rest zeroed). Reconstruction-critical information is forced into
    the early dims, so `z[:k]` is the best available k-dim code for every k —
    PCA's ordered-components property, with a fully nonlinear map. Late dims
    going quiet is dimensionality selection, not collapse. Diagnose with
    `truncation_curve(...)`.

    Jacobian orthogonality (`jac_ortho_lambda > 0`): penalizes the mean
    squared cosine between columns of the decoder Jacobian ∂x̂/∂z (each column
    is the feature-space direction latent i moves) on `jac_ortho_batch` rows
    per step, pushing different dims toward locally *disjoint* effects —
    curvilinear-orthogonal coordinates on the data manifold. Deliberately the
    Gram/cosine form rather than a Hessian-diagonal penalty: it does not force
    the decoder to be additive, so nonlinear interactions between latents stay
    representable. Evaluated at the codes the decoder actually consumed
    (truncated ones when nested dropout is on).

The adversary always consumes the *full* (untruncated) code, so `c`-leakage
is scrubbed from every dim regardless of truncation.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Optional, Sequence

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
class _FiLMConditioner(eqx.Module):
    """Shared non-linear embedding of the conditioning vector `c`.

    A small MLP (`c → embed_dim`, `layers` stages of Linear with GELU between)
    whose output `e` feeds every FiLM head. Being non-linear in `c`, it lets the
    per-layer γ/β bend with S/N and resolution rather than depend on them
    linearly. Computed once per source and shared by encoder + decoder.
    """
    layers: list

    def __init__(self, n_cond, embed_dim, layers, *, key):
        keys = jr.split(key, layers)
        mods, d = [], n_cond
        for i in range(layers):
            mods.append(eqx.nn.Linear(d, embed_dim, key=keys[i]))
            d = embed_dim
        self.layers = mods

    def __call__(self, c):
        h = c
        for i, lin in enumerate(self.layers):
            h = lin(h)
            if i < len(self.layers) - 1:
                h = jax.nn.gelu(h)
        return h


class FiLMBlock(eqx.Module):
    """Linear → FiLM(e) → GELU, where `e` is the conditioner's embedding of `c`.

    The FiLM head init is small so γ≈0, β≈0 → block ≈ Linear at start.
    """
    linear: eqx.nn.Linear
    film:   eqx.nn.Linear
    out_dim: int = eqx.field(static=True)

    def __init__(self, in_dim: int, out_dim: int, embed_dim: int, *, key):
        k1, k2 = jr.split(key, 2)
        self.linear = eqx.nn.Linear(in_dim, out_dim, key=k1)
        film = eqx.nn.Linear(embed_dim, 2 * out_dim, key=k2)
        # Shrink FiLM init so γ≈0, β≈0 → starts as identity modulation.
        film = eqx.tree_at(lambda l: l.weight, film, film.weight * 0.01)
        film = eqx.tree_at(lambda l: l.bias,   film, film.bias   * 0.0)
        self.film = film
        self.out_dim = out_dim

    def __call__(self, h, e):
        h = self.linear(h)
        gb = self.film(e)
        gamma, beta = gb[:self.out_dim], gb[self.out_dim:]
        h = (1.0 + gamma) * h + beta
        return jax.nn.gelu(h)


@jax.custom_vjp
def _grad_reverse(x, lam):
    """Identity in forward; in backward, multiplies the cotangent by `-lam`.

    Used to push the encoder to produce `z` that's hard to predict `c` from,
    while still letting the adversary head learn to be a good predictor.
    """
    return x

def _gr_fwd(x, lam):
    return x, lam

def _gr_bwd(res, g):
    lam = res
    return (-lam * g, jnp.zeros_like(lam))

_grad_reverse.defvjp(_gr_fwd, _gr_bwd)


class _AdvHead(eqx.Module):
    """Small MLP that predicts (a subset of) `c` from `z`. No FiLM; sees only `z`."""
    layers: list

    def __init__(self, latent: int, hidden: int, n_target: int, *, key):
        keys = jr.split(key, 3)
        self.layers = [
            eqx.nn.Linear(latent, hidden, key=keys[0]),
            eqx.nn.Linear(hidden, hidden, key=keys[1]),
            eqx.nn.Linear(hidden, n_target, key=keys[2]),
        ]

    def __call__(self, z):
        h = jax.nn.gelu(self.layers[0](z))
        h = jax.nn.gelu(self.layers[1](h))
        return self.layers[2](h)


class _Encoder(eqx.Module):
    blocks: list
    head: eqx.nn.Linear
    mask_token: jax.Array     # (F,) learned fill for randomly-masked (R) positions
    fill_token: jax.Array     # (F,) learned fill for intrinsically-masked (Q) positions
    n_features: int = eqx.field(static=True)

    def __init__(self, n_features, hidden, depth, latent, embed_dim, *, key):
        keys = jr.split(key, depth + 1)
        in_dim = 2 * n_features            # [x_filled | mask_code]
        blocks, d = [], in_dim
        for i in range(depth):
            blocks.append(FiLMBlock(d, hidden, embed_dim, key=keys[i]))
            d = hidden
        self.blocks = blocks
        self.head = eqx.nn.Linear(hidden, latent, key=keys[-1])
        # Start both tokens at 0 (≈ observed-mean for the ~zero-centred scaling),
        # then learn distinct fills for "predict me" (R) vs. "absent" (Q).
        self.mask_token = jnp.zeros(n_features)
        self.fill_token = jnp.zeros(n_features)
        self.n_features = n_features

    def __call__(self, x, m_orig, m_extra, e):
        # m_orig and m_extra are disjoint (R drawn from observed only), so the
        # three regions below never overlap and mask_code ∈ {0, 1, 2} exactly.
        visible = (1.0 - m_orig) * (1.0 - m_extra)
        x_filled = x * visible + self.fill_token * m_orig + self.mask_token * m_extra
        mask_code = m_orig + 2.0 * m_extra          # 0=observed, 1=Q, 2=R
        h = jnp.concatenate([x_filled, mask_code], axis=-1)
        for blk in self.blocks:
            h = blk(h, e)
        return self.head(h)


class _Decoder(eqx.Module):
    blocks: list
    head: eqx.nn.Linear

    def __init__(self, latent, hidden, depth, out_dim, embed_dim, *, key):
        keys = jr.split(key, depth + 1)
        blocks, d = [], latent
        for i in range(depth):
            blocks.append(FiLMBlock(d, hidden, embed_dim, key=keys[i]))
            d = hidden
        self.blocks = blocks
        self.head = eqx.nn.Linear(hidden, out_dim, key=keys[-1])

    def __call__(self, z, e):
        h = z
        for blk in self.blocks:
            h = blk(h, e)
        return self.head(h)


class _MAEModule(eqx.Module):
    conditioner: _FiLMConditioner
    encoder: _Encoder
    decoder: _Decoder
    adv_head: Optional[_AdvHead]
    adv_lambda: float = eqx.field(static=True)
    adv_warmup_steps: int = eqx.field(static=True)
    adv_target_dims: tuple = eqx.field(static=True)


# ---------------------------------------------------------------- jitted kernels
def _encode_single(model, x_i, m_orig_i, m_extra_i, c_i):
    e = model.conditioner(c_i)
    return model.encoder(x_i, m_orig_i, m_extra_i, e)


def _forward_single(model, x_i, m_orig_i, m_extra_i, c_i):
    e = model.conditioner(c_i)
    z = model.encoder(x_i, m_orig_i, m_extra_i, e)
    return model.decoder(z, e), z


def _draw_extra_mask(key, m_orig, mask_ratio):
    """Sample R ~ Bernoulli(mask_ratio) restricted to observed (¬Q) positions.

    Returns a {0,1} array disjoint from `m_orig` by construction.
    """
    observed = 1.0 - m_orig
    return ((jr.uniform(key, m_orig.shape) < mask_ratio) * (observed > 0)).astype(m_orig.dtype)


def _weighted_recon_loss(x_hat, x, m_orig, m_extra, masked_w, unmasked_w, quality_w):
    """Weighted MSE over reconstruction targets.

    Per position the weight is:
        Q  (m_orig=1)                          → quality_w  (imputed fill target)
        R  (m_extra=1, observed)               → masked_w   (held-out target)
        visible observed (¬Q, ¬R)              → unmasked_w (compression target)

    `R` is up-weighted by default (it carries the denoising signal); setting
    `unmasked_w=0` recovers the canonical held-out-only MAE objective.

    `quality_w` defaults to 0, so Q positions are normally *not* scored — the
    value was never measured, so the only "target" is whatever the catalog
    imputed into `x` (and the encoder saw `fill_token` there, not the value).
    Set it > 0 to fold those imputed values into the objective and measure the
    effect. Q and R are disjoint, so the three regions never overlap.
    """
    observed = 1.0 - m_orig
    w = observed * (unmasked_w * (1.0 - m_extra) + masked_w * m_extra) + quality_w * m_orig
    sq = (x_hat - x) ** 2
    return (sq * w).sum() / jnp.clip(w.sum(), 1.0, None)


def _draw_keep_mask(key, n_rows, latent_dim, p):
    """Per-row nested-dropout keep mask: k ~ Geometric(p), clipped to
    [1, latent_dim]; dims [0, k) are kept, the rest zeroed. E[k] ≈ 1/p, and
    rows whose draw exceeds latent_dim keep the full code.
    """
    u = jr.uniform(key, (n_rows,), minval=jnp.finfo(jnp.float32).tiny)
    k = 1 + jnp.floor(jnp.log(u) / jnp.log1p(-p)).astype(jnp.int32)
    k = jnp.clip(k, 1, latent_dim)
    return (jnp.arange(latent_dim)[None, :] < k[:, None]).astype(jnp.float32)


def _recon_and_z(model, x, m_orig, c, key, mask_ratio, masked_w, unmasked_w, quality_w, nd_p):
    """Shared kernel: artificial masking + forward + weighted recon loss.

    Returns (recon_loss, z, z_dec, e): the full code `z` (the adversary always
    sees every dim), the code the decoder actually consumed `z_dec`
    (nested-dropout-truncated when `nd_p > 0`, else `z` itself), and the
    conditioning embedding `e` (reused by the Jacobian-orthogonality penalty).
    """
    k_mask, k_nd = jr.split(key)
    m_extra = _draw_extra_mask(k_mask, m_orig, mask_ratio)
    e = jax.vmap(model.conditioner)(c)
    z = jax.vmap(model.encoder)(x, m_orig, m_extra, e)
    if nd_p > 0.0:
        keep = _draw_keep_mask(k_nd, z.shape[0], z.shape[1], nd_p)
        z_dec = z * keep
    else:
        z_dec = z
    x_hat = jax.vmap(model.decoder)(z_dec, e)
    recon = _weighted_recon_loss(x_hat, x, m_orig, m_extra, masked_w, unmasked_w, quality_w)
    return recon, z, z_dec, e


def _jac_ortho_loss(model, z, e):
    """Mean squared cosine between decoder-Jacobian columns.

    Column i of J = ∂x̂/∂z is the feature-space direction latent i moves at
    this point; the penalty drives distinct latents toward locally orthogonal
    (non-overlapping) effects without constraining each effect's shape or
    forbidding interactions. Cosines make it scale-invariant — it can't be
    gamed by shrinking the decoder's response — and near-dead dims (e.g.
    truncated away by nested dropout) contribute ≈ 0.
    """
    def per_row(z_i, e_i):
        J = jax.jacfwd(lambda zz: model.decoder(zz, e_i))(z_i)   # (F, L)
        G = J.T @ J                                              # (L, L)
        d = jnp.diag(G)
        cos = G / jnp.sqrt(jnp.clip(d[:, None] * d[None, :], 1e-12, None))
        L = G.shape[0]
        off = cos * (1.0 - jnp.eye(L))
        return (off ** 2).sum() / (L * (L - 1))
    return jax.vmap(per_row)(z, e).mean()


@eqx.filter_jit
def _loss_fn(model, x, m_orig, c, key, mask_ratio, masked_w, unmasked_w, quality_w,
             nd_p, jac_lambda, jac_batch, step):
    """Total loss = recon + adv (via gradient reversal) + jac_lambda · jac.

    `step` is a JAX scalar (current global step) used to ramp the adversarial
    weight linearly from 0 to `model.adv_lambda` over `model.adv_warmup_steps`.
    `nd_p`, `jac_lambda`, `jac_batch` arrive as Python scalars (static under
    filter_jit), so the disabled branches compile away entirely.
    Returns (total_loss, (recon_loss, adv_loss, jac_loss)) for `has_aux=True`.
    """
    recon, z, z_dec, e = _recon_and_z(
        model, x, m_orig, c, key, mask_ratio, masked_w, unmasked_w, quality_w, nd_p
    )
    total = recon

    if model.adv_head is not None:
        target_idx = jnp.array(model.adv_target_dims, dtype=jnp.int32)
        c_target = c[:, target_idx]
        warmup = jnp.maximum(model.adv_warmup_steps, 1).astype(jnp.float32)
        ramp = jnp.minimum(1.0, step.astype(jnp.float32) / warmup)
        lam = model.adv_lambda * ramp
        z_grl = _grad_reverse(z, lam)
        c_pred = jax.vmap(model.adv_head)(z_grl)
        adv = jnp.mean((c_pred - c_target) ** 2)
        total = total + adv
    else:
        adv = jnp.array(0.0)

    if jac_lambda > 0.0:
        n = min(int(jac_batch), z_dec.shape[0])
        jac = _jac_ortho_loss(model, z_dec[:n], e[:n])
        total = total + jac_lambda * jac
    else:
        jac = jnp.array(0.0)
    return total, (recon, adv, jac)


@eqx.filter_jit
def _val_loss_fn(model, x, m_orig, c, key, mask_ratio, masked_w, unmasked_w, quality_w, nd_p):
    """Validation: recon loss only (the metric we actually care about).

    Runs under the same nested-dropout regime as training so the number is
    comparable to `epoch_recon_loss` — which also means val losses are not
    comparable across runs with different `nested_dropout_p`.
    """
    recon, _, _, _ = _recon_and_z(
        model, x, m_orig, c, key, mask_ratio, masked_w, unmasked_w, quality_w, nd_p
    )
    return recon


@eqx.filter_jit
def _step_fn(model, opt_state, x, m, c, key, mask_ratio, masked_w, unmasked_w, quality_w,
             nd_p, jac_lambda, jac_batch, step, optimizer):
    (loss, aux), grads = eqx.filter_value_and_grad(_loss_fn, has_aux=True)(
        model, x, m, c, key, mask_ratio, masked_w, unmasked_w, quality_w,
        nd_p, jac_lambda, jac_batch, step
    )
    updates, opt_state = optimizer.update(grads, opt_state, model)
    model = eqx.apply_updates(model, updates)
    return model, opt_state, loss, aux


def _zeros_like_mask(m):
    return jnp.zeros_like(m)


@eqx.filter_jit
def _encode_fn(model, x, m, c):
    """Latent code using the catalog mask only (no extra masking)."""
    def per_row(x_i, m_i, c_i):
        return _encode_single(model, x_i, m_i, _zeros_like_mask(m_i), c_i)
    return jax.vmap(per_row)(x, m, c)


@eqx.filter_jit
def _decode_fn(model, z, c):
    """Decode latent codes `z` straight through the decoder (encoder skipped)."""
    def per_row(z_i, c_i):
        e = model.conditioner(c_i)
        return model.decoder(z_i, e)
    return jax.vmap(per_row)(z, c)


@eqx.filter_jit
def _trunc_err_fn(model, x, m, c, keep):
    """Per-source observed-MSE reconstructing from the truncated code z·keep.

    `keep` is a (latent_dim,) 0/1 vector shared by all rows; encoding uses the
    catalog mask only (no extra masking). Backs `truncation_curve`.
    """
    def per_row(x_i, m_i, c_i):
        e = model.conditioner(c_i)
        z = model.encoder(x_i, m_i, _zeros_like_mask(m_i), e)
        x_hat = model.decoder(z * keep, e)
        sq = (x_hat - x_i) ** 2
        w = 1.0 - m_i
        return (sq * w).sum() / jnp.clip(w.sum(), 1.0, None)
    return jax.vmap(per_row)(x, m, c)


@eqx.filter_jit
def _reconstruct_fn(model, x, m, c):
    def per_row(x_i, m_i, c_i):
        x_hat, _ = _forward_single(model, x_i, m_i, _zeros_like_mask(m_i), c_i)
        return x_hat
    return jax.vmap(per_row)(x, m, c)


@eqx.filter_jit
def _self_recon_err_fn(model, x, m, c, quality_w):
    """Per-source MSE over reconstruction targets (no extra masking).

    This is the headline anomaly score: it measures how well the 50→latent_dim
    code reconstructs every feature that was actually measured. Observed
    positions get weight 1; intrinsically-masked (Q) positions get weight
    `quality_w` (default 0, i.e. excluded — their only target is the catalog
    imputation). Set `quality_w>0` to also score the imputed Q positions.
    """
    def per_row(x_i, m_i, c_i):
        x_hat, _ = _forward_single(model, x_i, m_i, _zeros_like_mask(m_i), c_i)
        sq = (x_hat - x_i) ** 2
        w = (1.0 - m_i) + quality_w * m_i
        return (sq * w).sum() / jnp.clip(w.sum(), 1.0, None)
    return jax.vmap(per_row)(x, m, c)


@eqx.filter_jit
def _mc_recon_err_fn(model, x, m_orig, c, key, mask_ratio):
    """Per-source MSE on artificially hidden (R) positions, one MC draw.

    R is drawn from observed positions only, so every scored position has a
    real target.
    """
    m_extra = _draw_extra_mask(key, m_orig, mask_ratio)

    def per_row(x_i, m_orig_i, m_extra_i, c_i):
        x_hat, _ = _forward_single(model, x_i, m_orig_i, m_extra_i, c_i)
        sq = (x_hat - x_i) ** 2
        return (sq * m_extra_i).sum() / jnp.clip(m_extra_i.sum(), 1.0, None)

    return jax.vmap(per_row)(x, m_orig, m_extra, c)


@eqx.filter_jit
def _per_feature_loo_fn(model, x, m, c):
    """For each (source, feature j), reconstruct x_j after masking only j.

    Feature j is held out as an R position (learned mask token + R code); the
    catalog mask Q is left in place. If j is already a Q position it stays Q (no
    double-masking), and its reconstruction is the ordinary forward pass.

    Returns (N, F) reconstructed values at the masked diagonal.
    """
    F = x.shape[1]
    eye = jnp.eye(F, dtype=x.dtype)

    def per_row(x_i, m_orig_i, c_i):
        observed = 1.0 - m_orig_i

        def mask_j(j_onehot):
            m_extra = j_onehot * observed          # only newly hide it if observed
            x_hat, _ = _forward_single(model, x_i, m_orig_i, m_extra, c_i)
            return x_hat
        x_hats = jax.vmap(mask_j)(eye)             # (F, F)
        return jnp.diagonal(x_hats)                # (F,)

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
        cond_embed_dim: int = 16,
        cond_embed_layers: int = 2,
        # training
        mask_ratio: float = 0.5,
        masked_loss_weight: float = 1.0,
        unmasked_loss_weight: float = 0.3,
        quality_loss_weight: float = 0.0,
        batch_size: int = 8192,
        learning_rate: float = 1e-3,
        n_epochs: int = 10,
        grad_clip: float = 1.0,
        # latent structure regularizers (each off when 0)
        nested_dropout_p: float = 0.0,
        jac_ortho_lambda: float = 0.0,
        jac_ortho_batch: int = 128,
        # adversarial disentanglement (off when adv_lambda == 0)
        adv_lambda: float = 0.0,
        adv_target_dims: Optional[Sequence[int]] = None,
        adv_hidden: int = 64,
        adv_warmup_steps: int = 1000,
        seed: int = 0,
    ):
        if not (0.0 <= nested_dropout_p < 1.0):
            raise ValueError(f"nested_dropout_p must be in [0, 1), got {nested_dropout_p!r}")
        if jac_ortho_lambda < 0:
            raise ValueError(f"jac_ortho_lambda must be ≥ 0, got {jac_ortho_lambda!r}")
        if jac_ortho_lambda > 0 and latent_dim < 2:
            raise ValueError("jac_ortho_lambda > 0 needs latent_dim ≥ 2 (nothing to orthogonalize)")
        if jac_ortho_batch < 1:
            raise ValueError(f"jac_ortho_batch must be ≥ 1, got {jac_ortho_batch!r}")
        if adv_lambda < 0:
            raise ValueError(f"adv_lambda must be ≥ 0, got {adv_lambda!r}")
        if adv_warmup_steps < 0:
            raise ValueError(f"adv_warmup_steps must be ≥ 0, got {adv_warmup_steps!r}")
        if adv_target_dims is not None:
            for d in adv_target_dims:
                if not (0 <= int(d) < n_cond):
                    raise ValueError(
                        f"adv_target_dims contains {d!r}, must be in [0, {n_cond})"
                    )
        self.n_features = n_features
        self.n_cond = n_cond
        self.latent_dim = latent_dim
        self.hidden_dim = hidden_dim
        self.n_hidden_layers = n_hidden_layers
        self.cond_embed_dim = cond_embed_dim
        self.cond_embed_layers = cond_embed_layers
        self.mask_ratio = mask_ratio
        self.masked_loss_weight = masked_loss_weight
        self.unmasked_loss_weight = unmasked_loss_weight
        self.quality_loss_weight = float(quality_loss_weight)
        self.nested_dropout_p = float(nested_dropout_p)
        self.jac_ortho_lambda = float(jac_ortho_lambda)
        self.jac_ortho_batch = int(jac_ortho_batch)
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.n_epochs = n_epochs
        self.grad_clip = grad_clip
        self.adv_lambda = float(adv_lambda)
        self.adv_target_dims = (
            None if adv_target_dims is None else tuple(int(d) for d in adv_target_dims)
        )
        self.adv_hidden = int(adv_hidden)
        self.adv_warmup_steps = int(adv_warmup_steps)
        self.seed = seed

        self.history: dict = {}
        self._build()

    def _build(self):
        key = jr.PRNGKey(self.seed)
        kc, ke, kd, ka = jr.split(key, 4)
        conditioner = _FiLMConditioner(
            n_cond=self.n_cond,
            embed_dim=self.cond_embed_dim,
            layers=self.cond_embed_layers,
            key=kc,
        )
        encoder = _Encoder(
            n_features=self.n_features,        # input is [x_filled | mask_code], width 2F
            hidden=self.hidden_dim,
            depth=self.n_hidden_layers,
            latent=self.latent_dim,
            embed_dim=self.cond_embed_dim,
            key=ke,
        )
        decoder = _Decoder(
            latent=self.latent_dim,
            hidden=self.hidden_dim,
            depth=self.n_hidden_layers,
            out_dim=self.n_features,
            embed_dim=self.cond_embed_dim,
            key=kd,
        )

        if self.adv_lambda > 0.0:
            target_dims = (
                tuple(range(self.n_cond))
                if self.adv_target_dims is None
                else self.adv_target_dims
            )
            adv_head = _AdvHead(
                latent=self.latent_dim,
                hidden=self.adv_hidden,
                n_target=len(target_dims),
                key=ka,
            )
        else:
            target_dims = ()
            adv_head = None
        self._adv_target_dims_resolved = target_dims  # for inspection

        self.model = _MAEModule(
            conditioner=conditioner,
            encoder=encoder,
            decoder=decoder,
            adv_head=adv_head,
            adv_lambda=self.adv_lambda,
            adv_warmup_steps=self.adv_warmup_steps,
            adv_target_dims=target_dims,
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

        `mask` is the intrinsic / quality-flag mask Q (1 = imputed). The random
        pretext mask R is drawn internally each step from observed positions.

        Populates `self.history`:
            step_loss          per-step total training loss
            step_recon_loss    per-step reconstruction loss
            step_adv_loss      per-step adversarial loss (0 when adv off)
            step_jac_loss      per-step Jacobian-orthogonality penalty (0 when off)
            epoch_train_loss   mean total training loss per epoch
            epoch_recon_loss   mean reconstruction loss per epoch
            epoch_adv_loss     mean adversarial loss per epoch
            epoch_jac_loss     mean Jacobian-orthogonality penalty per epoch
            epoch_val_loss     mean validation (recon) loss per epoch (if val given)
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

        adv_on = self.model.adv_head is not None
        jac_on = self.jac_ortho_lambda > 0.0
        history = {
            "step_loss": [],
            "step_recon_loss": [],
            "step_adv_loss": [],
            "step_jac_loss": [],
            "epoch_train_loss": [],
            "epoch_recon_loss": [],
            "epoch_adv_loss": [],
            "epoch_jac_loss": [],
            "epoch_val_loss": [],
            "epoch_time": [],
            "best_val_loss": float("inf"),
            "best_epoch": -1,
        }
        has_val = val_x is not None and val_mask is not None and val_c is not None

        n_steps_per_epoch = (x.shape[0] + self.batch_size - 1) // self.batch_size
        use_tqdm = progress and tqdm is not None
        epoch_bar = tqdm(range(self.n_epochs), desc="epoch") if use_tqdm else range(self.n_epochs)

        global_step = 0
        for epoch in epoch_bar:
            t0 = time.time()

            train_losses = []
            recon_losses = []
            adv_losses = []
            jac_losses = []
            step_iter = self._iter_batches(x, mask, c, shuffle=True, rng=rng)
            if use_tqdm:
                step_iter = tqdm(step_iter, total=n_steps_per_epoch,
                                 desc=f"epoch {epoch}", leave=False)

            for step_i, (xb, mb, cb) in enumerate(step_iter):
                key, sub = jr.split(key)
                step_arr = jnp.asarray(global_step, dtype=jnp.float32)
                self.model, opt_state, loss, aux = _step_fn(
                    self.model, opt_state, xb, mb, cb, sub,
                    self.mask_ratio, self.masked_loss_weight, self.unmasked_loss_weight,
                    self.quality_loss_weight, self.nested_dropout_p,
                    self.jac_ortho_lambda, self.jac_ortho_batch, step_arr, optimizer,
                )
                global_step += 1
                lf = float(loss)
                rf = float(aux[0])
                af = float(aux[1])
                jf = float(aux[2])
                train_losses.append(lf)
                recon_losses.append(rf)
                adv_losses.append(af)
                jac_losses.append(jf)
                history["step_loss"].append(lf)
                history["step_recon_loss"].append(rf)
                history["step_adv_loss"].append(af)
                history["step_jac_loss"].append(jf)

                if use_tqdm:
                    recent = train_losses[-50:]
                    postfix = {"loss": f"{lf:.4f}", "mean": f"{np.mean(recent):.4f}"}
                    if adv_on:
                        postfix["recon"] = f"{rf:.4f}"
                        postfix["adv"] = f"{af:.4f}"
                    if jac_on:
                        postfix["jac"] = f"{jf:.4f}"
                    step_iter.set_postfix(**postfix)
                if log_every and step_i % log_every == 0:
                    extra = f" recon={rf:.4f} adv={af:.4f}" if adv_on else ""
                    if jac_on:
                        extra += f" jac={jf:.4f}"
                    print(f"  epoch {epoch} step {step_i:6d} loss={lf:.4f}{extra}")

            val_l = float("nan")
            if has_val:
                key, vsub = jr.split(key)
                val_losses = []
                for xb, mb, cb in self._iter_batches(val_x, val_mask, val_c,
                                                     shuffle=False, rng=rng):
                    vsub, kk = jr.split(vsub)
                    val_losses.append(float(_val_loss_fn(
                        self.model, xb, mb, cb, kk,
                        self.mask_ratio, self.masked_loss_weight, self.unmasked_loss_weight,
                        self.quality_loss_weight, self.nested_dropout_p,
                    )))
                val_l = float(np.mean(val_losses))
                history["epoch_val_loss"].append(val_l)
                if val_l < history["best_val_loss"]:
                    history["best_val_loss"] = val_l
                    history["best_epoch"] = epoch

            train_l = float(np.mean(train_losses))
            recon_l = float(np.mean(recon_losses))
            adv_l = float(np.mean(adv_losses))
            jac_l = float(np.mean(jac_losses))
            dt = time.time() - t0
            history["epoch_train_loss"].append(train_l)
            history["epoch_recon_loss"].append(recon_l)
            history["epoch_adv_loss"].append(adv_l)
            history["epoch_jac_loss"].append(jac_l)
            history["epoch_time"].append(dt)

            if use_tqdm:
                postfix = {"train": f"{train_l:.4f}", "dt": f"{dt:.1f}s"}
                if adv_on:
                    postfix["recon"] = f"{recon_l:.4f}"
                    postfix["adv"] = f"{adv_l:.4f}"
                if jac_on:
                    postfix["jac"] = f"{jac_l:.4f}"
                if has_val:
                    postfix["val"] = f"{val_l:.4f}"
                epoch_bar.set_postfix(**postfix)
            else:
                msg = f"epoch {epoch}: train={train_l:.4f}"
                if adv_on:
                    msg += f"  recon={recon_l:.4f}  adv={adv_l:.4f}"
                if jac_on:
                    msg += f"  jac={jac_l:.4f}"
                if has_val:
                    msg += f"  val={val_l:.4f}"
                msg += f"  ({dt:.1f}s)"
                print(msg)

            if ckpt_dir is not None:
                self.save(ckpt_dir / f"mae_epoch{epoch:03d}.eqx")

        self.history = history
        return self

    # --------------------------------------------------------------- baseline loss
    def loss_of(
        self,
        y_pred: np.ndarray,
        x: np.ndarray,
        mask: np.ndarray,
        *,
        mask_ratio: Optional[float] = None,
        n_draws: int = 10,
        seed: int = 0,
    ) -> float:
        """Training-objective reconstruction loss for an arbitrary `y_pred`.

        Mirrors `_weighted_recon_loss`: scores every observed position O = ¬Q,
        with held-out R positions weighted by `masked_loss_weight` and still-
        visible positions by `unmasked_loss_weight`, averaged over `n_draws`
        random R realizations. Numbers are directly comparable to
        `history['epoch_val_loss']` (which is recon-only).

        `y_pred` may be `(F,)` (broadcasts over rows — useful for median/mean
        baselines) or `(N, F)`.
        """
        x = np.asarray(x, dtype=np.float32)
        mask = np.asarray(mask, dtype=np.float32)
        if x.shape != mask.shape:
            raise ValueError(f"x {x.shape} and mask {mask.shape} must match")
        if x.shape[1] != self.n_features:
            raise ValueError(f"x has {x.shape[1]} features, expected {self.n_features}")

        y_pred = np.broadcast_to(np.asarray(y_pred, dtype=np.float32), x.shape)
        sq = (y_pred - x) ** 2

        mr = self.mask_ratio if mask_ratio is None else mask_ratio
        mw, uw, qw = self.masked_loss_weight, self.unmasked_loss_weight, self.quality_loss_weight
        rng = np.random.default_rng(seed)
        observed = 1.0 - mask
        losses = np.empty(n_draws, dtype=np.float64)
        for k in range(n_draws):
            m_extra = (rng.random(x.shape) < mr).astype(np.float32) * observed
            w = observed * (uw * (1.0 - m_extra) + mw * m_extra) + qw * mask
            denom = max(w.sum(), 1.0)
            losses[k] = (sq * w).sum() / denom
        return float(losses.mean())

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

    def decode(
        self,
        z: np.ndarray,
        c: np.ndarray,
        *,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> np.ndarray:
        """Reconstruct `x_hat` from latent codes `z` and conditioning `c`.

        The inverse-direction counterpart of `encode`: runs the decoder only
        (encoder skipped), so it works on arbitrary points in latent space —
        e.g. interpolations, samples, or `z` perturbed for counterfactuals.
        Conditioning `c` is FiLM-applied exactly as in the full forward pass.

        Returns `(N, n_features)`.
        """
        z = np.asarray(z, dtype=np.float32)
        c = np.asarray(c, dtype=np.float32)
        if z.ndim != 2 or z.shape[1] != self.latent_dim:
            raise ValueError(
                f"z must be (N, {self.latent_dim}), got {z.shape}"
            )
        if c.ndim != 2 or c.shape[1] != self.n_cond:
            raise ValueError(f"c must be (N, {self.n_cond}), got {c.shape}")
        if z.shape[0] != c.shape[0]:
            raise ValueError(
                f"z, c row-count mismatch: {z.shape[0]} / {c.shape[0]}"
            )

        out = np.empty((z.shape[0], self.n_features), dtype=np.float32)
        bs = batch_size or self.batch_size
        n_steps = (z.shape[0] + bs - 1) // bs
        it = range(0, z.shape[0], bs)
        if progress and tqdm is not None:
            it = tqdm(it, total=n_steps, desc="decode")
        for i in it:
            zb = jnp.asarray(z[i:i + bs])
            cb = jnp.asarray(c[i:i + bs])
            xh = _decode_fn(self.model, zb, cb)
            out[i:i + xh.shape[0]] = np.asarray(xh)
        return out

    def truncation_curve(
        self,
        x: np.ndarray,
        mask: np.ndarray,
        c: np.ndarray,
        *,
        ks: Optional[Sequence[int]] = None,
        batch_size: Optional[int] = None,
        progress: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Reconstruction error vs. latent truncation width k.

        For each k, encodes with the catalog mask only, zeroes `z[k:]`, decodes,
        and averages the per-source MSE over observed positions. The nested-
        dropout diagnostic: with `nested_dropout_p > 0` the curve shows how much
        each successive latent dim buys, and its knee is the effective
        dimensionality of the code. Without nested dropout the decoder never saw
        truncated codes, so expect the error to stay high until nearly all dims
        are kept — the flat-then-cliff shape is itself evidence the basis is
        arbitrary.

        Cost is one encode+decode pass per k; run on a subset (~1e5 rows).

        Returns `(ks, err)`: the truncation widths `(K,)` and the mean error
        at each `(K,)`.
        """
        ks = (np.arange(1, self.latent_dim + 1) if ks is None
              else np.asarray(list(ks), dtype=int))
        if ks.ndim != 1 or ks.size == 0 or ks.min() < 1 or ks.max() > self.latent_dim:
            raise ValueError(f"ks must be 1-D with values in [1, {self.latent_dim}]")
        bs = batch_size or self.batch_size
        rng = np.random.default_rng(0)
        err = np.empty(ks.size, dtype=np.float64)
        it = tqdm(ks, desc="truncation") if progress and tqdm is not None else ks
        for j, k in enumerate(it):
            keep = jnp.asarray((np.arange(self.latent_dim) < k).astype(np.float32))
            tot, cnt = 0.0, 0
            for xb, mb, cb in self._iter_batches(x, mask, c, shuffle=False,
                                                 rng=rng, batch_size=bs):
                eb = _trunc_err_fn(self.model, xb, mb, cb, keep)
                tot += float(eb.sum())
                cnt += eb.shape[0]
            err[j] = tot / max(cnt, 1)
        return ks, err

    def anomaly_score(
        self,
        x: np.ndarray,
        mask: np.ndarray,
        c: np.ndarray,
        *,
        mc_seeds: int = 0,
        mask_ratio: Optional[float] = None,
        quality_weight: Optional[float] = None,
        batch_size: Optional[int] = None,
        progress: bool = False,
        seed: int = 0,
    ) -> np.ndarray:
        """Per-source reconstruction-error score. Higher = more anomalous.

        Two modes:
          - `mc_seeds=0` (default): self-reconstruction MSE over the observed
            features, no extra masking. This is the headline score — it asks how
            well the 50→latent_dim code reproduces everything that was measured.
            Intrinsically-masked (Q) positions are weighted by `quality_weight`
            (defaults to `self.quality_loss_weight`, i.e. excluded unless set);
            pass an explicit value to override per call.
          - `mc_seeds>0`: averages over `mc_seeds` random masking draws at
            `mask_ratio` (defaults to the training value). MSE measured only on
            artificially hidden (R) positions, which are drawn from the observed
            set so every one has a real target. A held-out diagnostic; more
            robust to identity bypass on highly informative features. (Q never
            enters this mode — R is drawn from observed positions only.)

        Returns `(N,)` numpy array.
        """
        qw = self.quality_loss_weight if quality_weight is None else float(quality_weight)
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
                eb = _self_recon_err_fn(self.model, xb, mb, cb, qw)
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
                    acc = acc + _mc_recon_err_fn(self.model, xb, mb, cb, sub, mr)
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

        For each (source, feature j), masks only feature j (as an R position)
        and reconstructs it from the rest of the row. Returns:
            x_hat_loo : (N, F)  reconstructed value at the masked position
            sq_err    : (N, F)  (x - x_hat_loo)**2

        Note `sq_err[:, j]` is only meaningful where feature j was observed
        (`mask[:, j] == 0`); on intrinsically-masked positions there is no
        ground truth.

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
    _CONFIG_KEYS = (
        "n_features", "n_cond",
        "latent_dim", "hidden_dim", "n_hidden_layers",
        "cond_embed_dim", "cond_embed_layers",
        "mask_ratio", "masked_loss_weight", "unmasked_loss_weight",
        "quality_loss_weight",
        "batch_size", "learning_rate", "n_epochs", "grad_clip",
        "nested_dropout_p", "jac_ortho_lambda", "jac_ortho_batch",
        "adv_lambda", "adv_target_dims", "adv_hidden", "adv_warmup_steps",
        "seed",
    )

    def _config_dict(self) -> dict:
        cfg = {k: getattr(self, k) for k in self._CONFIG_KEYS}
        if cfg["adv_target_dims"] is not None:
            cfg["adv_target_dims"] = list(cfg["adv_target_dims"])
        return cfg

    @staticmethod
    def _config_path(path: str | Path) -> Path:
        path = Path(path)
        return path.with_suffix(path.suffix + ".json")

    def save(self, path: str | Path):
        """Serialize weights to `path` and architecture/config to `<path>.json`.

        The sidecar lets `MaskedAutoencoder.from_checkpoint(path)` rebuild the
        exact architecture without the caller knowing the hyperparameters.
        """
        path = Path(path)
        eqx.tree_serialise_leaves(str(path), self.model)
        self._config_path(path).write_text(json.dumps(self._config_dict(), indent=2))

    def load(self, path: str | Path):
        """Load weights into the current architecture (in place). Returns self.

        The current `MaskedAutoencoder` must have been built with hyperparameters
        matching the saved checkpoint. If a `<path>.json` sidecar exists it's
        cross-checked and a clear error is raised on mismatch. If you don't
        already have the architecture in hand, prefer `from_checkpoint(path)`.
        """
        cfg_path = self._config_path(path)
        if cfg_path.exists():
            saved = json.loads(cfg_path.read_text())
            current = self._config_dict()
            arch_keys = ("n_features", "n_cond", "latent_dim", "hidden_dim",
                         "n_hidden_layers", "cond_embed_dim", "cond_embed_layers",
                         "adv_lambda", "adv_target_dims", "adv_hidden")
            mismatches = {k: (current.get(k), saved.get(k))
                          for k in arch_keys
                          if current.get(k) != saved.get(k)}
            if mismatches:
                lines = [f"  {k}: in-memory={cur!r}, on-disk={dsk!r}"
                         for k, (cur, dsk) in mismatches.items()]
                raise ValueError(
                    "Architecture mismatch between in-memory model and saved "
                    f"checkpoint at {path}:\n" + "\n".join(lines) +
                    "\nUse MaskedAutoencoder.from_checkpoint(path) to rebuild "
                    "from the saved config, or pass the matching hyperparameters."
                )
        self.model = eqx.tree_deserialise_leaves(str(path), self.model)
        return self

    @classmethod
    def from_checkpoint(cls, path: str | Path) -> "MaskedAutoencoder":
        """Construct a `MaskedAutoencoder` from `<path>.json` and load weights.

        Requires the sidecar written by `save(...)`. For older checkpoints
        without a sidecar, build the model manually with the matching
        hyperparameters and call `.load(path)`.
        """
        cfg_path = cls._config_path(path)
        if not cfg_path.exists():
            raise FileNotFoundError(
                f"No config sidecar at {cfg_path}. This file is needed to "
                "reconstruct the architecture. For checkpoints saved before "
                "the sidecar existed, build MaskedAutoencoder(...) manually "
                "with the matching hyperparameters and call .load(path)."
            )
        cfg = json.loads(cfg_path.read_text())
        cfg = {k: v for k, v in cfg.items() if k in cls._CONFIG_KEYS}
        model = cls(**cfg)
        model.load(path)
        return model
