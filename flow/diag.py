import jax
import jax.numpy as jnp
import numpy as np
from functools import partial


def _unwrap(flow):
    """Accept either a flowjax distribution or a ConditionalFlow wrapper."""
    return getattr(flow, "flow", flow)


# ============================================================
# Version 1: Gradient-based attribution
# ============================================================
# Score: a_ij = x_ij * d log p(x_i | c_i) / d x_ij
# This is "feature i's contribution to the log-likelihood" in the
# input*gradient sense. Sign tells you whether the feature is
# pushing likelihood up or down; magnitude tells you how much.
# ============================================================

def make_gradient_attribution_fn(flow):
    """
    flow: a flowjax distribution (or ConditionalFlow wrapper) with .log_prob(x, condition=c)
    Returns a jit-compiled function mapping (x_batch, c_batch) -> (N, D) attributions.
    """
    flow = _unwrap(flow)

    def log_prob_single(x, c):
        # log_prob expects batched input; add/remove batch dim
        return flow.log_prob(x[None, :], condition=c[None, :])[0]

    # Gradient w.r.t. x for a single example
    grad_fn = jax.grad(log_prob_single, argnums=0)

    # Vectorize over the batch
    batched_grad = jax.vmap(grad_fn, in_axes=(0, 0))

    @jax.jit
    def attribution(x_batch, c_batch):
        grads = batched_grad(x_batch, c_batch)        # (N, D)
        return x_batch * grads                          # input * gradient

    return attribution


# ============================================================
# Version 2: Exact per-dimension decomposition via the latent
# ============================================================
# log p(x|c) = log p_z(z) + log|det J_{f}(x;c)|
#            = sum_i [-0.5 * z_i^2 - 0.5 log(2 pi)] + log|det J|
#
# The base-distribution term decomposes cleanly per *latent* dim.
# We attribute each latent's contribution back to input features
# via the Jacobian of the inverse map (or via input*grad on z_i).
#
# This gives you two arrays:
#   - latent_contrib: (N, D)  contribution of each z_i to log p_z
#   - feature_attrib: (N, D)  attribution of those contributions
#                              back to input features
# ============================================================

def make_latent_attribution_fn(flow):
    """
    Returns functions to:
      (a) compute per-latent log-density contributions (returns (z, contribs))
      (b) attribute them back to input features via Jacobian
    """
    flow = _unwrap(flow)

    def forward_single(x, c):
        # Map x -> z under the conditional flow.
        return flow.bijection.transform(x, condition=c)

    def latent_logp_contrib_single(x, c):
        z = forward_single(x, c)
        # Per-dim contribution to log p_z under standard normal base
        # -0.5 * z_i^2 - 0.5 * log(2*pi)
        contribs = -0.5 * z**2 - 0.5 * jnp.log(2.0 * jnp.pi)  # (D,)
        return z, contribs

    batched_latent_contrib = jax.jit(jax.vmap(latent_logp_contrib_single, in_axes=(0, 0)))

    # Jacobian of forward map: dz/dx, shape (D_z, D_x)
    jac_fn = jax.jit(jax.vmap(jax.jacfwd(forward_single, argnums=0), in_axes=(0, 0)))

    return batched_latent_contrib, jac_fn


def feature_attribution_from_latent(x_batch, c_batch, latent_contrib_fn, jac_fn):
    """
    Attribute each latent's log-density contribution back to input features.

    For each galaxy:
      latent_contrib[i]  = -0.5 z_i^2 - 0.5 log(2 pi)
      d(latent_contrib_i)/d(x_j) = -z_i * (dz_i/dx_j)
    We return per-feature attribution = sum_i d(latent_contrib_i)/d(x_j) * x_j
    (input*grad on the *base-distribution* term only; Jacobian-determinant term excluded).

    Returns:
      feature_attrib: (N, D) array, contribution of each input feature
                       to the base-distribution part of log p(x|c)
      latent_contrib: (N, D) per-latent log-density values
      jacobians:      (N, D, D) full Jacobian dz/dx for diagnostics
    """
    z_batch, latent_contrib = latent_contrib_fn(x_batch, c_batch)  # (N, D), (N, D)
    jacobians = jac_fn(x_batch, c_batch)                            # (N, D_z, D_x)

    # d(log p_z)/d(x_j) = sum_i (-z_i) * dz_i/dx_j
    dlogpz_dx = -jnp.einsum('ni,nij->nj', z_batch, jacobians)

    feature_attrib = x_batch * dlogpz_dx                            # input * grad
    return feature_attrib, latent_contrib, jacobians


# ============================================================
# Version 3: Per-feature conditional z-score (Monte Carlo)
# ============================================================
# For each source i and feature j:
#   mu_ij = E[x_j | c_i],  sd_ij = std[x_j | c_i]   (estimated by MC)
#   z_feat_ij = (x_ij - mu_ij) / sd_ij
#
# This answers "is feature j unusual on source i, given its observing
# conditions" without the input-magnitude bias of input*grad. Sign tells
# you whether x_j is above or below typical; magnitude is in conditional
# sigmas, so it is comparable across features.
#
# Cost: n_samples flow.sample calls per source. Use on a subset, not all
# 100M sources.
# ============================================================

def per_feature_anomaly_score(flow, x, c, *, n_samples=200, batch_size=64, seed=0):
    """
    Per-feature conditional z-score via Monte Carlo from p(x | c).

    Args:
        flow:        ConditionalFlow wrapper or flowjax distribution.
        x:           (N, D) features for the subset to inspect.
        c:           (N, n_cond) conditioning (already standardized like training).
        n_samples:   MC draws per source for estimating mu/sd of p(x_j | c_i).
        batch_size:  number of sources processed per JIT batch (memory control).
        seed:        PRNG seed.

    Returns:
        z_feat: (N, D) per-feature conditional z-scores.
        mu:     (N, D) conditional means E[x_j | c_i].
        sd:     (N, D) conditional stds  std[x_j | c_i].
    """
    flow = _unwrap(flow)
    x = np.asarray(x)
    c = np.asarray(c)
    N, D = x.shape

    def sample_one(key, c_i):
        return flow.sample(key, (n_samples,), condition=c_i)   # (K, D)

    batched_sample = jax.jit(jax.vmap(sample_one, in_axes=(0, 0)))

    rng = jax.random.PRNGKey(seed)
    mu = np.empty((N, D), dtype=np.float32)
    sd = np.empty((N, D), dtype=np.float32)

    for start in range(0, N, batch_size):
        end = min(start + batch_size, N)
        rng, sub = jax.random.split(rng)
        keys = jax.random.split(sub, end - start)
        samples = batched_sample(keys, jnp.asarray(c[start:end]))  # (B, K, D)
        mu[start:end] = np.asarray(samples.mean(axis=1))
        sd[start:end] = np.asarray(samples.std(axis=1))

    z_feat = ((x - mu) / np.where(sd > 0, sd, 1.0)).astype(np.float32)
    return z_feat, mu, sd


def report_top_anomalous_features(z_feat, x, mu, sd, feature_names,
                                  source_idx, top_k=5):
    """
    Print, for each source in `source_idx`, the top-k features ranked by
    |conditional z-score|, alongside the actual x value and the conditional
    typical range (mu ± sd).

    Args:
        z_feat:        (N, D) from per_feature_anomaly_score.
        x:             (N, D) features.
        mu, sd:        (N, D) from per_feature_anomaly_score.
        feature_names: list of D strings.
        source_idx:    iterable of row indices to inspect.
        top_k:         number of features to print per source.
    """
    name_w = max(8, max(len(n) for n in feature_names))
    for i in source_idx:
        order = np.argsort(-np.abs(z_feat[i]))[:top_k]
        print(f"\nsource {i}:")
        print(f"  {'feature':<{name_w}}  {'x':>9}  {'mu(c)':>9}  "
              f"{'sd(c)':>8}  {'z_feat':>8}")
        for j in order:
            print(f"  {feature_names[j]:<{name_w}}  "
                  f"{x[i, j]:>9.3f}  {mu[i, j]:>9.3f}  "
                  f"{sd[i, j]:>8.3f}  {z_feat[i, j]:>+8.2f}")
