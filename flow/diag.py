import jax
import jax.numpy as jnp
import numpy as np
from functools import partial

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
    flow: a flowjax distribution with .log_prob(x, condition=c)
    Returns a jit-compiled function mapping (x_batch, c_batch) -> (N, D) attributions.
    """
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
      (a) compute per-latent log-density contributions
      (b) attribute them back to input features via Jacobian
    """
    def forward_single(x, c):
        # Map x -> z under the conditional flow.
        # flowjax bijections expose a .transform method; adapt as needed
        # for your specific flow object (some use flow.bijection.transform).
        z = flow.bijection.transform(x, condition=c)
        return z

    def latent_logp_contrib_single(x, c):
        z = forward_single(x, c)
        # Per-dim contribution to log p_z under standard normal base
        # -0.5 * z_i^2 - 0.5 * log(2*pi)
        return -0.5 * z**2 - 0.5 * jnp.log(2.0 * jnp.pi)  # (D,)

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
    latent_contrib = latent_contrib_fn(x_batch, c_batch)  # (N, D)
    jacobians = jac_fn(x_batch, c_batch)                   # (N, D_z, D_x)

    # We need z to weight the rows of the Jacobian
    # Recompute z (cheap; or refactor to return it from latent_contrib_fn)
    def fwd(x, c):
        return flow.bijection.transform(x, condition=c)
    z_batch = jax.vmap(fwd)(x_batch, c_batch)              # (N, D)

    # d(log p_z)/d(x_j) = sum_i (-z_i) * dz_i/dx_j
    # Shape: (N, D_x) = einsum over i
    dlogpz_dx = -jnp.einsum('ni,nij->nj', z_batch, jacobians)

    feature_attrib = x_batch * dlogpz_dx                   # input * grad
    return feature_attrib, latent_contrib, jacobians

