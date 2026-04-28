import jax
import jax.numpy as jnp
import numpy as np
from functools import partial
import matplotlib.pyplot as plt


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

    rng = jax.random.key(seed)
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


# ============================================================
# Latent-space axes: late→early (or any labeled gradient)
# ============================================================
# Workflow:
#   z = model.transform(x, c)
#   u, s = fit_axis_in_z(z, label)                # direction in z
#   plot_axis_vs_label(s, label)                  # is the axis real?
#   xs, _ = model.walk_axis(seeds_x, seeds_c, u, alphas)
#   plot_axis_walk(model, seeds_x, seeds_c, u, alphas, feature_names)
#   v, cost = model.local_axis_in_x(x, c, u)      # per-source x-direction
#   plot_axis_loadings(u, v.mean(0), feature_names)
# ============================================================

def fit_axis_in_z(z, label, *, method="ridge", alpha=1.0):
    """Find a unit direction `u` in latent space that tracks `label`.

    Args:
        z:      (N, D) latent codes (from `model.transform`).
        label:  (N,) continuous (e.g. T-type) or binary (0/1).
        method: "ridge" for continuous labels, "lda" for binary.
        alpha:  ridge regularization (ignored for lda).

    Returns:
        u:     (D,) unit vector in z-space.
        score: (N,) projection `z @ u`.
    """
    z = np.asarray(z)
    label = np.asarray(label)
    if method == "ridge":
        from sklearn.linear_model import Ridge
        u = Ridge(alpha=alpha).fit(z, label).coef_.astype(np.float32)
    elif method == "lda":
        from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
        u = LinearDiscriminantAnalysis().fit(z, label).coef_[0].astype(np.float32)
    else:
        raise ValueError(f"unknown method {method!r}")
    n = float(np.linalg.norm(u))
    if n == 0:
        raise ValueError("fitted axis has zero norm")
    u = u / n
    return u, z @ u


def plot_axis_vs_label(score, label, ax=None, *, label_name="label", gridsize=50):
    """Hexbin of axis projection vs label, annotated with Pearson r."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    label = np.asarray(label)
    score = np.asarray(score)
    ax.hexbin(label, score, gridsize=gridsize, mincnt=1, cmap="viridis")
    r = float(np.corrcoef(label, score)[0, 1])
    ax.set_xlabel(label_name)
    ax.set_ylabel("z · u  (axis projection)")
    ax.set_title(f"axis vs {label_name}   (r = {r:.3f})")
    return ax


def plot_axis_walk(model, x_seeds, c_seeds, u, alphas, feature_names,
                   ax=None, *, kind="heatmap", center=True):
    """Walk seeds along `u` in z, plot mean walked features vs α.

    `kind="heatmap"` shows a feature × α image of mean walked x (recommended
    for ~50 features). `kind="lines"` overlays one line per feature.
    `center=True` subtracts the value at α≈0 so the panel shows change rather
    than absolute level.
    """
    xs, _ = model.walk_axis(x_seeds, c_seeds, u, alphas)        # (N, K, F)
    mean_x = xs.mean(axis=0)                                    # (K, F)
    alphas = np.asarray(alphas)
    if center:
        zero_idx = int(np.argmin(np.abs(alphas)))
        mean_x = mean_x - mean_x[zero_idx]

    if ax is None:
        _, ax = plt.subplots(figsize=(8, max(4, 0.18 * len(feature_names))))

    if kind == "heatmap":
        vmax = float(np.max(np.abs(mean_x))) or 1.0
        im = ax.imshow(
            mean_x.T, aspect="auto", cmap="RdBu_r",
            vmin=-vmax, vmax=vmax,
            extent=[alphas[0], alphas[-1], len(feature_names) - 0.5, -0.5],
        )
        ax.set_yticks(np.arange(len(feature_names)))
        ax.set_yticklabels(feature_names, fontsize=7)
        ax.set_xlabel("α  (steps along u in z)")
        ax.set_title("synthetic walk: Δ⟨x⟩ along axis" if center else "⟨x⟩ along walk")
        plt.colorbar(im, ax=ax, label=("Δ⟨x⟩" if center else "⟨x⟩"))
    elif kind == "lines":
        for j, name in enumerate(feature_names):
            ax.plot(alphas, mean_x[:, j], label=name, lw=1)
        ax.axhline(0, color="k", lw=0.5, alpha=0.5)
        ax.set_xlabel("α")
        ax.set_ylabel("⟨x⟩ along walk" + (" (centered)" if center else ""))
        if len(feature_names) <= 20:
            ax.legend(fontsize=7, ncol=2)
    else:
        raise ValueError(f"unknown kind {kind!r}")
    return ax


def plot_axis_loadings(u, v_mean, feature_names, axes=None, *, top_k=15):
    """Two panels: latent-space `u` (all dims) and top-|v_mean| feature loadings.

    `v_mean` is the mean over sources of `model.local_axis_in_x(...)[0]` —
    i.e. the average local feature-space direction corresponding to `u`.
    """
    u = np.asarray(u)
    v_mean = np.asarray(v_mean)
    top_k = min(top_k, len(feature_names))
    if axes is None:
        _, axes = plt.subplots(1, 2, figsize=(12, max(4, 0.3 * top_k)))

    axes[0].bar(np.arange(len(u)), u)
    axes[0].axhline(0, color="k", lw=0.5)
    axes[0].set_xlabel("latent dim")
    axes[0].set_ylabel("u_i")
    axes[0].set_title("axis direction in z")

    order = np.argsort(-np.abs(v_mean))[:top_k]
    axes[1].barh(np.arange(top_k), v_mean[order][::-1])
    axes[1].set_yticks(np.arange(top_k))
    axes[1].set_yticklabels([feature_names[j] for j in order][::-1], fontsize=8)
    axes[1].axvline(0, color="k", lw=0.5)
    axes[1].set_xlabel("⟨v_j⟩  (mean local x-direction)")
    axes[1].set_title(f"top {top_k} features driving the axis")
    return axes


def plot_logp_decomposition(log_pz, log_det, ax=None, *, label=None,
                            label_name="label", gridsize=60):
    """Scatter of `log p_z` vs `log|det J|` from `model.split_log_prob`.

    Without a label, hexbin density. With a label, per-source scatter colored
    by it — useful for spotting where late- vs early-type galaxies live in
    the (latent term, Jacobian term) plane.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    log_pz = np.asarray(log_pz)
    log_det = np.asarray(log_det)
    if label is None:
        ax.hexbin(log_pz, log_det, gridsize=gridsize, mincnt=1, cmap="viridis")
    else:
        sc = ax.scatter(log_pz, log_det, c=np.asarray(label), s=4,
                        cmap="coolwarm", alpha=0.6)
        plt.colorbar(sc, ax=ax, label=label_name)
    ax.set_xlabel("log p_z(z)   (latent term)")
    ax.set_ylabel("log|det J_{x→z}|   (Jacobian term)")
    ax.set_title("log p(x|c) decomposition")
    return ax


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
