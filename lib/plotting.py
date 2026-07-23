from matplotlib import pyplot as plt 
import numpy as np 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse, PathPatch
from matplotlib.path import Path
from scipy.stats import chi2, norm
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
 
 

def plot_cutouts(images, zp, pxscale, imsize=2, vmin=19, vmax=27):

    fig, axs = plt.subplots(1, len(images), figsize=(imsize*len(images), imsize))
    for ax, img in zip(axs, images):
        if img.dtype in [int, bool]:
            ax.imshow(img, origin='lower', cmap='gray', vmin=0, vmax=np.max(img)*0.8)
        else:
            ax.imshow(-2.5*np.log10(np.abs(img)/pxscale**2)+zp, origin='lower', cmap='gray', vmin=vmin, vmax=vmax) 
    
    for ax in axs:
        ax.axis('off')
    plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, top=1, bottom=0)
    return axs

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def make_smoothed_histogram(x, y, xlim, ylim, bins=50, sigma=1.5):
    """
    Build a smoothed 2D histogram of (x, y) over the given plot limits.

    Parameters
    ----------
    x, y : array-like
        Data arrays (NaNs/infs are masked out automatically).
    xlim, ylim : tuple
        (min, max) range for each axis — pass plt.xlim()/plt.ylim(), 
        or explicit bounds.
    bins : int
        Number of bins per axis.
    sigma : float
        Gaussian smoothing width in bin units.

    Returns
    -------
    xcenters, ycenters : 1D arrays
        Bin centers for each axis.
    H_smooth : 2D array, shape (nx, ny)
        Smoothed histogram, indexed [x_bin, y_bin].
    """
    x = np.asarray(x)
    y = np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]

    xedges = np.linspace(xlim[0], xlim[1], bins + 1)
    yedges = np.linspace(ylim[0], ylim[1], bins + 1)

    H, xedges, yedges = np.histogram2d(x, y, bins=[xedges, yedges])
    H_smooth = gaussian_filter(H, sigma=sigma)

    xcenters = 0.5 * (xedges[1:] + xedges[:-1])
    ycenters = 0.5 * (yedges[1:] + yedges[:-1])

    return xcenters, ycenters, H_smooth


def hdr_levels(H, fractions=(0.68, 0.95)):
    """
    Compute histogram-value levels enclosing the given cumulative 
    probability fractions (highest density region contours).
    """
    H_flat = np.sort(H.ravel())[::-1]  # descending
    cumsum = np.cumsum(H_flat)
    cumsum /= cumsum[-1]

    levels = []
    for f in fractions:
        idx = np.searchsorted(cumsum, f)
        idx = min(idx, len(H_flat) - 1)
        levels.append(H_flat[idx])

    levels = np.unique(levels)  # dedupe + sort ascending

    if len(levels) < 2:
        raise ValueError(
            f"Only {len(levels)} unique level(s) found for fractions={fractions}. "
            "Try more bins, less smoothing, or fewer/more separated fractions."
        )

    return levels

def plot_density_contour(x, y, xlim, ylim, bins=50, sigma=1.5,
                          fractions=(0.68, 0.95), colors='k', 
                          linewidths=1, ax=None, **contour_kwargs):
    """
    Convenience wrapper: build smoothed histogram + HDR levels, 
    draw contour on current axes.
    """
    ax = ax or plt.gca()
    
    xcenters, ycenters, H_smooth = make_smoothed_histogram(
        x, y, xlim, ylim, bins=bins, sigma=sigma
    )
    levels = hdr_levels(H_smooth, fractions=fractions)
    cs = ax.contour(xcenters, ycenters, H_smooth.T, levels=levels,
                      colors=colors, linewidths=linewidths, **contour_kwargs)
    return cs, levels




def plot_param_average_on_xygrid(x, y, z, xlim, ylim, gridsize=30, cmap='Spectral_r',
                                  reduce_C_function=np.mean, mincnt=3,
                                  alpha_range=(0.15, 1.0), alpha_clip_percentile=(0, 98),
                                  vmin=None, vmax=None, ax=None, 
                                  cbar_label=None, cbar_size='2%', cbar_pad=0,
                                  **hexbin_kwargs):
    """
    Hexbin-average z(x, y), with hex opacity scaled by point density,
    and a thin colorbar above the axes.
    """
    ax = ax or plt.gca()

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[mask], y[mask], z[mask]

    extent = (xlim[0], xlim[1], ylim[0], ylim[1])

    hb = ax.hexbin(x, y, C=z, reduce_C_function=reduce_C_function,
                    gridsize=gridsize, extent=extent, cmap=cmap,
                    mincnt=mincnt, vmin=vmin, vmax=vmax, **hexbin_kwargs)

    hb_counts = ax.hexbin(x, y, gridsize=gridsize, extent=extent, mincnt=mincnt)
    counts = hb_counts.get_array()
    hb_counts.remove()

    lo, hi = np.percentile(counts, alpha_clip_percentile)
    counts_clipped = np.clip(counts, lo, hi)
    if hi > lo:
        norm = (counts_clipped - lo) / (hi - lo)
    else:
        norm = np.ones_like(counts_clipped)
    alphas = alpha_range[0] + norm * (alpha_range[1] - alpha_range[0])

    hb.update_scalarmappable()
    facecolors = hb.get_facecolors()
    if facecolors.shape[0] != len(alphas):
        facecolors = hb.to_rgba(hb.get_array())
    facecolors[:, 3] = alphas
    hb.set_facecolors(facecolors)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # thin colorbar above the axes, aligned to axes width
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('top', size=cbar_size, pad=cbar_pad)
    cbar = plt.colorbar(hb, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cax.xaxis.set_label_position('top')
    if cbar_label:
        cbar.set_label(cbar_label)

    return hb, cbar



# ----------------------------------------------------------------------
# Gaussian mixture modelling
# ----------------------------------------------------------------------
def sigma_to_scale(n_sigma, sigma_mode="mahalanobis"):
    """
    Convert 'n sigma' into the Mahalanobis radius of the ellipse to draw.
 
    sigma_mode:
      'mahalanobis' : ellipse at Mahalanobis distance = n  (the usual "1-sigma
                      ellipse"; encloses 1 - exp(-n^2/2) of the probability,
                      i.e. 39% / 86% / 99% for n = 1, 2, 3).
      'coverage'    : ellipse enclosing the same probability as +/- n sigma in
                      1D (68.3% / 95.4% / 99.7%), i.e. chi2 with 2 dof.
    """
    n_sigma = np.atleast_1d(np.asarray(n_sigma, dtype=float))
    if sigma_mode == "mahalanobis":
        return n_sigma
    if sigma_mode == "coverage":
        p = 2.0 * norm.cdf(n_sigma) - 1.0
        return np.sqrt(chi2.ppf(p, df=2))
    raise ValueError("sigma_mode must be 'mahalanobis' or 'coverage'")
 
 
def cov_ellipse_params(cov):
    """Return (width, height, angle_deg) of the unit-Mahalanobis ellipse."""
    vals, vecs = np.linalg.eigh(cov)          # ascending eigenvalues
    vals = np.clip(vals, 0.0, None)
    major, minor = vecs[:, 1], vecs[:, 0]     # noqa: F841
    angle = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
    width = 2.0 * np.sqrt(vals[1])            # full axis length at 1 sigma
    height = 2.0 * np.sqrt(vals[0])
    return width, height, angle
 
 
def _ellipse_xy(mu, cov, s, npts=200):
    """Vertices of the ellipse at Mahalanobis radius s."""
    t = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=False)
    circ = np.column_stack([np.cos(t), np.sin(t)])
    try:
        L = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:                      # near-singular
        vals, vecs = np.linalg.eigh(cov)
        L = vecs @ np.diag(np.sqrt(np.clip(vals, 0, None)))
    return np.asarray(mu) + s * circ @ L.T
 
 
def _ring_patch(mu, cov, s_out, s_in=None, **kw):
    """
    Filled ellipse (s_in=None) or elliptical annulus between s_in and s_out.
 
    The annulus is a single patch with a hole, so overlapping translucent
    levels do NOT compound -- each shell shows exactly the alpha you asked for.
    """
    outer = _ellipse_xy(mu, cov, s_out)
    rings = [outer] if s_in is None else [outer, _ellipse_xy(mu, cov, s_in)[::-1]]
    verts, codes = [], []
    for r in rings:
        verts.append(np.vstack([r, r[:1]]))
        codes.append([Path.MOVETO] + [Path.LINETO] * (len(r) - 1) + [Path.CLOSEPOLY])
    return PathPatch(Path(np.vstack(verts), np.concatenate(codes)), **kw)
 
 
def component_mean_errors(gmm, X):
    """
    Uncertainty on each component mean: sigma_mu = sqrt(diag(Sigma_k) / N_k),
    with N_k = sum of the responsibilities (the effective number of members).
 
    Returns array of shape (n_components, 2).
    """
    resp = gmm.predict_proba(X)               # (n_samples, K)
    n_eff = resp.sum(axis=0)                  # (K,)
    covs = _full_covariances(gmm)
    var = np.array([np.diag(c) for c in covs])          # (K, 2)
    return np.sqrt(var / np.maximum(n_eff, 1.0)[:, None])
 
 
def _full_covariances(gmm):
    """Expand any covariance_type into full (K, 2, 2) matrices."""
    K, D = gmm.means_.shape
    ct, C = gmm.covariance_type, gmm.covariances_
    if ct == "full":
        return C
    if ct == "tied":
        return np.repeat(C[None], K, axis=0)
    if ct == "diag":
        return np.array([np.diag(c) for c in C])
    if ct == "spherical":
        return np.array([np.eye(D) * c for c in C])
    raise ValueError(f"unknown covariance_type {ct!r}")
 
 
# ----------------------------------------------------------------------
# fitting
# ----------------------------------------------------------------------
def fit_mixture(x, y, n_components, init_means=None, covariance_type="full", **gmm_kwargs):
    """
    Fit a Gaussian mixture seeded at the given means.
 
    init_means : (K, 2) array-like, e.g. [(mu_x0, mu_y0), (mu_x1, mu_y1), ...]
                 K is inferred from its length (3 for your case).
 
    Returns (gmm, X) where X is the (n, 2) array of finite points actually fit.
    """
    X = np.column_stack([np.ravel(np.asarray(x, float)),
                         np.ravel(np.asarray(y, float))])
    X = X[np.isfinite(X).all(axis=1)]
 
    if init_means:
        init_means = np.atleast_2d(np.asarray(init_means, float))
        
#     gmm = GaussianMixture(
#         n_components=n_components,
#         means_init=init_means,
#         covariance_type=covariance_type,
#         init_params='kmeans',
#         **gmm_kwargs,
#     )
    gmm = BayesianGaussianMixture(
        n_components=n_components,
        covariance_type=covariance_type,
        init_params='kmeans',
        **gmm_kwargs,
    )
    gmm.fit(X)
    return gmm, X
 
 
# ----------------------------------------------------------------------
# main entry point
# ----------------------------------------------------------------------
def gmm_plot(x, y, n_components, init_means=None, ax=None, n_sigma=(1, 2, 3), colors=None,
                 sigma_mode="mahalanobis", covariance_type="full",
                 err_scale=1.0, scatter=False,
                 alphas=None, alpha_max=0.55, alpha_min=0.10,
                 edge=True, lw=0.8, labels=None, gmm=None, **gmm_kwargs):
    """
    Fit a Gaussian mixture to (x, y) and plot it on `ax`.
 
    Parameters
    ----------
    x, y        : 1D arrays of data.
    init_means  : (K, 2) initial component centres.
    ax          : matplotlib axis (created if None).
    n_sigma     : scalar or sequence, e.g. 2 or (1, 2, 3) -> one contour each.
    colors      : list of K colours (one per component).
    sigma_mode  : 'mahalanobis' (default) or 'coverage'; see sigma_to_scale.
    err_scale   : multiply the mean-error cross arms (for visibility only).
    scatter     : if True, also scatter the points coloured by hard assignment.
    alphas      : opacity per level, ordered inner -> outer. Defaults to
                  linspace(alpha_max, alpha_min, len(n_sigma)), i.e. opacity
                  decreases with sigma.
    edge        : also stroke the ellipse outlines.
    gmm         : a pre-fit GaussianMixture; skips fitting if supplied.
 
    Returns
    -------
    gmm, ax
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
 
    if gmm is None:
        gmm, X = fit_mixture(x, y, n_components, init_means,
                             covariance_type=covariance_type, **gmm_kwargs)
    else:
        X = np.column_stack([np.ravel(x), np.ravel(y)])
        X = X[np.isfinite(X).all(axis=1)]
 
    K = gmm.n_components
    covs = _full_covariances(gmm)
    mu_err = component_mean_errors(gmm, X)
    scales = np.sort(sigma_to_scale(n_sigma, sigma_mode))   # ascending
 
    if alphas is None:
        alphas = np.linspace(alpha_max, alpha_min, len(scales))
    alphas = np.atleast_1d(alphas)
    if len(alphas) != len(scales):
        raise ValueError("alphas must have one entry per sigma level")
 
    if colors is None:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    colors = [colors[k % len(colors)] for k in range(K)]
 
    if scatter:
        lab = gmm.predict(X)
        for k in range(K):
            ax.scatter(*X[lab == k].T, s=6, color=colors[k],
                       alpha=0.25, lw=0, zorder=1)
 
    for k in range(K):
        mux, muy = gmm.means_[k]
        w0, h0, ang = cov_ellipse_params(covs[k])
 
        # n-sigma shells: filled, opacity decreasing outwards
        for j, s in enumerate(scales):
            s_in = scales[j - 1] if j > 0 else None      # hole -> no compounding
            ax.add_patch(_ring_patch(gmm.means_[k], covs[k], s, s_in,
                                     facecolor=colors[k], edgecolor="none",
                                     alpha=alphas[j], lw=0, zorder=2))
            if edge:
                ax.add_patch(Ellipse((mux, muy), s * w0, s * h0, angle=ang,
                                     facecolor="none", edgecolor=colors[k],
                                     lw=lw, alpha=min(1.0, alphas[j] + 0.25),
                                     zorder=3))
 
        # mean, as a cross whose arms are the uncertainty on the mean
        ex, ey = mu_err[k] * err_scale
        ax.errorbar(mux, muy, xerr=ex, yerr=ey, fmt="none",
                    ecolor=colors[k], elinewidth=2.0, capsize=0, zorder=5,
                    label=(labels[k] if labels is not None else
                           f"comp {k}  (w={gmm.weights_[k]:.2f})"))
        ax.plot(mux, muy, marker=".", ms=3, color=colors[k], zorder=6)
 
    ax.autoscale_view()
    return gmm, ax
 

def gmm_plot_pygmmi(gmms, ax=None, dims=(0, 1), nsig=1, ref=None,
             cov_average="arithmetic", show_replicates=False, scale_alpha=True,
             err_scale=1.0, color="C3", lw=2, label=None):
    """Plot the nsig-sigma contour of each GMM component.
 
    gmms : a single pygmmis.GMM, or a list of them (bootstrap replicates).
           With a list: components are matched across replicates, the contour is
           drawn from the AVERAGED covariance, and each mean gets an error bar
           from the bootstrap scatter.
 
    ref  : reference GMM defining component order (use your full-data fit).
 
    The two scales on the plot are different objects, do not conflate them:
      ellipse  -> the component's intrinsic width,  ~sqrt(Sigma_k)
      error bar-> the uncertainty on mu_k, smaller by ~sqrt(N_k)
 
    err_scale : multiply the error bars for visibility. At N ~ 2e4 the bootstrap
           error on mu is ~1/150 of the component width, i.e. sub-pixel. Blow it
           up and SAY SO in the caption -- do not silently mislead.
 
    Returns (ax, agg) where agg is the dict from aggregate_gmms.
    """
    from scipy.linalg import solve_triangular, logm, expm
    from scipy.optimize import linear_sum_assignment
    
    def match_components(gmms, ref=None):
        """Undo label switching across bootstrap replicates.

        Bootstrap fits return components in arbitrary order. Sorting by mean[:, 0]
        breaks as soon as two components overlap in x, so instead match each
        replicate to a reference fit by Hungarian assignment, costed by the
        Mahalanobis distance of the replicate means in the metric of the reference
        components (scale-free, unlike Euclidean).

        Returns (perms, ref), perms[b][a] = index in gmms[b] matching ref comp a.
        """
        if ref is None:
            ref = gmms[0]
        K = ref.K
        perms = []
        for g in gmms:
            cost = np.empty((K, K))
            for a in range(K):                         # rows: reference components
                L = np.linalg.cholesky(ref.covar[a])
                d = solve_triangular(L, (g.mean - ref.mean[a]).T, lower=True)
                cost[a] = (d ** 2).sum(0)              # squared Mahalanobis
            r, c = linear_sum_assignment(cost)
            perms.append(c[np.argsort(r)])
        return perms, ref


    def _mean_covar(covs, method="arithmetic"):
        """Average a stack of (B, D, D) SPD matrices.

        arithmetic : plug-in mean; unbiased for Sigma, but if replicate ellipses
                     differ in ORIENTATION it inflates the result ("determinant
                     swelling") -- the average of two crossed cigars is a disc.
        logeuclid  : expm(mean(logm)). Averages on the SPD manifold, no swelling.
                     Use when orientations scatter (weakly-constrained components).
        median     : elementwise; cheap and outlier-resistant.
        """
        if method == "arithmetic":
            return covs.mean(0)
        if method == "median":
            return np.median(covs, axis=0)
        if method == "logeuclid":
            return np.real(expm(np.mean([np.real(logm(C)) for C in covs], axis=0)))
        raise ValueError(method)



    def aggregate_gmms(gmms, ref=None, cov_average="arithmetic", q=(16, 50, 84)):
        """Collapse a bootstrap ensemble of pygmmis GMMs into mean parameters plus
        uncertainties. All arrays ordered like `ref`."""
        perms, ref = match_components(gmms, ref)
        K, D, B = ref.K, ref.mean.shape[1], len(gmms)

        amp = np.empty((B, K))
        mean = np.empty((B, K, D))
        covar = np.empty((B, K, D, D))
        for b, (g, p) in enumerate(zip(gmms, perms)):
            amp[b], mean[b], covar[b] = g.amp[p], g.mean[p], g.covar[p]

        lo, mid, hi = np.percentile(mean, q, axis=0)           # (K, D) each
        return dict(
            amp=amp.mean(0),
            amp_err=amp.std(0, ddof=1),
            mean=mean.mean(0),
            mean_err=np.stack([mid - lo, hi - mid]),           # (2, K, D), asymmetric
            mean_std=mean.std(0, ddof=1),
            covar=np.stack([_mean_covar(covar[:, k], cov_average) for k in range(K)]),
            samples=dict(amp=amp, mean=mean, covar=covar),
            n_boot=B,
        )


    def _ellipse(mu, C, s, **kw):
        vals, vecs = np.linalg.eigh(C)
        o = vals.argsort()[::-1]
        vals, vecs = vals[o], vecs[:, o]
        return Ellipse(mu, 2 * s * np.sqrt(vals[0]), 2 * s * np.sqrt(vals[1]),
                       angle=np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0])), **kw)



    single = not isinstance(gmms, (list, tuple))
    gmms = [gmms] if single else list(gmms)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    i, j = dims
    ij = np.ix_([i, j], [i, j])                    # correct marginal for D > 2
    nsigs = np.atleast_1d(nsig)
 
    agg = aggregate_gmms(gmms, ref=ref, cov_average=cov_average)
    K, amax = len(agg["amp"]), agg["amp"].max()
 
    for k in range(K):
        mu = agg["mean"][k, [i, j]]
        C = agg["covar"][k][ij]
        a = (0.25 + 0.75 * agg["amp"][k] / amax) if scale_alpha else 1.0
 
        if show_replicates and not single:         # honest about contour scatter
            for b in range(agg["n_boot"]):
                ax.add_patch(_ellipse(agg["samples"]["mean"][b, k, [i, j]],
                                      agg["samples"]["covar"][b, k][ij], nsigs[0],
                                      fc="none", ec=color, lw=1,
                                      alpha=0.18 * a, zorder=1))
 
        for s in nsigs:
            ax.add_patch(_ellipse(mu, C, s, fc="none", ec=color, lw=lw,
                                  alpha=a, zorder=3))
 
        if single:
            ax.plot(*mu, "+", color=color, ms=8, mew=2, alpha=a, zorder=4)
        else:
            label_plot = label if k == np.argmax(agg['amp']) else None
            err = err_scale * agg["mean_err"][:, k, [i, j]]   # [lo, hi] x [i, j]
            ax.errorbar(mu[0], mu[1], xerr=err[:, 0:1], yerr=err[:, 1:2],
                        fmt="+", color=color, ms=5, mew=1, elinewidth=1,
                        capsize=3, alpha=a, zorder=10, label=label_plot)

    return ax, agg