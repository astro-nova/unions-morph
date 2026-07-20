# Is the non-normal MAE latent distribution expected?

**Question.** The distribution of latent parameters predicted by the MAE is very
non-normal. X1 (the most important parameter) has sharp peaks and valleys.
Is this expected behaviour?

**Short answer.** The broad multimodality is expected and is arguably the point.
The narrow spike at 0 is a separate phenomenon and should not be accepted at
face value.

---

## 1. The multimodality is by design

There is no term anywhere in the objective that pushes `z` toward a Gaussian.
The loss is `recon + adv + jac_ortho` (`mae/model.py`) — no KL to a prior (that
would be a VAE), no base-density matching (that's `flow/`). A plain bottleneck
AE places `z` wherever reconstruction is cheapest, and if the galaxy population
is genuinely multimodal in morphology, the leading latent dims will be too.

If the "should be normal" intuition comes from the flow pipeline, that is the
right intuition for the wrong model. The README comparison table says it
directly:

| dimension     | `flow/`                            | `mae/`                        |
|---------------|------------------------------------|-------------------------------|
| latent space  | `(N, 51)` Gaussianized, invertible | `(N, latent_dim)` compressed  |

The flow Gaussianizes by construction. The MAE gives no distributional
guarantee whatsoever.

`mae/README.md` (lines 150–151) already calls this shot for nested dropout:

> Expect the leading dims to carry cluster identity — multimodal marginals
> there are a feature, not a bug.

A unimodal Gaussian X1 would mean there was no discoverable structure for
HDBSCAN to find.

## 2. The spike at zero is a different animal

Compare the scales in the histogram. The three broad modes (~-3, ~-1.2, ~+0.9)
are smooth and wide — those read as populations. The peak at 0 is *far*
narrower than any of them, and is flanked by deep depleted valleys on both
sides. That is not a cluster; it is a degenerate subpopulation collapsing onto
nearly identical `z`.

### Leading hypothesis: heavily-masked rows

From `mae/model.py:236`:

```python
x_filled = x * visible + self.fill_token * m_orig + self.mask_token * m_extra
```

A source where most of the 48 maskable features are intrinsically masked (`Q`)
gets `fill_token` in nearly every slot and `mask_code = 1` across the board. Its
encoder input then depends only on `c` via FiLM plus the handful of
always-observed features — so every such source lands in almost the same place.

In a UNIONS r-band catalog those are the failed Sérsic / failed CAS sources:
faint, small, poorly deblended. There are a lot of them, which matches a 370k
count needle. The flanking valleys fit too — those sources were pulled *out* of
the continuum onto the atom.

### The one-line check

```python
frac_masked = ds.mask.mean(axis=1)
spike = np.abs(z[:, 1]) < 0.05          # tune eps to the needle width
print(np.median(frac_masked[spike]), np.median(frac_masked[~spike]))
```

If the spike median is near 1 and the bulk near 0, confirmed. Also worth
checking `ds.c` for those rows — expect faint / low `sn_per_pixel`.

### Secondary candidate

If the above comes back flat: sentinel values in the 3 features with no stored
mask. `mae/data.py:250` scales as `(x - median) / iqr`, so any feature with an
atom at its median maps that atom to exactly 0.

### Ruled out

A late latent dim switched off by nested dropout would also give a zero spike,
but would not coexist with the broad trimodal continuum in the plot.

## 3. Why this matters downstream

If the masking cause is confirmed, it will distort clustering in a way that
looks like a result:

- **HDBSCAN** will find the atom as a huge ultra-dense cluster that means
  nothing but "no measurements available."
- **GMM** is worse — a component can collapse onto a near-zero-variance atom,
  sending the likelihood to infinity (singular covariance).

Recommended: cut on mask fraction before clustering rather than trying to
interpret that mode.

## 4. Open question

Is X1 the first latent dim, and is `nested_dropout_p > 0`? This changes the
read on whether the ordering interpretation applies.

---

## References

- `mae/README.md` — "Latent structure regularizers (optional)", comparison
  with `flow/`
- `mae/model.py:236` — encoder fill logic
- `mae/data.py:250` — feature scaling
