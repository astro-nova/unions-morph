# Masked Autoencoder pipeline

FiLM-conditioned masked autoencoder for the same UNIONS r-band morphology
catalog as `flow/`. Same dual goal — anomaly detection AND structure
discovery — approached through reconstruction error and a compressed
latent code instead of an exact density.

## Layout
```
mae/
  data.py     MAEDataset            HDF5 -> x, mask, c (mirrors FlowDataset)
  model.py    MaskedAutoencoder     one class: train / encode / reconstruct /
                                    anomaly_score / per_feature_anomaly /
                                    save / load
  run.py      example main
```

## Quick start (notebook)

```python
from mae import MAEDataset, MaskedAutoencoder

ds = MAEDataset(
    data_path="catalog.h5",
    feature_names=[...],          # 51 columns of x, in order
    masked_feature_names=[...],   # 48 features with stored masks; the rest get all-zero masks
)
x_tr, m_tr, c_tr = ds.train_arrays()
x_va, m_va, c_va = ds.val_arrays()

model = MaskedAutoencoder(n_features=ds.n_features, n_cond=ds.n_cond)
model.train(x_tr, m_tr, c_tr, val_x=x_va, val_mask=m_va, val_c=c_va)

z         = model.encode(x_va, m_va, c_va)         # (N, latent_dim)
x_hat     = model.reconstruct(x_va, m_va, c_va)    # (N, n_features)
score     = model.anomaly_score(x_va, m_va, c_va)  # (N,)
model.save("mae.eqx")
```

`MaskedAutoencoder.__init__` takes all model + training hyperparameters as
kwargs — built for use straight from a notebook.

## Design

- **Features `x`**: shape (N, 51), pre-scaled to ~(-10, 10).
- **Two kinds of mask, kept disentangled.** "Masking" changes two
  *independent* things — what the encoder is fed (input side) and what the
  loss scores (output side) — so we never fuse them into one bit:
  - **`Q` — intrinsic mask** (`mask` from the catalog, 1 = imputed): the
    value was never measured. Of the 51 features, 48 have stored masks; the
    other 3 get all-zero masks from `MAEDataset` so `mask` matches `x`'s
    column count. `Q` is **never** a loss target (no ground truth exists).
  - **`R` — random pretext mask**: positions artificially hidden each step
    so the latent must learn cross-feature structure. Drawn from the
    *observed* positions only, so `R ∩ Q = ∅` and every held-out position
    has a real target.
- **Conditioning `c`**: shape (N, 5) of standardized observing-condition
  vars (brightness, FWHM, sky median, sky σ, S/N). Always observed, so it
  sits outside the masking machinery entirely — never masked, never filled.
  Standardization stats are saved on the dataset (`ds.cond_mu`, `ds.cond_sd`).
- **Encoder**: input `[x_filled, mask_code]` of width `2 * n_features`, then
  `n_hidden_layers` of FiLM-modulated `Linear → GELU` to a tight bottleneck
  of size `latent_dim`. `x_filled` replaces masked entries with **distinct
  learned tokens** — a `mask_token` for `R` ("predict me") and a `fill_token`
  for `Q` ("absent"); both init at 0 ≈ observed-mean for the centred scaling,
  then learn. The single **bit-encoded** mask channel
  `mask_code = m_orig + 2·m_extra ∈ {0=observed, 1=Q, 2=R}` carries the state
  in one column instead of two (the types are mutually exclusive, so one
  channel suffices and the redundant all-zero dimension is dropped).
- **Decoder**: mirror — `latent_dim → hidden → ... → n_features`, also
  FiLM-modulated by `c`. Reconstructs all `n_features`.
- **FiLM**: each hidden activation is rescaled as
  `h ← (1 + γ(e)) ⊙ Linear(h) + β(e)`, applied at **every** encoder and
  decoder block. `e = conditioner(c)` is a small shared MLP that embeds `c`
  **non-linearly**, so γ/β can bend with S/N and resolution rather than
  depend on them linearly (a single linear FiLM head can't represent
  threshold-like noise/depth effects). The per-layer FiLM heads are
  initialized small so γ ≈ 0, β ≈ 0 — at step 0 the network is a plain MLP
  and conditioning is learned gradually. FiLM **modulates** the feature
  stream (γ can drive a unit toward 0, gating a channel when the conditioning
  says it is unreliable); cross-parameter degeneracies are still represented
  by the unconditioned `Linear` layers that mix channels between FiLM ops.
  The conditioning never travels the main connections as data, which is what
  concatenating `c` would do.

### Why FiLM rather than concatenating `c` to the input

Brightness influences nearly every measurement (S/N, PSF coupling,
deblending), so it needs to *modulate* the rest of the network's response,
not just add a constant offset. Concatenating `c` to the input layer
gives a fixed bias on each downstream activation; FiLM gives a per-channel
affine transform whose scale and shift depend non-linearly on `c`, at every
layer.
Practically: a galaxy that looks dim and a galaxy that looks bright should
end up at the same latent location if they're the same morphology — FiLM
lets the encoder "divide out" brightness instead of carrying it through.

### Adversarial disentanglement (optional)

FiLM lets the network *use* `c`, but doesn't stop the encoder from also
encoding `c`-correlated information into `z`. If you plot `z` colored by
e.g. r-magnitude and see a clear gradient, the encoder is leaking `c`
through, even though FiLM is doing its job. This matters when you want
`z` to capture *intrinsic* morphology rather than observational state.

Setting `adv_lambda > 0` adds an adversarial head — a small MLP that
predicts (a subset of) `c` from `z` — coupled to the encoder through a
gradient-reversal layer (GRL). The GRL is identity in forward and flips
the sign of the cotangent in backward, so a single MSE term

    adv_loss = ‖adv_head(GRL(z, λ)) − c[:, target_dims]‖²

simultaneously trains the head to be a *good* predictor and pushes the
encoder to make `z` *unpredictive* of the targeted condition vars. The
total objective is `recon_loss + adv_loss`.

`λ` is ramped linearly from 0 to `adv_lambda` over `adv_warmup_steps` so
the head becomes a competent predictor before its gradient meaningfully
shapes `z` — otherwise the encoder gets pulled around by an untrained
adversary. `adv_target_dims` chooses which columns of `c` to push out
(default: all of `c` — the standard "z ⊥ c" formulation). To target only
brightness, e.g. `adv_target_dims=(0,)`.

Cost: the head is a tiny MLP (~1–2K params vs ~150K+ in the encoder).
Per-step overhead is a few percent.

To verify the gradient is gone after training, regress `c` from `z`
(linear or MLP probe) on a held-out set — high R² means leakage remains;
near-zero R² for the targeted dims means the adversary did its job.

### Latent structure regularizers (optional)

A pure reconstruction loss pins down only the latent *subspace*, not its
basis: any invertible linear map of `z` can be absorbed by the decoder's
first layer at zero cost, so which direction each dim points is an accident
of initialization. Two independent knobs (both off by default) break that
symmetry so individual dims become meaningful. They compose with each other
and with the adversary; the adversarial head always sees the **full** code,
so `c`-leakage is scrubbed from every dim regardless of truncation.

**Nested dropout (`nested_dropout_p > 0`) — importance-ordered dims.**
Each training row draws `k ~ Geometric(p)` (clipped to `latent_dim`) and the
decoder sees only `z[:k]`, rest zeroed. Reconstruction-critical information
is forced into the early dims, so `z[:k]` is the best available k-dim code
for *every* k — PCA's ordered-components property, with a fully nonlinear
map (in the linear case this provably recovers PCA). Side effects you get
for free: pick the effective latent size *after* training from
`truncation_curve(...)` (its knee is the effective dimensionality), and
late dims going quiet is dimensionality selection, not collapse. Expect the
leading dims to carry cluster identity — multimodal marginals there are a
feature, not a bug. Guide: `E[k] ≈ 1/p`, so with `latent_dim=16`, `p=0.1`
truncates most rows (mean k ≈ 10, ~20% of rows keep the full code); `p → 0`
approaches plain training.

**Jacobian orthogonality (`jac_ortho_lambda > 0`) — non-overlapping
effects.** Penalizes the mean squared cosine between columns of the decoder
Jacobian `∂x̂/∂z` (column i = the feature-space direction latent i moves),
on `jac_ortho_batch` rows per step. Different dims are pushed to move
*disjoint* feature combinations at every point of the manifold —
curvilinear-orthogonal coordinates. This is deliberately the Gram/cosine
form and **not** a Hessian-diagonal penalty: it does not force the decoder
to be additive, so nonlinear interactions between latents (polar-coordinate
style) remain representable. Cosines make it scale-invariant (can't be
gamed by shrinking effects), and dims that nested dropout has killed
contribute ≈ 0, so the two regularizers don't fight. The penalty is
evaluated at the codes the decoder actually consumed (truncated ones when
nested dropout is on), which with ordering reads like Gram–Schmidt: each
newly activated dim adds an effect orthogonal to those before it.

After training, interpret axes with latent traversals: sweep one dim of `z`
across its population range, `decode(z, c)`, and watch which features
respond. If the population is clustered, traverse from cluster centroids
rather than the global median — per-cluster axis meanings are the honest
level of description for factors that only exist in some types.

### Training objective

The real objective is **compression — encode 50 features in a `latent_dim`
code and reconstruct them all** — with masking acting as regularization. So
every observed feature's reconstruction matters, masked or not.

Each batch we draw a random pretext mask
`m_extra ~ Bernoulli(mask_ratio)` **over the observed positions only**
(`m_orig = 0`), on top of the catalog mask `m_orig`. The encoder is fed the
learned tokens at the masked slots plus the bit-encoded `mask_code` channel.
The decoder reconstructs all 51 features.

The loss is an MSE over **every observed position** `O = ¬Q`, with the held-out
and still-visible regions weighted independently:

    loss = ( masked_w · Σ_{R} (x̂-x)²  +  unmasked_w · Σ_{visible∩O} (x̂-x)² ) / Σ weights

| region                          | `mask_code` | weight                 |
|---------------------------------|-------------|------------------------|
| `Q` intrinsically masked        | 1           | **0** (no ground truth)|
| `R` held-out (denoising signal) | 2           | `masked_loss_weight`   |
| visible & observed              | 0           | `unmasked_loss_weight` |

Both masked and unmasked features are scored. The held-out `R` positions carry
the MAE / denoising signal that forces the latent to learn cross-feature
structure, so they are **up-weighted by default** (`masked=1.0`,
`unmasked=0.3`); the visible positions keep the bottleneck honest as a
compressor of the full feature vector. The intrinsic mask `Q` is never scored —
it has no ground truth. Set `unmasked_loss_weight=0` for the strict held-out-only
MAE, or raise it toward `1.0` to weight plain reconstruction more.

The default **anomaly score** (`anomaly_score(...)`, `mc_seeds=0`) is the same
full-reconstruction error over *all* observed features — consistent with the
compression objective. The held-out MC variant (`mc_seeds>0`) is a secondary
diagnostic that scores only the artificially-hidden positions.

### Two ways to use the trained model

| use case            | call                         | what comes out                  |
|---------------------|------------------------------|---------------------------------|
| anomaly detection   | `anomaly_score(x, m, c)`     | per-source reconstruction error |
| pattern discovery   | `encode(x, m, c)`            | `(N, latent_dim)` latent code   |
| inspect a source    | `per_feature_anomaly(...)`   | leave-one-out recon per feature |
| effective dim       | `truncation_curve(x, m, c)`  | recon error vs. keeping `z[:k]` |

The latent code is a compressed, brightness-scrubbed embedding suitable
for HDBSCAN / GMM directly (same recommendation as for the flow's `z`:
operate in the embedding, not in `UMAP(embedding)`).

## Hyperparameters

All model and training knobs are constructor kwargs on `MaskedAutoencoder`.

### Architecture

| param | default | what it controls |
|-------|---------|------------------|
| `n_features` | — | dim of `x`. Match `ds.n_features`. |
| `n_cond` | — | dim of `c` (just the continuous conditioners; mask is a separate input). Match `ds.n_cond`. |
| `latent_dim` | 16 | bottleneck width. The single most important capacity knob — too small and distinct galaxy types collapse together; too large and the AE identity-bypasses, blunting the anomaly signal. |
| `hidden_dim` | 256 | width of every FiLM-MLP layer in encoder and decoder. |
| `n_hidden_layers` | 3 | depth of the encoder (decoder mirrors it). Each block is `Linear → FiLM(e) → GELU`. |
| `cond_embed_dim` | 16 | width of the shared conditioner that embeds `c` non-linearly into `e`, then feeds every FiLM head. Drives FiLM cost (≈ `2·n_hidden_layers · cond_embed_dim · hidden_dim`). |
| `cond_embed_layers` | 2 | depth of that conditioner MLP. `≥2` makes the embedding non-linear in `c`; `1` collapses it to a linear FiLM generator. |

- **`latent_dim`** — 8 is enough to separate the headline galaxy classes
  (early/late/irregular) but smears finer structure; 32 captures finer
  morphology + several artifact modes but reconstructs observed inputs
  almost trivially, which weakens `anomaly_score`. 16 is a reasonable
  default for ~50 features with multi-class structure.
- **`hidden_dim` / `n_hidden_layers`** — capacity to learn the conditional
  cross-feature manifold. Symptoms of too small: training loss flattens
  high; per-feature LOO recon is biased on common galaxy types.
  256×3 fits comfortably in 8 GB at batch 8192. Next step up is 512×3.

### Training

| param | default | what it controls |
|-------|---------|------------------|
| `mask_ratio` | 0.5 | fraction of *observed* features artificially hidden (`R`) each step. |
| `masked_loss_weight` | 1.0 | loss weight on the held-out `R` positions (the denoising target). |
| `unmasked_loss_weight` | 0.3 | loss weight on the still-visible observed positions. `R` is up-weighted relative to this; set to `0` for strict held-out-only MAE, or `1.0` to weight plain reconstruction equally. |
| `batch_size` | 8192 | examples per gradient step. |
| `learning_rate` | 1e-3 | Adam step size. |
| `n_epochs` | 10 | full passes over the training set. |
| `grad_clip` | 1.0 | max global grad norm. |
| `seed` | 0 | RNG for parameter init, shuffle order, and mask draws. |

- **`mask_ratio`** — vision MAE uses 0.75 because image patches are very
  redundant. Tabular morphology features are not nearly as redundant
  (Gini and M20 are correlated but Sersic n is not predictable from C/A/S
  alone). 0.4–0.6 is the useful range. Too high → loss flatlines high
  (model has no context). Too low → encoder learns to identity-pass
  visible features and the latent is uninformative. The ratio is over the
  *observed* positions, so mask budget is never spent on `Q` slots you drop
  from the loss anyway.
- **`masked_loss_weight` / `unmasked_loss_weight`** — the held-out–vs–visible
  balance. Default `1.0 / 0.3` up-weights the held-out `R` positions (the
  denoising signal) while still scoring every visible feature, since the goal
  is to reconstruct the whole vector through the bottleneck. `unmasked=0` is
  the canonical held-out-only MAE; `unmasked=1.0` weights plain reconstruction
  equally with the held-out signal. If `anomaly_score` separation weakens
  (identity bypass), lower `unmasked_loss_weight` toward the held-out emphasis.
- **`batch_size`** — larger = smoother gradients, fewer steps per epoch.
  Push up while it fits; 16k–32k is comfortable on 8 GB for this model.
- **`learning_rate`** — twitchy loss → drop to 5e-4. Slow convergence at
  large batch → bump to 2e-3.
- **`n_epochs`** — with ~100M sources, watch `epoch_val_loss` and stop
  when it flattens. 10 is a reasonable starting budget.

### Latent structure regularizers

Off by default. See "Latent structure regularizers (optional)" above.

| param | default | what it controls |
|-------|---------|------------------|
| `nested_dropout_p` | 0.0 | geometric truncation parameter; 0 disables. Each row keeps only `z[:k]`, `k ~ Geometric(p)`. `E[k] ≈ 1/p`: useful range ≈ `1/latent_dim` (gentle) to `2/latent_dim` (strong ordering). |
| `jac_ortho_lambda` | 0.0 | weight on the decoder-Jacobian orthogonality penalty; 0 disables (no `jacfwd` cost at all). The raw penalty is a mean squared cosine ∈ [0, 1] (random directions in 51-dim feature space sit near 1/51 ≈ 0.02), so start at 1.0 and move by ×3. |
| `jac_ortho_batch` | 128 | rows per step the penalty is evaluated on. Cost ≈ `latent_dim` decoder JVPs × this many rows — a few % overhead at defaults. |

- **`nested_dropout_p` too high** — early dims are overloaded and total
  recon degrades visibly; the truncation curve front-loads but the k =
  `latent_dim` error is clearly worse than the unregularized run. Lower `p`.
- **`nested_dropout_p` too low** — `truncation_curve` stays flat-high until
  large k (no ordering emerged). Raise `p`.
- **`jac_ortho_lambda`**: watch `history["epoch_jac_loss"]` — healthy runs
  decay and plateau low while recon stays near baseline. If recon stalls or
  gets noisy, the weight is too high. If `epoch_jac_loss` plateaus barely
  below its starting value, it's too low to matter.
- Validation loss is computed under the same truncation regime as training,
  so val numbers are **not comparable across runs with different
  `nested_dropout_p`** (use `loss_of`-style baselines or the truncation
  curve to compare).

### Adversarial disentanglement

Off by default. Turn on by setting `adv_lambda > 0`.

| param | default | what it controls |
|-------|---------|------------------|
| `adv_lambda` | 0.0 | weight on the adversarial loss after warmup. 0 disables the head entirely (no compute, no extra params). Start small (0.1–1.0) and raise if leakage persists. |
| `adv_target_dims` | `None` (= all of `c`) | tuple of column indices in `c` to push out of `z`. e.g. `(0,)` for rmag-only. |
| `adv_hidden` | 64 | hidden width of the adversary MLP (2 hidden layers). The head needs to be expressive enough to be a real predictor, but doesn't need to be deep. |
| `adv_warmup_steps` | 1000 | steps over which `λ` ramps linearly from 0 to `adv_lambda`. Set ~ 1 epoch's worth of steps. |

- **`adv_lambda` too high** — recon loss stalls, training is unstable, or
  `z` collapses (encoder over-prioritizes fooling the adversary). Drop
  `adv_lambda` or stretch `adv_warmup_steps`.
- **`adv_lambda` too low** — adv loss in `history["epoch_adv_loss"]` keeps
  falling indefinitely (head wins easily) and the rmag gradient in `z`
  doesn't soften. Raise `adv_lambda`.
- **Healthy run**: adv loss drops during warmup as the head learns, then
  rises and plateaus as the encoder fights back. Recon loss decreases
  monotonically (perhaps slightly worse than the no-adv baseline).

### Tuning by symptom

- Anomaly score has no separation between obvious artifacts and normal
  sources → `latent_dim` too large (identity bypass) **or** `mask_ratio`
  too low. Drop one or both.
- Latent space (`encode`) has no visible structure under HDBSCAN/GMM →
  `latent_dim` too small **or** training underfit (loss still falling).
  Bump `latent_dim`, then `hidden_dim`/`n_hidden_layers`, then epochs.
- Validation loss flat-but-high across epochs → underfit. Bigger model
  or smaller `mask_ratio`.
- `latent_dim` doubled and val-loss barely moved → bottleneck wasn't the
  binding constraint; the encoder/decoder MLPs are.
- `c`-leakage: cluster labels in `z` correlate strongly with brightness
  → FiLM gives the model a path to *use* `c` but does not force the
  encoder to factor `c` *out* of `z`. First try training longer and
  bumping `hidden_dim`. If the gradient persists and you want `z`
  intrinsic-only, enable adversarial disentanglement (`adv_lambda > 0`,
  `adv_target_dims=(0,)` for rmag-only or default-all-of-`c`). See
  "Adversarial disentanglement" above.
- Per-feature LOO error has one feature systematically off → that
  feature is the hardest to predict from the rest, which usually means
  it carries genuinely independent signal. Check whether it correlates
  with the catalog mask of another feature.

## Adapt to your HDF5

Override the path templates if needed:
```python
MAEDataset(
    ...,
    feature_path="features/{name}",
    mask_path="masks/{name}",
    cond_path="conditioning/{name}",
)
```
e.g. flat layout: `feature_path="{name}"`, `mask_path="{name}_mask"`.

## Training trackers

`model.train(...)` populates `model.history`:

| key                  | shape       | meaning                                            |
|----------------------|-------------|----------------------------------------------------|
| `step_loss`          | (n_steps,)  | total training loss per step (recon + adv + jac)   |
| `step_recon_loss`    | (n_steps,)  | reconstruction component per step                  |
| `step_adv_loss`      | (n_steps,)  | adversarial component per step (0 if `adv_lambda=0`) |
| `step_jac_loss`      | (n_steps,)  | Jacobian-orthogonality penalty per step (0 if off) |
| `epoch_train_loss`   | (n_epochs,) | mean total training loss per epoch                 |
| `epoch_recon_loss`   | (n_epochs,) | mean reconstruction loss per epoch                 |
| `epoch_adv_loss`     | (n_epochs,) | mean adversarial loss per epoch                    |
| `epoch_jac_loss`     | (n_epochs,) | mean Jacobian-orthogonality penalty per epoch      |
| `epoch_val_loss`     | (n_epochs,) | validation loss per epoch (recon only, if val given) |
| `epoch_time`         | (n_epochs,) | seconds per epoch                                  |
| `best_val_loss`      | float       | lowest val loss seen                               |
| `best_epoch`         | int         | epoch index of best val loss                       |

`progress=True` (default) shows nested tqdm bars (epoch + step). Pass
`checkpoint_dir=...` to write `mae_epoch{NNN}.eqx` after every epoch.

## Worked example: artifacts + clusters on UNIONS r-band

Goal: train one MAE on the bulk population, then use reconstruction error
to flag artifacts/mergers and the latent code to cluster the bulk into
galaxy types.

```python
import numpy as np
import matplotlib.pyplot as plt
from mae import MAEDataset, MaskedAutoencoder

# 1. Load. Same conventions as FlowDataset — standardizes c, fills in
#    zero-masks for the 3 features without stored masks.
ds = MAEDataset(
    data_path="catalogs/dataset.h5",
    feature_names=[
        "C", "A", "S", "G", "M20", "Gini",
        "rhalf_circ", "rhalf_ellip", "r20", "r80",
        "sersic_n", "sersic_ellip", "sersic_xc", "sersic_yc",
        # ... 51 total
    ],
    masked_feature_names=[
        "C", "A", "S", "G", "M20", "Gini",
        # ...
    ],
)
x_tr, m_tr, c_tr = ds.train_arrays()
x_va, m_va, c_va = ds.val_arrays()

# 2. Fit. Defaults are tuned for this catalog size; sweep from a notebook.
model = MaskedAutoencoder(
    n_features=ds.n_features,
    n_cond=ds.n_cond,
    latent_dim=16,
    hidden_dim=256,
    n_hidden_layers=3,
    mask_ratio=0.5,
    n_epochs=10,
    batch_size=16384,
    learning_rate=1e-3,
)
model.train(
    x_tr, m_tr, c_tr,
    val_x=x_va, val_mask=m_va, val_c=c_va,
    checkpoint_dir="runs/mae_v1",
)

# 3. Diagnostics: loss curves
h = model.history
fig, ax = plt.subplots(1, 2, figsize=(10, 3))
ax[0].plot(h["step_loss"]);                 ax[0].set_title("step loss")
ax[1].plot(h["epoch_train_loss"], label="train")
ax[1].plot(h["epoch_val_loss"],   label="val")
ax[1].legend(); ax[1].set_title("epoch loss")
print(f"best val loss {h['best_val_loss']:.4f} at epoch {h['best_epoch']}")

# 4. Anomaly score. MC mode is closer to the training objective.
score = model.anomaly_score(ds.x, ds.mask, ds.c, mc_seeds=8, progress=True)
thr = np.quantile(score, 0.999)
candidate_idx = np.where(score > thr)[0]

# 5. Latent embedding for HDBSCAN / GMM downstream.
z = model.encode(ds.x, ds.mask, ds.c, progress=True)   # (N, latent_dim)

# 6. Per-source diagnostic: which feature is the most surprising?
sub = candidate_idx[:200]
x_hat_loo, sq_err = model.per_feature_anomaly(
    ds.x[sub], ds.mask[sub], ds.c[sub]
)
worst_feat = sq_err.argmax(axis=1)   # which feature drives the anomaly

# 7. Persist. Save the model weights and the cond standardization stats so
#    the same MAE can score new tiles consistently. `save` writes
#    `mae.eqx` (weights) plus `mae.eqx.json` (architecture sidecar).
model.save("runs/mae_v1/mae.eqx")
np.savez("runs/mae_v1/cond_stats.npz",
         cond_mu=ds.cond_mu, cond_sd=ds.cond_sd)

# 8. Reload later. `from_checkpoint` reads the sidecar so you don't
#    have to remember n_features / latent_dim / etc.
later = MaskedAutoencoder.from_checkpoint("runs/mae_v1/mae.eqx")
assert np.allclose(
    later.anomaly_score(ds.x[:1000], ds.mask[:1000], ds.c[:1000]),
    model.anomaly_score(ds.x[:1000], ds.mask[:1000], ds.c[:1000]),
)
```

## Comparison with `flow/`

Both pipelines target the same dataset and the same dual goal (anomaly +
clustering). Choose by what you need:

| dimension                | `flow/`                                   | `mae/`                              |
|--------------------------|-------------------------------------------|-------------------------------------|
| anomaly score            | `-log p(x \| c)` — calibrated, per-source | reconstruction MSE — uncalibrated   |
| latent space             | `(N, 51)` Gaussianized, invertible        | `(N, latent_dim)` compressed        |
| handles missing values   | mask flags appended to `c`                | mask is a first-class encoder input |
| brightness conditioning  | concatenated into `c`                     | non-linear FiLM at every layer      |
| forward cost (per source)| 1 flow pass                               | 1 AE pass                           |
| per-feature score        | conditional MC z-score (`diag.py`)        | leave-one-out reconstruction        |

The flow gives you a probabilistic ranking; the MAE gives you a low-dim
embedding for free. They are complementary — the bulk population should
agree on what is anomalous, and disagreements between the two scorers
are themselves interesting.

## Notes / things to revisit

- `n_epochs=10` is a starting point. With ~100M sources you may converge
  in fewer; switch to a step budget if useful.
- I/O: `MAEDataset` loads everything to RAM (~27 GB for your data). Fine
  with 200 GB. For larger catalogs, swap to a chunked HDF5 reader.
- `anomaly_score(mc_seeds=0)` is the cheapest scoring mode but vulnerable
  to identity bypass on highly informative features. `mc_seeds=8` matches
  the training distribution and is the recommended default for downstream
  ranking.
- The per-feature LOO diagnostic costs F forward passes per source — fine
  for inspecting hundreds of candidates, not for the full 100M. For an
  always-on per-feature score, batch many sources at moderate
  `batch_size` and let JAX vmap absorb the cost.
- Latent `z` from `encode(...)` is what goes into HDBSCAN/GMM and the
  supervised anchoring step. Operate in `z` directly; UMAP only for
  visualizing pre-computed labels.

## Requirements

```
jax
jaxlib
equinox
optax
h5py
numpy
tqdm   # optional; only used for progress bars
```
