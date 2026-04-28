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
- **Mask**: per-feature binary (1 = imputed). Of the 51 features, 48 have
  stored masks; the remaining 3 had no bad values and get all-zero masks
  reconstructed by `MAEDataset` so `mask` has the same column count as `x`.
- **Conditioning `c`**: shape (N, 5) of standardized observing-condition
  vars (brightness, FWHM, sky median, sky σ, S/N). Standardization stats
  are saved on the dataset (`ds.cond_mu`, `ds.cond_sd`).
- **Encoder**: input `[x ⊙ (1-m), m]` of width `2 * n_features`, then
  `n_hidden_layers` of FiLM-modulated `Linear → GELU` to a tight bottleneck
  of size `latent_dim`.
- **Decoder**: mirror — `latent_dim → hidden → ... → n_features`, also
  FiLM-modulated by `c`.
- **FiLM**: each hidden activation is rescaled as
  `h ← (1 + γ(c)) ⊙ Linear(h) + β(c)` with `γ, β` linear in `c`. The FiLM
  head is initialized small so γ ≈ 0, β ≈ 0 — at step 0 the network is a
  plain MLP and conditioning is learned gradually.

### Why FiLM rather than concatenating `c` to the input

Brightness influences nearly every measurement (S/N, PSF coupling,
deblending), so it needs to *modulate* the rest of the network's response,
not just add a constant offset. Concatenating `c` to the input layer
gives a fixed bias on each downstream activation; FiLM gives a per-channel
affine transform whose scale and shift depend on `c`, at every layer.
Practically: a galaxy that looks dim and a galaxy that looks bright should
end up at the same latent location if they're the same morphology — FiLM
lets the encoder "divide out" brightness instead of carrying it through.

### Training objective

Each batch we draw an extra random mask `m_extra ~ Bernoulli(mask_ratio)`
on top of the catalog mask `m_orig`. The encoder sees zeros at all hidden
positions plus a binary `m_total = m_orig | m_extra` channel that flags
them. The decoder reconstructs all 51 features. Loss is MSE averaged over
positions that are *artificially hidden* and *originally observed*:

    loss = mean_{m_extra=1, m_orig=0} (x_hat - x)^2

i.e. the standard MAE / BERT-style denoising objective. Originally-imputed
positions are never used as ground truth (they aren't real). Artificially
masking forces the latent to capture cross-feature structure rather than
collapse to identity.

### Two ways to use the trained model

| use case            | call                         | what comes out                  |
|---------------------|------------------------------|---------------------------------|
| anomaly detection   | `anomaly_score(x, m, c)`     | per-source reconstruction error |
| pattern discovery   | `encode(x, m, c)`            | `(N, latent_dim)` latent code   |
| inspect a source    | `per_feature_anomaly(...)`   | leave-one-out recon per feature |

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
| `n_hidden_layers` | 3 | depth of the encoder (decoder mirrors it). Each block is `Linear → FiLM(c) → GELU`. |

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
| `mask_ratio` | 0.5 | fraction of features artificially hidden each step. The MAE objective is computed only on these positions. |
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
  visible features and the latent is uninformative.
- **`batch_size`** — larger = smoother gradients, fewer steps per epoch.
  Push up while it fits; 16k–32k is comfortable on 8 GB for this model.
- **`learning_rate`** — twitchy loss → drop to 5e-4. Slow convergence at
  large batch → bump to 2e-3.
- **`n_epochs`** — with ~100M sources, watch `epoch_val_loss` and stop
  when it flattens. 10 is a reasonable starting budget.

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
  → FiLM isn't separating out the conditioning. Train longer; if it
  persists, increase `hidden_dim` (FiLM head capacity scales with it).
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

| key                  | shape       | meaning                                  |
|----------------------|-------------|------------------------------------------|
| `step_loss`          | (n_steps,)  | training loss at every optimizer step    |
| `epoch_train_loss`   | (n_epochs,) | mean training loss per epoch             |
| `epoch_val_loss`     | (n_epochs,) | validation loss per epoch (if val given) |
| `epoch_time`         | (n_epochs,) | seconds per epoch                        |
| `best_val_loss`      | float       | lowest val loss seen                     |
| `best_epoch`         | int         | epoch index of best val loss             |

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
#    the same MAE can score new tiles consistently.
model.save("runs/mae_v1/mae.eqx")
np.savez("runs/mae_v1/cond_stats.npz",
         cond_mu=ds.cond_mu, cond_sd=ds.cond_sd)

# 8. Reload later (architecture must match — same n_features / n_cond /
#    latent_dim / hidden_dim / n_hidden_layers).
later = MaskedAutoencoder(n_features=ds.n_features, n_cond=ds.n_cond)
later.load("runs/mae_v1/mae.eqx")
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
| brightness conditioning  | concatenated into `c`                     | FiLM-modulated at every layer       |
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
