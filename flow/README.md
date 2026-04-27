# Conditional flow pipeline

Conditional rational-quadratic spline coupling flow on galaxy morphology
features, conditioned on observational state + per-feature mask flags.

## Layout
```
flow/
  data.py     FlowDataset    HDF5 -> x, mask, c (with mask reconstruction)
  model.py    ConditionalFlow  one class: train / evaluate / transform / save / load
  run.py      example main
```

## Quick start (notebook)

```python
from flow import FlowDataset, ConditionalFlow

ds = FlowDataset(
    data_path="catalog.h5",
    feature_names=[...],          # 51 columns of x, in order
    masked_feature_names=[...],   # 48 features with stored masks; the rest get all-zero masks
)
x_tr, c_tr = ds.train_arrays()
x_va, c_va = ds.val_arrays()

model = ConditionalFlow(n_features=ds.n_features, n_cond=ds.n_cond)
model.train(x_tr, c_tr, val_x=x_va, val_c=c_va)

log_prob = model.evaluate(x_va, c_va)   # per-source log p(x|c)
z        = model.transform(x_va, c_va)  # latent (Gaussianized)
model.save("flow.eqx")
```

`ConditionalFlow.__init__` takes all model + training hyperparameters as
kwargs — built for use straight from a notebook. Defaults match the old
`FlowConfig` (10 coupling layers, 8 spline knots, MLP 256x3, batch 8192,
lr 1e-3, 5 epochs).

## Design

- **Features `x`**: shape (N, 51), already in (-10, 10).
- **Masks**: per-feature binary (1 = imputed). Of the 51 features, 48 have
  stored masks; the remaining 3 had no bad values and get all-zero masks
  reconstructed by `FlowDataset` so `mask` has the same column count as `x`.
- **Conditioning `c`**: shape (N, 56) = `[standardized cond (5) | mask flags (51)]`.
  Standardization stats are saved on the dataset (`ds.cond_mu`, `ds.cond_sd`).
- **Model**: flowjax `coupling_flow` with `RationalQuadraticSpline(knots=8, interval=10)`,
  10 coupling layers, MLP conditioner 256x3.

## Adapt to your HDF5

Override the path templates if needed:
```python
FlowDataset(
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
`checkpoint_dir=...` to write `flow_epoch{NNN}.eqx` after every epoch.

## Worked example: anomaly hunt on UNIONS r-band

Goal: learn p(morphology | observing conditions) on the bulk population,
then flag the lowest-log-prob sources as anomaly candidates and embed the
cleaner sources for downstream clustering.

```python
import numpy as np
import matplotlib.pyplot as plt
from flow import FlowDataset, ConditionalFlow

# 1. Load. FlowDataset standardizes the 5 cond vars in place and
#    reconstructs zero-masks for the 3 features without stored masks.
ds = FlowDataset(
    data_path="catalogs/dataset.h5",
    feature_names=[
        "C", "A", "S", "G", "M20", "Gini",
        "rhalf_circ", "rhalf_ellip", "r20", "r80",
        "sersic_n", "sersic_ellip", "sersic_xc", "sersic_yc",
        # ... 51 total
    ],
    masked_feature_names=[
        # 48 of the above that have a stored mask in HDF5
        "C", "A", "S", "G", "M20", "Gini",
        # ...
    ],
)
x_tr, c_tr = ds.train_arrays()
x_va, c_va = ds.val_arrays()

# 2. Fit. n_features and n_cond come from the dataset; everything else
#    is a hyperparameter you can sweep from a notebook.
model = ConditionalFlow(
    n_features=ds.n_features,
    n_cond=ds.n_cond,
    n_epochs=10,
    batch_size=16384,
    learning_rate=5e-4,
)
model.train(
    x_tr, c_tr,
    val_x=x_va, val_c=c_va,
    checkpoint_dir="runs/r_band_v1",
)

# 3. Diagnostics: loss curves
h = model.history
fig, ax = plt.subplots(1, 2, figsize=(10, 3))
ax[0].plot(h["step_loss"]);                 ax[0].set_title("step loss")
ax[1].plot(h["epoch_train_loss"], label="train")
ax[1].plot(h["epoch_val_loss"],   label="val")
ax[1].legend(); ax[1].set_title("epoch loss")
print(f"best val loss {h['best_val_loss']:.4f} at epoch {h['best_epoch']}")

# 4. Anomaly score = -log p(x | c). Top 0.1% are candidates to inspect.
log_prob = model.evaluate(ds.x, ds.c, progress=True)
anomaly_score = -log_prob
thr = np.quantile(anomaly_score, 0.999)
candidate_idx = np.where(anomaly_score > thr)[0]

# 5. Latent embedding for HDBSCAN / UMAP downstream.
z = model.transform(ds.x, ds.c, progress=True)   # (N, 51), ~unit Gaussian

# 6. Persist. Save the model weights and the cond standardization stats so
#    the same flow can score new tiles consistently.
model.save("runs/r_band_v1/flow.eqx")
np.savez("runs/r_band_v1/cond_stats.npz",
         cond_mu=ds.cond_mu, cond_sd=ds.cond_sd)

# 7. Reload later (architecture must match -- same n_features / n_cond /
#    n_coupling_layers / n_bins / nn_width / nn_depth).
later = ConditionalFlow(n_features=ds.n_features, n_cond=ds.n_cond)
later.load("runs/r_band_v1/flow.eqx")
assert np.allclose(later.evaluate(ds.x[:1000], ds.c[:1000]),
                   log_prob[:1000])
```

## Notes / things to revisit

- `n_epochs=5` is a starting point. With ~100M sources you may converge in
  fewer epochs; switch to a step budget if useful.
- Batch size 8192 should fit comfortably in 8 GB; bump up if memory permits --
  larger batches give faster epoch wall-time on this kind of model.
- I/O: `FlowDataset` loads everything to RAM (~27 GB for your data). Fine
  with 200 GB. If you ever scale up, swap to an HDF5 chunked reader with
  prefetch.
- The flow learns p(x | c) including the imputed dimensions. Anomaly score
  for downstream use should subtract the contribution of masked dims --
  cheapest version: condition on the same mask but evaluate log-prob with
  imputed values resampled, average over a few draws. Worth a follow-up.
- Latent z from `model.transform(...)` is what goes into HDBSCAN/UMAP and the
  supervised anchoring step.

## Requirements

```
jax
jaxlib
flowjax
equinox
optax
h5py
numpy
tqdm   # optional; only used for progress bars
```
