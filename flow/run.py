"""Example end-to-end run. Edit feature lists / paths for your catalog."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .data import FlowDataset
from .model import ConditionalFlow


def main():
    out_dir = Path("./flow_runs/run0")
    out_dir.mkdir(parents=True, exist_ok=True)

    ds = FlowDataset(
        data_path="/path/to/catalog.h5",
        feature_names=[
            # 51 feature names, in the column order you want for x
        ],
        masked_feature_names=[
            # 48 of those that have a stored mask in HDF5
        ],
        # cond_var_names defaults to: r_mag, psf_fwhm, sky_median, sky_sigma, snr
    )
    print(f"x:    {ds.x.shape}")
    print(f"mask: {ds.mask.shape}  ({int(ds.mask.sum())} masked entries)")
    print(f"c:    {ds.c.shape}  (cond_dim = {ds.n_cond})")

    x_train, c_train = ds.train_arrays()
    x_val,   c_val   = ds.val_arrays()

    model = ConditionalFlow(
        n_features=ds.n_features,
        n_cond=ds.n_cond,
        n_coupling_layers=10,
        n_bins=8,
        spline_bound=10.0,
        nn_width=256,
        nn_depth=3,
        batch_size=8192,
        learning_rate=1e-3,
        n_epochs=5,
    )

    model.train(
        x_train, c_train,
        val_x=x_val, val_c=c_val,
        checkpoint_dir=out_dir,
    )

    val_lp = model.evaluate(x_val, c_val)
    np.save(out_dir / "val_log_prob.npy", val_lp)
    np.savez(out_dir / "cond_stats.npz",
             cond_mu=ds.cond_mu, cond_sd=ds.cond_sd)
    np.savez(out_dir / "history.npz", **{
        k: np.asarray(v) for k, v in model.history.items()
        if isinstance(v, list)
    })


if __name__ == "__main__":
    main()
