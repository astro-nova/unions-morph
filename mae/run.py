"""Example end-to-end run. Edit feature lists / paths for your catalog."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .data import MAEDataset
from .model import MaskedAutoencoder


def main():
    out_dir = Path("./mae_runs/run0")
    out_dir.mkdir(parents=True, exist_ok=True)

    ds = MAEDataset(
        data_path="/path/to/catalog.h5",
        feature_names=[
            # 51 feature names, in the column order you want for x
        ],
        masked_feature_names=[
            # 48 of those that have a stored mask in HDF5
        ],
        # cond_var_names defaults to: mag_ellip, fwhm, sky_median, sky_sigma, sn_per_pixel
    )
    print(f"x:    {ds.x.shape}")
    print(f"mask: {ds.mask.shape}  ({int(ds.mask.sum())} masked entries)")
    print(f"c:    {ds.c.shape}  (cond_dim = {ds.n_cond})")

    x_train, m_train, c_train = ds.train_arrays()
    x_val,   m_val,   c_val   = ds.val_arrays()

    model = MaskedAutoencoder(
        n_features=ds.n_features,
        n_cond=ds.n_cond,
        latent_dim=16,
        hidden_dim=256,
        n_hidden_layers=3,
        mask_ratio=0.5,
        batch_size=8192,
        learning_rate=1e-3,
        n_epochs=10,
    )

    model.train(
        x_train, m_train, c_train,
        val_x=x_val, val_mask=m_val, val_c=c_val,
        checkpoint_dir=out_dir,
    )

    val_score = model.anomaly_score(x_val, m_val, c_val)
    np.save(out_dir / "val_anomaly.npy", val_score)
    np.savez(out_dir / "cond_stats.npz",
             cond_mu=ds.cond_mu, cond_sd=ds.cond_sd)
    np.savez(out_dir / "history.npz", **{
        k: np.asarray(v) for k, v in model.history.items()
        if isinstance(v, list)
    })


if __name__ == "__main__":
    main()
