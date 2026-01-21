#!/usr/bin/env python
# submit_reduction.py - Python client-based submission
# Script to run `morph_parallel.py` in headless mode on CANFAR
from canfar.sessions import Session
from datetime import datetime

# Initialize session manager
session = Session()

# Set job parameters
job_name = f"unions-morphology-{datetime.now().strftime('%Y%m%d')}"
image = "images.canfar.net/skaha/astroflow:latest"
project="/arc/home/esazonova/unions-morph"
# data_path = f"{project}/data/{datetime.now().strftime('%Y%m%d')}"

imin = 146
imax = 166

# Or submit fixed job (guaranteed resources by specifying cores/ram)
job_ids = session.create(
    name=job_name,
    image=image,
    cores=16,
    ram=64,  # Having cores/ram makes it a fixed session
    cmd="python",
    args=f"{project}/scripts/morph_parallel_headless.py"
)

print(f"Submitted job(s): {job_ids}")
session.logs(job_ids, verbose=True)