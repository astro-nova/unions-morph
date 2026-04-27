#!/bin/bash
# =============================================================================
# Nibi sbatch: run scripts/morph_parallel_da.py on a single node / 90 cores
#
# Submit with:
#   sbatch --export=ALL,IMIN=0,IMAX=10 scripts/run_morph_parallel_da.sbatch
#
# One-time setup (on the LOGIN node, NOT inside a job):
#   PROJECT_DIR=/home/$USER/projects/def-taylor/$USER/unions-morph
#   module load StdEnv/2023 python/3.12.4 scipy-stack
#   module load cfitsio
#   virtualenv --no-download $PROJECT_DIR/venvs/morph
#   source $PROJECT_DIR/venvs/morph/bin/activate
#   pip install --no-index --upgrade pip
#   pip install --no-index -r $PROJECT_DIR/requirements.txt
#   pip install --no-index vos tqdm astropy scikit-image psutil
#   git clone https://github.com/astro-nova/statmorph-lsst $PROJECT_DIR/src/statmorph-lsst
#   pip install -e $PROJECT_DIR/src/statmorph-lsst     # editable, persists
#   deactivate
#   cadc-get-cert -u <your_cadc_username>              # writes ~/.ssl/cadcproxy.pem
#   # Proxy cert is valid 10 days; refresh on the login node when it expires.
# =============================================================================

#SBATCH --job-name=morph_da
#SBATCH --account=def-taylor
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=0
#SBATCH --time=1:00:00
#SBATCH --output=logs/morph_da_%j.out
#SBATCH --error=logs/morph_da_%j.err
#SBATCH --mail-user=liza.sazonova@uwaterloo.ca
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# -----------------------------------------------------------------------------
# 0. Inputs (pass --export=ALL,IMIN=...,IMAX=... to sbatch)
# -----------------------------------------------------------------------------
: "${IMIN:?IMIN must be set (e.g. --export=ALL,IMIN=0,IMAX=10)}"
: "${IMAX:?IMAX must be set}"

# Permanent locations under the def-taylor project allocation.
export PROJECT_DIR="${PROJECT_DIR:-/home/$USER/projects/def-taylor/$USER/unions-morph}"
VENV_DIR="${VENV_DIR:-$PROJECT_DIR/venvs/morph}"
STATMORPH_DIR="${STATMORPH_DIR:-$PROJECT_DIR/src/statmorph-lsst}"
CADC_CERT="${CADC_CERT:-/home/esazonov/.ssl/cadcproxy.pem}"

# Node-local fast disk for all intermediate I/O (cleaned at job end).
export WORK_DIR="$SLURM_TMPDIR/scratch"
echo "WORK_DIR=$WORK_DIR"
echo "SLURM_TMPDIR=$SLURM_TMPDIR"
ls -ld "$WORK_DIR" 2>&1 || true

# Final outputs land here (persistent, backed up).
export OUT_DIR="${OUT_DIR:-$PROJECT_DIR/catalogs}"
export LOG_DIR="${LOG_DIR:-$PROJECT_DIR/logs}"

echo "Job $SLURM_JOB_ID on $(hostname) at $(date)"
echo "IMIN=$IMIN  IMAX=$IMAX"
echo "PROJECT_DIR=$PROJECT_DIR"
echo "WORK_DIR=$WORK_DIR   (node-local scratch)"
echo "OUT_DIR=$OUT_DIR     (persistent)"

# -----------------------------------------------------------------------------
# 1. Preflight: required files must exist (fail fast, don't waste allocation)
# -----------------------------------------------------------------------------
[[ -d "$PROJECT_DIR"                                ]] || { echo "missing $PROJECT_DIR";        exit 1; }
[[ -d "$VENV_DIR"                                   ]] || { echo "missing venv $VENV_DIR (see one-time setup in header)"; exit 1; }
[[ -d "$STATMORPH_DIR"                              ]] || { echo "missing $STATMORPH_DIR (git clone statmorph-lsst in one-time setup)"; exit 1; }
[[ -f "$CADC_CERT"                                  ]] || { echo "missing CADC cert at $CADC_CERT; run 'cadc-get-cert -u <user>' on login node"; exit 1; }
[[ -f "$PROJECT_DIR/scripts/morph_parallel_da.py"   ]] || { echo "missing pipeline script";     exit 1; }
[[ -f "$PROJECT_DIR/catalogs/tiles_r.csv"           ]] || { echo "missing tile list";           exit 1; }

mkdir -p "$OUT_DIR" "$LOG_DIR" "$WORK_DIR" logs

# -----------------------------------------------------------------------------
# 2. Modules + persistent venv (built once on login node; only activated here)
# -----------------------------------------------------------------------------
module load StdEnv/2023 python/3.12.4 scipy-stack
source "$VENV_DIR/bin/activate"

# vos reads the proxy cert from this path
export VOSPACE_CERTFILE="$CADC_CERT"

# 90 worker processes × 1 BLAS thread each — no oversubscription
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# -----------------------------------------------------------------------------
# 3. Run. Python reads PROJECT_DIR / WORK_DIR / OUT_DIR / LOG_DIR from env,
#    so no hardcoded paths leak in. All per-tile FITS, weightmaps, segmaps,
#    and intermediate morph_batch_*.csv files are written to WORK_DIR, which
#    is $SLURM_TMPDIR — node-local NVMe, wiped at job end. Final per-tile
#    morph_tile_*.fits are written to OUT_DIR (project allocation).
# -----------------------------------------------------------------------------
cd "$PROJECT_DIR/scripts"
python -u morph_parallel_da.py --imin "$IMIN" --imax "$IMAX"

# -----------------------------------------------------------------------------
# 4. Belt-and-braces: if the Python script crashed mid-tile leaving a
#    morph_batch_*.csv in node-local scratch, copy it to OUT_DIR before the
#    node disk is wiped.
# -----------------------------------------------------------------------------
shopt -s nullglob
for f in "$WORK_DIR"/morph_batch_*.csv; do
  cp -v "$f" "$OUT_DIR/"
done

echo "Job $SLURM_JOB_ID done at $(date)"
