#!/bin/bash -l
  
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --job-name=submit
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:59:59

#!/usr/bin/env bash
###############################################################################
# batch_submit_all.sh
#
# Recursively submit every  submit_job.sh  found below a given directory.
# Usage examples
#   ./batch_submit_all.sh            # search from the current directory
#   ./batch_submit_all.sh  /path/to/workdir
#
# Exit codes:
#   0 – all jobs queued (or none found)
#   1 – sbatch not in $PATH
#   2 – supplied directory does not exist
###############################################################################

set -euo pipefail

# ---------- sanity checks ----------------------------------------------------
command -v sbatch >/dev/null 2>&1 || { echo "ERROR: sbatch not found in \$PATH"; exit 1; }

ROOT_DIR="${1:-.}"                         # default: current directory
[[ -d "$ROOT_DIR" ]] || { echo "ERROR: '$ROOT_DIR' is not a directory"; exit 2; }

# ---------- core loop --------------------------------------------------------
# Use -print0 / -d '' to cope with spaces, tabs, new-lines in file names.
find "$ROOT_DIR" -type f -name 'submit_job.sh' -print0 |
  while IFS= read -r -d '' jobfile; do
    jobdir="$(dirname "$jobfile")"
    jobname="$(basename "$jobfile")"

    echo "Submitting $jobfile …"
    (
      cd "$jobdir"
      sbatch "./$jobname"
    )
  done

echo "Done."
