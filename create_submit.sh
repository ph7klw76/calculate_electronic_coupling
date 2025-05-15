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
set -euo pipefail
IFS=$'\n\t'

# =============================================================================
# generate_and_submit.sh
#
# For each subfolder in ./output_gjf_files/, this script will:
#   1. Find all .gjf files
#   2. Sort them by string length (shortest ? longest)
#   3. Take the first three, strip “.gjf” to get the base names
#   4. Generate a Slurm batch script (submit_job.sh) with detailed comments
#   5. Make it executable and submit it with `sbatch`
# =============================================================================

BASE_DIR="./output_gjf_files/"

# Loop over every entry in BASE_DIR
for folder in "$BASE_DIR"/*/; do
  # Skip anything that isn’t a directory
  [ -d "$folder" ] || continue

  echo "Processing folder: $folder"

  # ---------------------------------------------------------------------------
  # 1) Gather all .gjf files directly under this folder
  #    - `find -maxdepth 1` ensures we don’t recurse deeper
  #    - `-printf '%f\n'` prints just the filename, not the full path
  # ---------------------------------------------------------------------------
  mapfile -t gjf_files < <(
    find "$folder" -maxdepth 1 -type f -name '*.gjf' -printf '%f\n' \
    | awk '{ printf "%d %s\n", length($0), $0 }' \
    | sort -n \
    | cut -d' ' -f2-
  )

  # If we found fewer than 3 .gjf files, skip this folder
  if [ "${#gjf_files[@]}" -lt 3 ]; then
    echo "Skipping (found <3 .gjf files)"
    continue
  fi

  # ---------------------------------------------------------------------------
  # 2) Take the first three filenames (shortest ? longest) and strip “.gjf”
  # ---------------------------------------------------------------------------
  file1="${gjf_files[0]%.gjf}"
  file2="${gjf_files[1]%.gjf}"
  file3="${gjf_files[2]%.gjf}"

  # Use the longest (third) base name as the Slurm job-name
  name="$file3"

  # Path to the batch script we’re about to generate
  submit_script="$folder/submit_job.sh"

  # ---------------------------------------------------------------------------
  # 3) Write out submit_job.sh with plenty of inline comments
  # ---------------------------------------------------------------------------
  cat > "$submit_script" <<EOF
#!/bin/bash -l

# -----------------------------------------------
# Slurm batch options
# -----------------------------------------------
#SBATCH --partition=cpu-epyc-genoa    # queue/partition name
#SBATCH --job-name=${name}            # job name shown in squeue
#SBATCH --output=${name}.out          # STDOUT ? ${name}.out
#SBATCH --error=${name}.err           # STDERR ? ${name}.err
#SBATCH --mem=40G                     # total RAM
#SBATCH --nodes=1                     # number of nodes
#SBATCH --ntasks=16                   # total tasks (MPI ranks)
#SBATCH --cpus-per-task=1             # threads per task
#SBATCH --qos=long                    # quality of service
#SBATCH --hint=nomultithread          # disable hyper-threading
#SBATCH --time=02-23:59:59            # max runtime (DD-HH:MM:SS)

# -----------------------------------------------
# Load Gaussian and source its environment
# -----------------------------------------------
module load gaussian/g09
source \$g09profile

# -----------------------------------------------
# 1) Run the first (shortest-named) geometry
# -----------------------------------------------
g09 <${file1}.gjf> ${file1}.log
sleep 5
# Rename the binary scratch output fort.7 for safekeeping
mv fort.7 fort.7-${file1}

# -----------------------------------------------
# 2) Run the second geometry
# -----------------------------------------------
g09 <${file2}.gjf> ${file2}.log
sleep 5
mv fort.7 fort.7-${file2}

# -----------------------------------------------
# 3) Run the third (longest-named) geometry
# -----------------------------------------------
g09 <${file3}.gjf> ${file3}.log

EOF

  # Make it executable
  chmod +x "$submit_script"

done

echo "All done."
