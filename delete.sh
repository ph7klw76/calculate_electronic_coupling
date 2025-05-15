#!/bin/bash -l
  
#SBATCH --partition=cpu-epyc-genoa
#SBATCH --job-name=delete
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --hint=nomultithread
#SBATCH --time=00-00:59:00

# Load dependencies and libraries
source /app/hpcx/2.17.1/hpcx-init.sh
hpcx_load
 

# Path to the file listing the directories to delete
LIST_FILE="file.txt"

# Check that the list file exists
if [[ ! -f "$LIST_FILE" ]]; then
  echo "Error: $LIST_FILE not found."
  exit 1
fi

# Read each line (directory name) and delete it
while IFS= read -r dir; do
  # Skip empty lines
  [[ -z "$dir" ]] && continue

  if [[ -d "$dir" ]]; then
    echo "Deleting directory: $dir"
    rm -rf -- "$dir"
  else
    echo "Directory not found, skipping: $dir"
  fi
done < "$LIST_FILE"
