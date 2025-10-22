#!/usr/bin/env bash
set -euo pipefail
# you need to vhange Omega in the gaussian file here

if [ $# -lt 4 ] || [ $# -gt 6 ]; then
  echo "Usage: $0 <gro-file> <index.ndx> <distances.txt> <target-distance> [tolerance] [nbo-value]"
  echo "Example: $0 output_whole.gro index.ndx nearest_molecule_distances.txt 0.40 0.40 0062100000"
  exit 1
fi

gro_file="$1"
index_file="$2"
distances_file="$3"
target="$4"
tol="${5:-0.0001}"
nbo="${6:-0043400000}"


# 1) Filter the distance file down to only those pairs whose AvgDist is within ±tol of target
filtered_list="filtered_pairs.txt"
awk -v tar="$target" -v tol="$tol" '
  BEGIN { FS = "[ \t]+" }
  NR>1 {
    avg = $3
    if (avg >= tar - tol && avg <= tar + tol) {
      print $1, $2, avg
    }
  }
' "$distances_file" > "$filtered_list"

echo "Found $(wc -l < "$filtered_list") pairs within $tol of $target nm ? processing…"

### Stage 1: extract individual PDBs with gmx_mpi trjconv
while read -r mol_i mol_j avg; do
    echo "? Mol $mol_i & $mol_j (avg=$avg nm)"

    output_file1="${mol_i}.pdb"
    output_file2="${mol_j}.pdb"

    # Only extract if the file doesn't already exist
    if [ ! -f "$output_file1" ]; then
        echo "  - extracting $output_file1"
        printf "r_%s\nq\n" "$mol_i" | gmx_mpi trjconv -f "$gro_file" -s "$gro_file" -n "$index_file" -o "$output_file1"
    else
        echo "  - skipping $output_file1 (already exists)"
    fi

    if [ ! -f "$output_file2" ]; then
        echo "  - extracting $output_file2"
        printf "r_%s\nq\n" "$mol_j" | gmx_mpi trjconv -f "$gro_file" -s "$gro_file" -n "$index_file" -o "$output_file2"
    else
        echo "  - skipping $output_file2 (already exists)"
    fi

done < "$filtered_list"

### Stage 2: build .xyz files
output_dir_xyz="output_xyz_files"
mkdir -p "$output_dir_xyz"

process_pdb() {
    awk '
    $1=="ATOM" {
      atom = substr($3,1,1)
      printf "%s %s %s %s\n", atom, $6, $7, $8
    }' "$1"
}

while read -r mol_i mol_j avg; do
    file1="${mol_i}.pdb"
    file2="${mol_j}.pdb"
    out_xyz="${output_dir_xyz}/${mol_i}_${mol_j}.xyz"

    echo "? building XYZ for ${mol_i}_${mol_j}"
    data1=$(process_pdb "$file1")
    data2=$(process_pdb "$file2")

    {
      total=$(( $(wc -l <<<"$data1") + $(wc -l <<<"$data2") ))
      echo "$total"
      echo "extracted"
      echo "$data1"
      echo "$data2"
    } > "$out_xyz"
done < "$filtered_list"

### Stage 3: build .gjf files
output_dir_gjf="output_gjf_files"
mkdir -p "$output_dir_gjf"
summary_file="summary.txt"
: > "$summary_file"

extract_unique_letters() {
    awk '
    $1=="ATOM" { print substr($3,1,1) }
    ' "$1" | sort -u | xargs
}

while read -r mol_i mol_j avg; do
    file1="${mol_i}.pdb"
    file2="${mol_j}.pdb"
    pair_dir="${output_dir_gjf}/${mol_i}_${mol_j}"
    mkdir -p "$pair_dir"

    # skip if PDBs missing
    if [[ ! -f "$file1" || ! -f "$file2" ]]; then
      echo "WARNING: Missing PDB for $mol_i or $mol_j" >> "$summary_file"
      continue
    fi

    data1=$(process_pdb "$file1")
    data2=$(process_pdb "$file2")
    letters1=$(extract_unique_letters "$file1")
    letters2=$(extract_unique_letters "$file2")
    combined_letters=$(printf "%s %s" "$letters1" "$letters2" | xargs -n1 | sort -u | xargs)

    pair_gjf="${pair_dir}/${mol_i}_${mol_j}.gjf"
    gjf1="${pair_dir}/${mol_i}.gjf"
    gjf2="${pair_dir}/${mol_j}.gjf"

    # Pair .gjf
    {
      echo "%nprocshared=16"
      echo "%mem=32GB"
      echo "# wb97x/gen guess=huckel nosymm pop=nboread IOp(3/107=${nbo},3/108=${nbo})"
      echo "# scf=(direct,nosymm)"
      echo ""
      echo "Pair_${mol_i}_${mol_j}"
      echo ""
      echo "0 1"
      echo "$data1"
      echo "$data2"
      echo ""
      echo "$combined_letters 0"
      echo "Def2SVP"
      echo "****"
      echo ""
      echo "\$NBO SAO=w53 FAO=W54 \$END"
    } > "$pair_gjf"

    # Single .gjf for mol_i
    {
      echo "%nprocshared=16"
      echo "%mem=32GB"
      echo "# wb97x/gen nosymm punch(MO) IOp(3/107=${nbo},3/108=${nbo})"
      echo "# scf=(direct,nosymm)"
      echo ""
      echo "Mol_${mol_i}"
      echo ""
      echo "0 1"
      echo "$data1"
      echo ""
      echo "$letters1 0"
      echo "Def2SVP"
      echo "****"
      echo ""
    } > "$gjf1"

    # Single .gjf for mol_j
    {
      echo "%nprocshared=16"
      echo "%mem=32GB"
      echo "# wb97x/gen nosymm punch(MO) IOp(3/107=${nbo},3/108=${nbo})"
      echo "# scf=(direct,nosymm)"
      echo ""
      echo "Mol_${mol_j}"
      echo ""
      echo "0 1"
      echo "$data2"
      echo ""
      echo "$letters2 0"
      echo "Def2SVP"
      echo "****"
      echo ""
    } > "$gjf2"

    # ...
done < "$filtered_list"

echo "All done. Summary in $summary_file."