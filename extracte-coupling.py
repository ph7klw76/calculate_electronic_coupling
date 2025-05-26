# %% 
# --- USER CONFIGURATION (edit these!) ---
input_dir    = './'          # directory containing output*_*.t
nn_file      = './nearest_neighbours_with_clusters-acceptor.txt'
output_file  = './e-coupling.txt'           # where to save your tab‐delimited output
search_key   = 't(181,181)='                     # substring to look for in each .t file

# %% 
# --- IMPORTS ---
import os
import re
import glob
import pandas as pd

# %% 
# --- FUNCTION DEFINITIONS ---
def extract_value_after_key(filepath, search_key):
    """
    Open the file and find the line containing `search_key`.
    Return the float after the search_key string, or None.
    """
    with open(filepath, 'r') as f:
        for line in f:
            if search_key in line:
                try:
                    return abs(float((line.split(search_key)[1].strip()).split('eV')[0]))
                except ValueError:
                    return None
    return None

def parse_filename(filename):
    """
    From 'output12_345.t' extract (12, 345), else (None, None).
    """
    m = re.match(r'^output(\d+)_(\d+)\.t$', filename)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None

# %% 
# --- MAIN LOGIC ---
# 1. Gather all .t files
pattern = os.path.join(input_dir, 'output*_*.t')
t_files = glob.glob(pattern)
if not t_files:
    raise FileNotFoundError(f"No files found matching {pattern!r}")

# 2. Extract values from each file
records = []
for path in t_files:
    fname = os.path.basename(path)
    mol_i, mol_j = parse_filename(fname)
    if mol_i is None:
        print(f"Skipping invalid filename: {fname}")
        continue

    val = extract_value_after_key(path, search_key)
    if val is None:
        print(f"Warning: '{search_key}' not found or unparsable in {fname}")
    records.append({
        'Mol_i': mol_i,
        'Mol_j_nearest': mol_j,
        'extracted_value': val
    })

df_vals = pd.DataFrame(records)
print(f"Extracted values from {len(df_vals)} files.")

# 3. Load the nearest‐neighbours lookup table
if not os.path.exists(nn_file):
    raise FileNotFoundError(f"Lookup file not found: {nn_file!r}")
df_nn = pd.read_csv(nn_file, sep=r'\s+', engine='python')

# 4. Merge to pull in Distance_nm
df_merged = pd.merge(
    df_vals,
    df_nn[['Mol_i', 'Mol_j_nearest', 'Distance_nm']],
    on=['Mol_i', 'Mol_j_nearest'],
    how='left'
)

# 5. Report any missing lookups
missing = df_merged[df_merged['Distance_nm'].isna()]
if not missing.empty:
    print("Warning: missing Distance_nm for these pairs:")
    print(missing[['Mol_i', 'Mol_j_nearest']].to_string(index=False))

# 6. Save to tab‐delimited file
df_merged.to_csv(output_file, sep='\t', index=False)
print(f"Saved merged results ({len(df_merged)} rows) to:\n  {output_file}")
