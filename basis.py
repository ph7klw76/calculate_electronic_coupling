import os
from collections import Counter, OrderedDict

# === USER SETTINGS ============================================================
# Path to your Gaussian input (.gjf/.com) file
gjf_path = "517.gjf"
# ==============================================================================
# def2‑SVP number of *spherical* atomic basis functions for elements of interest
DEF2SVP_BF = {
    "H": 5,   # 2s + 1p  → 2×1 + 1×3
    "B": 14,  # 3s + 2p + 1d
    "C": 14,
    "N": 14,
    "F": 14,
    "Si": 18, # 4s + 3p + 1d
    "P": 18,
    "S": 18,
    "Cl": 18,
    "Br": 32 # 5s + 4p + 3d  (no f‑functions in def2‑SVP column)
}

# --- helper to pull element list from coordinates ----------------------------
def parse_gjf_elements(file_text: str) -> list[str]:
    """Return a list of element symbols appearing in the coordinate block
    (after charge‑multiplicity line, until the next blank line)."""
    lines = file_text.splitlines()
    try:
        start = lines.index("0 1") + 1
    except ValueError:
        raise RuntimeError("Could not find the charge & multiplicity line (0 1).")
    elements = []
    for ln in lines[start:]:
        if not ln.strip():          # blank line terminates the coordinate section
            break
        sym = ln.strip().split()[0]
        elements.append(sym)
    return elements

# --- crunch the numbers -------------------------------------------------------
with open(gjf_path, "r") as fh:
    text = fh.read()

elems_in_file = parse_gjf_elements(text)
counts = Counter(e for e in elems_in_file if e in DEF2SVP_BF)

# per‑element basis‑function contributions
per_elem_bf = OrderedDict(
    (el, counts.get(el, 0) * DEF2SVP_BF[el]) for el in DEF2SVP_BF.keys()
)
total_bf = sum(per_elem_bf.values())

# --- print a concise report ---------------------------------------------------
print(f"File : {os.path.basename(gjf_path)}\n")
print("Element  Atoms  BF/atom  Total BF")
for el, tot in per_elem_bf.items():
    print(f"{el:>3} {counts.get(el,0):>6} {DEF2SVP_BF[el]:>8} {tot:>9}")
print("-" * 30)
print(f"TOTAL BASIS FUNCTIONS = {total_bf}")


