from collections import Counter, OrderedDict
import os, sys

# >>> USER SETTINGS <<<
gjf_path   = "your_file.gjf"       # any .gjf / .com file
basis_name = "def2-tzv"           # choose: sto-3g | def2-svp | def2-tzv
# ----------------------

BASIS_COUNTS = {
    "sto-3g": {  "H": 1,
                 "B": 5, "C": 5, "N": 5, "F": 5,
                 "Si": 9, "P": 9, "S": 9, "Cl": 9,
                 "Br": 18 },

    "def2-svp": { "H": 5,
                  "B":14, "C":14, "N":14, "F":14,
                  "Si":18, "P":18, "S":18, "Cl":18,
                  "Br":32 },

    "def2-tzv": { "H":14,
                  "B":30, "C":30, "N":30, "F":30,
                  "Si":34, "P":34, "S":34, "Cl":34,
                  "Br":48 }
}

def parse_elements(text):
    """Return element symbols from the Cartesian block."""
    lines = text.splitlines()
    # find charge–multiplicity line (two ints)
    start = next(i for i,l in enumerate(lines)
                 if len(l.split())>=2 and
                    all(t.lstrip('+-').isdigit() for t in l.split()[:2])) + 1
    elems = []
    for l in lines[start:]:
        if not l.strip(): break
        symbol = ''.join(filter(str.isalpha, l.split()[0]))
        elems.append(symbol)
    return elems

key = basis_name.lower()
if key not in BASIS_COUNTS:
    sys.exit(f"Unknown basis: {basis_name}")

with open(gjf_path) as f:
    atom_syms = parse_elements(f.read())

counts = Counter(s for s in atom_syms if s in BASIS_COUNTS[key])
totals = OrderedDict((e, counts[e]*BASIS_COUNTS[key][e])
                     for e in BASIS_COUNTS[key])

print(f"{os.path.basename(gjf_path)}  |  {basis_name.upper()}")
print("Elem  Atoms  BF/atom  Total")
for e, t in totals.items():
    print(f"{e:>3} {counts.get(e,0):>6} {BASIS_COUNTS[key][e]:>8} {t:>7}")
print("-"*34)
print("TOTAL BASIS FUNCTIONS =", sum(totals.values()))


