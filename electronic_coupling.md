# Understanding Electronic Coupling Between Molecules: Theory and Applications

Electronic coupling is a cornerstone concept in the field of molecular electronics, photochemistry, and biophysics. It quantifies the quantum mechanical interaction between electronic states of two molecules or molecular fragments and is essential for understanding processes like charge transfer, energy transfer, and electron transport. 

In this post, we delve into the theoretical framework of electronic coupling, discuss how it is calculated using molecular orbital data, and explore its applications in science and technology.

---

## What is Electronic Coupling?

Electronic coupling, often denoted as $V$, measures how strongly the electronic states of two molecular fragments interact. Mathematically, it is expressed as:

$$
V = \langle \psi_A \mid \hat{H} \mid \psi_B \rangle
$$

Here:
- $\psi_A$ and $\psi_B$: Molecular wavefunctions of fragment A and B, respectively.
- $\hat{H}$: The electronic Hamiltonian, governing interactions within the molecular system.

In simpler terms, $V$ describes how well electrons or excitons "couple" between two states or regions in a molecular system.

---

## Theoretical Framework

### 1. Fock and Overlap Matrices

Electronic coupling calculations rely on the Fock ($F$) and overlap ($S$) matrices derived from quantum chemical calculations:

- **Fock Matrix ($F$)**: Represents the total electronic interaction in the system, including Coulombic, exchange, and external field effects.
- **Overlap Matrix ($S$)**: Quantifies the spatial overlap of molecular orbitals, influencing the strength of electronic coupling.

The coupling strength is calculated as:

$$
V = \sum_{i,j} C_{Ai} F_{ij} C_{Bj}
$$

Where:
- $C_{Ai}$: Molecular orbital coefficient $i$ for fragment A.
- $C_{Bj}$: Molecular orbital coefficient $j$ for fragment B.
- $F_{ij}$: Element of the Fock matrix describing the interaction between orbitals $i$ and $j$.

---

### 2. Energy and Orbital Alignment

The magnitude of $V$ depends on the energy alignment of interacting orbitals. If the orbital energies $E_A$ and $E_B$ of fragments A and B are closely aligned, the coupling is typically stronger.

---

### 3. Two-State Model

A simplified two-state model is often used to calculate $V$:

$$
V = \frac{H_{AB} - \frac{1}{2}(E_A + E_B) S_{AB}}{1 - S_{AB}^2}
$$

Where:
- $H_{AB}$: The electronic coupling matrix element.
- $E_A$, $E_B$: Orbital energies of fragments A and B.
- $S_{AB}$: Overlap integral between the orbitals.

---

## How to Calculate Electronic Coupling?

### Input Requirements

To calculate electronic coupling between two fragments, the following inputs are needed:
- **Fock Matrix**: Represents electronic interactions across basis functions.
- **Overlap Matrix**: Describes the overlap between molecular orbitals.
- **Molecular Orbital Coefficients**: For both fragment A and fragment B.

These inputs are typically derived from quantum chemistry software like Gaussian, ORCA, or NWChem.

---

### Example Setup

#### Input Structure
- **InFile**: Contains the Fock and overlap matrices and molecular orbital coefficients for fragments A and B.

Example:
```plaintext
fock       filename
overlap    filename
molOrbA    filename
molOrbB    filename
```
paramFile: Specifies calculation parameters, such as the number of basis functions and orbital ranges.

```plaintext
nBasisFunctsA    number
nBasisFunctsB    number
orbitalsA        number number
orbitalsB        number number
```
## OutFilePrefix

Specifies the prefix for output files containing overlap, interaction energy, and electronic coupling values.

---

## Output

The calculation produces two key outputs:

1. **Total Interaction Energy**: Includes all interaction components.
2. **Electronic Coupling ($V$)**: The strength of interaction between specified molecular orbitals.

---

## Applications of Electronic Coupling

### 1. Organic Photovoltaics

In solar cells, electronic coupling controls the charge transfer between donor and acceptor molecules, influencing power conversion efficiency.

### 2. Photosynthesis and Light-Harvesting

Understanding coupling in excitonic energy transfer pathways sheds light on the efficiency of natural photosynthetic systems and synthetic light-harvesting complexes.

### 3. Molecular Electronics

Electronic coupling determines the conductivity of molecular wires and junctions, enabling the design of nanoscale electronic devices.

### 4. Catalysis

Coupling between reactant and catalyst orbitals influences reaction rates and selectivity.

---

## Challenges and Insights

While electronic coupling provides critical insights, its calculation is sensitive to:

1. **Choice of Basis Set**: Larger basis sets provide more accurate results but increase computational cost.
2. **Molecular Orientation**: Coupling is highly dependent on the relative spatial arrangement of fragments.
3. **Environmental Effects**: Solvent and polarization effects can significantly modulate coupling.

---
## Gaussian File

```plaintext
%mem=16GB
%nprocshared=8     
# gen nosymm punch(MO)   
# scf=(direct,nosymm) 

   carbazol

   0 1
   Cl        0.000000    0.000000    0.000000
   Cl        0.000000    0.000000    1.640000

Cl 0
 3-21G       
****
```

Cl 0: Specifies the element and its associated basis set.
3-21G: A minimal split-valence basis set applied to chlorine atoms.
****: Marks the end of the basis set definition.

The file (MO1.com) perform a single-point energy calculation, where the program computes the electronic energy and molecular orbitals for the given geometry (one molecule) without optimizing the structure.
Specific Outputs:
Energy of the molecule.
Molecular orbital coefficients (punched into a file, it create xxx.pun file).
Electron density and wavefunction information. 

The same can be done from MO2.com. Let say it looks like

```plaintext
%mem=16GB
%nprocshared=8    
# gen nosymm punch(MO)   
# scf=(direct,nosymm) 

   carbazol

   0 1
   F         2.500000    0.000000    0.000000
   F         2.500000    0.000000    1.340000

 F 0
 3-21G       
****
```


Note the geometries are extracted out from a pair of molecules (MO1 and MO2)
Let say the par is MO1_MO2.com

```plaintext
%mem=16GB           
%nprocshared=8    
#  gen guess=huckel nosymm pop=nboread 
# scf=(direct,nosymm)  

   Cb 1            

   0 1
Cl        0.000000    0.000000    0.000000
Cl        0.000000    0.000000    1.640000
F         2.500000    0.000000    0.000000
F         2.500000    0.000000    1.340000

 Cl F 0
  3-21G              
  ****

$NBO SAO=w53 FAO=W54 $END
```

## Keywords and Options

### 1. **`gen`**

Indicates that a general basis set will be provided in the input file (defined later in the file).

---

### 2. **`guess=huckel`**

Specifies that the initial guess for the molecular orbitals will use the Hückel method, which is particularly useful for:
- Initial wavefunction construction in systems with extended π-conjugation.
- Unusual bonding patterns.

---

### 3. **`nosymm`**

Disables the automatic symmetry detection and use of molecular symmetry in calculations. This ensures:
- All integrals and molecular orbitals are computed explicitly.
- Necessary adjustments for systems with low or broken symmetry.

---

### 4. **`pop=nboread`**

Requests population analysis using the Natural Bond Orbital (NBO) method. Additionally:
- Indicates that additional NBO options (such as SAO and FAO) will be provided in the input.

---

### 5. **`scf=(direct,nosymm)`**

- **`direct`**: Specifies that integrals needed for the self-consistent field (SCF) calculation are computed directly without being stored in memory. This is useful for:
  - Large systems.
  - Situations with limited disk space.

- **`nosymm`**: Ensures symmetry is ignored during the SCF calculation.

```plaintext
$NBO SAO=w53 FAO=W54 $END
```
---
## Additional Natural Bond Orbital (NBO) Analysis Options

- **`SAO=w53`**:  
  Requests specific analysis of Symmetric Atomic Orbitals (SAOs) with parameters written to the file `w53`.

- **`FAO=w54`**:  
  Requests additional analysis of Fragment Atomic Orbitals (FAOs) with parameters written to the file `w54`.

## What This Calculation Does

### 1. **Initial Guess**
- Uses the Hückel method to generate starting molecular orbitals.

---

### 2. **Self-Consistent Field (SCF) Calculation**
- Computes the electronic structure of the system using a direct SCF method without exploiting symmetry.

---

### 3. **Natural Bond Orbital (NBO) Analysis**
- Performs a detailed population analysis to evaluate bonding and electron distribution in terms of NBOs.
- Outputs specific information about:
  - Symmetric Atomic Orbitals (SAOs).
  - Fragment Atomic Orbitals (FAOs).

---

### 4. **Basis Set Application**
- Uses the **3-21G** basis set for all Cl and F atoms to approximate the molecular wavefunction.

---

### 5. **Output Files**
- **Wavefunction and Energy Information**: Standard Gaussian output file.
- **NBO Analysis Files**:
  - `w53`: Contains SAO-related data.
  - `w54`: Contains FAO-related data.

```python
import os
import numpy as np

HARTREE_TO_EV = 27.2114

def find_n_basis_functions(filename):
    """
    Search for the line containing 'NBasis=' in the given file
    and return the integer value following it.
    """
    filepath = os.path.join(os.getcwd(), filename)
    with open(filepath, 'r') as file:
        for line in file:
            if 'NBasis=' in line:
                return int(line.split()[1].strip())
    return None

def read_lower_triangular_matrix(filename, skip_lines=3):
    """
    Reads a lower-triangular matrix from a file, ignoring a given number of header lines.
    Returns a flat list of matrix elements.
    """
    elements = []
    with open(filename, 'r') as f:
        for idx, line in enumerate(f):
            if idx >= skip_lines:
                elements.extend(map(float, line.split()))
    return elements

def form_symmetric_matrix(elements, dim):
    """
    Given a 1D list of lower-triangular elements of a symmetric matrix,
    reconstruct the full symmetric (dim x dim) matrix.
    """
    mat = np.zeros((dim, dim))
    k = 0
    for i in range(dim):
        for j in range(i + 1):
            mat[i, j] = elements[k]
            k += 1
    # Symmetrize
    mat = mat + np.triu(mat.T, 1)
    return mat

def read_orbitals(filename, n_basis, coeffs_per_line=5, coeff_length=15):
    """
    Read molecular orbital coefficients from a formatted output file.
    The file contains 'Alpha' labels indicating the start of a new MO block.
    
    Parameters
    ----------
    filename : str
        Path to the MO coefficients file.
    n_basis : int
        Number of basis functions for this fragment.
    coeffs_per_line : int
        Number of coefficients on each line (for formatting).
    coeff_length : int
        Length of each coefficient field in the formatted file.
        
    Returns
    -------
    dict
        A dictionary keyed by orbital index with a list of float coefficients as values.
    """
    orbitals = {}
    orb_count = -1
    orb_flag = 0
    # Number of lines of coefficients for each orbital block
    n_orb_lines = (n_basis // coeffs_per_line) + 1

    with open(filename, 'r') as file:
        next(file)  # Skip the first line
        for line in file:
            cols = line.split()
            if orb_flag > 0:
                # Extract coefficients from fixed-width fields
                for i in range(coeffs_per_line):
                    coeff_str = line[i * coeff_length:(i + 1) * coeff_length].strip()
                    if coeff_str and 'Alpha' not in line:
                        # Convert Fortran D-notation (e.g., 1.234D+02) to E-notation (1.234e+02)
                        if 'D' in coeff_str:
                            base, power = coeff_str.split('D')
                            coeff_str = f"{float(base)}e{int(power)}"
                        orbitals[orb_count].append(float(coeff_str))
                orb_flag += 1
                # If we've read all lines for this orbital block, reset flag
                if orb_flag == n_orb_lines + 1:
                    orb_flag = 0

            # Check if this line indicates start of a new orbital block
            if len(cols) > 1 and cols[1] == 'Alpha':
                orb_count = int(cols[0])
                orbitals[orb_count] = []
                orb_flag = 1

    return orbitals

def parse_param_file(param_file):
    """
    Parse the parameter file, expecting:
      orbitalsA= start end
      orbitalsB= start end
      nBasisFunctsA= number
      nBasisFunctsB= number

    Returns: nBasisFunctsA, nBasisFunctsB, orbitalsA_begin, orbitalsA_end, orbitalsB_begin, orbitalsB_end
    """
    with open(param_file, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    
    # Expected format checking
    if len(lines) < 4:
        raise ValueError("Parameter file format is incorrect or incomplete.")

    # orbitalsA line
    if 'orbitalsA=' not in lines[0]:
        raise ValueError("Parameter file: Missing 'orbitalsA=' line.")
    orbitalsA_begin, orbitalsA_end = map(int, lines[0].split('=')[1].split())

    # orbitalsB line
    if 'orbitalsB=' not in lines[1]:
        raise ValueError("Parameter file: Missing 'orbitalsB=' line.")
    orbitalsB_begin, orbitalsB_end = map(int, lines[1].split('=')[1].split())

    # nBasisFunctsA line
    if 'nBasisFunctsA=' not in lines[2]:
        raise ValueError("Parameter file: Missing 'nBasisFunctsA=' line.")
    nBasisFunctsA = int(lines[2].split('=')[1])

    # nBasisFunctsB line
    if 'nBasisFunctsB=' not in lines[3]:
        raise ValueError("Parameter file: Missing 'nBasisFunctsB=' line.")
    nBasisFunctsB = int(lines[3].split('=')[1])

    return nBasisFunctsA, nBasisFunctsB, orbitalsA_begin, orbitalsA_end, orbitalsB_begin, orbitalsB_end

def parse_in_file(in_file):
    """
    Parse the input file expecting:
      Line 0: overlap filename
      Line 1: fock filename
      Line 2: molOrbA filename
      Line 3: molOrbB filename
    """
    with open(in_file, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]

    if len(lines) < 4:
        raise ValueError("Input file does not contain the required 4 lines.")

    return lines[0], lines[1], lines[2], lines[3]

def compute_coupling(orbCoeffsA, orbCoeffsB, O, F, nA, nB, orbA, orbB):
    """
    Compute the interaction (coupling) between MO orbA of fragment A and MO orbB of fragment B.
    
    Parameters
    ----------
    orbCoeffsA : dict
        Orbital coefficients for fragment A {orb_index: [coeffs]}
    orbCoeffsB : dict
        Orbital coefficients for fragment B {orb_index: [coeffs]}
    O : np.ndarray
        Overlap matrix (nA+nB x nA+nB)
    F : np.ndarray
        Fock matrix (nA+nB x nA+nB)
    nA : int
        Number of basis functions on fragment A
    nB : int
        Number of basis functions on fragment B
    orbA : int
        Orbital index for fragment A
    orbB : int
        Orbital index for fragment B
    
    Returns
    -------
    float
        Coupling t(orbA, orbB) in eV
    """
    # S(A,A), E(A,A)
    Saa = sum(orbCoeffsA[orbA][i]*O[i,j]*orbCoeffsA[orbA][j]
              for i in range(nA) for j in range(nA))
    Eaa = sum(orbCoeffsA[orbA][i]*F[i,j]*orbCoeffsA[orbA][j]
              for i in range(nA) for j in range(nA))
    Eaaev = Eaa * HARTREE_TO_EV

    # S(B,B), E(B,B)
    Sbb = sum(orbCoeffsB[orbB][i]*O[i+nA, j+nA]*orbCoeffsB[orbB][j]
              for i in range(nB) for j in range(nB))
    Ebb = sum(orbCoeffsB[orbB][i]*F[i+nA, j+nA]*orbCoeffsB[orbB][j]
              for i in range(nB) for j in range(nB))
    Ebbev = Ebb * HARTREE_TO_EV

    # S(A,B), E(A,B)
    Sab = sum(orbCoeffsA[orbA][i]*O[j+nA, i]*orbCoeffsB[orbB][j]
              for i in range(nA) for j in range(nB))
    Eab = sum(orbCoeffsA[orbA][i]*F[j+nA, i]*orbCoeffsB[orbB][j]
              for i in range(nA) for j in range(nB))
    Eabev = Eab * HARTREE_TO_EV

    # Compute coupling
    tABev = (Eabev - ((Eaaev+Ebbev)*Sab/2)) / (1 - Sab*Sab)

    return Saa, Eaa, Eaaev, Sbb, Ebb, Ebbev, Sab, Eab, Eabev, tABev

def coupling(argv):
    """
    Main function to compute coupling between MOs of two fragments.
    The arguments are expected to be a list:
      argv[1]: inFile name
      argv[2]: paramFile name
      argv[3]: outFile prefix
    """
    in_file = os.path.join(os.getcwd(), argv[1])
    param_file = os.path.join(os.getcwd(), argv[2])
    out_prefix = argv[3]

    out_full = os.path.join(os.getcwd(), out_prefix + '.out')
    out_short = os.path.join(os.getcwd(), out_prefix + '.t')

    # Parse input & parameter files
    overlap_file, fock_file, molOrbFileA, molOrbFileB = parse_in_file(in_file)
    nBasisFunctsA, nBasisFunctsB, orbitalsA_begin, orbitalsA_end, orbitalsB_begin, orbitalsB_end = parse_param_file(param_file)

    # Construct full file paths
    overlap_file = os.path.join(os.getcwd(), overlap_file)
    fock_file = os.path.join(os.getcwd(), fock_file)
    molOrbFileA = os.path.join(os.getcwd(), molOrbFileA)
    molOrbFileB = os.path.join(os.getcwd(), molOrbFileB)

    # Dimensions and matrix element count
    nA, nB = nBasisFunctsA, nBasisFunctsB
    nAB = nA + nB
    nMatrixElems = (nAB * (nAB + 1)) // 2
    print("No. of elements in overlap/Fock matrix =", nMatrixElems)

    # Read matrices
    fock_elements = read_lower_triangular_matrix(fock_file, skip_lines=3)
    F = form_symmetric_matrix(fock_elements, nAB)

    overlap_elements = read_lower_triangular_matrix(overlap_file, skip_lines=3)
    O = form_symmetric_matrix(overlap_elements, nAB)

    # Read orbital coefficients
    orbCoeffsA = read_orbitals(molOrbFileA, nBasisFunctsA)
    orbCoeffsB = read_orbitals(molOrbFileB, nBasisFunctsB)

    # Compute couplings and write results
    with open(out_full, 'w') as fOut, open(out_short, 'w') as fOut2:
        for orbA in range(orbitalsA_begin, orbitalsA_end + 1):
            for orbB in range(orbitalsB_begin, orbitalsB_end + 1):
                Saa, Eaa, Eaaev, Sbb, Ebb, Ebbev, Sab, Eab, Eabev, tABev = compute_coupling(
                    orbCoeffsA, orbCoeffsB, O, F, nA, nB, orbA, orbB
                )

                # Detailed output
                print("  ****************************************", file=fOut)
                print("", file=fOut)
                print(f"Interaction between MO {orbA} on fragment A and MO {orbB} on fragment B", file=fOut)
                print("", file=fOut)
                print(f"S(A,A) = {Saa} E(A,A) = {Eaa} au ({Eaaev} eV)", file=fOut)
                print(f"S(B,B) = {Sbb} E(B,B) = {Ebb} au ({Ebbev} eV)", file=fOut)
                print("", file=fOut)
                print(f"S(A,B) = {Sab} E(A,B) = {Eab} au ({Eabev} eV)", file=fOut)
                print("", file=fOut)
                print(f"t({orbA},{orbB})= {tABev} eV", file=fOut)
                print("", file=fOut)
                print("", file=fOut)

                # Short output (just the coupling)
                print(f"t({orbA},{orbB})= {tABev} eV", file=fOut2)


coupling(['a', 'inFile.in', 'paramFile.txt', 'output'])


```

# Python Code for Computing Electronic Coupling Between Molecular Orbitals

The Python code computes the electronic coupling between molecular orbitals (MOs) of two molecular fragments, referred to as fragment A and fragment B, using quantum chemistry data such as the Fock matrix, overlap matrix, and MO coefficients. This guide provides step-by-step instructions for running the code with detailed explanations of the workflow and input/output requirements.

---

## Purpose

The code calculates the following:

1. **Self-overlap ($S$)**:  
   Overlap integrals for molecular orbitals within fragments.

2. **Self-energy ($E$)**:  
   Energies of the molecular orbitals on each fragment.

3. **Cross-overlap ($S(A,B)$)**:  
   Overlap integral between molecular orbitals of the two fragments.

4. **Electronic Coupling ($t_{AB}$)**:  
   Quantum mechanical interaction between MOs of fragment A and fragment B.

# Input Files and Parameters

## 1. InFile: Input File Descriptions
The `InFile` specifies the file paths for:
- **Fock Matrix (`fockFile`)**: Contains electronic interaction information.
- **Overlap Matrix (`overlapFile`)**: Measures the spatial overlap of molecular orbitals.
- **Molecular Orbital Coefficients (`molOrbFileA`, `molOrbFileB`)**: For fragments A and B.

---

## 2. paramFile: Parameters for the Calculation
The `paramFile` specifies:
- **Number of basis functions**:
  - `nBasisFunctsA`: For fragment A.
  - `nBasisFunctsB`: For fragment B.
- **Molecular orbital ranges**:
  - `orbitalsA_begin` to `orbitalsA_end` (for fragment A).
  - `orbitalsB_begin` to `orbitalsB_end` (for fragment B).

paramFile looks like:

``` plaintext
 BasisFunctions1=	26   
 BasisFunctions2=	18   
 orbitalsA=	16 19
 orbitalsB=	9 10
```
there are 2 flourines at 3-21G there are 9 basis function with 2 flourines you have 18.
# Basis Functions in 3-21G for Fluorine (F)

## Calculation of Total Basis Functions

### 1. **Core Orbitals**
- The core orbitals for fluorine include the **1s orbital**:
  - **1s orbital**: Represented by 1 contracted function.
- **Total core functions**:
  $$
  1
  $$

---

### 2. **Valence Orbitals**
- The valence shell for fluorine consists of the **2s** and **2p** orbitals:
  - **2s orbital**: Split into 2 functions (1 inner + 1 outer).
  - **2p orbitals**: Each $p$ orbital ($p_x$, $p_y$, $p_z$) is split into 2 functions (1 inner + 1 outer).

- **Total valence functions**:
  $$
  2 + 3 \times 2 = 8
  $$

---

### 3. **Total Basis Functions**
- **Core functions**:
  $$
  1
  $$
  
- **Valence functions**:
  $$
  8
  $$
- **Total basis functions for F (3-21G)**:
  $$
  1 + 8 = 9
  $$

The number can be obtained from fort.7 automatically generated by the gaussian code
---

## 3. outFilePrefix: Output File Prefix
Two output files are generated:
1. **Detailed Results (`filename.out`)**: Includes overlaps, energies, and coupling information.
2. **Summary (`filename.t`)**: Contains only the coupling values.

---

# Workflow

## 1. Helper Function: `findnBasisFuncts`
Reads the basis set information (`nBasisFunctsA`, `nBasisFunctsB`) from a file:
- Locates the line containing `NBasis=`.
- Extracts and returns the integer value of the number of basis functions.

---

## 2. Main Function: `coupling(argv)`

### **Step 1: File Handling**
- Resolves paths for the `InFile`, `paramFile`, and output files.
- Reads the input files and loads the Fock, overlap, and orbital coefficient files.

you should have the xxx.in file where the pun files and 53, 54 files are generated by the gaussian
``` plaintext
Cl2-F2.FILE.53
Cl2-F2.FILE.54
Cl2-F2-MO1.com.pun
Cl2-F2-MO2.com.pun
```

### **Step 2: Reading the InFile**
Reads file names for:
- Fock matrix (`fockFile`).
- Overlap matrix (`overlapFile`).
- Molecular orbital coefficients for fragments A (`molOrbFileA`) and B (`molOrbFileB`).

### **Step 3: Reading the paramFile**
Extracts:
- Orbital ranges:
  - `orbitalsA_begin` to `orbitalsA_end` (for fragment A).
  - `orbitalsB_begin` to `orbitalsB_end` (for fragment B).
- Number of basis functions for A (`nBasisFunctsA`) and B (`nBasisFunctsB`) using `findnBasisFuncts`.

### **Step 4: Reading Matrices**
- **Fock Matrix (`F`)** and **Overlap Matrix (`O`)**:
  - Both are stored in a lower triangular format in their respective files.
  - The code reconstructs the full symmetric matrices $F$ and $O$ from these values.

### **Step 5: Reading Molecular Orbitals**
- Molecular orbital coefficients for fragments A and B are read from their respective files.
- Coefficients are processed to handle scientific notation in the form `D` (common in quantum chemistry outputs).

---

## Step 6: Coupling Calculation

For each pair of molecular orbitals $(\text{orbA}, \text{orbB})$:

1. **Self-overlap**:
   
$$
   S(X, X) = \sum_{i, j} C_X[i] \cdot O[i, j] \cdot C_X[j]
$$

3. **Self-energy**:
   
$$
   E(X, X) = \sum_{i, j} C_X[i] \cdot F[i, j] \cdot C_X[j]
$$
   - Converted from Hartree to eV using:
     
$$
     1 \, \text{Hartree} = 27.2114 \, \text{eV}
$$

4. **Cross-term**:
   - **Overlap**:
     
$$
     S(A, B) = \sum_{i, j} C_A[i] \cdot O[p, q] \cdot C_B[j]
$$
   - **Interaction**:
     
$$
     E(A, B) = \sum_{i, j} C_A[i] \cdot F[p, q] \cdot C_B[j]
$$
   - Index $p, q$ accounts for offsets due to the basis functions of A.

4. **Electronic Coupling**:
   
$$
   t_{AB} = \frac{E(A, B) - \frac{1}{2} \big(E(A, A) + E(B, B)\big) S(A, B)}{1 - S(A, B)^2}
$$

---

### **Step 7: Writing Results**
Outputs:
1. **Detailed Results**: Overlap and energy terms for each orbital pair $(\text{orbA}, \text{orbB})$ in `filename.out`.
2. **Coupling Value**: $t_{AB}$ in `filename.t`.

---

# Key Equations

### **1. Self-energy **($E(X, X)$):
$$
E(X, X) = \sum_{i, j} C_X[i] \cdot F[i, j] \cdot C_X[j]
$$

### **2. Overlap Integral ($S(X, X)$)**:
$$
S(X, X) = \sum_{i, j} C_X[i] \cdot O[i, j] \cdot C_X[j]
$$

### **3. Electronic Coupling ($t_{AB}$)**:
$$
t_{AB} = \frac{E(A, B) - \frac{1}{2} \big(E(A, A) + E(B, B)\big) S(A, B)}{1 - S(A, B)^2}
$$

## Conclusion

Electronic coupling is a vital parameter for understanding and optimizing molecular interactions in applications ranging from photovoltaics to molecular electronics. By leveraging Fock and overlap matrices, quantum chemistry tools enable precise calculation of $V$, guiding the design of advanced materials and devices.

```
![image](https://github.com/user-attachments/assets/081cd0b7-09f1-436a-b6e0-6e64caf79015)

for hole transport

get the energy from the neutral state and get the energy of an electron removed. the difference will be lamda 1
optimized the hole state and get the energy and then remove the hole and get the energy. the difference is lamda 2
the total will be lamda1 and lamda 2

