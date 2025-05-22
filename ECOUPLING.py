# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 20:41:49 2024

@author: User
"""

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
