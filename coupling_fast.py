# -*- coding: utf-8 -*-
"""
ECOUPLING.py  —  Import-safe, vectorised & **parallel-JIT-compiled** module for
molecular-orbital coupling calculations.

Main fixes in this revision (v1.1)
----------------------------------
1. **Robust MO-coefficient reader** – correct ceiling division for the number of
   coefficient lines and *hard-truncate* each coefficient list to `n_basis`
   elements so mismatches cannot propagate to the BLAS dot products.
2. **Defensive slicing** – in the driver we now slice the coefficient vectors
   explicitly (`cA = coeffs[:nA]`) so even malformed input cannot upset the
   algebra.
3. **Safer dot3 implementation** – rewrite as `np.dot(v1, M @ v2)` which never
   attempts an invalid vector × matrix multiplication order.
4. Still `@njit(parallel=True, fastmath=True)` compatible; keeps working both
   as **script** and **importable module**.
"""

from __future__ import annotations

import sys
from typing import List, Dict, Tuple

import numpy as np

# ───────────────────────── Optional Numba wrapper ───────────────────────────────
try:
    from numba import njit  # type: ignore

    def _jit(func):  # noqa: D401
        return njit(parallel=True, fastmath=True)(func)  # compile if available

except ModuleNotFoundError:  # no-Numba fallback

    def _jit(func):  # noqa: D401
        return func

# ─────────────────────────────── Constants ──────────────────────────────────────
HARTREE_TO_EV = 27.2114

# ───────────────────────────── I/O Utilities ────────────────────────────────────

def _ceildiv(a: int, b: int) -> int:
    """Small helper: ceiling division."""
    return (a + b - 1) // b


def read_lower_triangular_matrix(filename: str, *, skip_lines: int = 3) -> np.ndarray:
    elems: List[float] = []
    with open(filename, "r") as fh:
        for idx, line in enumerate(fh):
            if idx < skip_lines:
                continue
            elems.extend(map(float, line.split()))
    return np.asarray(elems, dtype=np.float64)


def form_symmetric_matrix(elements: np.ndarray, dim: int) -> np.ndarray:
    mat = np.zeros((dim, dim), dtype=np.float64)
    k = 0
    for i in range(dim):
        for j in range(i + 1):
            mat[i, j] = elements[k]
            k += 1
    return mat + np.triu(mat.T, 1)


def read_orbitals(
    filename: str,
    n_basis: int,
    *,
    coeffs_per_line: int = 5,
    coeff_length: int = 15,
) -> Dict[int, List[float]]:
    """Parse Gaussian-style MO coefficient block and **guarantee** len(coeffs)==n_basis."""

    orbitals: Dict[int, List[float]] = {}
    orb_count = -1
    orb_flag = 0
    n_orb_lines = _ceildiv(n_basis, coeffs_per_line)  # corrected formula

    with open(filename, "r") as fh:
        next(fh)  # skip header line
        for line in fh:
            cols = line.split()
            if orb_flag:
                for i in range(coeffs_per_line):
                    seg = line[i * coeff_length : (i + 1) * coeff_length].strip()
                    if seg and "Alpha" not in seg:
                        if "D" in seg:  # Fortran D-notation
                            base, power = seg.split("D")
                            seg = f"{float(base)}e{int(power)}"
                        orbitals[orb_count].append(float(seg))
                orb_flag += 1
                if orb_flag == n_orb_lines + 1:
                    # truncate in case of over-read & reset flag
                    orbitals[orb_count] = orbitals[orb_count][:n_basis]
                    orb_flag = 0

            if len(cols) > 1 and cols[1] == "Alpha":
                orb_count = int(cols[0])
                orbitals[orb_count] = []
                orb_flag = 1

    # Final safety pass – ensure every list exactly n_basis long
    for k, v in orbitals.items():
        if len(v) < n_basis:
            raise ValueError(
                f"Orbital {k}: expected {n_basis} coeffs, got {len(v)} – file truncated?"
            )
        orbitals[k] = v[:n_basis]
    return orbitals


def parse_param_file(param_file: str) -> Tuple[int, int, int, int, int, int]:
    with open(param_file, "r") as fh:
        lines = [l.strip() for l in fh if l.strip()]
    if len(lines) < 4:
        raise ValueError("Parameter file incomplete.")
    orbA_beg, orbA_end = map(int, lines[0].split("=")[1].split())
    orbB_beg, orbB_end = map(int, lines[1].split("=")[1].split())
    nA = int(lines[2].split("=")[1])
    nB = int(lines[3].split("=")[1])
    return nA, nB, orbA_beg, orbA_end, orbB_beg, orbB_end


def parse_in_file(in_file: str) -> Tuple[str, str, str, str]:
    with open(in_file, "r") as fh:
        lines = [l.strip() for l in fh if l.strip()]
    if len(lines) < 4:
        raise ValueError("Input file must have four non-blank lines.")
    return tuple(lines)  # type: ignore[return-value]

# ─────────────────────────── Numerical core ─────────────────────────────────────

def _slice_blocks(O: np.ndarray, F: np.ndarray, nA: int):
    O_AA, O_BB, O_AB = O[:nA, :nA], O[nA:, nA:], O[:nA, nA:]
    F_AA, F_BB, F_AB = F[:nA, :nA], F[nA:, nA:], F[:nA, nA:]
    return O_AA, O_BB, O_AB, F_AA, F_BB, F_AB


@_jit
def _dot3(v1: np.ndarray, M: np.ndarray, v2: np.ndarray) -> float:  # pragma: no cover
    """Compute v1ᵀ M v2 safely (matrix right-multiplies **first**)."""
    return float(np.dot(v1, M @ v2))


@_jit
def compute_coupling_vec(
    cA: np.ndarray,
    cB: np.ndarray,
    O_AA: np.ndarray,
    O_BB: np.ndarray,
    O_AB: np.ndarray,
    F_AA: np.ndarray,
    F_BB: np.ndarray,
    F_AB: np.ndarray,
) -> Tuple[float, float, float, float, float, float, float]:  # pragma: no cover
    Saa = _dot3(cA, O_AA, cA)
    Eaa = _dot3(cA, F_AA, cA)
    Sbb = _dot3(cB, O_BB, cB)
    Ebb = _dot3(cB, F_BB, cB)
    Sab = _dot3(cA, O_AB, cB)
    Eab = _dot3(cA, F_AB, cB)
    t_ev = (Eab * HARTREE_TO_EV - 0.5 * (Eaa + Ebb) * HARTREE_TO_EV * Sab) / (
        1.0 - Sab * Sab
    )
    return Saa, Eaa, Sbb, Ebb, Sab, Eab, t_ev

# ─────────────────────────── Driver routine ─────────────────────────────────────

def coupling(argv: List[str]) -> None:
    if len(argv) != 4:
        raise SystemExit("Usage: ... inFile paramFile outPrefix")

    _, in_file, param_file, out_prefix = argv
    overlap_file, fock_file, molOrbFileA, molOrbFileB = parse_in_file(in_file)
    nA, nB, orbA_beg, orbA_end, orbB_beg, orbB_end = parse_param_file(param_file)

    nAB = nA + nB
    print(f"[ECOUPLING] Basis sizes: nA={nA} nB={nB}  →  full dim {nAB}")

    F = form_symmetric_matrix(read_lower_triangular_matrix(fock_file), nAB)
    O = form_symmetric_matrix(read_lower_triangular_matrix(overlap_file), nAB)

    orbCoeffsA = read_orbitals(molOrbFileA, nA)
    orbCoeffsB = read_orbitals(molOrbFileB, nB)

    blocks = _slice_blocks(O, F, nA)
    out_full = out_prefix + ".out"
    out_short = out_prefix + ".t"

    with open(out_full, "w") as fh_full, open(out_short, "w") as fh_short:
        for orbA in range(orbA_beg, orbA_end + 1):
            cA = np.asarray(orbCoeffsA[orbA][:nA], dtype=np.float64)
            for orbB in range(orbB_beg, orbB_end + 1):
                cB = np.asarray(orbCoeffsB[orbB][:nB], dtype=np.float64)
                Saa, Eaa, Sbb, Ebb, Sab, Eab, t_ev = compute_coupling_vec(
                    cA, cB, *blocks
                )
                # detailed log
                print("  ****************************************", file=fh_full)
                print(
                    f"Interaction between MO {orbA} on A and MO {orbB} on B", file=fh_full
                )
                print(file=fh_full)
                print(
                    f"S(A,A)={Saa:.8f} E(A,A)={Eaa:.8f} au ({Eaa * HARTREE_TO_EV:.8f} eV)",
                    file=fh_full,
                )
                print(
                    f"S(B,B)={Sbb:.8f} E(B,B)={Ebb:.8f} au ({Ebb * HARTREE_TO_EV:.8f} eV)",
                    file=fh_full,
                )
                print(file=fh_full)
                print(
                    f"S(A,B)={Sab:.8f} E(A,B)={Eab:.8f} au ({Eab * HARTREE_TO_EV:.8f} eV)",
                    file=fh_full,
                )
                print(file=fh_full)
                print(f"t({orbA},{orbB}) = {t_ev:.8f} eV", file=fh_full)
                print(file=fh_full)
                # short log
                print(f"t({orbA},{orbB}) = {t_ev:.8f} eV", file=fh_short)

    print("[ECOUPLING] Done →", out_full, out_short)

# ───────────────────────────── Script entry-point ────────────────────────────────

if __name__ == "__main__":
    coupling(sys.argv)
