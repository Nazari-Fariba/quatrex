#!/usr/bin/env python
"""Spinor test inputs for Quatrex (interleaved basis, index = 2*orbital + spin).

Reads and writes with Quatrex's own load_hdf5_dict / save_hdf5_dict.

Modes
  spinor    M -> kron(M, I2). With --zeeman J (--direction n) adds J (sigma . n) on every
            orbital at R = (0,0,0). With --split S the orbitals with index < S get
            direction n, the others --direction-right (same J). S counts spatial orbitals
            of the input matrix (before the spin expansion). Hamiltonian only.
  shift     M -> M + shift * identity at R = 0 (spinless reference). With --split S the
            orbitals with index >= S get --shift-right instead.
  potential repeats every value of a per-orbital potential.npy twice (one per spin).
  selftest  checks the algebra without HDF5.

Examples
  python spinor_tools.py spinor ham.h5 out.h5 --zeeman 0.05 --direction 0 0 1 --split 1300 --direction-right 0 0 -1
  python spinor_tools.py shift  ham.h5 out.h5 --shift 0.05 --shift-right -0.05 --split 1300
  python spinor_tools.py potential potential.npy out_potential.npy
"""
import argparse
import re

import numpy as np
from scipy import sparse


def sigma_n(direction):
    nx, ny, nz = np.asarray(direction, float) / np.linalg.norm(direction)
    return np.array([[nz, nx - 1j * ny], [nx + 1j * ny, -nz]])


def zeeman_term(n, j, direction, split=None, direction_right=None):
    left = np.ones(n, bool) if split is None else np.arange(n) < split
    dr = direction if direction_right is None else direction_right
    return (sparse.kron(sparse.diags(left.astype(float)), j * sigma_n(direction), format="coo")
            + sparse.kron(sparse.diags((~left).astype(float)), j * sigma_n(dr), format="coo"))


def spinor_expand(m, zeeman=0.0, direction=(0, 0, 1), split=None, direction_right=None):
    dense = not sparse.issparse(m)
    a = sparse.coo_matrix(m)
    out = sparse.kron(a, sparse.identity(2), format="coo")
    if zeeman != 0.0:
        out = out + zeeman_term(a.shape[0], zeeman, direction, split, direction_right)
    return out.toarray() if dense else sparse.coo_matrix(out)


def shift_diag(m, shift, shift_right=None, split=None):
    dense = not sparse.issparse(m)
    n = m.shape[0]
    d = np.full(n, shift, dtype=complex)
    if split is not None and shift_right is not None:
        d[split:] = shift_right
    a = sparse.coo_matrix(m) + sparse.diags(d)
    return a.toarray() if dense else sparse.coo_matrix(a)


def selftest():
    rng = np.random.default_rng(1)
    n, j, s = 8, 0.37, 5
    h = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h = h + h.conj().T
    # 1. uniform Zeeman: eigenvalues e +/- J for any direction
    e0 = np.linalg.eigvalsh(h)
    for d in ((0, 0, 1), (1, 2, 2)):
        e = np.sort(np.linalg.eigvalsh(spinor_expand(h, j, d)))
        print("uniform", d, np.allclose(e, np.sort(np.concatenate([e0 + j, e0 - j]))))
    # 2. decoupled halves with different directions: e_left +/- J and e_right +/- J
    hd = h.copy(); hd[:s, s:] = 0; hd[s:, :s] = 0
    el, er = np.linalg.eigvalsh(hd[:s, :s]), np.linalg.eigvalsh(hd[s:, s:])
    e = np.sort(np.linalg.eigvalsh(spinor_expand(hd, j, (0, 0, 1), s, (1, 2, 2))))
    ref = np.sort(np.concatenate([el + j, el - j, er + j, er - j]))
    print("decoupled halves, n_right arbitrary:", np.allclose(e, ref))
    # 3. antiparallel: spin blocks equal spinless chains with opposite shifts on the two halves
    ap = spinor_expand(h, j, (0, 0, 1), s, (0, 0, -1))
    up = shift_diag(h, j, -j, s)
    dn = shift_diag(h, -j, j, s)
    print("antiparallel up block == H + (+J, -J):", np.allclose(ap[0::2, 0::2], up),
          "| down block == H + (-J, +J):", np.allclose(ap[1::2, 1::2], dn),
          "| spin-flip:", np.abs(ap[0::2, 1::2]).max())
    print("hermitian:", np.abs(ap - ap.conj().T).max() < 1e-12)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["spinor", "shift", "potential", "selftest"])
    p.add_argument("src", nargs="?")
    p.add_argument("dst", nargs="?")
    p.add_argument("--zeeman", type=float, default=0.0)
    p.add_argument("--direction", type=float, nargs=3, default=[0, 0, 1])
    p.add_argument("--direction-right", type=float, nargs=3, default=None)
    p.add_argument("--split", type=int, default=None)
    p.add_argument("--shift", type=float, default=0.0)
    p.add_argument("--shift-right", type=float, default=None)
    a = p.parse_args()
    if a.mode == "selftest":
        selftest()
        raise SystemExit
    if a.mode == "potential":
        v = np.load(a.src)
        np.save(a.dst, np.repeat(v, 2))
        print(f"wrote {a.dst}: {v.shape} -> {np.load(a.dst).shape}")
        raise SystemExit
    from qttools.utils.hdf5_utils import load_hdf5_dict, save_hdf5_dict

    data = load_hdf5_dict(a.src)
    out = {}
    for key, m in data.items():
        at_origin = tuple(int(x) for x in re.findall(r"-?\d+", key)) == (0, 0, 0)
        if a.mode == "spinor":
            out[key] = spinor_expand(m, a.zeeman if at_origin else 0.0, a.direction,
                                     a.split, a.direction_right)
        else:
            out[key] = shift_diag(m, a.shift, a.shift_right, a.split) if at_origin else m
    save_hdf5_dict(a.dst, out)
    print(f"wrote {a.dst}: {len(out)} keys, shape {next(iter(out.values())).shape}")
