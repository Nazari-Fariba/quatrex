#!/usr/bin/env python
"""Does H(R) need the division by ndegen(R)? Compare with the DFT bands in <seed>.eig.

Inside the frozen window [dis_froz_min, dis_froz_max] the Wannier-interpolated
band energies at the k-points of the mp_grid must equal the DFT energies.

Usage:  python check_degeneracy.py Fe_hr.dat Fe.eig Fe.win
Needs only numpy.
"""
import re
import sys

import numpy as np


def read_hr(path):
    lines = open(path).read().split("\n")
    nw, nr = int(lines[1]), int(lines[2])
    nl = -(-nr // 15)
    ndeg = np.array(" ".join(lines[3:3 + nl]).split(), dtype=int)
    body = np.array([l.split() for l in lines[3 + nl:3 + nl + nr * nw * nw]], float)
    R = body[:, :3].astype(int).reshape(nr, nw, nw, 3)[:, 0, 0]
    val = (body[:, 5] + 1j * body[:, 6]).reshape(nr, nw, nw)
    return R, ndeg, np.transpose(val, (0, 2, 1))


def read_win(path):
    text = open(path).read()
    def num(key):
        m = re.search(rf"^\s*{key}\s*=\s*([-+0-9.deDE]+)", text, re.M | re.I)
        return float(m.group(1).lower().replace("d", "e"))
    block = re.search(r"begin\s+kpoints(.*?)end\s+kpoints", text, re.S | re.I).group(1)
    ks = np.array([[float(x.lower().replace("d", "e")) for x in l.split()[:3]]
                   for l in block.strip().split("\n") if l.strip()])
    return num("dis_froz_min"), num("dis_froz_max"), ks


def read_eig(path, nk):
    d = np.loadtxt(path)
    nb = int(d[:, 0].max())
    return d[:, 2].reshape(nk, nb)       # lines are ordered: band fastest, k slowest


if __name__ == "__main__":
    hr, eig, win = sys.argv[1:4]
    R, ndeg, H = read_hr(hr)
    fmin, fmax, ks = read_win(win)
    E = read_eig(eig, len(ks))
    print(f"{len(R)} R vectors, {len(ks)} k-points, {E.shape[1]} DFT bands, "
          f"frozen window [{fmin}, {fmax}] eV")
    variants = {"divided by ndegen": H / ndeg[:, None, None], "not divided": H}
    for name, Hs in variants.items():
        worst, count = 0.0, 0
        for k, e in zip(ks, E):
            frozen = np.sort(e[(e >= fmin) & (e <= fmax)])
            Hk = np.einsum("r,rab->ab", np.exp(2j * np.pi * (R @ k)), Hs)
            w = np.linalg.eigvalsh((Hk + Hk.conj().T) / 2)
            worst = max(worst, np.abs(w[:len(frozen)] - frozen).max())
            count += len(frozen)
        print(f"{name:18s}: max |Wannier - DFT| over {count} frozen band energies = {worst:.4f} eV")