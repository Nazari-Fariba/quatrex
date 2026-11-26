#!/usr/bin/env python
"""Convert a Wannier90 <seed>_hr.dat into a Quatrex hamiltonian.h5.

One dense (num_wann x num_wann) complex128 matrix per hopping vector, stored
under the key '[R1,R2,R3]' with the attribute format='ndarray', which is what
qttools.utils.hdf5_utils.load_hdf5_dict expects.

By default H(R) is divided by the Wigner-Seitz degeneracy ndegen(R), as
Wannier90 does when it interpolates (H(k) = sum_R exp(2 pi i k.R) H(R) / ndegen(R)).
Use --no-degeneracy to skip the division.

Usage:
    python w90_hr_to_quatrex.py Fe_hr.dat hamiltonian.h5 [--no-degeneracy]
"""
import sys

import numpy as np


def read_hr(path):
    lines = open(path).read().split("\n")
    nw, nr = int(lines[1]), int(lines[2])
    nl = -(-nr // 15)                       # 15 degeneracies per line
    ndeg = np.array(" ".join(lines[3:3 + nl]).split(), dtype=int)
    body = np.array([l.split() for l in lines[3 + nl:3 + nl + nr * nw * nw]], float)
    R = body[:, :3].astype(int).reshape(nr, nw, nw, 3)[:, 0, 0]
    val = (body[:, 5] + 1j * body[:, 6]).reshape(nr, nw, nw)   # [r, n (slow), m (fast)]
    H = np.transpose(val, (0, 2, 1))                            # H[r, row m, col n]
    return R, ndeg, H


def write_quatrex_h5(path, R, H):
    data = {f"[{r[0]},{r[1]},{r[2]}]": m.astype(np.complex128) for r, m in zip(R, H)}
    try:
        from qttools.utils.hdf5_utils import save_hdf5_dict

        save_hdf5_dict(path, data)           # exactly the format Quatrex writes
    except ImportError:
        import h5py

        with h5py.File(path, "w") as f:
            for key, m in data.items():
                f.create_dataset(key, data=m).attrs["format"] = "ndarray"


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    R, ndeg, H = read_hr(src)
    if "--no-degeneracy" not in sys.argv[3:]:
        H = H / ndeg[:, None, None]
    write_quatrex_h5(dst, R, H)
    print(f"wrote {dst}: {len(R)} hopping vectors, {H.shape[1]} x {H.shape[2]} each, "
          f"degeneracy {'ignored' if '--no-degeneracy' in sys.argv[3:] else 'applied'}")