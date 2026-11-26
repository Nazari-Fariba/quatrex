#!/usr/bin/env python
"""Compare the spin-valve cases.  Usage: python compare_valve.py [WORKDIR]"""
import glob
import sys

import numpy as np

d = (sys.argv[1] if len(sys.argv) > 1 else "/tmp/valve") + "/"


def T(c):
    f = sorted(glob.glob(f"{d}{c}/**/transmission_lr.npy", recursive=True))
    assert f, f"no transmission_lr.npy for {c}"
    t = np.load(f[0])
    return t.sum(axis=0) if t.ndim > 1 else t


up, dn, P, AP, X, Y, D45, D135 = (T(c) for c in ("upAP", "dnAP", "P", "AP", "X", "Y", "D45", "D135"))
print("T scale:", np.abs(P).max())
print("1. AP vs T(up chain) + T(down chain)   :", np.abs(AP - (up + dn)).max(), "(exact reference)")
print("   P vs AP differ (valve effect)       :", np.abs(P - AP).max(), "(must be clearly > 0)")
print("2. X vs Y, right direction x or y       :", np.abs(X - Y).max(), "(azimuthal independence)")
i = int(np.argmax(np.abs(P - AP)))
print(f"energy index {i} with the largest valve effect: T(0)={P[i]:.5g}  T(45)={D45[i]:.5g}  "
      f"T(90)={X[i]:.5g}  T(135)={D135[i]:.5g}  T(180)={AP[i]:.5g}")
print("   (T(0) and T(180) are the exact limits; the three angles in between should lie between them smoothly)")
