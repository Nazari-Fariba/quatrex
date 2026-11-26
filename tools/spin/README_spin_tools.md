# Spin and Wannier90 test tools for Quatrex

Scripts used to test a spinor basis in the QTBM and to prepare a Wannier90 Hamiltonian.
They were developed on the CP2K carbon-chain example and on bcc Fe. They are test helpers,
not part of the Quatrex package. Paths inside the shell scripts assume the Quatrex checkout
is at `$HOME/quatrex` and the carbon-chain example inputs are unchanged.

| File | Purpose |
|---|---|
| `spinor_tools.py` | Builds spinor test inputs (interleaved basis, index = 2 x orbital + spin): `spinor` (kron(M, I2) plus a Zeeman term J sigma.n, optionally with a different direction in the right half), `shift` (spinless reference with a diagonal shift), `potential` (repeats a per-orbital potential for the two spins), `selftest` (algebra check without HDF5). Reads and writes with `qttools.utils.hdf5_utils`. |
| `valve_cases.sh` | Builds the spin-valve cases on the carbon chain (parallel, antiparallel, 45, 90, 135 degrees, and two spinless references) and a run script. `MODE=devspin` writes configs for the upstream `dev-spin` branch. |
| `compare_valve.py` | Compares the valve results: antiparallel against the sum of two spinless chains (exact reference), azimuthal independence, angle scan. |
| `devspin_cases.sh` | Builds two test cases for upstream `dev-spin` (`spin = true`). |
| `w90_hr_to_quatrex.py` | Converts a Wannier90 `*_hr.dat` into a Quatrex `hamiltonian.h5` (one dense matrix per hopping vector, with the `format` attribute). Divides H(R) by the degeneracy ndegen(R); `--no-degeneracy` skips this. |
| `check_degeneracy.py` | Checks the degeneracy treatment against the DFT bands in `*.eig` inside the frozen window of the disentanglement. Needs only numpy. |

Results obtained with these tools and open points are in the report
"Non-collinear transmission in Quatrex: status and next steps".

Notes

- Run one case at a time; the spinor carbon-chain cases use about 3.5 GB of memory.
- The case builders delete and recreate their case directories: start each builder once.
- Spin-orbit coupling has not been tested with these tools yet.
