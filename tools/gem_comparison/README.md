# GEM, MATLAB and Python comparison (no-flow CCE)

C1-nC8, z = 0.923 / 0.077, 68 F, jL3 (FWI-only). Three EOS sets: bulk, EF2 54 nm, B1 11.25 nm.

| File | What it is |
| --- | --- |
| `CCEV_C_jL3.dat` | GEM 2025.20 deck: 225 x 3 closed blocks, one row per EOS set, column i at pressure i |
| `pressures_all.csv` | the 225 block pressures (psia): 100 psi from 2000 to 6500, 25 psi from 3300 to 3900, 10 psi within 150 psi and 1 psi within 20 psi of each dew point |
| `python_<set>.csv` | Python PR reference at the same pressures, using only the numbers in the GEM deck |
| `compare_cce.py` | reads `python_`, `matlab_` and `gem_<set>.csv` from one folder and writes `cce_comparison.png` and `comparison_summary.csv` |

Steps:

1. Run `GEM_MATLAB_comparison.m` from the repository root. It writes `results/gem_comparison/matlab_<set>.csv` and `matlab_dewpoints.csv`.
2. Run `CCEV_C_jL3.dat` in GEM and convert the time-zero grid output to `gem_<set>.csv` with the same columns.
3. Copy the `python_<set>.csv` files next to them and run `python -I tools/gem_comparison/compare_cce.py results/gem_comparison`.

Expected from the Python reference (dew point, last two-phase / first single-phase block, maximum liquid saturation):

| Set | Pdew (psia) | Blocks | Max SL |
| --- | --- | --- | --- |
| bulk | 4065.6 | 4065 / 4066 | 0.212 at 3550 psia |
| EF2 54 nm | 4205.1 | 4205 / 4206 | 0.200 at 3575 psia |
| B1 11.25 nm | 4164.7 | 4164 / 4165 | 0.191 at 3500 psia |
