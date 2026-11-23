# GEM against MATLAB (no-flow CCE)

C1-nC8, z = 0.923 / 0.077, 68 F, jL3 (FWI-only). Three EOS sets: bulk, EF2 54 nm, B1 11.25 nm.

| File | What it is |
| --- | --- |
| `CCEV_C_jL3.dat` | GEM 2025.20 deck: 225 x 3 closed blocks, one row per EOS set, column i at pressure i |
| `pressures_all.csv` | the 225 block pressures (psia): 100 psi from 2000 to 6500, 25 psi from 3300 to 3900, 10 psi within 150 psi and 1 psi within 20 psi of each dew point |
| `read_gem_cce.py` | reads the GEM `.out` grid printout after the first time step and writes `gem_<set>.csv` |
| `compare_gem_matlab.py` | compares `gem_<set>.csv` with `matlab_<set>.csv`; writes `gem_matlab_summary.csv` and `gem_matlab_cce.png` |
| `python_<set>.csv` | Python PR reference at the deck pressures, using only the numbers in the GEM deck (kept for reference) |

Steps:

1. Run `GEM_MATLAB_comparison.m` from the repository root. It writes `results/gem_comparison/matlab_<set>.csv` and `matlab_dewpoints.csv`.
2. Run `CCEV_C_jL3.dat` in GEM and copy `CCEV_C_jL3.out` to `results/gem_comparison/`.
3. `python -I tools/gem_comparison/read_gem_cce.py results/gem_comparison/CCEV_C_jL3.out`
4. `python -I tools/gem_comparison/compare_gem_matlab.py results/gem_comparison`

Use the printout after the first time step, not time zero. At time zero GEM prints its
initialisation, which labels blocks between about 4000 psia and each dew point
single-phase. The first step flashes them at fixed moles and volume and raises their
pressure by up to 76 psi. A closed block keeps the feed composition, so it stays a point
on the same CCE at GEM's final pressure; the comparison interpolates MATLAB to it.

Result (23 Nov 2026), GEM 2025.20 against MATLAB:

| Set | MATLAB Pdew (psia) | GEM last two-phase / first single-phase (psia) |
| --- | --- | --- |
| bulk | 4065.64 | 4065.3 / 4066.0 |
| EF2 54 nm | 4205.17 | 4204.5 / 4205.1 |
| B1 11.25 nm | 4164.70 | 4164.3 / 4165.0 |

Largest differences away from the dew points: liquid saturation 0.0006, phase mass
densities 0.012 %, nC8 mole fractions 0.0005 (GEM prints 3 decimals). The EF2 bracket
ends 0.07 psi below the MATLAB dew point because a liquid saturation below 0.0005 prints
as zero.
