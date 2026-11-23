"""Read the time-zero grid printout of the no-flow GEM deck CCEV_C_jL3 and write
gem_<set>.csv (P_psia, phases, SL, rhoL_kgm3, rhoV_kgm3, x_nC8, y_nC8, VISO_cp, VISG_cp).

Row J of the deck is EOS set J (1 bulk, 2 EF2 54 nm, 3 B1 11.25 nm); column I is one
block. The printout at time `want` (days, as GEM prints it) is read; the default "1" is
the state after the first time step. At time 0 GEM prints its initialisation, which
labels blocks between about 4000 psia and the dew point single-phase; the first step
flashes them at fixed moles and volume and moves their pressure by up to 76 psi, so the
pressures in the output are GEM's final block pressures, not the list in the deck. GEM prints saturation and
compositions to 3 decimals, so differences below 0.0005 are below its print precision.
Run:  python -I read_gem_cce.py CCEV_C_jL3.out [out_folder] [time] [prefix]
"""
import os, re, sys, csv

SETS = {1: "bulk", 2: "EF2_54", 3: "B1_11.25"}
LBFT3 = 16.018463          # lb/ft3 -> kg/m3
TITLES = {"Oil Saturation": "SO", "Gas Saturation": "SG", "Pressure  (psia)": "P",
          "Oil Mass Density (lb/cuft)": "DENO", "Gas Mass Density (lb/cuft)": "DENG",
          "Oil Viscosity, cp": "VISO", "Gas Viscosity, cp": "VISG",
          "C8 Oil Composition, Mole Fraction": "X8", "C8 Gas Composition, Mole Fraction": "Y8"}


def read(path, want="1"):
    lines = open(path, errors="replace").read().replace("\r", "").split("\n")
    tabs, cur, cols, seen_time = {}, None, [], None
    for k, ln in enumerate(lines):
        m = re.match(r"\s*Time = (\S+)\s+\*{10}", ln)
        if m:
            if m.group(1) != want:
                cur = None
                continue
            seen_time = m.group(1)
            title = lines[k + 1].strip()
            cur = TITLES.get(lines[k + 1].strip()) if title else None
            if cur and cur in tabs:
                cur = None                              # keep the first copy
            if cur:
                tabs[cur] = {}
            continue
        if cur is None:
            continue
        h = re.match(r"\s*I =\s*(.*)", ln)
        if h:
            cols = [int(c) for c in h.group(1).split()]
            continue
        r = re.match(r"\s*J=\s*(\d+)\s+(.*)", ln)
        if r and cols:
            vals = r.group(2).split()
            for i, v in zip(cols, vals):
                tabs[cur][(i, int(r.group(1)))] = float(v)
    return tabs


def main(path, outdir, want="1", prefix="gem"):
    t = read(path, want)
    missing = [k for k in TITLES.values() if k not in t]
    if missing:
        sys.exit(f"tables not found: {missing}")
    ni = max(i for i, _ in t["P"])
    for j, tag in SETS.items():
        with open(os.path.join(outdir, f"{prefix}_{tag}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["P_psia", "phases", "SL", "rhoL_kgm3", "rhoV_kgm3", "x_nC8", "y_nC8", "VISO_cp", "VISG_cp"])
            rows = []
            for i in range(1, ni + 1):
                so = t["SO"][(i, j)]
                two = so > 0
                rows.append([t["P"][(i, j)], 2 if two else 1, so,
                             t["DENO"][(i, j)] * LBFT3 if two else t["DENG"][(i, j)] * LBFT3,
                             t["DENG"][(i, j)] * LBFT3,
                             t["X8"][(i, j)] if two else t["Y8"][(i, j)], t["Y8"][(i, j)],
                             t["VISO"][(i, j)], t["VISG"][(i, j)]])
            rows.sort()
            w.writerows(rows)
    print(f"read {ni} columns x {len(SETS)} sets")


if __name__ == "__main__":
    a = sys.argv
    out = a[2] if len(a) > 2 else os.path.dirname(os.path.abspath(a[1]))
    main(a[1], out, a[3] if len(a) > 3 else "1", a[4] if len(a) > 4 else "gem")
