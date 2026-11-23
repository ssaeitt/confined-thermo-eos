"""Compare GEM, MATLAB and the Python PR reference block by block.

Inputs in one folder, per set tag (bulk, EF2_54, B1_11.25):
  python_<tag>.csv   Python reference (this script's folder)
  matlab_<tag>.csv   from GEM_MATLAB_comparison.m
  gem_<tag>.csv      GEM blocks, same columns (P_psia, phases, SL, rhoL_kgm3, rhoV_kgm3, x_nC8, y_nC8)
Any file that is missing is skipped. Writes comparison_summary.csv and cce_comparison.png.
Run:  python -I compare_cce.py <folder>
"""
import csv, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAGS = [("bulk", "Bulk"), ("EF2_54", "EF2 54 nm"), ("B1_11.25", "B1 11.25 nm")]
COLS = ["SL", "rhoL_kgm3", "rhoV_kgm3", "x_nC8", "y_nC8"]


def load(path):
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        r = list(csv.DictReader(f))
    d = {k: np.array([float(x[k]) if x[k] not in ("", "NaN", "nan") else np.nan for x in r])
         for k in r[0] if k in ["P_psia", "phases"] + COLS}
    o = np.argsort(d["P_psia"])
    return {k: v[o] for k, v in d.items()}


def bracket(d):
    two = d["P_psia"][d["phases"] == 2]
    one = d["P_psia"][d["phases"] == 1]
    return (two.max() if two.size else np.nan, one.min() if one.size else np.nan)


def main(folder):
    rows = []
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), sharey=True)
    for ax, (tag, name) in zip(axs, TAGS):
        ref = load(os.path.join(folder, f"python_{tag}.csv"))
        for src, style in [("python", dict(color="#8A93A3", lw=4, label="Python reference")),
                           ("matlab", dict(color="#1F5FA8", lw=1.8, label="MATLAB")),
                           ("gem", dict(color="#B5521F", ls="none", marker="o", ms=3.5, label="GEM blocks"))]:
            d = load(os.path.join(folder, f"{src}_{tag}.csv"))
            if d is None:
                continue
            ax.plot(d["P_psia"], d["SL"], **style)
            lo, hi = bracket(d)
            row = {"set": tag, "source": src, "last_two_phase_psia": lo, "first_single_phase_psia": hi,
                   "MLDO": np.nanmax(d["SL"]), "P_MLDO_psia": d["P_psia"][np.nanargmax(d["SL"])]}
            if src != "python" and ref is not None:
                common = np.intersect1d(ref["P_psia"], d["P_psia"])
                ir = np.searchsorted(ref["P_psia"], common)
                idd = np.searchsorted(d["P_psia"], common)
                both2 = (ref["phases"][ir] == 2) & (d["phases"][idd] == 2)
                row["phase_mismatches"] = int(np.sum(ref["phases"][ir] != d["phases"][idd]))
                for c in COLS:
                    a, b = ref[c][ir], d[c][idd]
                    m = both2 if c in ("rhoL_kgm3", "x_nC8", "y_nC8") else np.isfinite(a) & np.isfinite(b)
                    row["maxabs_d_" + c] = float(np.nanmax(np.abs(a[m] - b[m]))) if m.any() else np.nan
            rows.append(row)
        ax.set_title(name)
        ax.set_xlabel("Pressure (psia)")
        ax.set_xlim(2000, 4500)
    axs[0].set_ylabel("Liquid saturation, CCE (fraction of volume)")
    axs[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "cce_comparison.png"), dpi=150)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("set", "source"), k))
    with open(os.path.join(folder, "comparison_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
