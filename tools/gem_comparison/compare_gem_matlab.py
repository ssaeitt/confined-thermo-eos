"""Compare GEM blocks with the MATLAB confined EOS (CCE of C1-nC8, z = 0.923).

Reads from one folder:
  gem_<set>.csv      from read_gem_cce.py (state after the first time step)
  matlab_<set>.csv   from GEM_MATLAB_comparison.m
  matlab_dewpoints.csv
GEM's block pressures move during the first time step, so MATLAB is interpolated to
them, using two-phase MATLAB rows only for two-phase quantities. GEM prints saturations
and mole fractions to 3 decimals, so differences up to about 0.0005 are print precision.
Writes gem_matlab_summary.csv and gem_matlab_cce.png.
Run:  python -I compare_gem_matlab.py <folder>
"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SETS = [("bulk", "Bulk", "#6B6A66"), ("EF2_54", "EF2 54 nm", "#1F5FA8"), ("B1_11.25", "B1 11.25 nm", "#B5651D")]


def load(path):
    rows = list(csv.DictReader(open(path)))
    d = {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if k != "converged"}
    o = np.argsort(d["P_psia"])
    return {k: v[o] for k, v in d.items()}


def main(folder):
    dew = {r["set"]: float(r["Pdew_psia"]) for r in csv.DictReader(open(os.path.join(folder, "matlab_dewpoints.csv")))}
    out = []
    fig, axs = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [1.5, 1]})
    for tag, name, col in SETS:
        g = load(os.path.join(folder, f"gem_{tag}.csv"))
        m = load(os.path.join(folder, f"matlab_{tag}.csv"))
        Pd, P = dew[tag], g["P_psia"]
        m2 = m["phases"] == 2
        two = (g["phases"] == 2) & (P <= m["P_psia"][m2].max())
        ip2 = lambda c: np.interp(P[two], m["P_psia"][m2], m[c][m2])
        one = g["phases"] == 1
        away = np.abs(P - Pd) > 1.0
        mSL = np.interp(P, m["P_psia"], m["SL"])
        out.append({
            "set": tag, "Pdew_matlab_psia": Pd,
            "gem_last_two_phase_psia": P[g["phases"] == 2].max(),
            "gem_first_single_phase_psia": P[one].min(),
            "phase_mismatch_away_from_dew": int(np.sum((g["phases"] != np.where(P < Pd, 2, 1)) & away)),
            "max_abs_dSL": np.max(np.abs(g["SL"] - mSL)[away]),
            "max_rel_drhoL_pct": 100 * np.max(np.abs(g["rhoL_kgm3"][two] / ip2("rhoL_kgm3") - 1)),
            "max_rel_drhoV_pct": 100 * max(np.max(np.abs(g["rhoV_kgm3"][two] / ip2("rhoV_kgm3") - 1)),
                                           np.max(np.abs(g["rhoV_kgm3"][one] / np.interp(P[one], m["P_psia"], m["rhoV_kgm3"]) - 1))),
            "max_abs_dx_nC8": np.max(np.abs(g["x_nC8"][two] - ip2("x_nC8"))),
            "max_abs_dy_nC8": np.max(np.abs(g["y_nC8"][two] - ip2("y_nC8"))),
            "max_SL_gem": g["SL"].max(), "max_SL_matlab": m["SL"].max()})
        for a in axs:
            a.plot(m["P_psia"], m["SL"], color=col, lw=2.4, label=f"{name}, MATLAB")
            a.plot(P, g["SL"], "o", ms=4.5, mfc="white", mec=col, mew=1.4, label=f"{name}, GEM blocks")
        axs[1].axvline(Pd, color=col, ls=":", lw=1.5)
    axs[0].set_xlim(2000, 4400); axs[0].set_ylim(0, 0.23)
    axs[0].set_xlabel("Pressure (psia)"); axs[0].set_ylabel("Liquid saturation (fraction of volume)")
    axs[1].set_xlim(4040, 4215); axs[1].set_ylim(0, 0.10); axs[1].set_xlabel("Pressure (psia)")
    axs[1].set_title("Near the dew points (dotted: MATLAB dew point)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(os.path.join(folder, "gem_matlab_cce.png"), dpi=150)
    with open(os.path.join(folder, "gem_matlab_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    for r in out:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
