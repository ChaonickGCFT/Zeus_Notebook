"""
One-off: regenerate v0.2 assets (heatmap, chi2 table, diffractive CSV, MCMC trace).
Run from repo root:  python notebooks/_v02_generate_outputs.py
"""
from __future__ import annotations

import csv
import gzip
import io
import re
import sys
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "v0.2_assets"
OUT.mkdir(parents=True, exist_ok=True)

# Only table Q² above this reuse the largest table Q² ≤ cut (here: 20k reuses 12k).
FREEZE_Q2_CUT = 12_000.0

HEPDATA_2010 = "https://www.hepdata.net/download/submission/ins836107/1/csv"
DIFR = "https://www.hepdata.net/download/submission/ins447269/1/csv"


def parse_pct(tok: str) -> float:
    tok = tok.strip().strip('"').rstrip("%")
    return float(tok) / 100.0


def load_ins836107() -> pd.DataFrame:
    import tarfile

    raw = urlopen(HEPDATA_2010).read()
    sig_rows, f2_rows = [], []
    with tarfile.open(fileobj=gzip.GzipFile(fileobj=io.BytesIO(raw)), mode="r|") as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith(".csv"):
                continue
            text = tar.extractfile(member).read().decode("utf-8", errors="replace")
            q2 = None
            mode = None
            for line in text.splitlines():
                line = line.strip()
                if not line:
                    continue
                m = re.match(r"#:\s*Q\*\*2\s*\[GeV\^2\]\s*,\s*([0-9.eE+-]+)", line)
                if m:
                    q2 = float(m.group(1))
                    mode = None
                    continue
                if line.startswith("#"):
                    continue
                ul = line.upper()
                if ul.startswith("X,") and "SIG" in ul and "REDUCED" in ul.replace(" ", ""):
                    mode = "sig"
                    continue
                if ul.startswith("X,") and "F2" in ul and "SIG" not in ul:
                    mode = "f2"
                    continue
                if mode == "sig" and q2 is not None:
                    parts = next(csv.reader([line]))
                    parts = [p.strip() for p in parts]
                    if len(parts) < 6:
                        continue
                    try:
                        x, y, ecm = float(parts[0]), float(parts[1]), float(parts[2])
                        sig = float(parts[3])
                    except ValueError:
                        continue
                    rel_parts = []
                    i = 4
                    while i + 1 < len(parts):
                        try:
                            dp = abs(parse_pct(parts[i]))
                            dm = abs(parse_pct(parts[i + 1]))
                            rel_parts.append(0.5 * (dp + dm))
                        except ValueError:
                            break
                        i += 2
                    stat_rel = rel_parts[0] if rel_parts else 0.0
                    comb_rel = float(np.sqrt(np.sum(np.square(rel_parts)))) if rel_parts else stat_rel
                    sig_rows.append(
                        {
                            "Q2": q2,
                            "x": x,
                            "y": y,
                            "Ecm": ecm,
                            "sigma_r": sig,
                            "sigma_r_stat_rel": stat_rel,
                            "sigma_r_comb_rel": comb_rel,
                        }
                    )
                elif mode == "f2" and q2 is not None:
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) < 4:
                        continue
                    try:
                        x, y, ecm = float(parts[0]), float(parts[1]), float(parts[2])
                    except ValueError:
                        continue
                    f2tok = parts[3].strip()
                    if f2tok in ("-", "", "nan"):
                        continue
                    try:
                        f2 = float(f2tok)
                    except ValueError:
                        continue
                    norm_rel = 0.005
                    if len(parts) >= 6 and "%" in parts[4]:
                        try:
                            norm_rel = max(norm_rel, abs(parse_pct(parts[4])))
                        except ValueError:
                            pass
                    f2_rows.append(
                        {
                            "Q2": q2,
                            "x": x,
                            "y": y,
                            "Ecm": ecm,
                            "F2": f2,
                            "F2_norm_rel": norm_rel,
                        }
                    )
    sig_df = pd.DataFrame(sig_rows)
    f2_df = pd.DataFrame(f2_rows)
    keys = ["Q2", "x", "y", "Ecm"]
    sig_df = sig_df.drop_duplicates(subset=keys, keep="first")
    f2_df = f2_df.drop_duplicates(subset=keys, keep="first")
    merge_cols = keys + ["sigma_r", "sigma_r_stat_rel", "sigma_r_comb_rel"]
    merged = f2_df.merge(sig_df[merge_cols], on=keys, how="left")
    floor = merged["F2_norm_rel"].clip(lower=0.005)
    merged["F2_rel_err_stat"] = merged["sigma_r_stat_rel"].fillna(floor)
    merged["F2_rel_err"] = merged["sigma_r_comb_rel"].fillna(floor)
    df = merged[np.isfinite(merged["F2"]) & (merged["F2"] > 0)].copy()
    df["logx"] = np.log(df["x"])
    df["logQ2"] = np.log(df["Q2"])
    df["logF2"] = np.log(df["F2"])
    return df


def fit_global_quad(d: pd.DataFrame):
    lx, lq = d["logx"].values, d["logQ2"].values
    X = np.column_stack([np.ones(len(d)), lx, lq, lx**2, lq**2, lx * lq])
    coef, _, _, _ = np.linalg.lstsq(X, d["logF2"], rcond=None)
    return coef


def log_f2_global_quad(x, Q2, coef):
    c0, c1, c2, c3, c4, c5 = coef
    lx = np.log(x)
    lq = np.log(Q2)
    return c0 + c1 * lx + c2 * lq + c3 * lx**2 + c4 * lq**2 + c5 * lx * lq


def fit_per_q2_surfaces(d: pd.DataFrame):
    out = {}
    for q2, g in d.groupby("Q2"):
        n = len(g)
        if n < 2:
            out[q2] = ("none", None)
        elif n < 4:
            X = np.column_stack([np.ones(n), g["logx"].values])
            c, _, _, _ = np.linalg.lstsq(X, g["logF2"].values, rcond=None)
            out[q2] = ("line", c)
        else:
            lx = g["logx"].values
            X = np.column_stack([np.ones(n), lx, lx**2])
            c, _, _, _ = np.linalg.lstsq(X, g["logF2"].values, rcond=None)
            out[q2] = ("quad", c)
    return out


def freeze_per_q2_above(per_q2: dict, d: pd.DataFrame, q_cut: float = 500.0) -> dict:
    uq = np.sort(d["Q2"].unique())
    le = uq[uq <= q_cut]
    q_ref = float(le.max()) if len(le) else float(uq.min())
    ref = None
    for k, v in per_q2.items():
        if np.isclose(float(k), q_ref, rtol=0.0, atol=1e-6):
            ref = v
            break
    if ref is None or ref[1] is None:
        return per_q2
    rk, rc = ref[0], np.array(ref[1], dtype=float, copy=True)
    out = {}
    for k, v in per_q2.items():
        fk = float(k)
        if fk > q_cut:
            out[k] = (rk, rc.copy())
        else:
            if v[1] is None:
                out[k] = v
            else:
                out[k] = (v[0], np.array(v[1], dtype=float, copy=True))
    return out


def log_f2_per_q2(x, Q2, per_q2, coef_quad):
    x = np.asarray(x, dtype=float)
    Q2 = np.asarray(Q2, dtype=float)
    out = np.empty_like(x, dtype=float)
    xr, qr, orv = x.ravel(), Q2.ravel(), out.ravel()
    for i in range(xr.size):
        xi, qi = xr[i], qr[i]
        ent = per_q2.get(qi)
        if ent is None:
            ent = per_q2.get(float(qi))
        kind, c = ent if ent is not None else ("none", None)
        if kind == "quad" and c is not None:
            lx = np.log(xi)
            orv[i] = c[0] + c[1] * lx + c[2] * lx**2
        elif kind == "line" and c is not None:
            orv[i] = c[0] + c[1] * np.log(xi)
        else:
            orv[i] = log_f2_global_quad(xi, qi, coef_quad)
    return out


def load_diffractive(url: str) -> pd.DataFrame:
    import tarfile

    raw = urlopen(url).read()
    rows = []
    with tarfile.open(fileobj=gzip.GzipFile(fileobj=io.BytesIO(raw)), mode="r|") as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith(".csv"):
                continue
            text = tar.extractfile(member).read().decode("utf-8", errors="replace")
            beta = np.nan
            q2 = np.nan
            mode = False
            for line in text.splitlines():
                line = line.strip()
                if not line:
                    continue
                mb = re.search(r"BETA.*,\s*([0-9.eE+-]+)", line)
                if line.startswith("#:") and mb:
                    beta = float(mb.group(1))
                m = re.match(r"#:\s*Q\*\*2\s*\[GeV\^2\]\s*,\s*([0-9.eE+-]+)", line)
                if m:
                    q2 = float(m.group(1))
                if line.startswith("#:"):
                    continue
                ul = line.upper()
                if ul.startswith("X,") and "F2" in ul:
                    mode = True
                    continue
                if mode:
                    parts = next(csv.reader([line]))
                    parts = [p.strip() for p in parts]
                    if len(parts) < 4:
                        continue
                    try:
                        xv = float(parts[0])
                    except ValueError:
                        continue
                    if parts[1].strip() in ("-", "", "nan"):
                        continue
                    try:
                        yv = float(parts[1])
                    except ValueError:
                        continue
                    rels = []
                    j = 2
                    while j + 1 < len(parts):
                        try:
                            rels.append(
                                0.5 * (abs(parse_pct(parts[j])) + abs(parse_pct(parts[j + 1])))
                            )
                        except ValueError:
                            break
                        j += 2
                    comb = float(np.sqrt(np.sum(np.square(rels)))) if rels else 0.03
                    rows.append(
                        {
                            "Q2": q2,
                            "beta": beta,
                            "x": xv,
                            "xF2D3": yv,
                            "rel_err": comb,
                            "table": member.name.split("/")[-1],
                        }
                    )
    return pd.DataFrame(rows)


def main():
    df = load_ins836107()
    coef_quad = fit_global_quad(df)
    per_q2 = freeze_per_q2_above(fit_per_q2_surfaces(df), df, FREEZE_Q2_CUT)
    xv = df["x"].values
    Q2v = df["Q2"].values
    f2v = df["F2"].values
    sig_rel = df["F2_rel_err"].values
    log_perq = log_f2_per_q2(xv, Q2v, per_q2, coef_quad)
    pred = np.exp(np.clip(log_perq, -80.0, 80.0))
    res_rel = (f2v - pred) / f2v
    sig_f2 = np.maximum(sig_rel * f2v, 1e-15)
    chi2 = ((f2v - pred) / sig_f2) ** 2

    print("Global |rel res| percentiles (freeze-Q2):", np.percentile(np.abs(res_rel), [50, 90, 95, 99]))
    print("Fraction |rel res| < 7%:", float(np.mean(np.abs(res_rel) < 0.07)))
    uq = np.sort(df["Q2"].unique())
    q_ref = float(uq[uq <= FREEZE_Q2_CUT].max())
    print(f"Freeze: table Q2 > {FREEZE_Q2_CUT:g} reuse coeffs from table Q2 = {q_ref:g} GeV^2")

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    hb = ax.hexbin(
        np.log10(Q2v),
        np.log10(xv),
        C=res_rel,
        reduce_C_function=np.mean,
        gridsize=28,
        cmap="RdBu_r",
        vmin=-0.07,
        vmax=0.07,
        mincnt=1,
    )
    ax.set_xlabel(r"$\log_{10}(Q^2/\mathrm{GeV}^2)$")
    ax.set_ylabel(r"$\log_{10}(x)$")
    ax.set_title(
        rf"Mean rel. residual (per-$Q^2$; $Q^2>{FREEZE_Q2_CUT:g}$ reuse ${q_ref:g}$ GeV$^2$ template)"
    )
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label("mean rel. residual")
    fig.tight_layout()
    fig.savefig(OUT / "residual_heatmap_v02.png", dpi=160)
    plt.close()

    dfv = df.copy()
    dfv["chi2"] = chi2
    dfv["abs_rel_res"] = np.abs(res_rel)
    per_q2_tbl = (
        dfv.groupby("Q2", sort=True)
        .agg(
            n=("chi2", "size"),
            sum_chi2=("chi2", "sum"),
            mean_abs_rel_res=("abs_rel_res", "mean"),
            frac_lt_7pct=("abs_rel_res", lambda s: float(np.mean(s < 0.07))),
        )
        .reset_index()
        .sort_values("sum_chi2", ascending=False)
    )
    per_q2_tbl.to_csv(OUT / "chi2_per_q2.csv", index=False)
    print("Top 10 table Q2 by sum(chi2):")
    print(per_q2_tbl.head(10).to_string(index=False))

    try:
        bins = np.geomspace(dfv["Q2"].min(), dfv["Q2"].max(), 18)
        dfv["Q2_bin"] = pd.cut(dfv["Q2"], bins=bins)
        top = (
            dfv.groupby("Q2_bin", observed=False)["chi2"]
            .sum()
            .sort_values(ascending=False)
            .head(12)
        )
    except Exception:
        dfv["Q2_bin"] = pd.qcut(dfv["Q2"], q=12, duplicates="drop")
        top = (
            dfv.groupby("Q2_bin", observed=False)["chi2"]
            .sum()
            .sort_values(ascending=False)
            .head(12)
        )
    top.to_csv(OUT / "chi2_top_bins.csv", header=["sum_chi2"])

    ddf = load_diffractive(DIFR)
    ddf = ddf[np.isfinite(ddf["xF2D3"]) & (ddf["xF2D3"] > 0)].copy()

    def fit_grp(g: pd.DataFrame) -> pd.Series:
        if len(g) < 2:
            return pd.Series(np.nan, index=g.index)
        lx = np.log(g["x"].values)
        ly = np.log(np.maximum(g["xF2D3"].values, 1e-30))
        X = np.column_stack([np.ones(len(g)), lx])
        c, _, _, _ = np.linalg.lstsq(X, ly, rcond=None)
        predl = c[0] + c[1] * lx
        return pd.Series(np.exp(predl), index=g.index)

    ddf["model_gcft_template"] = np.nan
    for (_, _), g in ddf.groupby(["Q2", "beta"]):
        ddf.loc[g.index, "model_gcft_template"] = fit_grp(g)
    ddf["delta_rel_model_minus_data"] = (ddf["model_gcft_template"] - ddf["xF2D3"]) / ddf["xF2D3"]
    ddf.to_csv(OUT / "sigma_diff_gcft_vs_data.csv", index=False)

    # --- lightweight MCMC on eta (same as notebook MAP setup) ---
    from scipy import optimize as spo

    sigma_eta = 0.04
    mu_log10_x0 = -3.0
    w_xion = 1.15
    x0_map = 10.0**mu_log10_x0
    phi = np.exp(-0.5 * (np.log(xv / x0_map) / w_xion) ** 2)

    def neg_log_post(eta):
        pred2 = np.exp(np.clip(log_perq + eta * phi, -80.0, 80.0))
        chi = np.sum(((f2v - pred2) / sig_f2) ** 2)
        return 0.5 * chi + 0.5 * (eta / sigma_eta) ** 2

    def logpost(eta):
        return -neg_log_post(eta)

    rng = np.random.default_rng(0)
    n = 18_000
    eta_tr = np.zeros(n)
    cur = 0.0
    lp_cur = logpost(cur)
    acc = 0
    for i in range(n):
        prop = cur + rng.normal(0, 0.015)
        if abs(prop) > 0.35:
            eta_tr[i] = cur
            continue
        lp_p = logpost(prop)
        if np.log(rng.random()) < lp_p - lp_cur:
            cur = prop
            lp_cur = lp_p
            acc += 1
        eta_tr[i] = cur
    np.savez_compressed(OUT / "mcmc_eta_trace.npz", eta=eta_tr, accept_rate=acc / n)
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(eta_tr, lw=0.4, alpha=0.8)
    ax.set_xlabel("step")
    ax.set_ylabel(r"$\eta$")
    ax.set_title("Metropolis–Hastings trace (xion amplitude, diagonal σ)")
    fig.tight_layout()
    fig.savefig(OUT / "mcmc_eta_trace.png", dpi=140)
    plt.close()

    print("Wrote:", list(OUT.iterdir()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
