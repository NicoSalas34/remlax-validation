"""Figures of a run_cascade.R test (sum framework).

    python3 simulation/cogwas/fig_sumtest.py <run_dir> <out_dir> [--rep 1]

fig_recovery   BLUP accuracy and log relative error of the variances by effect,
               over the replicates, for both null models (as in Salas et al. 2026)
fig_blups      simulated against predicted values, one replicate
fig_manhattan  GWAS of one replicate (interaction null model)
fig_cogwas     strongest co-GWAS pairs of one replicate, focal against neighbour position
"""
import glob
import os
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.expanduser("~/.claude/skills/figure-style"))
from figure_style import apply_figure_style, panel_letter  # noqa: E402

MODELS = {"additive": ("Additive null model", "#8c8c8c"),
          "interaction": ("Null model with DGE x IGE interaction", "#1f5fa8")}
EFFECTS = ["DGE", "IGE intra", "IGE inter", "GxG intra", "GxG inter"]
VAR_OF = {"var_DGE": "DGE", "var_IGE_intra": "IGE intra", "var_IGE_inter": "IGE inter",
          "var_GxG_intra": "GxG intra", "var_GxG_inter": "GxG inter", "var_residual": "Residual"}
TICK = {"DGE": "DGE", "IGE intra": "IGE\nwheat", "IGE inter": "IGE\nalfalfa",
        "GxG intra": "DGE x IGE\nwheat", "GxG inter": "DGE x IGE\nalfalfa", "Residual": "Residual"}


def load(run, key):
    files = sorted(glob.glob(os.path.join(run, f"*_{key}.csv")))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame()


def boxes(ax, df, cats, value, models, floor=None):
    """Boxes with raw points, one pair of boxes per category. Values below
    `floor` (estimates at the zero boundary on a log scale) are drawn as open
    triangles at the floor."""
    w = 0.34
    for j, m in enumerate(models):
        col = MODELS[m][1]
        for i, c in enumerate(cats):
            v = df.loc[(df.variant == m) & (df.cat == c), value].dropna().to_numpy()
            if not len(v):
                continue
            if floor is not None:
                v = np.maximum(v, floor - 1)
            x = i + (j - (len(models) - 1) / 2) * w
            ax.boxplot(np.maximum(v, floor) if floor is not None else v, positions=[x], widths=w * 0.8, showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor=mpl.colors.to_rgba(col, 0.25), edgecolor=col),
                       medianprops=dict(color=col, lw=1.4), whiskerprops=dict(color=col),
                       capprops=dict(color=col))
            rng = np.random.default_rng(i * 10 + j)
            jx = x + rng.uniform(-w * 0.25, w * 0.25, len(v))
            low = v < floor if floor is not None else np.zeros(len(v), bool)
            ax.scatter(jx[~low], v[~low], s=7, color=col, zorder=3, lw=0)
            if low.any():
                ax.scatter(jx[low], np.full(low.sum(), floor), s=16, marker="v", facecolors="none",
                           edgecolors=col, lw=0.8, zorder=3, clip_on=False)
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([TICK[c] for c in cats])
    ax.set_xlim(-0.6, len(cats) - 0.4)


def fig_recovery(run, out):
    acc = load(run, "acc").rename(columns={"effect": "cat"})
    vc = load(run, "vc")
    nrep = acc.rep.nunique()
    v = vc[vc.parameter.isin(VAR_OF)].copy()
    v["cat"] = v.parameter.map(VAR_OF)
    v["lre"] = np.log(v.estimated / v.simulated)
    rg = vc[vc.parameter == "r_g"].copy()
    rg["cat"] = "DGE"
    rg["err"] = rg.estimated - rg.simulated
    models = [m for m in MODELS if m in set(acc.variant)]
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[6, 1.3], hspace=0.55, wspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    boxes(ax1, acc, EFFECTS, "accuracy", models)
    ax1.set_ylabel("BLUP accuracy\n(correlation with simulated)")
    ax1.set_title("Additive effects are predicted well, interaction effects poorly", loc="left")
    ax1.set_ylim(min(-0.05, acc.accuracy.min() - 0.05), 1.02)
    ax1.axhline(0, color="0.6", lw=0.6, zorder=0)
    ax2 = fig.add_subplot(gs[1, 0])
    cats = EFFECTS + ["Residual"]
    boxes(ax2, v, cats, "lre", models, floor=-4)
    ax2.set_ylim(-4.3, 3)
    ax2.axhline(0, color="0.3", lw=0.8, ls="--", zorder=0)
    ax2.set_ylabel("log(estimated / simulated)\nvariance")
    ax2.set_title("Additive variances are recovered, DGE x IGE variances are not", loc="left")
    ax3 = fig.add_subplot(gs[1, 1])
    boxes(ax3, rg, ["DGE"], "err", models)
    ax3.set_xticklabels(["r(DGE, IGE)"])
    ax3.axhline(0, color="0.3", lw=0.8, ls="--", zorder=0)
    ax3.set_ylabel("estimated - simulated")
    handles = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(MODELS[m][1], 0.25), edgecolor=MODELS[m][1],
                                 label=MODELS[m][0]) for m in models]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.99), ncol=1, frameon=False)
    for ax, l in zip([ax1, ax2, ax3], "abc"):
        panel_letter(ax, l, dx=-0.12 if ax is not ax3 else -0.6)
    fig.text(0.01, 0.005, f"{nrep} replicates of scenario S2; points are replicates; open triangles: variance estimated at zero", fontsize=6, color="0.35")
    save(fig, out, "fig_recovery")


def fig_blups(run, out, rep):
    bl = load(run, "blup_lines"); bf = load(run, "blup_fam"); gx = load(run, "gxg")
    tag = sorted(bl.tag.unique())[rep - 1]
    bl = bl[(bl.tag == tag) & (bl.variant == "interaction")]
    bf = bf[(bf.tag == tag) & (bf.variant == "interaction")]
    gx = gx[gx.tag == tag] if len(gx) else gx
    panels = [("DGE", bl.dge, bl.dge_hat, "lines"), ("IGE wheat", bl.ige, bl.ige_hat, "lines"),
              ("IGE alfalfa", bf.inter, bf.inter_hat, "families")]
    if len(gx):
        for e, lab in [("GxG intra", "DGE x IGE wheat"), ("GxG inter", "DGE x IGE alfalfa")]:
            g = gx[gx.effect == e]
            panels.append((lab, g.truth, g.blup, "pairs"))
    fig, axes = plt.subplots(1, len(panels), figsize=(7.2, 2.3))
    for ax, (lab, x, y, unit), l in zip(axes, panels, "abcde"):
        big = len(x) > 2000
        ax.scatter(x, y, s=1.5 if big else 6, alpha=0.15 if big else 0.7, color="#1f5fa8", lw=0, rasterized=big)
        r = np.corrcoef(x, y)[0, 1]
        nn = f"{len(x) / 1000:.1f}k" if len(x) >= 1000 else str(len(x))
        ax.set_title(f"{lab}\nr = {r:.2f}\n{nn} {unit}", loc="left")
        ax.set_xlabel("simulated")
        ax.axhline(0, color="0.8", lw=0.5, zorder=0); ax.axvline(0, color="0.8", lw=0.5, zorder=0)
        panel_letter(ax, l, dx=-0.25, dy=1.35)
    axes[0].set_ylabel("BLUP")
    fig.tight_layout(w_pad=0.8)
    save(fig, out, "fig_blups")


def genome_x(df, chr_col="chr", pos_col="pos", lengths=None):
    lengths = lengths if lengths is not None else df.groupby(chr_col)[pos_col].max()
    off = lengths.cumsum() - lengths
    return df[pos_col].to_numpy() + off.reindex(df[chr_col]).to_numpy(), off, lengths


def manhattan(ax, d, col, loci, eff, title, lengths):
    x, off, L = genome_x(d, lengths=lengths)
    p = d[col].to_numpy(); ok = np.isfinite(p)
    lp = -np.log10(p)
    cols = np.where(d.chr.to_numpy() % 2 == 0, "#7b9cc4", "#2b4d75")
    ax.scatter(x[ok], lp[ok], s=2, c=cols[ok], lw=0, rasterized=True)
    thr = -np.log10(0.05 / ok.sum())
    ax.axhline(thr, color="#c0392b", lw=0.7, ls="--")
    lo = loci[loci.effect == eff]
    if len(lo):
        lx = lo.pos.to_numpy() + off.reindex(lo.chr).to_numpy()
        ax.scatter(lx, np.full(len(lx), -0.06 * max(thr, np.nanmax(lp[ok]))), marker="^", s=18,
                   color="#e08214", clip_on=False, zorder=4)
    mids = off + L / 2
    ax.set_xticks(mids.to_numpy()); ax.set_xticklabels(mids.index.astype(str))
    ax.set_xlim(0, (off + L).max())
    ax.set_ylim(0, max(thr, np.nanmax(lp[ok])) * 1.1)
    ax.set_title(title, loc="left")
    ax.set_ylabel("-log10 p")


def fig_manhattan(run, out, rep):
    sc = load(run, "scan"); loci_all = load(run, "loci")
    tag = sorted(load(run, "truth").tag.unique())[rep - 1]
    files = sorted(glob.glob(os.path.join(run, f"{tag}_scan.csv")))
    sc = pd.read_csv(files[0]); sc = sc[sc.variant == "interaction"]
    loci = pd.read_csv(os.path.join(run, f"{tag}_loci.csv"))
    w = sc[sc.species == "wheat"]; a = sc[sc.species == "alfalfa"]
    Lw = w.groupby("chr").pos.max(); La = a.groupby("chr").pos.max()
    fig, axes = plt.subplots(4, 1, figsize=(7.2, 6.4))
    manhattan(axes[0], w, "p_dir", loci, "major direct", "Direct effect, wheat chip", Lw)
    manhattan(axes[1], w, "p_ind", loci, "major indirect intra", "Indirect effect on wheat neighbours, wheat chip", Lw)
    manhattan(axes[2], w, "p_sim", loci, "similarity", "Similarity given direct and indirect effects, wheat chip", Lw)
    manhattan(axes[3], a, "p_ind", loci, "major indirect inter", "Indirect effect of alfalfa on wheat, alfalfa family chip", La)
    axes[3].set_xlabel("chromosome")
    for ax, l in zip(axes, "abcd"):
        panel_letter(ax, l, dx=-0.08, dy=1.08)
    fig.text(0.01, 0.005, "Triangles: simulated QTL (off chip). Dashed line: Bonferroni 5 %.",
             fontsize=6, color="0.35")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    save(fig, out, "fig_manhattan")


def fig_cogwas(run, out, rep):
    tag = sorted(load(run, "truth").tag.unique())[rep - 1]
    cg = pd.read_csv(os.path.join(run, f"{tag}_cg.csv")); cg = cg[cg.variant == "interaction"]
    sc = pd.read_csv(os.path.join(run, f"{tag}_scan.csv"))
    loci = pd.read_csv(os.path.join(run, f"{tag}_loci.csv"))
    Lw = sc[sc.species == "wheat"].groupby("chr").pos.max()
    La = sc[sc.species == "alfalfa"].groupby("chr").pos.max()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    spec = [("intra", Lw, "wheat neighbour", [("pair intra focal", "pair intra neighbour", "pair"),
                                              ("similarity", "similarity", "similarity")]),
            ("inter", La, "alfalfa neighbour", [("pair inter focal", "pair inter neighbour", "pair")])]
    vmax = -np.log10(cg.p.min())
    for ax, (kind, Ln, nlab, causal), l in zip(axes, spec, "ab"):
        d = cg[cg.kind == kind].sort_values("p", ascending=False)
        fx, offw, _ = genome_x(d, "focal_chr", "focal_pos", Lw)
        nx, offn, _ = genome_x(d, "neigh_chr", "neigh_pos", Ln)
        lp = -np.log10(d.p.to_numpy())
        sig = d.p.to_numpy() < d.threshold.to_numpy()
        im = ax.scatter(fx, nx, c=lp, s=np.where(sig, 14, 2), cmap="viridis", vmin=lp.min(), vmax=vmax, lw=0)
        for fe, ne, lab in causal:
            f = loci[loci.effect == fe]; nn = loci[loci.effect == ne]
            for k in range(len(f)):
                cx = f.pos.iloc[k] + offw[f.chr.iloc[k]]
                ny = nn.pos.iloc[min(k, len(nn) - 1)] + offn[nn.chr.iloc[min(k, len(nn) - 1)]]
                ax.scatter(cx, ny, s=90, facecolors="none", edgecolors="#e08214", lw=1.2, zorder=5)
        for o in offw.to_numpy()[1:]:
            ax.axvline(o, color="0.88", lw=0.5, zorder=0)
        for o in offn.to_numpy()[1:]:
            ax.axhline(o, color="0.88", lw=0.5, zorder=0)
        ax.set_xlim(0, (offw + Lw).max()); ax.set_ylim(0, (offn + Ln).max())
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("focal wheat SNP (genome position)")
        ax.set_ylabel(f"{nlab} SNP (genome position)")
        n_sig = int(sig.sum())
        ax.set_title(f"{'Within wheat' if kind == 'intra' else 'Wheat x alfalfa'}: {n_sig} significant pairs",
                     loc="left")
        panel_letter(ax, l, dx=-0.08)
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("-log10 p")
    fig.text(0.01, 0.005, "5,000 strongest pairs per co-GWAS; grey lines: chromosomes; large dots pass Bonferroni; "
             "orange circles: simulated pairs (diagonal: similarity).", fontsize=6, color="0.35")
    save(fig, out, "fig_cogwas")


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300)
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    over = [(a.get_text()[:30], b.get_text()[:30]) for i, (a, ba) in enumerate(texts)
            for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    print(name, "text overlaps:", over if over else "none")
    plt.close(fig)


if __name__ == "__main__":
    run, out = sys.argv[1], sys.argv[2]
    rep = int(sys.argv[sys.argv.index("--rep") + 1]) if "--rep" in sys.argv else 1
    apply_figure_style()
    fig_recovery(run, out)
    fig_blups(run, out, rep)
    fig_manhattan(run, out, rep)
    fig_cogwas(run, out, rep)
