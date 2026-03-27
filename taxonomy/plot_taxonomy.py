#!/usr/bin/env python3
"""
plot_taxonomy.py

  1. Cladograma circular  -- Dominio > Reino > Filo
     Filos rotulados sao marcados com simbolos numerados (sem texto no circulo)
     Contagem de PDBs mantida no anel colorido de cada dominio
     Painel de legenda lateral explica elementos e lista os filos

  2. Treemap              -- Dominio > Reino
     Texto condicional por tamanho de celula; legenda fora do grafico

Saidas:
  taxonomy/taxonomy_cladogram.png
  taxonomy/taxonomy_treemap.png
  data_exploration/imgs/taxonomy_cladogram.png  (copia para o README)
  data_exploration/imgs/taxonomy_treemap.png    (copia para o README)
"""

import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squarify

SCRIPT_DIR = Path(__file__).parent
OUT_DIR    = SCRIPT_DIR
IMG_DIR    = SCRIPT_DIR.parent / "data_exploration" / "imgs"

DOMAIN_ORDER = ["Eukaryota", "Bacteria", "Archaea", "Viruses", "Unclassified"]
DOMAIN_COLORS = {
    "Eukaryota":    "#5B9BD5",   # azul empoeirado
    "Bacteria":     "#E07B54",   # terracota suave
    "Archaea":      "#5FAD74",   # verde sage
    "Viruses":      "#9B6BB5",   # lavanda
    "Unclassified": "#9DAAB2",   # cinza azulado
}

# Cladograma -- radii
R_PHYLUM     = 1.00
R_KINGDOM    = 0.60
R_DOMAIN     = 0.26
R_BAND_INNER = 1.10
R_BAND_OUTER = 1.18
GAP          = 0.030   # gap angular entre dominios (fracao de 2pi)
LABEL_MIN    = 50      # filos com n_pdbs >= LABEL_MIN recebem simbolo numerado


# ---------------------------------------------------------------------------
# Dados
# ---------------------------------------------------------------------------
def load_phyla() -> pd.DataFrame:
    summary = pd.read_csv(OUT_DIR / "taxonomy_summary.csv")
    lineage = pd.read_csv(OUT_DIR / "taxonomy_lineage.csv")
    merged = summary.merge(
        lineage[["ncbi_taxonomy_id", "domain", "kingdom", "phylum"]],
        on="ncbi_taxonomy_id", how="left",
    )
    for col in ("domain", "kingdom", "phylum"):
        merged[col] = merged[col].fillna("Unclassified")
    phyla = merged.groupby(["domain", "kingdom", "phylum"])["n_pdbs"].sum().reset_index()
    dom_ord = {d: i for i, d in enumerate(DOMAIN_ORDER)}
    phyla["_d"] = phyla["domain"].map(dom_ord).fillna(99)
    phyla = (phyla.sort_values(["_d", "kingdom", "n_pdbs"], ascending=[True, True, False])
             .drop(columns="_d").reset_index(drop=True))
    return phyla


# ---------------------------------------------------------------------------
# Cladograma
# ---------------------------------------------------------------------------
def _arc(ax, r, a1, a2, color, lw, alpha=1.0, zorder=3):
    if a1 > a2:
        a1, a2 = a2, a1
    n = max(4, int((a2 - a1) / 0.008))
    theta = np.linspace(a1, a2, n)
    ax.plot(theta, np.full(n, r), color=color, lw=lw,
            alpha=alpha, solid_capstyle="round", zorder=zorder)


def assign_angles(phyla: pd.DataFrame) -> np.ndarray:
    n_domains = phyla["domain"].nunique()
    total_gap = GAP * 2 * np.pi * n_domains
    angle_per_leaf = (2 * np.pi - total_gap) / len(phyla)
    angles = np.empty(len(phyla))
    cur, prev_dom = 0.0, None
    for i, row in phyla.iterrows():
        if row["domain"] != prev_dom:
            if prev_dom is not None:
                cur += GAP * 2 * np.pi
            prev_dom = row["domain"]
        angles[i] = cur
        cur += angle_per_leaf
    return angles


def _draw_legend(ax_leg, phyla, phylum_nums):
    """Painel lateral: explicacao dos elementos + lista de filos numerados."""
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(0, 1)
    ax_leg.axis("off")

    y   = 0.98
    dh  = 0.036   # espaco entre linhas
    xi  = 0.05    # x do icone
    xt  = 0.22    # x do texto
    GRAY = "#444444"

    def title(text):
        nonlocal y
        ax_leg.text(xi, y, text, fontsize=16, fontweight="bold",
                    va="top", color="#111111")
        y -= dh * 1.7

    def row_marker(label, marker, ms, fc, ec="white"):
        nonlocal y
        ax_leg.plot(xi + 0.07, y, marker, ms=ms, color=fc,
                    markeredgewidth=0.6, markeredgecolor=ec, zorder=5,
                    transform=ax_leg.transData, clip_on=False)
        ax_leg.text(xt, y, label, va="center", fontsize=14, color=GRAY)
        y -= dh

    def row_line(label, lw, color):
        nonlocal y
        ax_leg.plot([xi, xi + 0.13], [y, y], lw=lw, color=color,
                    solid_capstyle="round")
        ax_leg.text(xt, y, label, va="center", fontsize=14, color=GRAY)
        y -= dh

    def row_band(label, color):
        nonlocal y
        ax_leg.add_patch(mpatches.FancyBboxPatch(
            (xi, y - 0.007), 0.13, 0.015,
            boxstyle="square,pad=0", fc=color, ec="none", alpha=0.85,
            transform=ax_leg.transData, clip_on=False))
        ax_leg.text(xt, y, label, va="center", fontsize=14, color=GRAY)
        y -= dh

    def row_symbol(label, color):
        nonlocal y
        ax_leg.text(xi + 0.07, y, "①", ha="center", va="center",
                    fontsize=12, color="white", fontweight="bold",
                    bbox=dict(boxstyle="circle,pad=0.20", fc=color,
                              ec="white", lw=0.5),
                    transform=ax_leg.transData, clip_on=False)
        ax_leg.text(xt, y, label, va="center", fontsize=14, color=GRAY)
        y -= dh

    def sep():
        nonlocal y
        y -= dh * 0.4
        ax_leg.axhline(y, xmin=0.02, xmax=0.98, lw=0.5, color="#dddddd")
        y -= dh * 0.6

    # ---- Elementos ----
    title("Elements")
    row_marker("Phylum (leaf node)",            "o",  5,  "#888888")
    row_marker("Kingdom (internal node)",        "o",  9,  "#888888")
    row_marker("Domain (internal node)",         "D", 12,  "#888888")
    row_marker("Root (common ancestor)",         "o", 12,  "#111111", ec="#111111")
    row_line  ("Thin branch: phylum → kingdom",  0.7, "#aaaaaa")
    row_line  ("Branch: kingdom → domain → root", 2.2, "#888888")
    row_line  ("Arc: sibling grouping at each level", 1.6, "#2196F3")
    row_band  ("Domain ring (n = total PDBs in domain)", "#F4511E")
    row_symbol("Numbered phylum (see table below)", "#8E24AA")

    sep()

    # ---- Dominios ----
    title("Domains of Life")
    for dom in DOMAIN_ORDER:
        if dom not in phyla["domain"].values:
            continue
        color = DOMAIN_COLORS[dom]
        n = phyla[phyla["domain"] == dom]["n_pdbs"].sum()
        ax_leg.add_patch(mpatches.FancyBboxPatch(
            (xi, y - 0.007), 0.12, 0.015,
            boxstyle="square,pad=0", fc=color, ec="none",
            transform=ax_leg.transData, clip_on=False))
        ax_leg.text(xt, y, f"{dom}  ({n:,} PDBs)",
                    va="center", fontsize=14, color=color, fontweight="bold")
        y -= dh

    sep()

    # ---- Filos numerados (duas colunas) ----
    title(f"Labeled Phyla  (≥{LABEL_MIN} PDBs)")

    sorted_nums = sorted(phylum_nums.items(), key=lambda x: x[1])
    n_total = len(sorted_nums)
    split   = (n_total + 1) // 2

    # Coluna 1: itens 1..split   Coluna 2: itens split+1..n_total
    col_cfg = [
        (sorted_nums[:split],       xi,       xt,       y),
        (sorted_nums[split:], xi + 0.50, xt + 0.50, y),
    ]

    y_after = y
    for items, xi_c, xt_c, y_c in col_cfg:
        yy = y_c
        for idx, num in items:
            row = phyla.loc[idx]
            color = DOMAIN_COLORS.get(row["domain"], "#9E9E9E")
            ax_leg.text(xi_c + 0.07, yy, str(num),
                        ha="center", va="center", fontsize=11,
                        color="white", fontweight="bold",
                        bbox=dict(boxstyle="circle,pad=0.24",
                                  fc=color, ec="white", lw=0.5),
                        transform=ax_leg.transData, clip_on=False)
            name = row["phylum"]
            if len(name) > 18:
                name = name[:17] + "…"
            ax_leg.text(xt_c, yy,
                        f"{name}  ({int(row['n_pdbs']):,})",
                        va="center", fontsize=13, color=color)
            yy -= dh
        y_after = min(y_after, yy)

    # Atualiza y global para baixo das duas colunas
    y = y_after


def plot_cladogram(phyla: pd.DataFrame, angles: np.ndarray) -> None:
    # Figura two-panel: cladograma (esq) + legenda (dir)
    fig = plt.figure(figsize=(32, 20))
    gs  = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.0], wspace=0.0)
    ax      = fig.add_subplot(gs[0], projection="polar")
    ax_leg  = fig.add_subplot(gs[1])

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    kingdom_groups: dict = {}
    domain_groups:  dict = {}
    for i, row in phyla.iterrows():
        kingdom_groups.setdefault((row["domain"], row["kingdom"]), []).append(i)
        domain_groups.setdefault(row["domain"], []).append(i)

    # Mapeamento filo -> numero (top filos por n_pdbs)
    labeled = phyla[phyla["n_pdbs"] >= LABEL_MIN].sort_values("n_pdbs", ascending=False)
    phylum_nums = {idx: num + 1 for num, idx in enumerate(labeled.index)}

    # ---- ramos filo -> reino ----
    for (dom, _), idxs in kingdom_groups.items():
        color = DOMAIN_COLORS.get(dom, "#9E9E9E")
        angs  = angles[idxs]
        a_min, a_max = angs.min(), angs.max()
        for a in angs:
            ax.plot([a, a], [R_KINGDOM, R_PHYLUM],
                    color=color, lw=0.7, alpha=0.45, zorder=2)
        _arc(ax, R_KINGDOM, a_min, a_max, color, lw=1.4)
        ax.plot((a_min + a_max) / 2, R_KINGDOM, "o", color=color, ms=5.5,
                markeredgewidth=0.5, markeredgecolor="white", zorder=6)

    # ---- ramos reino -> dominio ----
    for dom, idxs in domain_groups.items():
        color = DOMAIN_COLORS.get(dom, "#9E9E9E")
        angs  = angles[idxs]
        a_min, a_max = angs.min(), angs.max()
        for (d2, _), kidxs in kingdom_groups.items():
            if d2 != dom:
                continue
            k_angs = angles[kidxs]
            k_ang  = (k_angs.min() + k_angs.max()) / 2
            ax.plot([k_ang, k_ang], [R_DOMAIN, R_KINGDOM],
                    color=color, lw=1.5, alpha=0.7, zorder=2)
        _arc(ax, R_DOMAIN, a_min, a_max, color, lw=2.5)
        d_ang = (a_min + a_max) / 2
        ax.plot(d_ang, R_DOMAIN, "D", color=color, ms=10,
                markeredgewidth=0.8, markeredgecolor="white", zorder=7)
        ax.plot([d_ang, d_ang], [0.04, R_DOMAIN], color=color, lw=2.5, zorder=2)

    # ---- raiz ----
    ax.plot(0, 0, "o", color="#111111", ms=12, zorder=10)

    # ---- folhas: ponto simples ou simbolo numerado ----
    for i, row in phyla.iterrows():
        color = DOMAIN_COLORS.get(row["domain"], "#9E9E9E")
        a = angles[i]
        if i in phylum_nums:
            ax.text(a, R_PHYLUM, str(phylum_nums[i]),
                    ha="center", va="center", fontsize=12,
                    color="white", fontweight="bold", zorder=8,
                    bbox=dict(boxstyle="circle,pad=0.25",
                              fc=color, ec="white", lw=0.8))
        else:
            ax.plot(a, R_PHYLUM, "o", color=color, ms=3.5,
                    markeredgewidth=0.4, markeredgecolor="white", zorder=5)

    # ---- anel de dominio com contagem em branco ----
    for dom, idxs in domain_groups.items():
        color = DOMAIN_COLORS.get(dom, "#9E9E9E")
        angs  = angles[idxs]
        a_min, a_max = angs.min(), angs.max()
        n     = max(4, int((a_max - a_min) / 0.008))
        theta = np.linspace(a_min, a_max, n)
        ax.fill_between(theta, R_BAND_INNER, R_BAND_OUTER,
                        color=color, alpha=0.88, zorder=4)
        # contagem de PDBs centrada no anel (texto branco)
        d_ang     = (a_min + a_max) / 2
        n_pdbs_dom = phyla.loc[idxs, "n_pdbs"].sum()
        rot = -np.degrees(d_ang)
        if np.pi / 2 < d_ang < 3 * np.pi / 2:
            rot += 180
        ax.text(d_ang, (R_BAND_INNER + R_BAND_OUTER) / 2,
                f"{n_pdbs_dom:,}",
                ha="center", va="center", rotation=rot,
                rotation_mode="anchor",
                fontsize=14, fontweight="bold", color="white", zorder=5)

    ax.set_ylim(0, 1.35)
    ax.axis("off")
    ax.set_title(
        "Taxonomic Diversity of S40 Non-redundant CATH Structures\n"
        "NCBI hierarchy: Domain  >  Kingdom  >  Phylum",
        fontsize=20, pad=18, color="#333333",
    )

    # ---- painel de legenda ----
    _draw_legend(ax_leg, phyla, phylum_nums)

    out = OUT_DIR / "taxonomy_cladogram.png"
    plt.savefig(out, dpi=700, bbox_inches="tight", facecolor="white")
    plt.close()
    shutil.copy(out, IMG_DIR / out.name)
    print(f"Cladograma salvo: {out}")


# ---------------------------------------------------------------------------
# Treemap
# ---------------------------------------------------------------------------
def plot_treemap(phyla: pd.DataFrame) -> None:
    dom_ord = {d: i for i, d in enumerate(DOMAIN_ORDER)}
    kingdoms = (phyla.groupby(["domain", "kingdom"])["n_pdbs"].sum().reset_index())
    kingdoms["_d"] = kingdoms["domain"].map(dom_ord).fillna(99)
    kingdoms = (kingdoms.sort_values(["_d", "n_pdbs"], ascending=[True, False])
                .drop(columns="_d").reset_index(drop=True))

    W, H = 13.0, 8.5
    sizes  = kingdoms["n_pdbs"].values.astype(float)
    colors = [mpl.colors.to_rgba(DOMAIN_COLORS.get(d, "#9E9E9E"), alpha=0.88)
              for d in kingdoms["domain"]]
    norm  = squarify.normalize_sizes(sizes, W, H)
    rects = squarify.squarify(norm, x=0, y=0, dx=W, dy=H)

    fig, ax = plt.subplots(figsize=(W + 2.8, H))

    for rect, row, color in zip(rects, kingdoms.itertuples(), colors):
        dx, dy = rect["dx"], rect["dy"]
        ax.add_patch(mpatches.FancyBboxPatch(
            (rect["x"], rect["y"]), dx, dy,
            boxstyle="square,pad=0",
            facecolor=color, edgecolor="white", linewidth=1.8,
        ))
        cx, cy = rect["x"] + dx / 2, rect["y"] + dy / 2
        count  = f"({row.n_pdbs:,})"
        if dx >= 1.2 and dy >= 0.65:
            fs = min(10, max(6.5, dx * 2.0))
            ax.text(cx, cy, f"{row.kingdom}\n{count}",
                    ha="center", va="center", fontsize=fs,
                    color="white", fontweight="bold", multialignment="center",
                    linespacing=1.4)
        elif dx >= 0.5 and dy >= 0.35:
            ax.text(cx, cy, count, ha="center", va="center",
                    fontsize=6.5, color="white", fontweight="bold")

    # Legenda fora do plot
    handles = [mpatches.Patch(color=DOMAIN_COLORS[d], label=d)
               for d in DOMAIN_ORDER if d in kingdoms["domain"].values]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=True, framealpha=0.95, fontsize=11,
              title="Taxon Domains", title_fontsize=11, edgecolor="#cccccc")

    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_title(
        "Distribution of S40 Non-redundant CATH Structures by Kingdom\n"
        "(area proportional to number of PDBs)",
        fontsize=13, pad=12,
    )
    ax.axis("off")

    out = OUT_DIR / "taxonomy_treemap.png"
    plt.savefig(out, dpi=700, bbox_inches="tight", facecolor="white")
    plt.close()
    shutil.copy(out, IMG_DIR / out.name)
    print(f"Treemap salvo: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    phyla  = load_phyla()
    labeled = (phyla["n_pdbs"] >= LABEL_MIN).sum()
    print(f"Dominios: {phyla['domain'].nunique()}  |  "
          f"Reinos: {phyla['kingdom'].nunique()}  |  "
          f"Filos: {phyla['phylum'].nunique()}  |  "
          f"Filos numerados: {labeled}")
    angles = assign_angles(phyla)
    plot_cladogram(phyla, angles)
    plot_treemap(phyla)


if __name__ == "__main__":
    main()
