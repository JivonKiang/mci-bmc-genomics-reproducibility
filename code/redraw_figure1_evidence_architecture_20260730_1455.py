from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PKG = Path(os.environ.get("MN_SUBMISSION_PACKAGE", "."))
MAIN = PKG / "03_Figures"
MERGED = MAIN / "merged_source"
OUT_NAME = "Figure1_study_design_nature_rebuild"

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 9,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
    }
)


def add_box(ax, x, y, w, h, title, body, face, edge, title_size=11.7, body_size=8.2):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.7,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h * 0.72,
        title,
        ha="center",
        va="center",
        fontsize=title_size,
        fontweight="bold",
        color=edge,
        linespacing=1.1,
    )
    ax.text(
        x + w / 2,
        y + h * 0.34,
        body,
        ha="center",
        va="center",
        fontsize=body_size,
        color="#272727",
        linespacing=1.2,
    )


def arrow(ax, start, end, rad=0.0, lw=1.35):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=lw,
            color="#566573",
            shrinkA=4,
            shrinkB=5,
            connectionstyle=f"arc3,rad={rad}",
        )
    )


def draw():
    fig, ax = plt.subplots(figsize=(12, 9.4))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    blue_face, blue_edge = "#EAF2F8", "#2C3E50"
    teal_face, teal_edge = "#E8F8F5", "#117A65"
    peach_face, peach_edge = "#FDEBD0", "#B9770E"
    violet_face, violet_edge = "#F5EEF8", "#7D3C98"
    red_face, red_edge = "#FDEDEC", "#922B21"
    neutral_face, neutral_edge = "#F4F6F7", "#566573"

    ax.text(
        6,
        9.65,
        "MCI-centred evidence architecture",
        ha="center",
        va="center",
        fontsize=17,
        fontweight="bold",
        color="#172B3A",
    )

    add_box(
        ax,
        0.45,
        8.55,
        11.1,
        0.75,
        "Primary question and evidence hierarchy",
        "Can a de novo blood programme mark the MCI stage, remain stable within development batches and inform progression or mechanism without overclaiming?",
        blue_face,
        blue_edge,
        title_size=12.2,
        body_size=8.8,
    )

    x = [0.45, 4.25, 8.05]
    w = 3.5
    ys = [6.75, 5.10, 3.45]
    h = [1.15, 1.25, 1.25]
    centers = [v + w / 2 for v in x]

    # Column 1: the primary stage-programme path.
    add_box(
        ax,
        x[0],
        ys[0],
        w,
        h[0],
        "Stage discovery and cohort audit",
        "GSE63060 n = 155; GSE63061 n = 142\nHC, MCI and AD labels retained\nage and sex adjusted contrasts",
        blue_face,
        blue_edge,
    )
    add_box(
        ax,
        x[0],
        ys[1],
        w,
        h[1],
        "Candidate funnel",
        "95-gene registry\n40 MCI-peak + 55 MCI-trough\n83 direction-consistent in GSE63061\n12-gene lock: 6 peak + 6 trough",
        teal_face,
        teal_edge,
    )
    add_box(
        ax,
        x[0],
        ys[2],
        w,
        h[2],
        "Locked score output",
        "MCI peak and AD trough score\nwithin-system stability in AddNeuroMed\nfixed weights; no clinical threshold",
        teal_face,
        teal_edge,
    )

    # Column 2: external and progression tests.
    add_box(
        ax,
        x[1],
        ys[0],
        w,
        h[0],
        "Independent and progression cohorts",
        "GSE249477 n = 62 blood RNA-seq\nGSE282742 86 unique subjects\n52 baseline P/S and 11 paired P-MCI/AD",
        peach_face,
        peach_edge,
        title_size=11.2,
    )
    add_box(
        ax,
        x[1],
        ys[1],
        w,
        h[1],
        "Label and mapping audit",
        "Alias-aware mapping: 90/95 registry\n11/12 locked-panel entries\nsubject and paired-subject units preserved\nno endpoint relabelling",
        peach_face,
        peach_edge,
        title_size=11.2,
    )
    add_box(
        ax,
        x[1],
        ys[2],
        w,
        h[2],
        "External and progression results",
        "GSE249477: 10/12 mapped; directions reversed\nGSE282742: no FDR-positive feature\npaired P-MCI to AD score P = 0.558",
        red_face,
        red_edge,
        title_size=11.2,
    )

    # Column 3: orthogonal biology and mechanism probes.
    add_box(
        ax,
        x[2],
        ys[0],
        w,
        h[0],
        "Cellular and orthogonal layers",
        "GSE134578: 22,428 cells, 13 samples\nGSE285831 brain LCM RNA\nGSE136243 and GSE150693 context",
        violet_face,
        violet_edge,
        title_size=11.2,
    )
    add_box(
        ax,
        x[2],
        ys[1],
        w,
        h[1],
        "Single-cell analysis stack",
        "UMAP and cell composition\nCellChat source-target summaries\nDoRothEA TF activity\ndonor-aware pseudobulk and statistics",
        violet_face,
        violet_edge,
        title_size=11.2,
    )
    add_box(
        ax,
        x[2],
        ys[2],
        w,
        h[2],
        "Mechanism probes",
        "2,804 MCI CD8-state cells from two donors\nvirtual perturbation: no biological-response gate\n103 cis-eQTL rows; MR, network and structure remain exploratory",
        violet_face,
        violet_edge,
        title_size=11.2,
    )

    add_box(
        ax,
        0.75,
        0.70,
        10.5,
        1.25,
        "Evidence synthesis and interpretation boundary",
        "Supported: development-system stability.  Not supported: cross-cohort transportability.\nUnresolved: clinical progression and causal mechanism; no clinical biomarker or therapeutic target is claimed.",
        neutral_face,
        neutral_edge,
        title_size=12.0,
        body_size=9.0,
    )

    # Edges terminate at box boundaries, so arrows never cross text.
    for cx in centers:
        arrow(ax, (6, 8.55), (cx, 7.90), rad=(cx - 6) * 0.002)
    for i, cx in enumerate(centers):
        arrow(ax, (cx, ys[0]), (cx, ys[1] + h[1]))
        arrow(ax, (cx, ys[1]), (cx, ys[2] + h[2]))
        arrow(ax, (cx, ys[2]), (6, 1.95), rad=(cx - 6) * 0.012)

    ax.text(2.20, 8.28, "PRIMARY STAGE PROGRAMME", ha="center", va="center", fontsize=7.4, fontweight="bold", color=blue_edge)
    ax.text(6.00, 8.28, "EXTERNAL AND PROGRESSION TESTS", ha="center", va="center", fontsize=7.4, fontweight="bold", color=peach_edge)
    ax.text(9.80, 8.28, "ORTHOGONAL CONTEXT", ha="center", va="center", fontsize=7.4, fontweight="bold", color=violet_edge)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    for out_dir in (MAIN, MERGED):
        out_dir.mkdir(parents=True, exist_ok=True)
        base = out_dir / OUT_NAME
        fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
        fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")

    plt.close(fig)
    print(MERGED / f"{OUT_NAME}.svg")


if __name__ == "__main__":
    draw()
