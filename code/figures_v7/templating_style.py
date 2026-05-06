"""shared paper style for v7 figures.

every panel script in code/figures_v7/ should:
    from templating_style import apply_paper_style, OKABE_ITO
    apply_paper_style()

panels target 85 mm x 60 mm physical size (figsize ~ (3.346, 2.362) inches).
fonts are 8pt for labels, 7pt for ticks; pdf.fonttype=42 keeps text editable
in inkscape for downstream composite assembly.
"""

import seaborn as sns
import matplotlib.pyplot as plt


def apply_paper_style():
    """apply paper-grade rcParams and seaborn context globally."""
    sns.set_context("paper")
    sns.set_style("ticks")
    plt.rcParams.update({
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 1.0,
        "lines.linewidth": 1.5,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,   # editable text in pdf (truetype, not type-3)
        "ps.fonttype": 42,
    })


# okabe-ito palette: 8 colorblind-safe categorical colors
OKABE_ITO = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "skyblue": "#56B4E9",
    "black": "#000000",
}
