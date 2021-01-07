from matplotlib import rcParams
from collections import OrderedDict

# ---------------------------------------------------
# Color sets
# ---------------------------------------------------
# Standard tableau 20 set
tableau = OrderedDict(
    [
        ("blue", "#0173B2"),
        ("orange", "#DE8F05"),
        ("green", "#029E73"),
        ("red", "#D55E00"),
        ("purple", "#CC78BC"),
        ("brown", "#CA9161"),
        ("pink", "#FBAFE4"),
        ("grey", "#949494"),
        ("yellow", "#ECE133"),
        ("turquoise", "#56B4E9"),
    ]
)

fontsize = 20
nearly_black = "#161616"
light_grey = "#EEEEEE"
lighter_grey = "#F5F5F5"
white = "#ffffff"
grey = "#7F7F7F"

master_formatting = {
    "axes.formatter.limits": (-3, 3),
    "xtick.major.pad": 7,
    "ytick.major.pad": 7,
    "ytick.color": nearly_black,
    "xtick.color": nearly_black,
    "axes.labelcolor": nearly_black,
    "axes.spines.bottom": False,
    "axes.spines.left": False,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.axisbelow": True,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "font.size": fontsize,
    "font.family": "serif",
    # 'font.serif' : 'Source Serif Pro',
    "text.usetex": True,
    "savefig.bbox": "tight",
    "axes.facecolor": lighter_grey,
    "axes.labelpad": 10.0,
    "axes.labelsize": fontsize * 0.8,
    "axes.titlepad": 30,
    "axes.titlesize": fontsize,
    "axes.grid": False,
    "grid.color": white,
    "lines.markersize": 7.0,
    "lines.scale_dashes": False,
    "xtick.labelsize": fontsize * 0.8,
    "ytick.labelsize": fontsize * 0.8,
    "legend.fontsize": fontsize * 0.8,
    "lines.linewidth": 4,
}

for k, v in master_formatting.items():
    rcParams[k] = v

color_cycle = tableau.values()
try:
    from matplotlib import cycler

    rcParams["axes.prop_cycle"] = cycler(color=color_cycle)
except Exception:
    raise
