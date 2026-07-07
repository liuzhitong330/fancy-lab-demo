"""
Two SVG figures for the Fancy lab demo:
1. myelin_heatmap.svg  — Wnt / OPC / OL / vascular genes across 7 brain regions
                         (corpus callosum excluded — ISH zero for fiber tract)
2. opc_vascular.svg    — scatter: OPC index vs vascular index per region,
                         illustrating the vascular scaffold model
"""
import csv, math
from pathlib import Path

OUT = Path(__file__).parent

# Drop Corpus_callosum (all zeros) — fiber tract has no cell bodies for ISH
STRUCTURES = ["Isocortex","Hippocampus","Cerebellum","Striatum",
              "Hypothalamus","Olfact_bulb","Brainstem"]
STRUCT_LABELS = ["Isocortex","Hippocampus","Cerebellum","Striatum",
                 "Hypothalamus","Olf. bulb","Brainstem"]

GROUPS = [
    ("Wnt pathway",      ["Axin2","Ctnnb1","Fzd1","Tcf7l2","Wnt3a"],          "#e67e22"),
    ("OPC progenitors",  ["Pdgfra","Cspg4"],                                    "#2980b9"),
    ("OL lineage / myelin", ["Sox10","Mbp","Plp1","Mog","Cnp"],                "#27ae60"),
    ("Vascular",         ["Pecam1","Cdh5","Cldn5","Pdgfrb"],                   "#8e44ad"),
]
ALL_GENE_ORDER = [g for _, gs, _ in GROUPS for g in gs]
GENE_COLOR = {g: c for _, gs, c in GROUPS for g in gs}

rows_by_gene = {}
with open(OUT / "expression.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows_by_gene[r["gene"]] = [float(r[s]) for s in STRUCTURES]


def row_norm(vals):
    mn, mx = min(vals), max(vals)
    return [(v - mn) / (mx - mn) if mx > mn else 0.0 for v in vals]


def green_red(v):
    """0 → cool green/white, 1 → warm red (myelination color scheme)"""
    if v <= 0.5:
        t = v / 0.5
        r = int(235 + (255 - 235) * t)
        g = int(245 + (255 - 245) * t)
        b = int(235 + (255 - 235) * t)
        # low: pale green (#ebf5eb) → white
        r = int(235 + (255 - 235) * t)
        g = int(245 + (255 - 245) * t)
        b = int(235 + (255 - 235) * t)
    else:
        t = (v - 0.5) / 0.5
        r = 255
        g = int(255 + (70 - 255) * t)
        b = int(255 + (70 - 255) * t)
    return f"rgb({r},{g},{b})"


# ── Figure 1: Heatmap ─────────────────────────────────────────────────────────
CELL_W = 70
CELL_H = 26
LEFT   = 90
TOP    = 105
n_genes   = len(ALL_GENE_ORDER)
n_structs = len(STRUCTURES)
W = LEFT + n_structs * CELL_W + 70
H = TOP + n_genes * CELL_H + 100

BS_IDX = STRUCTURES.index("Brainstem")
HIPP_IDX = STRUCTURES.index("Hippocampus")

# Column highlights
col_bg = ""
for si, color, label in [(BS_IDX, "#fff8f0", "most myelinated"),
                          (HIPP_IDX, "#f0f8ff", "OPC-rich")]:
    x = LEFT + si * CELL_W
    col_bg += (f'<rect x="{x}" y="{TOP-54}" width="{CELL_W}" height="{n_genes*CELL_H+58}" '
               f'fill="{color}" opacity="0.7"/>')

cells = ""
group_dividers = ""
gi_abs = 0
for grp_label, gs, grp_color in GROUPS:
    for gene in gs:
        if gene not in rows_by_gene:
            gi_abs += 1; continue
        vals = rows_by_gene[gene]
        nv   = row_norm(vals)
        for si, v in enumerate(nv):
            x = LEFT + si * CELL_W
            y = TOP + gi_abs * CELL_H
            color = green_red(v)
            cells += (f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" '
                      f'fill="{color}" stroke="white" stroke-width="1.5"/>')
            raw = vals[si]
            if raw > 1.5 or (v > 0.82 and raw > 0.1):
                fc = "white" if v > 0.9 else "#333"
                cells += (f'<text x="{x+CELL_W/2:.0f}" y="{y+CELL_H/2+4:.0f}" '
                          f'text-anchor="middle" font-size="8" fill="{fc}">{raw:.1f}</text>')
        gi_abs += 1
    group_dividers += (f'<line x1="{LEFT}" y1="{TOP+gi_abs*CELL_H}" '
                       f'x2="{LEFT+n_structs*CELL_W}" y2="{TOP+gi_abs*CELL_H}" '
                       f'stroke="#bbb" stroke-width="1.5"/>')

# Column headers
col_hdrs = ""
for si, lbl in enumerate(STRUCT_LABELS):
    x = LEFT + si * CELL_W + CELL_W // 2
    fc = "#c05000" if si == BS_IDX else ("#1a5c8a" if si == HIPP_IDX else "#444")
    fw = "bold" if si in (BS_IDX, HIPP_IDX) else "normal"
    col_hdrs += (f'<text x="{x}" y="{TOP-12}" text-anchor="middle" font-size="10" '
                 f'fill="{fc}" font-weight="{fw}">{lbl}</text>')

# Column sublabels
col_hdrs += (f'<text x="{LEFT+BS_IDX*CELL_W+CELL_W//2}" y="{TOP-1}" text-anchor="middle" '
             f'font-size="8.5" fill="#c05000" font-style="italic">↑ most myelinated</text>')
col_hdrs += (f'<text x="{LEFT+HIPP_IDX*CELL_W+CELL_W//2}" y="{TOP-1}" text-anchor="middle" '
             f'font-size="8.5" fill="#1a5c8a" font-style="italic">↑ OPC-rich</text>')

# Row labels
row_lbls = ""
gi_abs = 0
for _, gs, grp_color in GROUPS:
    for gene in gs:
        y = TOP + gi_abs * CELL_H + CELL_H // 2 + 4
        row_lbls += (f'<text x="{LEFT-5}" y="{y}" text-anchor="end" font-size="10.5" '
                     f'fill="{GENE_COLOR[gene]}" font-weight="600" '
                     f'font-family="ui-monospace,Menlo,monospace">{gene}</text>')
        gi_abs += 1

# Colorbar
cb_x = LEFT + n_structs * CELL_W + 10
cb_h = n_genes * CELL_H
colorbar = (f'<defs><linearGradient id="cb1" x1="0" y1="1" x2="0" y2="0">'
            f'<stop offset="0%" stop-color="{green_red(0)}"/>'
            f'<stop offset="50%" stop-color="{green_red(0.5)}"/>'
            f'<stop offset="100%" stop-color="{green_red(1)}"/>'
            f'</linearGradient></defs>'
            f'<rect x="{cb_x}" y="{TOP}" width="13" height="{cb_h}" '
            f'fill="url(#cb1)" stroke="#ccc" stroke-width="0.5"/>'
            f'<text x="{cb_x+6}" y="{TOP-4}" text-anchor="middle" font-size="8" fill="#555">high</text>'
            f'<text x="{cb_x+6}" y="{TOP+cb_h+10}" text-anchor="middle" font-size="8" fill="#555">low</text>')

svg1 = f"""<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{W//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Wnt Pathway, OPC/OL Lineage, and Vascular Gene Expression Across Brain Regions
  </text>
  <text x="{W//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas · ISH expression energy (row-normalized per gene)
  </text>
  <text x="{W//2}" y="53" text-anchor="middle" font-size="10" fill="#555">
    Note: corpus callosum excluded — ISH is zero in fiber tracts (no cell bodies)
  </text>
  {col_bg}{colorbar}{col_hdrs}{cells}{group_dividers}{row_lbls}
</svg>"""

with open(OUT / "myelin_heatmap.svg", "w") as f:
    f.write(svg1)
print("Wrote myelin_heatmap.svg")


# ── Figure 2: OPC vs Vascular scatter ────────────────────────────────────────
# OPC index = mean(Pdgfra, Cspg4) normalized
# Vascular index = mean(Pecam1, Cdh5, Cldn5)
# Each dot = one brain region

opc_genes  = ["Pdgfra", "Cspg4"]
vasc_genes = ["Pecam1", "Cdh5", "Cldn5"]
myelin_genes = ["Mbp", "Plp1", "Mog"]

def mean_vals(gene_list, struct_idx):
    total = 0.0
    n = 0
    for g in gene_list:
        if g in rows_by_gene:
            total += rows_by_gene[g][struct_idx]
            n += 1
    return total / n if n else 0.0

opc_scores   = [mean_vals(opc_genes,    i) for i in range(len(STRUCTURES))]
vasc_scores  = [mean_vals(vasc_genes,   i) for i in range(len(STRUCTURES))]
myelin_scores = [mean_vals(myelin_genes, i) for i in range(len(STRUCTURES))]

SW, SH = 580, 460
PAD_L = 65
PAD_B = 70
PAD_T = 85
AW = SW - PAD_L - 40
AH = SH - PAD_B - PAD_T

max_opc  = max(opc_scores)  * 1.15
max_vasc = max(vasc_scores) * 1.15
max_myel = max(myelin_scores)

def sx2(v): return PAD_L + v / max_vasc * AW
def sy2(v): return PAD_T + AH - v / max_opc  * AH

REGION_COLORS = {
    "Isocortex":    "#2980b9",
    "Hippocampus":  "#c0392b",
    "Cerebellum":   "#27ae60",
    "Striatum":     "#e67e22",
    "Hypothalamus": "#8e44ad",
    "Olfact_bulb":  "#16a085",
    "Brainstem":    "#7f8c8d",
}

dots = ""
for i, reg in enumerate(STRUCTURES):
    cx = sx2(vasc_scores[i])
    cy = sy2(opc_scores[i])
    # dot size by myelin score (normalized)
    r = 7 + myelin_scores[i] / max_myel * 12
    col = REGION_COLORS[reg]
    dots += f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{col}" opacity="0.85"/>'
    # label
    lx = cx + r + 4
    ly = cy + 4
    # adjust label for brainstem (bottom right)
    if reg == "Brainstem":
        lx = cx - r - 4
        ly = cy - r - 4
        anchor = "end"
    else:
        anchor = "start"
    dots += (f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="10" fill="{col}" '
             f'text-anchor="{anchor}" font-weight="500">{reg.replace("_"," ")}</text>')

# Axes + ticks
yticks = ""
xticks = ""
for t in [0, 0.5, 1.0, 1.5, 2.0, 2.5]:
    if t > max_opc: break
    yy = sy2(t)
    yticks += (f'<line x1="{PAD_L-4}" y1="{yy:.1f}" x2="{PAD_L+AW}" y2="{yy:.1f}" '
               f'stroke="#f0f0f0" stroke-width="1"/>'
               f'<text x="{PAD_L-7}" y="{yy+3:.1f}" text-anchor="end" font-size="9" fill="#666">{t:.1f}</text>')
for t in [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    if t > max_vasc: break
    xx = sx2(t)
    xticks += (f'<line x1="{xx:.1f}" y1="{PAD_T}" x2="{xx:.1f}" y2="{PAD_T+AH}" '
               f'stroke="#f0f0f0" stroke-width="1"/>'
               f'<text x="{xx:.1f}" y="{PAD_T+AH+14}" text-anchor="middle" font-size="9" fill="#666">{t:.1f}</text>')

# Size legend
size_legend = (
    f'<circle cx="{PAD_L+AW-50}" cy="{PAD_T+30}" r="7" fill="#ccc" opacity="0.7"/>'
    f'<circle cx="{PAD_L+AW-20}" cy="{PAD_T+30}" r="19" fill="#ccc" opacity="0.7"/>'
    f'<text x="{PAD_L+AW-50}" y="{PAD_T+50}" text-anchor="middle" font-size="8.5" fill="#666">low myelin</text>'
    f'<text x="{PAD_L+AW-20}" y="{PAD_T+58}" text-anchor="middle" font-size="8.5" fill="#666">high myelin</text>'
)

svg2 = f"""<svg viewBox="0 0 {SW} {SH}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{SW//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    OPC Marker vs Vascular Marker Co-Distribution
  </text>
  <text x="{SW//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Each dot = one brain region · dot size = myelin gene expression level
  </text>
  <text x="{SW//2}" y="53" text-anchor="middle" font-size="10" fill="#555">
    OPC index = mean(Pdgfra, Cspg4) · Vascular index = mean(Pecam1, Cdh5, Cldn5)
  </text>
  <line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+AH}" stroke="#aaa" stroke-width="1.5"/>
  <line x1="{PAD_L}" y1="{PAD_T+AH}" x2="{PAD_L+AW}" y2="{PAD_T+AH}" stroke="#aaa" stroke-width="1.5"/>
  {yticks}{xticks}{dots}{size_legend}
  <text x="{PAD_L+AW/2:.0f}" y="{SH-12}" text-anchor="middle" font-size="10.5" fill="#444">
    Vascular index (mean Pecam1 + Cdh5 + Cldn5 expression energy)
  </text>
  <text transform="rotate(-90)" x="-{PAD_T+AH//2}" y="16" text-anchor="middle" font-size="10.5" fill="#444">
    OPC index (mean Pdgfra + Cspg4)
  </text>
</svg>"""

with open(OUT / "opc_vascular.svg", "w") as f:
    f.write(svg2)
print("Wrote opc_vascular.svg")
