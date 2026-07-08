"""
Pull Allen Mouse Brain Atlas ISH expression for three gene groups:
  1. Wnt pathway genes  — known inhibitors of OPC differentiation
  2. OPC / OL lineage   — progenitor markers, mature myelin markers
  3. Vascular markers   — endothelial + pericyte, for OPC-vessel co-distribution

Structures chosen to span gray/white matter gradients:
  Corpus callosum (CC) 776 — pure white matter
  Isocortex 315            — cortical gray matter + myelin layers
  Hippocampus 1089
  Cerebellum 512           — heavy myelination
  Striatum 672             — mixed, active remyelination model
  Hypothalamus 1097
  Olfactory bulb 507
  Brainstem 313            — heavily myelinated, frequent MS lesion site
"""
import requests, csv, time
from pathlib import Path

OUT  = Path(__file__).parent
BASE = "http://api.brain-map.org/api/v2/data/query.json"

GENES = {
    "wnt": ["Axin2", "Ctnnb1", "Fzd1", "Tcf7l2", "Wnt3a"],
    "opc_ol": ["Pdgfra", "Cspg4", "Sox10", "Mbp", "Plp1", "Mog", "Cnp"],
    "vascular": ["Pecam1", "Cdh5", "Cldn5", "Pdgfrb"],
}
ALL_GENES = [g for grp in GENES.values() for g in grp]

STRUCTURES = {
    "Corpus_callosum": 776,
    "Isocortex":       315,
    "Hippocampus":     1089,
    "Cerebellum":      512,
    "Striatum":        672,
    "Hypothalamus":    1097,
    "Olfact_bulb":     507,
    "Midbrain":        313,
}

def get_dataset(gene):
    url = (f"{BASE}?criteria=model::SectionDataSet,"
           f"rma::criteria,genes[acronym$eq'{gene}'],"
           f"[failed$eqfalse],"
           f"products[abbreviation$eqMouse]&num_rows=5")
    r = requests.get(url, timeout=20).json()
    if not r.get("msg"):
        return None
    for row in r["msg"]:
        if row.get("plane_of_section_id") == 2:
            return row["id"]
    return r["msg"][0]["id"] if r["msg"] else None

def get_expression(dataset_id, struct_ids):
    id_str = ",".join(str(i) for i in struct_ids)
    url = (f"{BASE}?criteria=model::StructureUnionize,"
           f"rma::criteria,[section_data_set_id$eq{dataset_id}],"
           f"[structure_id$in{id_str}]&num_rows=50")
    r = requests.get(url, timeout=20).json()
    return {row["structure_id"]: row.get("expression_energy", 0.0) or 0.0
            for row in r.get("msg", [])}

struct_ids = list(STRUCTURES.values())
id_to_name = {v: k for k, v in STRUCTURES.items()}
gene_to_grp = {g: grp for grp, gs in GENES.items() for g in gs}

rows = []
for gene in ALL_GENES:
    print(f"  {gene}...", end=" ", flush=True)
    ds = get_dataset(gene)
    if ds is None:
        print("no dataset"); continue
    expr = get_expression(ds, struct_ids)
    row = {"gene": gene, "group": gene_to_grp[gene], "dataset_id": ds}
    for sid, name in id_to_name.items():
        row[name] = round(expr.get(sid, 0.0), 4)
    rows.append(row)
    print(f"ds={ds} ✓")
    time.sleep(0.3)

fieldnames = ["gene","group","dataset_id"] + list(STRUCTURES.keys())
with open(OUT / "expression.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"\nSaved expression.tsv ({len(rows)} genes)")
for r in rows:
    vals = {k: r[k] for k in STRUCTURES}
    top  = max(vals, key=vals.get)
    print(f"  {r['gene']:12s}  peak={top} ({vals[top]:.3f})")
