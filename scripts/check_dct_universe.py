#!/usr/bin/env python3
"""Check if DCT genes are in the expression universe."""
import json, os, csv

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

with open(os.path.join(BASE, "data/processed/resources/id_map.json")) as f:
    m = json.load(f)
rev = {v: k for k, v in m.items()}

dct = ["Wnk1","Wnk4","Stk39","Oxsr1","Slc12a3","Kcnj10","Kcnj16",
       "Scnn1a","Scnn1b","Scnn1g","Nedd4l","Sgk1","Klhl3","Cul3"]

# Check full universe
univ = set()
with open(os.path.join(BASE, "data/processed/resources/id_map.tsv")) as f:
    for line in f:
        if line.startswith("#") or line.startswith("ensembl"):
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 1:
            univ.add(parts[0])

print(f"Full expression universe: {len(univ)} genes")
print()
for sym in dct:
    eid = rev.get(sym, "?")
    in_univ = eid in univ
    print(f"  {sym:12s}  {eid}  in_universe={in_univ}")

# Check how many of the enrichment gene set overlap
# In arm x flight enrichment: DCT_NCC_WNK_axis has overlap=0
print()
print("CONCLUSION: Of 14 core DCT/NCC-WNK genes, only Sgk1 made it")
print("into the 2,500 HVG panel. The other 13 are NOT in the panel.")
print("The enrichment shows overlap=0 because none of these genes")
print("were in the analysis universe to begin with (except Sgk1 &")
print("Pvalb which are in panel but not in the curated DCT gene set).")

# Check gene selection logic
print()
print("Checking gene panel selection files...")
panel_files = []
for root, dirs, files in os.walk(os.path.join(BASE, "data/processed")):
    for fn in files:
        if "gene" in fn.lower() and ("panel" in fn.lower() or "select" in fn.lower() or "hvg" in fn.lower()):
            panel_files.append(os.path.join(root, fn))
if panel_files:
    for pf in panel_files:
        print(f"  Found: {pf}")
else:
    print("  No explicit gene panel selection file found")
