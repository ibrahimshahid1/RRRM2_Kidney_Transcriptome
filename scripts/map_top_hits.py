#!/usr/bin/env python3
"""Map top Ensembl IDs from Feb 19 regression results to gene symbols,
and check DCT/NCC-WNK gene status in results."""

import json
import csv
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN = os.path.join(BASE, "data/results/run_20260219_115338_2500g")

# Load ID map
with open(os.path.join(BASE, "data/processed/resources/id_map.json")) as f:
    m = json.load(f)

print(f"Loaded {len(m)} gene mappings\n")

# TASK 1: Map top regression hits
groups = {
    "FLIGHT MAIN EFFECT": [
        "ENSMUSG00000025735","ENSMUSG00000087881","ENSMUSG00002076852",
        "ENSMUSG00000094932","ENSMUSG00000074527","ENSMUSG00000064604",
        "ENSMUSG00000015962","ENSMUSG00000092746","ENSMUSG00000119662",
        "ENSMUSG00000073125","ENSMUSG00000019768","ENSMUSG00000057836",
        "ENSMUSG00000118839","ENSMUSG00000065208","ENSMUSG00000054555",
        "ENSMUSG00000077611","ENSMUSG00000067825","ENSMUSG00000080538",
        "ENSMUSG00000107383","ENSMUSG00000033763",
    ],
    "AGE x FLIGHT": [
        "ENSMUSG00000078896","ENSMUSG00000065208","ENSMUSG00000087881",
        "ENSMUSG00000082810","ENSMUSG00000025479","ENSMUSG00000069307",
        "ENSMUSG00000076617","ENSMUSG00000073640","ENSMUSG00000084017",
        "ENSMUSG00002076852","ENSMUSG00000034919","ENSMUSG00000065105",
        "ENSMUSG00000006574","ENSMUSG00000019232","ENSMUSG00000080538",
        "ENSMUSG00000017724","ENSMUSG00000010205","ENSMUSG00000034855",
        "ENSMUSG00000039037","ENSMUSG00000064427",
    ],
    "ARM x FLIGHT": [
        "ENSMUSG00000006154","ENSMUSG00000084391","ENSMUSG00000044468",
        "ENSMUSG00000028553","ENSMUSG00000019894","ENSMUSG00000030283",
        "ENSMUSG00000087881","ENSMUSG00000092746","ENSMUSG00000113441",
        "ENSMUSG00000081059","ENSMUSG00000029804","ENSMUSG00000064604",
        "ENSMUSG00000028989","ENSMUSG00000095675","ENSMUSG00000036814",
        "ENSMUSG00000013611","ENSMUSG00000039542","ENSMUSG00000096596",
        "ENSMUSG00000021207","ENSMUSG00000015962",
    ],
    "3-WAY (Age x Arm x Flight)": [
        "ENSMUSG00000031380","ENSMUSG00000081059","ENSMUSG00000083012",
        "ENSMUSG00000021135","ENSMUSG00000068037","ENSMUSG00000084017",
        "ENSMUSG00000097827","ENSMUSG00000044645","ENSMUSG00000036278",
        "ENSMUSG00000028970","ENSMUSG00000028553","ENSMUSG00000094932",
        "ENSMUSG00000026676","ENSMUSG00000082810","ENSMUSG00000100426",
        "ENSMUSG00000102204","ENSMUSG00000084391","ENSMUSG00000026205",
        "ENSMUSG00002076304","ENSMUSG00000022304",
    ],
}

for label, ids in groups.items():
    print(f"=== {label} (top 20 q_BH) ===")
    for i, g in enumerate(ids, 1):
        sym = m.get(g, "???")
        print(f"  {i:2d}. {g} -> {sym}")
    print()

# TASK 2: Find DCT/NCC-WNK genes in the map
dct_symbols = [
    "Wnk1","Wnk4","Stk39","Oxsr1","Slc12a3","Kcnj10","Kcnj16",
    "Scnn1a","Scnn1b","Scnn1g","Nedd4l","Sgk1","Klhl3","Cul3",
]
# Also add some related markers
extra_dct = ["Trpm6","Calb1","Pvalb","Trpv5","Slc8a1","Atp2b1","Atp2b4"]

# Reverse map: symbol -> ensembl
rev = {}
for eid, sym in m.items():
    rev[sym] = eid

print("=== DCT/NCC-WNK PATHWAY GENES ===")
dct_ensembl = {}
for sym in dct_symbols:
    eid = rev.get(sym, "NOT IN MAP")
    dct_ensembl[sym] = eid
    print(f"  {sym:12s} -> {eid}")

print("\n=== ADDITIONAL DCT MARKERS ===")
for sym in extra_dct:
    eid = rev.get(sym, "NOT IN MAP")
    dct_ensembl[sym] = eid
    print(f"  {sym:12s} -> {eid}")

# TASK 3: Check DCT genes in regression results
reg_files = [
    "gene_flight_effect.tsv",
    "gene_age_flight_interaction.tsv",
    "gene_arm_flight_interaction.tsv",
    "gene_age_arm_flight_3way.tsv",
]

# Collect all valid DCT Ensembl IDs
dct_ids = {v: k for k, v in dct_ensembl.items() if v != "NOT IN MAP"}

print("\n=== DCT/NCC-WNK GENES IN REGRESSION RESULTS ===")
for fname in reg_files:
    fpath = os.path.join(RUN, "phase6_regression", fname)
    print(f"\n--- {fname} ---")
    found = False
    with open(fpath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["gene"] in dct_ids:
                sym = dct_ids[row["gene"]]
                found = True
                print(f"  {sym:12s} | edges={row['n_edges']:>4s} | beta={row['median_beta']:>8s} | t={row['median_t']:>8s} | p={row['p_stouffer']:>12s} | q={row['q_BH']:>12s}")
    if not found:
        print("  (none found)")

# TASK 4: Check DCT genes in rewiring results
rew_files = [
    "ISS_T_OLD_FLT_minus_GND_rewiring_agg.tsv",
    "ISS_T_YNG_FLT_minus_GND_rewiring_agg.tsv",
    "LAR_OLD_FLT_minus_GND_rewiring_agg.tsv",
    "LAR_YNG_FLT_minus_GND_rewiring_agg.tsv",
]

print("\n=== DCT/NCC-WNK GENES IN REWIRING RESULTS ===")
for fname in rew_files:
    fpath = os.path.join(RUN, "phase3_rewiring", fname)
    print(f"\n--- {fname.replace('_rewiring_agg.tsv','')} ---")
    found = False
    with open(fpath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        all_rows = list(reader)
        total = len(all_rows)
        for i, row in enumerate(all_rows):
            if row["gene"] in dct_ids:
                sym = dct_ids[row["gene"]]
                rank = i + 1
                pctile = 100.0 * (1 - rank / total)
                found = True
                print(f"  {sym:12s} | rewiring={row['rewiring_mean']:>8s} | std={row['rewiring_std']:>8s} | rank={rank:>5d}/{total} ({pctile:.1f}%ile)")
    if not found:
        print("  (none found)")

# TASK 5: Check permutation pvals for DCT genes
perm_files = [
    "ISS_T_OLD_FLT_minus_GND_perm_pvals.tsv",
    "ISS_T_YNG_FLT_minus_GND_perm_pvals.tsv",
    "LAR_OLD_FLT_minus_GND_perm_pvals.tsv",
    "LAR_YNG_FLT_minus_GND_perm_pvals.tsv",
]

print("\n=== DCT/NCC-WNK GENES IN PERMUTATION RESULTS ===")
for fname in perm_files:
    fpath = os.path.join(RUN, "phase6_uncertainty", fname)
    print(f"\n--- {fname.replace('_perm_pvals.tsv','')} ---")
    found = False
    with open(fpath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["gene"] in dct_ids:
                sym = dct_ids[row["gene"]]
                found = True
                print(f"  {sym:12s} | obs={row['rewiring_abs_obs']:>10s} | p_perm={row['p_perm']:>8s} | q_BH={row['q_BH']:>8s}")
    if not found:
        print("  (none found)")

# TASK 6: Summary of enrichment gene set overlap
print("\n=== DCT_NCC_WNK_axis GENE SET IN ENRICHMENT ===")
print("Note: In ALL enrichment files, DCT_NCC_WNK_axis shows:")
print("  - overlap = 0 genes (in gene_flight_effect, age_flight, 3-way regression enrichment)")
print("  - hits = 0 in ALL rewiring-based grounding files")
print("  - This means NO DCT/NCC-WNK genes overlapped with significant gene lists")
