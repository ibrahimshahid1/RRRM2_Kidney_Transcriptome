import pandas as pd
import os

# CONFIG
GENE_LIST = "data/processed/networks/phase2/phase2_genes.txt"

# The "Distal Nephron / NCC-WNK" Hit List
pathway_genes = {
    "ENSMUSG00000031766": "Slc12a3 (NCC)",  # The Transporter
    "ENSMUSG00000035112": "Wnk4",           # The Regulator
    "ENSMUSG00000045962": "Wnk1",           # The Isoform Switch (KS-WNK1 vs L-WNK1)
    "ENSMUSG00000019970": "Sgk1",           # The Master Kinase (Aldosterone/Stress)
    "ENSMUSG00000031431": "Tsc22d3 (GILZ)", # The Stress Mediator
    "ENSMUSG00000027030": "Stk39 (SPAK)",   # The Direct Phosphorylator of NCC
    "ENSMUSG00000036737": "Oxsr1 (OSR1)",   # The Backup Phosphorylator
    "ENSMUSG00000014164": "Klhl3",          # The Ubiquitin Ligase (Degrades WNK4)
    "ENSMUSG00000001786": "Fbxo7",          # WNK Ubiquitination (SCF complex component)
    "ENSMUSG00000020893": "Per1",           # Circadian Clock (Regulates NCC transcription)
}

# EXECUTION
if not os.path.exists(GENE_LIST):
    print(f"Error: Could not find {GENE_LIST}")
    exit()

print(f"Checking {len(pathway_genes)} pathway genes in your Network Universe...")

# Load your Phase 2 Gene List
with open(GENE_LIST, 'r') as f:
    network_genes = set(line.strip() for line in f)

print(f"\n{'Gene Symbol':<15} | {'ID':<20} | {'Status'}")
print("-" * 50)

hits = 0
for ens_id, symbol in pathway_genes.items():
    if ens_id in network_genes:
        status = "[+] IN NETWORK"
        hits += 1
    else:
        status = "[-] Not Selected (Low Variance)"
    print(f"{symbol:<15} | {ens_id:<20} | {status}")

print("-" * 50)
print(f"Total Pathway Hits: {hits}/{len(pathway_genes)}")

if "ENSMUSG00000002250" in network_genes:
    print("\n[CRITICAL SUCCESS] Sgk1 is in the network.")
    print("Since Sgk1 is the upstream regulator, its presence allows you to model")
    print("the *regulatory pressure* on NCC, even if NCC itself is stable.")