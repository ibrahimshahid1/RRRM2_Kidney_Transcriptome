import numpy as np
import pandas as pd
import os
import networkx as nx
from node2vec import Node2Vec  # pip install node2vec
from scipy.linalg import orthogonal_procrustes

# ==========================================
# CONFIG
# ==========================================
BASE_DIR = "data"
# Input 1: The Group-Mean Networks you generated
NETWORK_PATH = os.path.join(BASE_DIR, "processed/networks/phase1/group_mean_edges.npz")
# Input 2: The Anchor Files you just created
ANCHOR_DIR = os.path.join(BASE_DIR, "processed/networks/phase1")
# Input 3: The Gene Universe
GENE_LIST = os.path.join(BASE_DIR, "processed/networks/phase1/phase1_genes.txt")

OUT_DIR = os.path.join(BASE_DIR, "results/embeddings")
os.makedirs(OUT_DIR, exist_ok=True)

# Node2Vec Hyperparameters (PDF Section 4.5)
DIMENSIONS = 128
WALK_LENGTH = 80
NUM_WALKS = 20   # Use 20 for quick test, 200 for final paper
P = 0.25         # Return parameter (local structure)
Q = 4.0          # In-out parameter (global structure)
WORKERS = 4 

# ==========================================
# HELPER FUNCTIONS
# ==========================================

def load_data():
    print("Loading gene list...")
    with open(GENE_LIST, 'r') as f:
        genes = [line.strip() for line in f]
    
    print("Loading network bundle...")
    data = np.load(NETWORK_PATH)
    # The .npz contains the edges. We need to reconstruct the matrix or graph.
    # Assuming the previous script saved correlation matrices or similar.
    # Let's inspect the keys.
    networks = {k: data[k] for k in data.files}
    return genes, networks

def embed_network(adj_matrix, gene_names):
    """
    Runs node2vec on a weighted adjacency matrix.
    """
    # 1. Prepare Graph
    # Take absolute values (we care about connection strength, not sign for embedding)
    w_matrix = np.abs(adj_matrix)
    np.fill_diagonal(w_matrix, 0)
    
    # Thresholding for speed: Keep only strong edges (e.g., top 10% or > 0.1)
    # Node2vec is slow on fully connected graphs.
    threshold = np.percentile(w_matrix, 90) # Top 10% edges
    w_matrix[w_matrix < threshold] = 0
    
    G = nx.from_numpy_array(w_matrix)
    
    # 2. Run Node2Vec
    n2v = Node2Vec(G, dimensions=DIMENSIONS, walk_length=WALK_LENGTH, 
                   num_walks=NUM_WALKS, p=P, q=Q, workers=WORKERS, quiet=True)
    model = n2v.fit(window=10, min_count=1, batch_words=4)
    
    # 3. Extract Embeddings in correct order
    # Node2Vec might shuffle node IDs, so we map back to gene_names indices
    emb_ordered = np.zeros((len(gene_names), DIMENSIONS))
    for i in range(len(gene_names)):
        if str(i) in model.wv:
            emb_ordered[i, :] = model.wv[str(i)]
    
    return emb_ordered

def align_embeddings(base_emb, target_emb, anchors, gene_names):
    """
    Procrustes Alignment: Rotates 'target' to match 'base' using anchors.
    """
    # Find indices of anchor genes
    anchor_indices = [i for i, g in enumerate(gene_names) if g in anchors]
    
    if len(anchor_indices) < 10:
        print(f"WARNING: Only {len(anchor_indices)} anchors found. Alignment may be unstable.")
        return target_emb
        
    # Extract anchor subsets
    A = base_emb[anchor_indices]
    B = target_emb[anchor_indices]
    
    # Compute Rotation Matrix R
    # Minimize || A - B*R ||
    R, scale = orthogonal_procrustes(B, A)
    
    # Apply Rotation to ALL genes in target
    target_aligned = np.dot(target_emb, R)
    return target_aligned

# ==========================================
# EXECUTION
# ==========================================

genes, networks = load_data()
print(f"Loaded {len(networks)} networks (groups) for {len(genes)} genes.")

# Define the Spaceflight Comparisons (Aim 2)
comparisons = [
    # Reference (Ground)          Target (Flight)
    ("YNG|ISS-T|GC",            "YNG|ISS-T|FLT"),
    ("OLD|ISS-T|GC",            "OLD|ISS-T|FLT"),
    # Controls
    ("YNG|LAR|GC",              "YNG|LAR|FLT"),
]

results = []

for ref_name, target_name in comparisons:
    print(f"\n--- Processing: {ref_name} vs {target_name} ---")
    
    if ref_name not in networks or target_name not in networks:
        print("Skipping: Network not found in .npz bundle")
        continue

    # 1. Embed Reference
    print(f"Embedding Reference ({ref_name})...")
    emb_ref = embed_network(networks[ref_name], genes)
    
    # 2. Embed Target
    print(f"Embedding Target ({target_name})...")
    emb_target = embed_network(networks[target_name], genes)
    
    # 3. Load Anchors
    # Filename format: anchors_YNG_ISS-T_GC_vs_YNG_ISS-T_FLT_k150.tsv
    safe_ref = ref_name.replace("|", "_")
    safe_tar = target_name.replace("|", "_")
    anchor_file = f"anchors_{safe_ref}_vs_{safe_tar}_k150.tsv"
    anchor_path = os.path.join(ANCHOR_DIR, anchor_file)
    
    if os.path.exists(anchor_path):
        print(f"Loading anchors from: {anchor_file}")
        anchors_df = pd.read_csv(anchor_path, sep='\t')
        anchor_genes = anchors_df['gene'].tolist()
    else:
        print(f"ERROR: Anchor file missing: {anchor_path}")
        continue

    # 4. Align
    print("Aligning flight network to ground network...")
    emb_target_aligned = align_embeddings(emb_ref, emb_target, anchor_genes, genes)
    
    # 5. Calculate Rewiring (Cosine Distance)
    # Cosine Sim = dot(u, v) / (|u| |v|)
    # Normalize vectors first
    norm_ref = emb_ref / np.linalg.norm(emb_ref, axis=1, keepdims=True)
    norm_tar = emb_target_aligned / np.linalg.norm(emb_target_aligned, axis=1, keepdims=True)
    
    cosine_sim = np.sum(norm_ref * norm_tar, axis=1)
    rewiring = 1 - cosine_sim
    
    # 6. Save top hits
    df_res = pd.DataFrame({
        'gene': genes,
        'rewiring_score': rewiring,
        'cosine_sim': cosine_sim,
        'comparison': f"{target_name}_vs_{ref_name}"
    })
    
    # Filter out anchors (they shouldn't be top hits)
    df_res = df_res[~df_res['gene'].isin(anchor_genes)]
    
    # Sort
    df_res = df_res.sort_values('rewiring_score', ascending=False)
    
    out_file = os.path.join(OUT_DIR, f"rewiring_{safe_tar}_vs_{safe_ref}.csv")
    df_res.to_csv(out_file, index=False)
    print(f"Saved results to: {out_file}")
    
    print("Top 5 Rewired Genes:")
    print(df_res.head(5)[['gene', 'rewiring_score']])

print("\nDONE. The top genes in the CSVs are your 'Silent Shifters'.")