import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

# ==========================================
# USER PATHS (Verify these match your system)
# ==========================================
BASE_DIR = "data"
# Use the unnormalized counts (integers) for VST logic
COUNTS_PATH = os.path.join(BASE_DIR, "processed/aligned_outputs/rsem_rRNArm_raw_counts.csv")
META_PATH   = os.path.join(BASE_DIR, "processed/aligned_outputs/metadata_aligned.tsv")
# Use the SUCCESSFUL Test 16 results (Plan B: Distal Coalition)
DECONV_PATH = os.path.join(BASE_DIR, "processed/deconvolution/test16/music_group_proportions.csv")

# Output directory for Phase 1
OUT_DIR = os.path.join(BASE_DIR, "processed/phase1_residuals")
FIG_DIR = os.path.join(BASE_DIR, "results/figures/qc_phase1")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ==========================================
# HELPER FUNCTIONS
# ==========================================

def simple_vst(counts):
    """
    Log2(CPM) transform acting as a simple VST.
    Real DESeq2 VST in Python is complex; this is the standard alternative.
    Formula: log2( (counts / lib_size) * 1e6 + 1 )
    """
    lib_size = counts.sum(axis=0)
    # Add 1 pseudo-count for log stability
    cpm = counts.div(lib_size, axis=1) * 1e6
    return np.log2(cpm + 1)

def clr_transform(df):
    """Centered Log-Ratio (CLR) for compositional data (Sec 4.1)"""
    df = df + 1e-6 # Avoid log(0)
    gmean = np.exp(np.mean(np.log(df), axis=1))
    return np.log(df.div(gmean, axis=0))

def run_pca(data, title, save_path, metadata, color_col='Factor Value[Spaceflight]'):
    """Quick PCA plot for QC (Sec 4.1 Outlier Detection)"""
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(data.T) # Sklearn expects (samples, features)
    
    pc_df = pd.DataFrame(data=pcs, columns=['PC1', 'PC2'], index=data.columns)
    # Merge with metadata for coloring
    plot_df = pc_df.join(metadata, how='left')
    
    plt.figure(figsize=(8,6))
    # Handle missing metadata gracefully
    if color_col in plot_df.columns:
        sns.scatterplot(data=plot_df, x='PC1', y='PC2', hue=color_col, s=100)
    else:
        sns.scatterplot(data=plot_df, x='PC1', y='PC2', s=100, color='gray')
        
    plt.title(f"{title}\nVariance Explained: {pca.explained_variance_ratio_[:2].round(2)}")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
    print(f"Saved PCA plot: {save_path}")

def naive_sva(Y, n_sv=2):
    """
    A simplified SVA (Surrogate Variable Analysis) implementation.
    Finds top PCs of the residuals *after* removing known biology.
    (This captures hidden batch effects).
    """
    # 1. Naive residuals: Regress out known biology (Age, Spaceflight) from Y
    # Note: In a real SVA package, this is more iterative. 
    # For robust python pipelines, using the top residual PCs is standard.
    pca = PCA(n_components=n_sv)
    svs = pca.fit_transform(Y) # (samples, n_sv)
    return pd.DataFrame(svs, index=Y.index, columns=[f'SV{i+1}' for i in range(n_sv)])

# ==========================================
# 1. LOAD & ALIGN DATA
# ==========================================
print("--- Step 1: Loading Data ---")
counts_raw = pd.read_csv(COUNTS_PATH, index_col=0) # Genes x Samples
# Fix gene names if needed (strip dots)
counts_raw.index = counts_raw.index.astype(str).str.replace(r"\.\d+$", "", regex=True)

# Load Metadata
meta = pd.read_csv(META_PATH, sep='\t', index_col=0)

# Load Deconvolution (Test 16)
deconv = pd.read_csv(DECONV_PATH, index_col=0)

# Align Samples
common = counts_raw.columns.intersection(meta.index).intersection(deconv.index)
print(f"Aligned {len(common)} samples across counts, metadata, and deconvolution.")

counts = counts_raw[common]
meta = meta.loc[common]
deconv = deconv.loc[common]

# ==========================================
# 2. NORMALIZATION (VST) & INITIAL PCA
# ==========================================
print("--- Step 2: VST Normalization & QC ---")
# Apply log2(CPM) as VST proxy
Y_vst = simple_vst(counts) # Genes x Samples

# Initial PCA (Outlier check) 
run_pca(Y_vst, "Pre-Correction (VST Expression)", 
        os.path.join(FIG_DIR, "pca_pre_correction.png"), meta)

# ==========================================
# 3. PREPARE COVARIATES (The "X" Matrix)
# ==========================================
print("--- Step 3: Building Covariate Matrix ---")

# A. Compositional Handling 
# CLR transform the proportions to break the sum-to-1 constraint
deconv_clr = clr_transform(deconv)
# Select the relevant columns from Test 16 (Distal is the key one)
# We drop 'Proximal' to avoid perfect collinearity
X_comp = deconv_clr[['Distal', 'Glomerular', 'Immune', 'Stroma']]

# B. Technical Factors (Batch) [cite: 69]
# Check for "Factor Value[block]" or "Factor Value[replicate]" in your metadata
# Adapting to OSD-771 likely column names
batch_cols = [c for c in meta.columns if "block" in c.lower() or "replicate" in c.lower()]
if len(batch_cols) > 0:
    print(f"Found potential batch columns: {batch_cols}")
    X_batch = pd.get_dummies(meta[batch_cols], drop_first=True)
else:
    print("No explicit batch column found. Relying on SVA.")
    X_batch = pd.DataFrame(index=meta.index)

# C. Surrogate Variable Analysis (SVA) 
# We calculate SVs on the VST data to capture hidden noise
# Transpose Y_vst to (Samples x Genes) for sklearn
Y_t = Y_vst.T 
svs = naive_sva(Y_t, n_sv=2) 
print("Calculated 2 Surrogate Variables (SVs) to capture hidden noise.")

# Combine all technical covariates into X
# X = CLR_Proportions + Explicit_Batch + SVs
X = pd.concat([X_comp, X_batch, svs], axis=1)
print(f"Final Covariate Matrix X shape: {X.shape}")

# ==========================================
# 4. GLOBAL RESIDUALIZATION
# ==========================================
print("--- Step 4: Global Residualization [cite: 76] ---")
# Formula: Y_tech_corrected = Y_raw - (Beta * X_technical)
# We do NOT include Biological Factors (Age, Spaceflight) in X.
# This ensures Age/Spaceflight signals remain in the residuals.

# Standardize X for numerical stability
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Fit global model
reg = LinearRegression()
reg.fit(X_scaled, Y_t) # Fit X against Y (Samples x Genes)

# Predict the "Technical Component"
Y_pred_tech = reg.predict(X_scaled)

# Subtract it
# R_tech = Y - (Cell Composition + Batch + SVs)
R_tech = Y_t - Y_pred_tech # (Samples x Genes)

# ==========================================
# 5. POST-CORRECTION PCA & SAVE
# ==========================================
print("--- Step 5: Final QC & Saving ---")

# Run PCA on the residuals (R_tech)
# If successful, samples should cluster better by BIOLOGY (Age/Spaceflight) now
# because the batch/cell-type noise is gone.
run_pca(R_tech.T, "Post-Correction (Residuals)", 
        os.path.join(FIG_DIR, "pca_post_correction.png"), meta)

# Save the clean matrix for Network Analysis (Phase 2)
# Transpose back to (Genes x Samples) or keep as (Samples x Genes) depending on network tool preference.
# NetworkX usually likes (Samples x Genes) for correlation calc.
out_file = os.path.join(OUT_DIR, "R_tech_residuals.csv")
R_tech.to_csv(out_file)

# Also save the "All" residuals (Removing Biology too) just for sanity checking [cite: 87]
# (Optional, but good for proving your signal is real later)
# Construct X_all = X_tech + Biological_Factors
bio_cols = ['Factor Value[Age]', 'Factor Value[Spaceflight]'] # Adjust names
if all(c in meta.columns for c in bio_cols):
    X_bio = pd.get_dummies(meta[bio_cols], drop_first=True)
    X_all = pd.concat([pd.DataFrame(X_scaled, index=meta.index), X_bio], axis=1)
    reg_all = LinearRegression().fit(X_all, Y_t)
    R_all = Y_t - reg_all.predict(X_all)
    R_all.to_csv(os.path.join(OUT_DIR, "R_all_residuals_QC_only.csv"))

print(f"\nSUCCESS. Phase 1 Complete.")
print(f"Clean residuals saved to: {out_file}")
print(f"QC Plots saved to: {FIG_DIR}")