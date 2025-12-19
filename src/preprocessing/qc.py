# src/preprocessing/qc.py
import pandas as pd
import numpy as np
from pybiomart import Server
import time

def filter_protein_coding(counts_df, retries=5):
    """
    ADAPTED FROM OLD SCRIPT: 
    Filters the expression matrix to keep only protein-coding genes using Ensembl BioMart.
    """
    print("Querying Ensembl BioMart for protein-coding genes...")
    
    # Try connecting to the server with retries (as in your old script)
    server = None
    for i in range(retries):
        try:
            server = Server(host='http://www.ensembl.org', use_cache=False)
            dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['mmusculus_gene_ensembl']
            
            # Query for gene IDs and types
            # Note: We assume your index is Ensembl IDs. If it's Symbols, change 'ensembl_gene_id' to 'external_gene_name'
            gene_info = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype'])
            
            # Filter for protein_coding
            coding_genes = gene_info[gene_info['Gene type'] == 'protein_coding']['Gene stable ID']
            
            # Filter the dataframe
            # Assuming counts_df index matches 'Gene stable ID'
            overlap = counts_df.index.intersection(coding_genes)
            filtered_df = counts_df.loc[overlap]
            
            print(f"BioMart Filter: Reduced from {counts_df.shape[0]} to {filtered_df.shape[0]} protein-coding genes.")
            return filtered_df
            
        except Exception as e:
            print(f"BioMart query failed (Attempt {i+1}/{retries}). Error: {e}")
            time.sleep(5)
            
    print(" WARNING: BioMart filtering failed after all retries. Returning original data.")
    return counts_df

# ... (Keep your existing detect_outliers function here) ...