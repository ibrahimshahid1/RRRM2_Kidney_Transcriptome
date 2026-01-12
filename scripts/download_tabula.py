import cellxgene_census
import tiledbsoma as soma

with cellxgene_census.open_soma() as census:
    # Pick organism and datasetp
    exp = census["rna"]["mouse"].experiment("Tabula Muris Senis")

    # Read into an AnnData
    adata = cellxgene_census.get_anndata(
        census,
        organism="mouse",
        census_version="latest",
        var_value_filter="gene_id == gene_id",
    )

print(adata)
adata.write_h5ad("tabula_muris_senis_kidney.h5ad")
