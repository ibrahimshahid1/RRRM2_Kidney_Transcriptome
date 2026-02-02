# In R
suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Mm.eg.db)
})

genes <- read.delim("data/results/phase6_uncertainty/ISS_T_YNG_FLT_minus_GC_perm_pvals.tsv")$gene
map <- AnnotationDbi::select(org.Mm.eg.db, keys = genes, keytype = "ENSEMBL", columns = c("SYMBOL"))
map <- map[!duplicated(map$ENSEMBL), ]
write.table(map, "data/results/phase6_uncertainty/ensembl_to_symbol.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE
)
