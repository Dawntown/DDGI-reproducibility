library(MuDataSeurat)

setwd("raw")

seurat_HEK293FT_obj = readRDS("GSE269596_CPA_HEK293FT.Rds")
seurat_K562_obj = readRDS("GSE269596_CPA_K562.Rds")



DefaultAssay(seurat_HEK293FT_obj) <- "RNA"
seurat_HEK293FT_obj = JoinLayers(seurat_HEK293FT_obj)
DefaultAssay(seurat_K562_obj) <- "RNA"
seurat_K562_obj = JoinLayers(seurat_K562_obj)

WriteH5AD(seurat_HEK293FT_obj, "GSE269596_CPA_HEK293FT.h5ad", assay = "RNA")
WriteH5AD(seurat_K562_obj, "GSE269596_CPA_K562.h5ad", assay = "RNA")

