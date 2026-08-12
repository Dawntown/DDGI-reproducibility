library(MuDataSeurat)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

setwd("raw")

resting_cell_obj = LoadH5Seurat("HuTcellsCRISPRaPerturbSeq_Resting.h5Seurat")
stimulated_cell_obj = LoadH5Seurat("HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat")


Assays(resting_cell_obj)


DefaultAssay(resting_cell_obj) <- "RNA"
DefaultAssay(stimulated_cell_obj) <- "RNA"

WriteH5AD(resting_cell_obj, "HuTcellsCRISPRaPerturbSeq_Resting.h5ad", assay = "RNA")
WriteH5AD(stimulated_cell_obj, "HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5ad", assay = "RNA")

