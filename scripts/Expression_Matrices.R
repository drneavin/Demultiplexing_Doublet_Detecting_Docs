library(Seurat)
library(tidyverse)
library(data.table)


dir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/"
outdir <- file.path(dir, "Data4externalCollaborators/2021_09_Kleeman/output/")
dir.create(outdir)

### Read in Data ###
fibro_seurat <- readRDS(file.path(dir,"/Fibroblast_eQTLs/output/ClusterCharacterization/2019-09-13/FibroblastSeuratObjectCombat.rds"))
iPSC_seurat
