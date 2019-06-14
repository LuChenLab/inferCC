# scripts to random select cells using ENN from UBL package


args = commandArgs(trailingOnly = T)

rds = args[1]

rds = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Basal/Basal_cells.rds"

load_package <- function() {
    library("UBL")
    library("Seurat")
}

suppressPackageStartupMessages(load_package())

root_dir = paste(dirname(rds), "module", sep = "/")
dir.create(root_dir, showWarnings = F)
setwd(root_dir)

if (!file.exists("../paga/normalized_counts.csv.gz")) {
    stop("Cannot find paga/normalized_counts.csv.gz")
}


# load rds
obj <- readRDS(rds)

# read normalized counts
r = gzfile("../paga/normalized_counts.csv.gz")
normalized_counts = read.csv(r, row.names = 1)


