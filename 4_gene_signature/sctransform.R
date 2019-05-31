# script to make sctransformed gene counts

seurat3_path = "/mnt/raid62/Lung_cancer_10x/scripts/seurat3"

# install.packages("Seurat", lib = seurat3_path)

library(Seurat, lib.loc = seurat3_path)
library(sctransform)

args = commandArgs(trailingOnly = TRUE)

# args = c(getwd())

r = gzfile(paste(args[1], "raw.csv.gz", sep = "/"))
counts = read.csv(r, row.names = 1)


r = gzfile(paste(args[1], "meta.csv.gz", sep = "/"))
meta = read.csv(r, row.names = 1)


obj <- CreateSeuratObject(counts = counts, meta.data = meta)
obj <- SCTransform(object = obj, verbose = TRUE)

data = as.matrix(GetAssayData(obj, "data"))
colnames(data) = gsub("\\.", "-", colnames(data), perl=F)
colnames(data) = gsub("^X", "", colnames(data), perl=T)

w = gzfile(paste(args[1], "normalized_counts.csv.gz", sep = "/"), "w+")
write.csv(counts, w)
close(w)

