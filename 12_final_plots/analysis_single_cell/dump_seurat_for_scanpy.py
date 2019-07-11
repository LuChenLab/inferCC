#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.07.11
"""
import os
import sys
from glob import glob
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool

RSCRIPT = """
library(Seurat)

dump_seurat_for_scanpy <- function(object, path) {
    
    dir.create(path, showWarnings = F)
    
    gz = gzfile(paste(path, "raw.csv.gz", sep = "/"), "w+")
    write.csv(as.matrix(object@raw.data), file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "scale.csv.gz", sep = "/"), "w+")
    write.csv(object@scale.data, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "meta.csv.gz", sep = "/"), "w+")
    write.csv(object@meta.data, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "pca.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$pca@cell.embeddings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "pca_gene.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$pca@gene.loadings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "harmony.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$harmony@cell.embeddings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "harnomy_gene.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$harmony@gene.loadings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "tsne.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$tsne@cell.embeddings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "umap.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$umap@cell.embeddings, file = gz)
    close(gz)
}
"""


def call(path):
    u"""
    :param path: path to seurat obj
    """
    scripts = RSCRIPT + "\nobj <- readRDS(\"{0}\")\nsetwd(\"{1}\")\ndump_seurat_for_scanpy(obj, \"scanpy\")\n".format(
        os.path.abspath(path),
        os.path.dirname(os.path.abspath(path))
    )

    temp = path + "_temp.R"
    
    with open(temp, "w+") as w:
        w.write(scripts)
        
    check_call("Rscript " + temp, shell = True)
    
    os.remove(temp)


def main(input_dir, n_jobs=20):
    u"""
    main
    """
    input_dir = os.path.abspath(input_dir)
    
    files = glob(os.path.join(input_dir, "*/*/seurat.rds"))
    
    with Pool(n_jobs) as p:
        p.map(call, files)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
