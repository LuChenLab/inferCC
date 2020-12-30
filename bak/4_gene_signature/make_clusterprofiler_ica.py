#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Perform clusterProfiler to genes from ICA
"""
import os
from glob import glob
import json
from multiprocessing import Pool
from subprocess import  check_call, CalledProcessError


def make_clusterprofiler(args):
    u"""

    :return:
    """
    genes, output = args

    outdir = os.path.dirname(output)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    temp_genes = output + ".gene_tmp"
    with open(temp_genes, "w+") as w:
        w.write("\n".join(genes))

    rscript = """
    genes = "{0}"
    output = "{1}"
    """.format(temp_genes, output)

    rscript += """
    library(clusterProfiler)
    library(openxlsx)
    library(DOSE)
    library(org.Hs.eg.db) 
    
    genes = read.table(genes)
    
    eg <- bitr(unique(genes[, 1]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    
    
    kk <- enrichKEGG(gene     = eg$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)
    kk <- as.data.frame(kk)
    
        
    res = NULL
    for(i in c("BP", "CC", "MF")) {
        ego <- enrichGO(gene      = eg$ENTREZID,
                        keyType       = 'ENTREZID',
                        OrgDb         = org.Hs.eg.db,
                        ont           = i,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)
        
        if(is.null(ego)) {
            return(ego)
        }
        ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
        
        ego = as.data.frame(ego)
        
        if (nrow(ego) > 0) {
            ego$ONTOLOGY = i
        }
        
        if (is.null(res)) {
            res = ego
        } else {
            res = rbind(res, ego)
        }
    }
    
    do <- enrichDO(gene     = eg$ENTREZID,
                   ont           = "DO",
                   pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize     = 5,
                   maxGSSize     = 500,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
    
    do = as.data.frame(do)
    
    wb = createWorkbook()

    addWorksheet(wb, "KEGG")
    writeData(wb, 1, kk)
    
    addWorksheet(wb, "GO")
    writeData(wb, 2, res)
    
    addWorksheet(wb, "DOSE")
    writeData(wb, 3, do)
    
    saveWorkbook(wb, file = output, overwrite = T)
    """

    temp_r = output + ".r"
    with open(temp_r, "w+") as w:
        w.write(rscript)

    try:
        check_call("Rscript {0}".format(temp_r), shell=True)
    except CalledProcessError as err:
        print(err)

    os.remove(temp_r)
    os.remove(temp_genes)


def main(path, output, n_jobs=20):
    files = glob(os.path.join(path, "Monocytes/gene_module_sctransform/RENN_ICA_data.json"))

    output = os.path.abspath(output)
    if not os.path.exists(output):
        os.makedirs(output)

    tasks = []
    for i in files:
        cell_name = os.path.basename(os.path.dirname(os.path.dirname(i)))

        with open(i) as r:
            data = json.load(r)

        for cluster, genes in data.items():
            if len(genes) >= 3:
                temp_output = os.path.join(output, "{0}_{1}.xlsx".format(cell_name, cluster))
                tasks.append([genes, temp_output])

    with Pool(n_jobs) as p:
        p.map(make_clusterprofiler, tasks)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)