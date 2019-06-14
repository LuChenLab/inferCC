#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Create at 2019.05.28

Using mfuzz to find out the gene modules of different cells
"""

import os
import glob
import sys

from subprocess import check_call, CalledProcessError
from multiprocessing import Pool

from argparse import ArgumentParser, ArgumentError


def command_line(argv):
    parser = ArgumentParser()

    parser.add_argument(
        "-i",
        help="Path to input directory",
        type=str,
        required=True
    )

    parser.add_argument(
        "-p",
        help="CPU numbers (default: %(default)s)",
        type=int,
        default=1
    )

    if len(argv) <= 1:
        parser.print_help()
        exit(0)
    else:
        try:
            return parser.parse_args(argv[1:])
        except ArgumentError as err:
            print(err)
            exit(err.__hash__())


def _mfuzz_(input_directory: str):
    u"""
    call mfuzz to calculate
    :param input_directory: path to input directory
    :return:
    """
    input_directory = os.path.abspath(input_directory)

    if "Granulocyte_ADC" not in input_directory:
        return

    rds = ""
    for i in glob.glob(os.path.join(input_directory, "*.rds")):
        if "monocle" not in i and "slingshot" not in i and "wgcna" not in i.lower():
            rds = i
            break

    if not rds:
        print("{0} not find any rds".format(input_directory))
        return

    rscript = """
    rds = "{0}"
    root.dir = "{1}"
    output = "mfuzz_gene_module"
    heatmap_output = paste(output, "heatmap", sep = "/")
    """.format(
        rds,
        input_directory
    )

    rscript += """
    setwd(root.dir)
    dir.create(output, FALSE)
    dir.create(heatmap_output, FALSE)
    set.seed(1)
    
    
    library(ggplot2)
    library(Mfuzz)
    library(Seurat)
    library(openxlsx)
    
    options(stringsAsFactors = F)
    
    obj <- readRDS(rds)
    r = gzfile("paga/normalized_counts.csv.gz")
    data = read.csv(r, row.names = 1)
    
    stage_spec = read.xlsx("annotation_results_by_stage.xlsx")
    stage_spec = stage_spec[stage_spec$p_val_adj < 0.05 & stage_spec$avg_logFC > 0, ]
    
    if (file.exists("annotation_results_by_cluster.xlsx")) {
        cluster_spec = read.xlsx("annotation_results_by_cluster.xlsx")
        cluster_spec = cluster_spec[cluster_spec$p_val_adj < 0.05 & cluster_spec$avg_logFC > 0, ]
        
        colnames(data) = gsub("\\\\.", "-", colnames(data), perl=F)
        colnames(data) = gsub("^X", "", colnames(data), perl=T)
        
        
        # plot(density(log10(rowSums(data) + 1)))
        # data = data[log10(rowSums(data) + 1) / ncol(data) > 0.001, ]
        
        # check stage markers
        p <- DotPlot(obj, genes.plot=unique(stage_spec$gene), cols.use = c("lightgrey", "blue"),
          col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
          scale.by = "radius", scale.min = NA, scale.max = NA, group.by="res.0.6",
          plot.legend = FALSE, do.return = TRUE, x.lab.rot = T
        ) + coord_flip()
        
        ggsave(paste(output, "dotplot_stage.png", sep = "/"), width = 12, height = 24, dpi = 300, units = "in")
        
        
        # find out the cluster modules
        temp_data = data[unique(intersect(cluster_spec$gene, rownames(data))),]
        expr = ExpressionSet(
            as.matrix(temp_data)
        )
          
        cl = mfuzz(expr, c=16, m=1.25)
        res = as.data.frame(cl$cluster)
        write.table(res, file = paste(output, "cluster_data.txt"), quote=F, col.name=F, sep = "\t")
         
        for (i in unique(cl$cluster)) {
            temp_genes = names(cl$cluster[cl$cluster == i])
            
            p <- DoHeatmap(
                obj, 
                group.by = "res.0.6", 
                genes.use = temp_genes, 
                cex.col=0,
                slim.col.label = TRUE, 
                remove.key = TRUE,
                do.plot = F
            )
            
            
            height = length(temp_genes) / 8
            
            if(height < 5) {
                height = 5
            } else if(height > 40){
                height = 40
            }
            
            ggsave(
                paste(heatmap_output, paste0("cluster_", i, ".png"), sep = "/"), 
                width = 6,
                height = height,
                dpi = 600, 
                units = "in",
                limitsize = F
            )
        }
        
        # remove cluster-specific genes in stage markers
        cluster_spec = cluster_spec[cluster_spec$p_val_adj < 0.01 & cluster_spec$avg_logFC > 0.6, ]
        stage_spec = stage_spec[!stage_spec$gene %in% cluster_spec$gene, ]
        
        p <- DotPlot(obj, genes.plot=unique(stage_spec$gene), cols.use = c("lightgrey", "blue"),
          col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
          scale.by = "radius", scale.min = NA, scale.max = NA, group.by="res.0.6",
          plot.legend = FALSE, do.return = TRUE, x.lab.rot = T
        ) + coord_flip()
        
        ggsave(paste(output, "dotplot_stage_filtered.png", sep = "/"), width = 12, height = 24, dpi = 300, units = "in") 
    }
    
    # find out stage modules
    temp_data = data[unique(stage_spec$gene),]
    expr = ExpressionSet(
        as.matrix(temp_data)
    )
      
    cl = mfuzz(expr, c=16, m=1.25)
     
    for (i in unique(cl$cluster)) {
        temp_genes = names(cl$cluster[cl$cluster == i])
        
        p <- DoHeatmap(
            obj, 
            group.by = "Stage", 
            genes.use = temp_genes, 
            cex.col=0,
            slim.col.label = TRUE, 
            remove.key = TRUE,
            do.plot = F
        )
        
        
        height = length(temp_genes) / 8
        
        if(height < 5) {
            height = 5
        } else if(height > 40){
            height = 40
        }
        
        ggsave(
            paste(heatmap_output, paste0("stage_", i, ".png"), sep = "/"), 
            width = 6,
            height = height,
            dpi = 600, 
            units = "in",
            limitsize = F
        )
    }
    
    res = as.data.frame(cl$cluster)
    write.table(res, file = paste(output, "stage_data.txt"), quote=F, col.name=F, sep = "\t")
    """

    temp_script = os.path.join(input_directory, "temp.R")
    with open(temp_script, "w+") as w:
        w.write(rscript)

    try:
        with open(os.devnull) as w:
            check_call("Rscript {0}".format(temp_script), shell=True, stderr=w, stdout=w)
        os.remove(temp_script)
    except CalledProcessError as err:
        print(err)


def main(argv):
    u"""
    :param argv: sys.argv
    :return:
    """
    args = command_line(argv)

    files = glob.glob(os.path.join(args.i, "*_SCC"))

    files += glob.glob(os.path.join(args.i, "*_ADC"))

    with Pool(args.p) as p:
        p.map(_mfuzz_, files)

    pass


if __name__ == '__main__':
    main(sys.argv)
