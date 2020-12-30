#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.26
"""
from glob import glob
import pandas as pd
import os
from tqdm import tqdm
from subprocess import check_call
from multiprocessing import Pool
import sys


def read_genes(path):
    data = pd.read_excel(path, sheet_name = 1)
    genes = set()
    for i in data["Genes"]:
        for j in i.split("|"):
            genes.add(j)
    return genes


def call(cmd):
    rds, cell, output = cmd
    
    Rscript = """
    library(reshape2)
    library(dplyr)

    # args = commandArgs(trailingOnly = T)
    obj <- readRDS("{0}")

    genes = read.table("/mnt/raid62/Lung_cancer_10x/tmod_stage/SYSU/gene.txt")

    data = melt(as.matrix(obj@raw.data[genes[,1],]))

    data = merge(data, obj@meta.data[, c("Stage", "Cells")], by.x = "Var2", by.y = "Cells")

    temp = data %>% select(Var1, value, Stage) %>% group_by(Stage, Var1) %>% mutate(m
    = mean(value)) %>% select(Var1, Stage, m) %>% unique()

    temp = as.data.frame(temp)

    temp$cell = "{1}"

    write.csv(temp, file = "{2}")
    """.format(rds, cell, output)
    
    temp = f"{output}.R"
    with open(temp, "w+") as w:
        w.write(Rscript)
        
    check_call(f"Rscript {temp}", shell = True)
    
    os.remove(temp)
    
def main():
    output = os.path.abspath(sys.argv[2])
    if not os.path.exists(output):
        os.makedirs(output, mode=0o777, exist_ok=False)
    
    tasks = []
    for i in glob(os.path.join(os.path.abspath(sys.argv[1]), "*_SCC")):
        cell = os.path.basename(i)
        for j in glob(os.path.join(i, "*.rds")):
            if "monocle" not in j and "slingshot" not in j and "wgcna" not in j.lower():
                tasks.append([
                    j, cell, os.path.join(output, f"{cell}.csv")
                ])

    for i in glob(os.path.join(sys.argv[1], "*_ADC")):
        cell = os.path.basename(i)
        for j in glob(os.path.join(i, "*.rds")):
            if "monocle" not in j and "slingshot" not in j and "wgcna" not in j.lower():
                tasks.append([
                    j, cell, os.path.join(output, f"{cell}.csv")
                ]) 
                
    with Pool(20) as p:
        p.map(call, tasks)
     
        
if __name__ == '__main__':
    main()