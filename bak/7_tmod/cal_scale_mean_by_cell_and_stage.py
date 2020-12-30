#!/usr/bin/env python
#-*- utf-8 -*-
u"""
@since 2019.06.17
This scripts is used to collect the scale data for radar plots
"""
import os
from glob import glob
from multiprocessing import Pool
from subprocess import check_call
import sys


def call(args):
    u"""
    call Rscript
    :param args: path to rds; cell name and output
    """
    rds, cell, output = args
    
    Rscript = """
    library(dplyr)
    library(reshape2)
    
    obj <- readRDS("{0}")
    
    data <- melt(as.matrix(obj@scale.data))
    colnames(data) <- c("gene", "Cells", "value")
    
    data <- merge(data, obj@meta.data[, c("Cells", "Stage")], by = "Cells")
    data$Cell = "{1}"
    
    res = as.data.frame(data %>% group_by(Cells, Stage) %>% mutate(m = mean(value)) %>% select(gene, Cells, Stage, m) %>% unique())
    
    write.csv(res, file = "{2}") 
    
    """.format(
        rds, cell, output
    )
    
    temp = "{0}.R".format(output)
    with open(temp, "w+") as w:
        w.write(Rscript)
        
    check_call("Rscript {0}".format(temp), shell = True)
    
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
                
    with Pool(10) as p:
        p.map(call, tasks)
    

if __name__ == '__main__':
    main()
        
    