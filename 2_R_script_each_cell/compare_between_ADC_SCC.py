#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.21
"""
import os
from glob import glob
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool


__dir__ = os.path.dirname(os.path.abspath(__file__))
RSCRIPT = os.path.join(__dir__, "compare_between_ADC_SCC.R")
RSCRIPT1 = os.path.join(__dir__, "make_scatter_between_clusters.R")


CELLS = [
    "Alveolar_II",
    "B_cells",
    "Basal",
    "Exhaust_T",
    "CD4",
    "Dendritic",
    "Endothelial",
    "Fibroblasts",
    "Goblet",
    "Granulocyte",
    "Mast",
    "Monocytes",
    "Neuroendocine",
    "NK",
    "T_cells",
    "Treg"
]

ALL_CELL = {
   "Alveolar_II",
    "B_cells",
    "Basal",
    "CD4",
    "Ciliated",
    "Club",
    "Dendritic",
    "Endothelial",
    "Erythoroid_precursor",
    "Exhaust_T",
    "Fibroblasts",
    "Goblet",
    "Granulocyte",
    "Mast",
    "Monocytes",
    "Neuroendocine",
    "NK",
    "T_cells",
    "Treg" 
}


def call(args):
    u"""
    Call R
    :param args: path to rds and path to output directory
    """

    with open(os.devnull) as w:
        try:
            check_call("Rscript {0} {1} {2}".format(args[0], args[1], args[2]), shell=True, stdout = w, stderr = w)
        except CalledProcessError as err:
            print(err)

        
def main(input_dir, output_dir, n_jobs=10, all_cell=False):
    u"""
    Main
    :param input_dir path to directory contains all results
    :param output_dir path to output directory
    :param n_jobs: run how many tasks at same time
    """
    
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    tasks = []
    
    if not all_cell:    
        for i in CELLS:
            path = os.path.join(input_dir, i)
            
            rds = glob(os.path.join(path, "*.rds"))
            
            for j in rds:
                if "monocle" not in j and "slingshot" not in j and "wgcna" not in j.lower():
                    tasks.append([RSCRIPT, j, os.path.join(output_dir, i)])
    
    else:
        for i in ALL_CELL:
            path = os.path.join(input_dir, i)
            rds = glob(os.path.join(path, "*.rds"))
            
            for j in rds:
                if "monocle" not in j and "slingshot" not in j and "wgcna" not in j.lower():
                    tasks.append([RSCRIPT1, j, os.path.join(output_dir, i)])
                    
            
            rds = glob(os.path.join(path + "_SCC", "*.rds"))
            for j in rds:
                if "monocle" not in j and "slingshot" not in j and "wgcna" not in j.lower():
                    tasks.append([RSCRIPT1, j, os.path.join(output_dir, i + "_SCC")])
            
            rds = glob(os.path.join(path + "_ADC", "*.rds"))
            for j in rds:
                if "monocle" not in j and "slingshot" not in j and "wgcna" not in j.lower():
                    tasks.append([RSCRIPT1, j, os.path.join(output_dir, i + "_ADC")])
                
    with Pool(n_jobs) as p:
        p.map(call, tasks)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
