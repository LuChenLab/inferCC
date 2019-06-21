#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.21
"""
import os
from glob import glob
from subprocess import check_call
from multiprocessing import Pool


__dir__ = os.path.dirname(os.path.abspath(__file__))
RSCRIPT = os.path.join(__dir__, "compare_between_ADC_SCC.R")


CELLS = [
    "Alveolar_II",
    "B_cells",
    "Basal",
    "Exhaust_T",
    "CD4",
    "Dendritic",
    "Endothelial",
    "Fibroblast",
    "Goblet",
    "Granulocyte",
    "Mast",
    "Monocytes",
    "Neuroendocine",
    "NK",
    "T_cells",
    "T_reg"
]


def call(args):
    u"""
    Call R
    :param args: path to rds and path to output directory
    """

    with open(os.devnull) as w:
        check_call("Rscript {0} {1} {2}".format(RSCRIPT, args[0], args[1]), shell=True, stdout = w, stderr = w)
        
def main(input_dir, output_dir, n_jobs=10):
    u"""
    Main
    :param input_dir path to directory contains all results
    :param output_dir path to output directory
    :param n_jobs: run how many tasks at same time
    """
    
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    tasks = []
    
    for i in CELLS:
        path = os.path.join(input_dir, i)
        
        rds = glob(os.path.join(path, "*.rds"))
        
        for j in rds:
            if "monocle" not in j and "slingshot" not in j and "wgcna" not in j.lower():
                tasks.append([j, os.path.join(output_dir, i)])
                
    with Pool(n_jobs) as p:
        p.map(call, tasks)
    

if __name__ == '__main__':
    from fire import Fire
    Fire(main)
