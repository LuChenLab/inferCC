#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.12

This script is used to call makt_tmod_mfuzz_heatmap.R
"""

import os
import re

import pandas as pd

from glob import glob
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool


__dir__ = os.path.abspath(os.path.dirname(__file__))

CELLS_NAME = {
    "Alveolar_II": "APII",
    "B_cells": "B",
    "Basal": "Basal",
    "CD4": "CD4",
    "Ciliated": "Cili",
    "Club": "Club",
    "Dendritic": "Dend",
    "Endothelial": "Endp",
    "Erythroid_precursor": "Eryt",
    "Fibroblasts": "Fibr",
    "Goblet": "Goblet",
    "Granulocyte": "Granu",
    "Mast": "Mast",
    "Monocytes": "Mono",
    "NK": "NK",
    "T_cells": "T",
    "Treg": "Treg",
    "Neuroendocrine": "Neuro",
    "Exhaust_T": "Exha"
}

def call(cmd):
    u"""
    Wrapper of check_call
    :param command to run
    """
    with open(os.devnull, "w+") as w:
        try:
            check_call(cmd, stdout = w, stderr = w, shell = True)
        except CalledProcessError as err:
            print(err)


def main(path, overlap_file, output, n_jobs = 10):
    u"""
    main
    :param path: input directory
    :param overlap_file: path to tmod and mfuzz overlap file
    :param output: ooutput directory
    """
    
    # cells = set()
    # with open(overlap_file) as r:
    #     for line in r:
    #         cells.add(line.split()[0])
    data = pd.read_excel(overlap_file)
    cells = data["Cell_name"].unique()           
    
    tasks = []
    for cell in cells:
        clean_cell = re.sub("_(ADC|SCC)(_bak)?$", "", cell)
        for rds in glob(os.path.join(path, clean_cell, "*.rds")):
            if "monocle" not in rds and "slingshot" not in rds and "wgcna" not in rds.lower():
                tasks.append(
                    "Rscript {0} {1} {2} {3} {4} {5}".format(
                        os.path.join(__dir__, "make_tmod_mfuzz_heatmap.R"),
                        os.path.abspath(overlap_file),
                        cell, 
                        os.path.abspath(rds),
                        os.path.abspath(output),
                        CELLS_NAME[clean_cell]
                    )
                )
                
                # tasks.append(
                #     "Rscript {0} {1} {2} {3} {4}".format(
                #         os.path.join(__dir__, "make_radar_plots.R"), 
                #         os.path.abspath(overlap_file), 
                #         os.path.abspath(rds),
                #         os.path.join(os.path.abspath(output), cell, CELLS_NAME[clean_cell]),
                #         cell
                #     )
                # )
    
    with Pool(n_jobs) as p:
        p.map(call, tasks)
        

if __name__ == '__main__':
    from fire import Fire
    Fire(main)        
