#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.07.01

Used to make different plots
"""

import os
from glob import glob

from subprocess import check_call, CalledProcessError
from multiprocessing import Pool

import pandas as pd



__dir__ = os.path.dirname(os.path.abspath(__file__))


CELLS = {
    "Pulmonary alveolar type II cells": "Alveolar_II",
    "B cells": "B_cells",
    "Basal cells": "Basal",
    "CD4+ T cells": "CD4",
    "Ciliated cells": "Ciliated",
    "Club cells": "Club",
    "Dendritic": "Dendritic",
    "Endothelial cells": "Endothelial",
    "Erythroid-like and erythroid precursor cells": "Erythroid_precursor",
    "Exhaust T": "Exhaust_T",
    "Fibroblasts": "Fibroblasts",
    "Airway goblet cells": "Goblet",
    "Granulocyte": "Granulocyte",
    "Mast": "Mast",
    "Monocytes": "Monocytes",
    "Neuroendocrine cell": "Neuroendocrine",
    "NK": "NK",
    "T cells": "T_cells",
    "Treg": "Treg"
}


def __call__(cmd):
    u"""
    Call R
    :param cmd
    """
    print(cmd)
    with open(os.devnull, "w+") as w:
        try:
            check_call(cmd, shell=True, stderr=w, stdout=w)
        except CalledProcessError as err:
            print(err)


def __get_rds__(rds):
    for i in rds:
        if "monocle" not in i and "slingshot" not in i and "wgcna" not in i.lower():
            return i

def runr(input_dir, metacell_dir, output_dir, n_jobs = 20):
    u"""
    main
    """
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    metacell_dir = os.path.abspath(metacell_dir)
    
    tasks = []
    for i in CELLS.values():
        
        rds = __get_rds__(glob(os.path.join(input_dir, i, "*.rds")))
        
        # if os.path.exists(os.path.join(metacell_dir, i)):
        #     tasks.append(
        #         "Rscript {0} {1} {2} {3} {4} {5}".format(
        #             os.path.join(__dir__, "make_meta_cell_heatmap.R"),
        #             rds,
        #             os.path.join(metacell_dir, i, "mat.{0}.Rda".format(i)),
        #             os.path.join(metacell_dir, i, "mc.{0}_mc_f.Rda".format(i)), 
        #             os.path.join(metacell_dir, i, "gstat.{0}.Rda".format(i)),
        #             os.path.join(output_dir, i)
        #         )
        #     )

        # rds = args[1]
        # scmat = args[2]
        # mc = args[3]
        # gstat = args[4]
        # output_dir = args[5]
        
        for j in ["ADC", "SCC"]:
            temp = "{0}_{1}".format(i, j)
            
            if not os.path.exists(os.path.join(metacell_dir, temp)):
                continue
            
            rds = __get_rds__(glob(os.path.join(input_dir, temp, "*.rds"))) 
        
            tasks.append(
                "Rscript {0} {1} {2} {3} {4} {5}".format(
                    os.path.join(__dir__, "make_meta_cell_heatmap.R"),
                    rds,
                    os.path.join(metacell_dir, temp, "mat.{0}.Rda".format(temp)),
                    os.path.join(metacell_dir, temp, "mc.{0}_mc_f.Rda".format(temp)), 
                    os.path.join(metacell_dir, temp, "gstat.{0}.Rda".format(temp)),
                    os.path.join(output_dir, temp)
                )
            )
            
    with Pool(n_jobs) as p:
        p.map(__call__, tasks)


def __convert__(args):
    u"""
    convert img to pdf
    :param path input directory
    """
    path, output_dir = args
    if not os.listdir(path):
        return
    
    pngs = glob(os.path.join(path, "*.png"))

    temp_md = os.path.join(path, "temp.md")
        
    with open(temp_md, "w+") as w:
        w.write("# {0}\n\n".format(os.path.basename(path).replace("_", "\_")))
        
        for i in sorted(pngs):
            w.write("### {}\n\n".format(os.path.basename(i).replace("_", "\_")))
            w.write("![]({})\n\n".format(i))
                    
    
    check_call(
        "python {0} -md {1} -pdf {2}".format(
            os.path.join(os.path.dirname(__dir__), "mdconverter/converter.py"),
            temp_md,
            os.path.join(output_dir, os.path.basename(path) + ".pdf")
        ),
        shell=True
    )
    
    os.remove(temp_md)


def pdf(input_dir, output_dir, n_jobs = 20):
    u"""
    make pdf
    :param input_dir
    """
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    tasks = [os.path.join(input_dir, x) for x in os.listdir(input_dir)]
    tasks = [[x, output_dir] for x in tasks if os.path.isdir(x)]
    
    with Pool(n_jobs) as p:
        p.map(__convert__, tasks)

        
if __name__ == '__main__':
    from fire import Fire
    Fire()


