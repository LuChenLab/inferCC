#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.28

Just used to call r scripts
"""
import os
import re
from glob import glob

from subprocess import check_call, CalledProcessError
from multiprocessing import Pool


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



def runr(input_dir, output_dir, n_jobs = 30):
    u"""
    main function
    :param input_dir
    :param output_dir
    :param n_jobs
    """
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    
    tasks = []
    for i in glob(os.path.join(input_dir, "*/*.rds")):
        if "monocle" not in i and "slingshot" not in i and "wgcna" not in i.lower():
            cell_name = os.path.basename(os.path.dirname(i))
            
            tasks.append(
                "Rscript {0} {1} {2}".format(
                    os.path.join(__dir__, "cluster.R"),
                    i,
                    os.path.join(output_dir, cell_name)
                )
            )

    with Pool(n_jobs) as p:
        p.map(__call__, sorted(tasks))



def __convert__(args):
    u"""
    convert img to pdf
    :param path input directory
    """
    path, output_dir = args
    if not os.listdir(path):
        return
    
    pngs = []
    for i in glob(os.path.join(path, "*.png")):
        if re.search(r"\d+_\w+.png", os.path.basename(i)):
            pngs.append(int(os.path.basename(i).split("_")[0]))
    pngs = set(pngs)

    temp_md = os.path.join(path, "temp.md")
        
    with open(temp_md, "w+") as w:
        w.write("# {0}\n\n".format(os.path.basename(path)))
        
        for i in sorted(pngs):
            w.write("### {}\n\n".format(i))
            for j in ["gene.png", "volca.png"]:
                file = os.path.join(path, "{0}_{1}".format(i, j))
                if os.path.exists(file):
                    w.write("![]({})\n\n".format(file))
                    
    
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


if __name__ == "__main__":
    from fire import Fire
    Fire()
