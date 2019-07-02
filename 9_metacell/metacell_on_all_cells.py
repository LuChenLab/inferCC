#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.28

Just used to call r scripts
"""
import os
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



def runr(input_dir, output_dir, n_jobs = 20):
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
            
            __call__(
                "Rscript {0} {1} {2} {3}".format(
                    os.path.join(__dir__, "metacell_on_all_cells.R"),
                    i,
                    os.path.join(output_dir, cell_name), 
                    "res.0.6"
                )
            )


if __name__ == "__main__":
    from fire import Fire
    Fire()
