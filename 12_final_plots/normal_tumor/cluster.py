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
    # print(cmd)
    with open(os.devnull, "w+") as w:
        try:
            check_call(cmd, shell=True, stderr=w, stdout=w)
        except CalledProcessError as err:
            print(err)



def runr(input_dir, n_jobs = 30):
    u"""
    main function
    :param input_dir
    :param output_dir
    :param n_jobs
    """
    input_dir = os.path.abspath(input_dir)
    
    tasks = []
    for i in glob(os.path.join(input_dir, "*/seurat.rds")):
        
        tasks.append(
            "Rscript {0} {1} {2}".format(
                os.path.join(__dir__, "cluster.R"),
                i,
                os.path.join(os.path.dirname(i), "clusters")
            )
        )

    with Pool(n_jobs) as p:
        p.map(__call__, sorted(tasks))


if __name__ == "__main__":
    from fire import Fire
    Fire(runr)
