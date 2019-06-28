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
    # print(cmd)
    with open(os.devnull, "w+") as w:
        try:
            check_call(cmd, shell=True, stdout=w, stderr=w)
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
    for i in CELLS.values():
        temp_dir = os.path.join(input_dir, i)
        for j in glob(os.path.join(temp_dir, "*.rds")):
            if "monocle" not in j and "slingshot" not in j and "wgcna" not in j.lower():
                tasks.append(
                    "Rscript {0} {1} {2}".format(
                        os.path.join(__dir__, "normal_cancer.R"),
                        j,
                        os.path.join(output_dir, i)
                    )
                )

    with Pool(n_jobs) as p:
        p.map(__call__, tasks)



def __convert__(path):
    u"""
    convert img to pdf
    :param path input directory
    """
    if not os.listdir(path):
        return
    
    disease = ["ADC", "SCC"]
    pngs = {
        "Top10 variable genes": "gene.png",
        "KEGG": "kegg.png",
        "GO": "go.png",
        "DO": "do.png"
    }

    temp_md = os.path.join(path, "temp.md")
    
    with open(temp_md, "w+") as w:
        for i in disease:
            w.write("# {0}\n\n".format(i))
            
            for title, j in pngs.items():
                file = os.path.join(path, "{0}_{1}".format(i, j))
                
                if os.path.exists(file):
                    w.write("### {}\n\n".format(title))
                    w.write("![]({})\n\n".format(file))
                    
    
    check_call(
        "python {0} -md {1} -pdf {2}".format(
            os.path.join(os.path.dirname(__dir__), "mdconverter/converter.py"),
            temp_md,
            os.path.join(os.path.dirname(path), os.path.basename(path) + ".pdf")
        ),
        shell=True
    )
    
    os.remove(temp_md)


def pdf(input_dir, n_jobs = 20):
    u"""
    make pdf
    :param input_dir
    """
    input_dir = os.path.abspath(input_dir)
    tasks = [os.path.join(input_dir, x) for x in os.listdir(input_dir)]
    tasks = [x for x in tasks if os.path.isdir(x)]
    
    with Pool(n_jobs) as p:
        p.map(__convert__, tasks)
        

if __name__ == "__main__":
    from fire import Fire
    Fire()


        