#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.28

Just used to call r scripts
"""
import os

from subprocess import check_call, CalledProcessError
from multiprocessing import Pool


__dir__ = os.path.dirname(os.path.abspath(__file__))


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
    for parent, dirs, files in os.walk(input_dir):
        for f in files:
            if f == "seurat.rds":
                cell_name = parent.replace(input_dir, "").strip("/")
               
                tasks.append(
                    "Rscript {0} {1} {2} 1".format(
                        os.path.join(__dir__, "metacell_on_all_cells.R"),
                        os.path.join(parent, f),
                        os.path.join(output_dir, cell_name)
                    )
                )

    with Pool(n_jobs) as p:
        p.map(__call__, sorted(tasks))
    

if __name__ == "__main__":
    from fire import Fire
    Fire(runr)
    