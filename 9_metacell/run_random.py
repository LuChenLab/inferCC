#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.07.01

Just used to call r scripts
"""
import os
from glob import glob

from subprocess import check_call, CalledProcessError
from multiprocessing import Pool


__dir__ = os.path.dirname(os.path.abspath(__file__))


def __call__(cmd):
    u"""
    Call R
    :param cmd
    """
    print(cmd)

    try:
        check_call(cmd, shell=True)
    except CalledProcessError as err:
        print(err)



def runr(input_dir, output_dir):
    u"""
    main function
    :param input_dir
    :param output_dir
    :param n_jobs
    """
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for i in glob(os.path.join(input_dir, "metacell_random_selected_*.rds")):
        print(i)
        if "metacell_random_selected_NM.rds" in i:
            continue
        
        cell_name = os.path.basename(i)
        cell_name = cell_name.split(".")[0]
        
        __call__(
            "Rscript {0} {1} {2} {3}".format(
                os.path.join(__dir__, "metacell_on_all_cells.R"),
                i,
                output_dir,
                cell_name
            )
        )


if __name__ == "__main__":
    from fire import Fire
    Fire(runr)