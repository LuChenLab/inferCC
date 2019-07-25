#!/usr/bin/env python3
# -*- conding: utf-8 -*-
u"""
@since 2019.07.24
"""
import os
from glob import glob

from subprocess import check_call, CalledProcessError

from multiprocessing import Pool


__dir__ = os.path.abspath(os.path.dirname(__file__))


def call(cmd):
    with open(os.devnull, "w+") as w:
        try:
            check_call(cmd, stdout=w, stderr=w, shell=True)
        except CalledProcessError as err:
            print(err)


def main(input_dir, metacell_dir, meta, n_jobs = 20):
    u"""
    Main function
    :param input_dir
    :param metacell_dir
    :param n_jobs: n_jobs
    """
    
    input_dir = os.path.abspath(input_dir)
    metacell_dir = os.path.abspath(metacell_dir)
    
    rds = glob(os.path.join(input_dir, "*/seurat.rds"))
    
    tasks = []
    
    __rscripts__ = os.path.join(__dir__, "correlation_with_metainfo.R")
    
    for i in rds:
        cell_name = os.path.basename(os.path.dirname(i))
        
        tasks.append(
            f"Rscript {__rscripts__} {i} {os.path.join(metacell_dir, cell_name)} {os.path.join(os.path.dirname(i), 'correlation')} {os.path.abspath(meta)}"
        )
    
    with Pool(n_jobs) as p:
        p.map(call, tasks)
    

if __name__ == '__main__':
    from fire import Fire
    Fire(main)
