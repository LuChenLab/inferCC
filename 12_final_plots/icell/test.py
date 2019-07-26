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


def main(input_dir, output_dir, n_jobs = 20):
    u"""
    Main function
    :param input_dir
    :param metacell_dir
    :param n_jobs: n_jobs
    """
    
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    
    rds = glob(os.path.join(input_dir, "LUAD/*/seurat.rds"))
    
    tasks = []
    
    __rscripts__ = os.path.join(__dir__, "test.R")
    
    for i in rds:
        cell_name = os.path.basename(os.path.dirname(i))
        
        j = i.replace('LUAD', 'LUSC')
        
        if not os.path.exists(j):
            continue
        
        tasks.append(
            f"Rscript {__rscripts__} {i} {j} {os.path.join(output_dir, cell_name)}"
        )
    
    with Pool(n_jobs) as p:
        p.map(call, tasks)
    

if __name__ == '__main__':
    from fire import Fire
    Fire(main)