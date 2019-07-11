#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.07.09
"""
import os

from glob import glob
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool


__dir__ = os.path.abspath(os.path.dirname(__file__))



def call(cmd):
    u"""
    run rscript
    """
    with open(os.devnull) as w:
        try:
            check_call(cmd, shell = True, stdout = w, stderr = w)
        except CalledProcessError as err:
            print(err)
            
            
def main(input_dir, n_jobs = 15):
    u"""
    Main
    """
    input_dir = os.path.abspath(input_dir)
    
    tasks = []
    for i in glob(os.path.join(input_dir, "*/markers_by_stage.xlsx")):
        target_dir = os.path.dirname(i)
        tasks.append(
            "Rscript {0} {1} {2}".format(
                os.path.join(__dir__, "run_mfuzz_on_stage.R"),
                target_dir,
                os.path.join(target_dir, "mfuzz")
            )
        )
    print(tasks) 
    with Pool(n_jobs) as p:
        p.map(call, tasks)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)

