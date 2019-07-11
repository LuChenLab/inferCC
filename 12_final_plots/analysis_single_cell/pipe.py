#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.07.07
"""
import os
from glob import glob
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool

__dir__ = os.path.dirname(os.path.abspath(__file__))




def call(cmd):
    u"""
    
    """
    with open(os.devnull) as w:
        try:
            check_call(cmd, shell = True, stdout = w, stderr = w)
        except CalledProcessError as err:
            print(err)
            
def main(input_dir, output_dir, n_jobs = 20):
    u"""
    basic_script = args[1]
    input_rds = args[2]
    image_dir = args[3]
    trajectory = args[4]
    """
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)

    tasks = []
    for i in os.listdir(input_dir):
        if i.endswith("rds"):
            tasks.append(
                "Rscript {0} {1} {2} {3} DO".format(
                    os.path.join(__dir__, "pipe_seurat_monocle_slingshot.R"),
                    os.path.join(__dir__, "basic.R"),
                    os.path.join(input_dir, i),
                    os.path.join(output_dir, i.split(".")[0].replace("+", ""))
                )
            )
            
                
    with Pool(n_jobs) as p:
        p.map(call, tasks)
        

if __name__ == '__main__':
    from fire import Fire
    Fire(main)

