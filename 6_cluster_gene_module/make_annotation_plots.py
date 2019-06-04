#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.03
@author Zhang Yiming

Just call Rscript in batch
"""
import os
from glob import glob
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool
from fire import Fire


__dir__ = os.path.join(os.path.abspath(__file__))


def call(cmd):
    try:
        with open(os.devnull, "w+") as w:
            check_call(cmd, shell=True, stdout=w, stderr=w)
    except CalledProcessError as err:
        print(err)


def main(path, n_jobs=15):
    files = glob(os.path.join(path, "*/cluster_gene_module/results.xlsx"))

    tasks = []
    for i in files:
        tasks.append(f"Rscript {os.path.join(__dir__, 'annotation_plots.R')} {os.path.dirname(i)}")

    with Pool(n_jobs) as p:
        p.map(call, tasks)


if __name__ == '__main__':
    Fire(main)
