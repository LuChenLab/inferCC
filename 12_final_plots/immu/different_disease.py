#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.07.22
"""

import os
from glob import glob
from subprocess import check_call, CalledProcessError

from multiprocessing import Pool


__dir__ = os.path.abspath(os.path.dirname(__file__))


def call(cmd):
    u"""
    Call
    """

    with open(os.devnull, "w+") as w:
        try:
            check_call(cmd, stdout = w, stderr = w, shell = True)
        except CalledProcessError as err:
            print(err)


def main(input_dir, n_jobs = 10):
    u"""
    Main function
    """

    input_dir = os.path.abspath(input_dir)

    tasks = []

    for i in glob(os.path.join(input_dir, "*/seurat.rds")):
        cmd = f"Rscript {os.path.join(__dir__, 'different_disease.R')} {i} {os.path.join(os.path.dirname(i), 'disease')}"
        tasks.append(cmd)

    with Pool(n_jobs) as p:
        p.map(call, tasks)


if __name__ == '__main__':
    from fire import Fire

    Fire(main)
