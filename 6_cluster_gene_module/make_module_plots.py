#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.03
"""

import os
from glob import glob
from multiprocessing import Pool
from subprocess import check_call, CalledProcessError

from fire import Fire


__dir__ = os.path.dirname(os.path.abspath(__file__))


def call(cmd):
    try:
        with open(os.devnull, "w+") as w:
            check_call(cmd, shell=True, stdout=w, stderr=w)
    except CalledProcessError as err:
        print(err)


def main(path, top_n=50, n_jobs=20):
    u"""
    Main function
    :param path: path to input
    :param output:
    :return:
    """
    files = glob(os.path.join(path, "*_SCC/annotation_results_by_cluster.xlsx"))
    files += glob(os.path.join(path, "*_ADC/annotation_results_by_cluster.xlsx"))

    tasks = []

    for i in files:
        print(i)
        dirname = os.path.abspath(os.path.dirname(i))

        rds = glob(os.path.join(dirname, "*.rds"))

        for j in rds:
            if "monocle" in j or "slingshot" in j or "wgcna" in j.lower():
                continue
            tasks.append(f"Rscript {os.path.join(__dir__, 'make_heatmap_dotplot.R')} {dirname} {j} {top_n} \"{'res.0.6' if 'Basal_SCC' not in j else 'res.0.4' }\"")

    with Pool(n_jobs) as p:
        p.map(call, tasks)


if __name__ == '__main__':
    Fire(main)
