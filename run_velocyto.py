#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by Zhang at 2019.01.17
"""

import os
import sys
import argparse as ap
import subprocess as sp
from progressbar import ProgressBar


def call(cmd):
    u"""

    :param cmd:
    :return:
    """
    print(" ".join(cmd))
    try:
        sp.check_call(" ".join(cmd), shell=True)
    except sp.CalledProcessError:
        call(cmd)


def main(args):
    u"""
    main
    :param args:
    :return:
    """
    gene_gtf = os.path.abspath(args.gtf)
    rmsk = os.path.abspath(args.rmsk)
    indir = os.path.abspath(args.input)



    cmds = [
        "velocyto",
        "run10x",
        "-@ 20",
        "-m",
        rmsk,
        None,
        gene_gtf
    ]

    files = [os.path.join(indir, i) for i in os.listdir(indir) if not i.endswith("sh")]

    with ProgressBar(max_value=len(files), redirect_stdout=True) as bar:
        for idx, f in enumerate(files):
            bar.update(idx)
            cmds[-2] = f
            call(cmds)


if __name__ == '__main__':
    parser = ap.ArgumentParser("Run valocyto")
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input directory, cellranger output",
        type=str,
        required=True
    )
    parser.add_argument(
        "-g",
        "-gtf",
        help="Path to cellranger prepared genes.gtf",
        type=str,
        required=True
    )
    parser.add_argument(
        "-r",
        "--rmsk",
        type=str,
        required=True
    )

    if len(sys) <= 1:
        parser.print_help()
    else:
        args = parser.parse_args(sys.argv[1:])
        main(args)
