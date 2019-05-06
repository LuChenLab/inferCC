#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu at 2019.01.17
"""

import os
import sys
import argparse as ap
import subprocess as sp


def call(cmd):
    u"""

    :param cmd:
    :return:
    """
    print(cmd)
    try:
        with open(os.devnull, "w+") as w:
            sp.check_call(cmd, shell=True)
    except sp.CalledProcessError:
        # print(cmd)
        pass


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

    files = []
    for i in os.listdir(indir):
        sample = os.path.join(indir, i)

        if os.path.isfile(sample):
            continue
        files.append(sample)

    for f in files:
        cmds[-2] = f

        call(" ".join(cmds))


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
        "--gtf",
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

    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        args = parser.parse_args(sys.argv[1:])
        main(args)
