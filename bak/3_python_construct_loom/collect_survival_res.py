#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@Created at 2019.05.10
This script is used to collect all pvalue and median survival days calculated by TCGAanlyze_survival (TCGAbiolinks)
"""
import os
import sys
from glob import glob

from fire import Fire
import pandas as pd
from tqdm import tqdm



def main(input_dir: str, output_file: str) -> None:
    u"""
    main function

    this is used to collect all the data calculated by TCGAanalyze_survival from TCGAbiolinks package

    :param input_dir: path to input dir
    :param output_file: path to output file
    :return None
    """
    # check file existance
    input_dir = os.path.abspath(input_dir)
    output_file = os.path.abspath(output_file)

    assert os.path.exists(input_dir), "{0} not found".format(input_dir)

    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # find all files
    files = glob(os.path.join(input_dir, "*.txt"))

    data = []
    for i in tqdm(files):
        temp_data = {}
        with open(i) as r:
            for line in r:
                lines = line.split()
                temp_data[lines[0]] = lines[1]
        data.append(temp_data)

    data = pd.DataFrame(data)

    data.to_csv(output_file)


if __name__ == '__main__':
    Fire(main)

    main("Survival_plots/LUSC", "LUSC_survival.csv")
