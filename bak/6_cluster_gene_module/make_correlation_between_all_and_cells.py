#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.03
@author Zhang Yiming

Scripts to make correlation density plots between total cells and specific cells
"""
import os
from glob import glob
from multiprocessing import Pool

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from matplotlib import pyplot as plt
import seaborn as sns
from tqdm import tqdm


N_ROWS = 100


def read_data(path: list) -> dict:

    if path[0] == "All":
        data = pd.read_csv(path[1], index_col=0, nrows=N_ROWS, sep="\t")
    else:
        data = pd.read_csv(path[1], index_col=0, nrows=N_ROWS)

    return {path[0]: data}


def calculate(args):
    res, key = args
    temp_r_values, temp_label = [], []
    for gene in tqdm(set(res["All"].index) & set(res[key].index)):
        temp = np.corrcoef(res["All"].loc[gene, res[key].columns], res[key].loc[gene, :])[0, 1]

        if not np.isnan(temp):
            temp_r_values.append(temp)
            temp_label.append(key)

    return temp_r_values, temp_label


def main(total_counts: str, input_dir: str, output: str, n_jobs=15):
    u"""
    make plots
    :param total_counts:
    :param input_dir:
    :param output:
    :return:
    """
    files = glob(os.path.join(os.path.abspath(input_dir), "*/paga/normalized_counts.csv.gz"))

    tasks = [["All", os.path.abspath(total_counts)]]

    for x in files:
        if "SCC" in x or "ADC" not in x or "Tumor" not in x or "P01" not in x or "Survival" not in x or "SMOTE" not in x:
            tasks.append([
                os.path.basename(os.path.dirname(os.path.dirname(x))),
                x
            ])

    res = {}
    with Pool(n_jobs) as p:
        temp_res = p.map(read_data, tasks)

    p.close()
    p.join()

    for i in temp_res:
        res.update(i)

    temp_r_values = []
    temp_label = []

    tasks = [[res, x] for x in res.keys() if x != "All"]
    with Pool(n_jobs) as p:
        temp = p.map(calculate, tasks)
    p.close()
    p.join()

    for temp_res in temp:
        temp_r_values += temp_res[0]
        temp_label += temp_res[1]

    # for key in tqdm(keys):
    #     for gene in tqdm(set(res["All"].index) & set(res[key].index)):
    #         temp = np.corrcoef(res["All"].loc[gene, res[key].columns], res[key].loc[gene, :])[0, 1]
    #
    #         if not np.isnan(temp):
    #            temp_r_values.append(temp)
    #            temp_label.append(key)

    print(len(temp_r_values), len(temp_label))
    r_values = pd.DataFrame([temp_label, temp_r_values]).transpose()

    print(r_values.head())

    r_values.columns = ["cells", "r_value"]

    r_values.to_csv(f"{output}.csv")

    fig, ax = plt.subplots(figsize=(8, 6))

    sns.violinplot(x="cells", y="r_value", hue="cells", data=r_values, ax=ax)
    plt.tight_layout()
    plt.savefig(output, dpi=600)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)




