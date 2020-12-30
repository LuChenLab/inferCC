#!/usr/bin/env python3
# -*- coding:utf- 8 -*-
u"""
@since 2019.06.10

Just used for pre-processing and testings different conditions
"""
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import seaborn as sns

import xml.etree.ElementTree as ET
import json
import pandas as pd
import os
from glob import glob
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool
from matplotlib_venn import venn2
from tqdm import tqdm


__dir__ = os.path.dirname(os.path.abspath(__file__))


CELLS_NAME = {
    "Alveolar_II": "APII",
    "B_cells": "B",
    "Basal": "Basal",
    "CD4": "CD4",
    "Ciliated": "Cili",
    "Club": "Club",
    "Dendritic": "Dend",
    "Endothelial": "Endp",
    "Erythroid_precursor": "Eryt",
    "Fibroblasts": "Fibr",
    "Goblet": "Goblet",
    "Granulocyte": "Granu",
    "Mast": "Mast",
    "Monocytes": "Mono",
    "NK": "NK",
    "T_cells": "T",
    "Treg": "Treg",
    "Neuroendocrine": "Neuro",
    "Exhaust_T": "Exha"
}


def iter_all_clusters(path):
    u"""
    read the xlsx file, and yield the information that needed
    :param path:
    :return:
    """
    data = pd.read_excel(path)

    for _, row in data.iterrows():
        yield [
            "{0}_{1}".format(row[0], row[1]),
            "stage_{0}_data.json".format(row[2]),
            row[3],
            row[4],
            "{0}.{1}.{2}.M{3}.{4}".format(row[1], CELLS_NAME[row[0]], row[4], row[2], row[3])
        ]


def call(cmd):
    with open(os.devnull, "w+") as w:
        try:
            check_call(cmd, shell=True, stderr=w, stdout=w)
        except CalledProcessError as err:
            print(err)


def main(path, xlsx, msig, output, n_jobs=20):
    u"""
    main function
    :param path:
    :param xlsx:
    :param output:
    :param n_jobs
    :return:
    """
    if not os.path.exists(output):
        os.makedirs(output)

    total_gene = set()

    tasks = []
    output_files = []
    for tmod_dir, j_data, ica_id, stage, outfile in iter_all_clusters(xlsx):
        outfile = os.path.join(output, "{0}.xlsx".format(outfile))

        input_data = os.path.join(path, tmod_dir, "tmod_stage", j_data)
        for rds in glob(os.path.join(path, tmod_dir, "*.rds")):
            if "monocle" not in rds and "slingshot" not in rds and "wgcna" not in rds.lower():

                with open(input_data) as r:
                    input_data = json.load(r)

                try:
                    genes = input_data[ica_id]
                except KeyError:
                    genes = input_data[str(ica_id)]

                genes = [x.strip() for x in genes]

                total_gene = total_gene | set(genes)

                temp = "{0}.txt".format(outfile)
                with open(temp, "w+") as w:
                    w.write("\n".join(genes))

                # temp_freq = []
                # for key, value in res.items():
                #     temp_freq.append(len(set(genes) & set(value)) / len(value))

                # temp_freq = pd.DataFrame(temp_freq)
                # temp_freq.columns = ["value"]
                # temp_freq["label"] = [os.path.basename(outfile).replace(".xlsx", "") for _ in range(temp_freq.shape[0])]
                # temp_freq["cell"] = [os.path.basename(outfile).split(".")[1] for _ in range(temp_freq.shape[0])]
                #
                # if percentage is None:
                #     percentage = temp_freq
                # else:
                #     percentage = pd.concat([percentage, temp_freq])

                tasks.append("Rscript {0} {1} {2} {3} {4}".format(
                    os.path.join(__dir__, "make_tmod.R"),
                    os.path.abspath(temp),
                    os.path.abspath(msig),
                    os.path.abspath(rds),
                    os.path.abspath(outfile)
                ))

                output_files.append(os.path.abspath(outfile))

    # fig, ax = plt.subplots()
    # sns.violinplot(
    #     x="label",
    #     y="value",
    #     hue="cell",
    #     data=percentage,
    #     ax=ax,
    #     dodge=False
    # )
    #
    # ax.tick_params(axis='x', rotation=90)
    # plt.yscale("log")
    # plt.legend(loc=9, ncol=4)
    # plt.tight_layout()
    # plt.savefig(os.path.join(output, "stats.png"), dpi=600)

    # with Pool(n_jobs) as p:
    #     p.map(call, tasks)

    print("merging")
    msd, logfc, pval, pc = [], [], [], []
    for i in tqdm(output_files):
        label = os.path.basename(i).replace(".xlsx", "")
        temp = pd.read_excel(i, sheet_name=1)
        temp["module"] = [label for _ in range(temp.shape[0])]

        if temp.shape[0] > 0:
            msd.append(temp)

        temp = pd.read_excel(i, sheet_name=2)
        temp["module"] = [label for _ in range(temp.shape[0])]
        if temp.shape[0] > 0:
            logfc.append(temp)

        temp = pd.read_excel(i, sheet_name=3)
        temp["module"] = [label for _ in range(temp.shape[0])]
        if temp.shape[0] > 0:
            pval.append(temp)

        temp = pd.read_excel(i, sheet_name=4)
        temp["module"] = [label for _ in range(temp.shape[0])]
        if temp.shape[0] > 0:
            pc.append(temp)

    writer = pd.ExcelWriter(os.path.join(output, "tmod.xlsx"), engine='xlsxwriter')

    pd.concat(msd, sort=False).to_excel(writer, sheet_name="MSD", index=False)
    pd.concat(logfc, sort=False).to_excel(writer, sheet_name="logFC", index=False)
    pd.concat(pval, sort=False).to_excel(writer, sheet_name="p_adj_val", index=False)
    pd.concat(pc, sort=False).to_excel(writer, sheet_name="PC", index=False)

    writer.save()

    [os.remove(x) for x in output_files]


if __name__ == '__main__':
    from fire import Fire
    Fire(main)

