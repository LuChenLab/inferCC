#!/usr/bin/env python3
#-*-coding:utf-8 -*-
u"""

"""
import os
import json
import logging
import sys
from statistics import median, mean
from multiprocessing import Pool
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
# from scipy.stats import wilcoxon
from anndata import AnnData
from imblearn.under_sampling import EditedNearestNeighbours, RepeatedEditedNearestNeighbours
from sklearn.metrics import calinski_harabasz_score

logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s] - [%(name)s] - [%(levelname)s] - %(message)s')
logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s] - [%(name)s] - [%(levelname)s] - %(message)s')
logger = logging.getLogger(__name__)


__date__ = "2019.05.15"
__author__ = "ygidtu"
__email__ = "ygidtu@gmail.com"


def calculate_tsi(array) -> int:
    u"""
    calculate the tissue specific index of input array
    """

    if len(array) == 0:
        return 0
        
    if len(array) == 1:
        return 1

    array = array / max(array)
    return sum(array) / (len(array) - 1)


def calculate_tsi_of_dataframe(args) -> dict:
    u"""
    calculate the tissue specific index based on input dataframe
    :param data:
    :param key: the group key
    :return:
    """
    key, data = args

    temp_res = {}
    for idx, row in data.iterrows():
        tsi = calculate_tsi(row)

        temp_res[idx] = tsi

    return {key: temp_res}


def calculate_tsi_of_dict(args) -> dict:
    u"""
    calculate the tissue specific index based on input dict
    :param args:
        - gene: str.
        - data: dict.
    :return:
    """

    gene, data = args

    return {gene: calculate_tsi(list(data[gene].values()))}


def calculate_expression_by_group(args):
    u"""
    :param args: combinations of
        - data: AnnData
        - groups: str, column in data.obs
        - key: str, gene name
        - method: median or mean
        - log_method: log2 or log10
    :return:
    """
    key, data, groups, method, log_method, expressed_perc = args

    idents = data.obs[groups].unique()

    res = {}
    expressed_perc_by_group = []
    for i in idents:
        temp_idx = data.obs.loc[data.obs[groups] == i, :]

        temp_value = data.to_df().loc[temp_idx.index, key]

        # _, p = wilcoxon(temp_value)

        # print(sum(temp_value > 0) / len(temp_value), expressed_perc)

        expressed_perc_by_group.append(sum(temp_value > 0) / len(temp_value))

        res[i] = log_method(method(temp_value) + 1)

    if any([x > expressed_perc for x in expressed_perc_by_group]):
        return {key: res}

    return None


class GeneSignature:
    u"""
    class to handle the gene signature issues
    """
    def __init__(
            self,
            data: AnnData,
            group: str,
            expression_method="mean",
            expression_log_method="log2",
            tsi_low: float = 0.15,
            tsi_high: float = 0.85,
            log_fc: float = 1.25,
            expressed_perc: float = 0.1,
            # p_value: float = 0.05,
            processes: int = 10
    ):
        u"""
        init this class
        :param data: AnnData,
        :param group: str, column in data.obs
        :param tsi_low: float, the low threshold of tsi
        :param tsi_high: float, the high threshold of tsi
        :param expression_method: should be mean or median
        :param log_fc: the fold change threshold of mean gene expression
                        in first highly expressed group against second highly expressed group
        :param expressed_perc: the threshold of expressed percentage each group
        :param p_value: threshold of p value
        :param processes: how many cpu to use
        """

        if expression_method not in ("mean", "median"):
            raise ValueError("expression methods should be mean or median")

        methods = {
            "median": median,
            "mean": mean,
            "log2": np.log2,
            "log10": np.log10
        }

        self.data = data
        self.group = group
        self.tsi_low = tsi_low
        self.tsi_high = tsi_high
        self.log_fc = log_fc
        self.expression_method = methods[expression_method]
        self.expression_log_method = methods[expression_log_method]
        self.processes = processes
        # self.p_value = p_value
        self.expressed_perc = expressed_perc
        self.genes_expression = {}
        self.genes_most_expressed_group = {}
        self.genes_conserved_each_groups = {}
        self.genes_differential_expressed_between_groups = {}

    def __highly_variable_genes__(self):
        u"""
        find highly variable genes between groups
        :return:
        """
        logger.info("Find highly variable genes")

        args = []

        data = self.data.to_df()

        data = data.iloc[:, list(data.sum(axis=1) > 0)]

        for i in data.columns:
            args.append(
                [
                    i, self.data,
                    self.group, self.expression_method,
                    self.expression_log_method,
                    self.expressed_perc,
                    # self.p_value
                ]
            )

        with Pool(self.processes) as p:
            res = p.map(calculate_expression_by_group, args)

        [self.genes_expression.update(x) for x in res if x is not None]

        for gene, cells in self.genes_expression.items():
            temp = [[k, v] for k, v in cells.items()]

            if len(temp) > 0:
                temp = sorted(temp, key=lambda x: x[1], reverse=True)

                self.genes_most_expressed_group[gene] = temp[0][0]

    def __get_conserved_genes_inside_idents__(self):
        u"""
        calculate the tsi of each ident
        :return:
        """
        logger.info("Find conserved genes inside each group")

        idents = self.data.obs[self.group].unique()

        target_genes = {}
        for gene, expression in self.genes_expression.items():

            expression = sorted(expression.values(), reverse=True)

            if len(expression) == 0:
                continue
            elif len(expression) == 1:
                temp = 99
            elif expression[1] == 0:
                temp = expression[0] / 0.01
            else:
                temp = expression[0] / mean(expression)

            if temp > self.log_fc:
                most_expressed_cell = self.genes_most_expressed_group[gene]
                temp_genes = target_genes.get(most_expressed_cell, [])
                temp_genes.append(gene)
                target_genes[most_expressed_cell] = temp_genes

        for k in target_genes.values():
            logger.debug(len(k))

        tasks = []
        for i in idents:
            temp_idx = self.data.obs.loc[self.data.obs[self.group] == i, :]

            try:
                temp_expression = self.data.to_df().loc[temp_idx.index, target_genes[i]].transpose()

                tasks.append([i, temp_expression, ])
            except KeyError:
                continue

        with Pool(self.processes) as p:
            temp_res = p.map(calculate_tsi_of_dataframe, tasks)

        [self.genes_conserved_each_groups.update(i) for i in temp_res]

    def __get_gene_variable_between_groups__(self):
        u"""
        get genes that differential expressed between different groups
        :return:
        """
        logger.info("Find highly variable genes between groups")

        tasks = []
        for value in self.genes_conserved_each_groups.values():
            for gene, tsi in value.items():
                if tsi < self.tsi_low:
                    tasks.append([gene, self.genes_expression, ])

        with Pool(self.processes) as p:
            temp_res = p.map(calculate_tsi_of_dict, tasks)

        [self.genes_differential_expressed_between_groups.update(i) for i in temp_res]

    def dump(self, output: str):
        u"""

        :return:
        """
        self.__highly_variable_genes__()
        self.__get_conserved_genes_inside_idents__()
        self.__get_gene_variable_between_groups__()

        if not os.path.exists(output):
            os.makedirs(output)

        with open(os.path.join(output, "genes_conserved.json"), "w+") as w:
            json.dump(self.genes_conserved_each_groups, w, indent=4)

        with open(os.path.join(output, "genes_expression.json"), "w+") as w:
            json.dump(self.genes_expression, w, indent=4)

        with open(os.path.join(output, "genes_variable.json"), "w+") as w:
            json.dump(self.genes_differential_expressed_between_groups, w, indent=4)

        groups = sorted(self.genes_conserved_each_groups.keys())

        with open(os.path.join(output, "final.txt"), "w+") as w:
            w.write("group\tgene\ttsi_in_group\ttsi_between_group\t{0}\n".format("\t".join(groups)))

            for gene, tsi in self.genes_differential_expressed_between_groups.items():
                if tsi > self.tsi_high:

                    group = self.genes_most_expressed_group[gene]

                    tsi_in_group = self.genes_conserved_each_groups[group][gene]

                    if tsi_in_group > self.tsi_low:
                        continue

                    temp = [
                        group, gene,
                        str(tsi_in_group),
                        str(tsi)
                    ]

                    for i in groups:
                        try:
                            temp.append(str(self.genes_expression[gene][i]))
                        except KeyError:
                            temp.append("NA")

                    w.write("\t".join(temp) + "\n")


def load_from_csv(input_dir: str) -> (AnnData, AnnData):
    u"""
    load data from csv files
    :param input_dir:
    :return:
    """
    logger.info("read")
    mtx = pd.read_csv(os.path.join(input_dir, "normalized_counts.csv.gz"), index_col=0, engine="c")
    meta = pd.read_csv(os.path.join(input_dir, "meta.csv.gz"), index_col=0, engine="c")
    meta = meta.loc[meta.index, :]

    mtx = mtx.transpose()

    data = AnnData(mtx, obs=meta)
    data.obs = meta

    logger.info("enn")
    enn = EditedNearestNeighbours(n_jobs=10, return_indices=True)

    mtx_enn, group_enn, idx_enn = enn.fit_resample(mtx, meta["Stage"])

    data_enn = AnnData(mtx.iloc[list(idx_enn), :], meta.iloc[idx_enn, :])

    data_enn.obs = meta.iloc[idx_enn, :]

    logger.info("Repeated enn")
    renn = RepeatedEditedNearestNeighbours(n_jobs=10, return_indices=True)

    mtx_renn, group_renn, idx_renn = renn.fit_resample(mtx, meta["Stage"])

    data_renn = AnnData(mtx.iloc[list(idx_renn), :], meta.iloc[idx_renn, :])

    data_renn.obs = meta.iloc[idx_renn, :]

    return data, data_enn, data_renn


def make_dotplot(coord: pd.DataFrame, data: AnnData, filename: str):
    u"""
    make dotplot of tsne
    :param coord: pd.DataFrame
    :param data:
    :param filename:
    :return:
    """

    colors = {
        "Benign": "#3E6596",
        "I": "#EC7A21",
        "II": "#D73F47",
        "III": "#65A9A3",
        "IV": "#4A933E",
        "ADC": "#063373",
        "LCC": "#FECC1B",
        "Normal": "#73BDFF",
        "Other": "#DA6906",
        "SCC": "#A0C807",
        "LELC": "#6A0019",
        "CPI": "#C5000B",
        "LBT": "#0084D1",
        "Bronichi": "#FB290F",
        "Lymph node": "#488F17",
        "Tumor": "#FC8210"
    }

    coord = coord.loc[data.obs.index, :]
    coord["Stage"] = data.obs.loc[coord.index, "Stage"]

    fig, ax = plt.subplots()

    groups = coord["Stage"].unique()

    for group in groups:
        temp = coord.loc[coord["Stage"] == group, :]
        ax.scatter(x=temp["tSNE_1"], y=temp["tSNE_2"], c=colors[group], label=group, s=0.5)
    ax.legend(loc="best")

    plt.tight_layout()
    plt.savefig(filename, dpi=600)


if __name__ == '__main__':
    # data, data_enn, data_renn = load_from_csv(sys.argv[1])
    #
    # coord = pd.read_csv(os.path.join(sys.argv[1], "tsne.csv.gz"), index_col=0)
    #
    # make_dotplot(coord, data, os.path.join(sys.argv[2], "tsne.png"))
    # make_dotplot(coord, data_enn, os.path.join(sys.argv[2], "tsne_enn.png"))
    # make_dotplot(coord, data_renn, os.path.join(sys.argv[2], "tsne_renn.png"))
    #
    # GeneSignature(data, "Stage").dump(sys.argv[2])
    # GeneSignature(data_enn, "Stage").dump(os.path.join(sys.argv[2], "enn"))
    # GeneSignature(data_renn, "Stage").dump(os.path.join(sys.argv[2], "renn"))

    data = pd.read_csv("paga/normalized_counts.csv.gz", index_col=0)
    temp = data.to_dict()

    x, y = [], []
    for key, values in temp.items():

        key = key.split("_")[0]

        for i in values.values():
            x.append(key)
            y.append(i)

    temp = pd.DataFrame({"x": x, "y": y})

    import seaborn as sns
    sns.violinplot(x="x", y="y", hue="x", data=temp)
