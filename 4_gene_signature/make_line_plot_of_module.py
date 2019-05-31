#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
make line plot and pie chart of module gene expression
"""

import os
import logging
import json
import matplotlib
from glob import glob
from itertools import combinations
import math
matplotlib.use("Agg")
import scanpy as sc

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import pandas as pd
from subprocess import check_call, CalledProcessError
from scipy.stats import mannwhitneyu, zscore
from matplotlib_venn import venn2, venn3
from multiprocessing import Pool
from tqdm import tqdm


logging.basicConfig(level=logging.INFO, format='[%(asctime)s] - [%(name)s] - [%(levelname)s] - %(message)s')
logger = logging.getLogger(__name__)


__date__ = "2019.05.21"
__author__ = "ygidtu"
__email__ = "ygidtu@gmail.com"
__dir__ = os.path.abspath(os.path.dirname(__file__))


COLORS = {
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


def load_gene_expression(
        counts: str,
        meta: str,
        genes_use=None,
        group_by=None,
        groups=None
):
    u"""
    load expression data by meta data
    :param counts: path to counts csv
    :param meta: path to meta data
    :param genes_use: which genes to use
    :param group_by: column in meta
    :param groups: str, items of group_by, eg: first,second
    :return: counts and two list -> groups and data
    """
    counts = pd.read_csv(counts, index_col=0, engine="c")
    meta = pd.read_csv(meta, index_col=0, engine="c")

    if genes_use:
        counts = counts.iloc[set(genes_use) & set(counts.index), :]

    res = []
    if group_by and groups:
        groups = groups.split(",")

        for group in groups:
            temp_meta = meta.iloc[meta[groups] == group, :]
            temp_data = counts.loc[:, temp_meta.index]

            res.append(temp_data)
    else:
        res.append(res)

    return groups, res


class GeneExpression(object):
    u"""
    read gene expression and meta info
    """

    def __init__(self, counts: str, meta: str, group_by="res.0.6"):
        u"""
        read data
        :param counts:
        :param meta:
        :param group_by:
        """

        self.counts = pd.read_csv(counts, index_col=0, engine="c")
        self.meta = pd.read_csv(meta, index_col=0, engine="c")
        self.group_by = group_by

    def get_gene_exp(self, gene):
        u"""
        get gene expression
        :param gene:
        :param group:
        :return:
        """
        groups = self.meta[self.group_by].unique()

        if gene not in self.counts.index:
            return None

        res = {}
        for group in groups:
            idx = self.meta.loc[self.meta[self.group_by] == group, :]

            res[group] = list(self.counts.loc[gene, idx.index])

        return res

    def get_gene_zscore(self, gene):
        u"""

        :param gene:
        :return:
        """
        res = {}
        for key, value in self.get_gene_exp(gene).items():
            res[key] = zscore(value)

        return res


def _make_line_plot_(args):
    # data, clt, ax = args
    # fig, ax = plt.subplots(figsize=(12, 6))
    # plt.tight_layout()

    data_dict, out_dir, key, columns = args

    data = pd.DataFrame()

    for c1, c2 in zip(data_dict, columns):
        data[c2] = c1

    del data_dict

    # print(data.head())
    # make plots
    # data["zscore"] = np.log10(data["zscore"] + 1)
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.violinplot(x="cluster", y="zscore", data=data, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("zscore")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "{0}_violin.png".format(key)), dpi=600)

    fig, ax = plt.subplots(figsize=(12, 6))
    data["value"] = np.log10(data["value"] + 1)
    sns.violinplot(x="cluster", y="value", data=data, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("log10(normalized counts + 1)")
    # plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "{0}_violin_value.png".format(key)), dpi=600)

    # clusters = sorted(data["cluster"].unique())
    # nrow = math.ceil(len(clusters) / ncol)

    # fig = plt.figure(figsize=(6 * ncol, 4 * nrow))
    # gs = gridspec.GridSpec(nrows=nrow, ncols=ncol, figure=fig, wspace=.1,  hspace=.3)

    # max_y = max(data["value"])

    # tasks = []
    # for c_idx, clt in enumerate(clusters):
    #     curr_row = c_idx // ncol
    #     curr_col = c_idx % ncol
    #
    #     curr_ax = plt.subplot(gs[curr_row, curr_col])
    #     # tasks.append([temp, clt, curr_ax])
    #     temp = data.loc[data["cluster"] == clt, :]
    #
    #     sns.lineplot(x="idx", y="value", hue="gene", data=temp, ax=curr_ax, legend=False)
    #     curr_ax.set_title(clt)
    #     curr_ax.spines['top'].set_color('none')
    #     curr_ax.spines['right'].set_color('none')
    #     curr_ax.set_ylabel('')
    #     curr_ax.set_xlabel('')
    #     curr_ax.set_ylim(bottom=0, top=max_y)
    #
    #     plt.xticks([])

    data = data.loc[:, ["gene", "cluster", "zscore"]]
    data = data.groupby(["gene", "cluster"]).mean().reset_index()

    height = 6 + data["gene"].unique().shape[0] // 6
    if height > 45:
        height = 45
    fig, ax = plt.subplots(figsize=(12, height))
    sns.lineplot(x="cluster", y="zscore", hue="gene", ax=ax, data=data)
    ax.set_xticks(data["cluster"].unique())
    ax.set_xlabel("")
    ax.legend(loc=9, ncol=8, bbox_to_anchor=(0.5, -0.05))
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "{0}_lineplot.png".format(key)), dpi=600)
    plt.close(fig)


def _read_module_(target_dir, method):
    u"""
    read module from json and txt file
    :param target_dir:
    :param key:
    :return:
    """
    data = {}
    out_files = glob(os.path.join(target_dir, "{0}_*_data.json".format(method)))

    # read data from json files
    for j in out_files:
        key = os.path.basename(j).split("_")[1]

        with open(j) as r:
            temp = json.load(r)

        for k, v in temp.items():
            if len(v) >= 3:
                data["{0}_{1}".format(key, k)] = v

    # read data from txt files
    out_files = glob(os.path.join(target_dir, "{0}_*_data.txt".format(method)))

    for j in out_files:
        logger.info(j)

        key = os.path.basename(j).split("_")[1]

        temp_data = {}
        with open(j) as r:
            for line in r:
                lines = line.split()
                temp = temp_data.get(lines[1], [])
                temp.append(lines[0])
                temp_data[lines[1]] = temp

        for cluster, genes in temp_data.items():
            data["{0}_{1}".format(key, cluster)] = genes

    return data


def plot(counts, meta, group_by, target_dir, n_jobs=10):
    #
    data_labels = ["Imbalanced", "ENN", "RENN"]

    counts = GeneExpression(counts=counts, meta=meta, group_by=group_by)

    for i in data_labels:
        # if i in ["Imbalanced", "ENN",]:
        #     continue

        out_dir = os.path.join(target_dir, "violin/{0}".format(i))

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # format data
        data = _read_module_(target_dir, i)

        tasks = []
        # key -> method_cluster
        # value -> genes
        # cluster -> gene module cluster id
        # temp -> gene expression {cell cluster id: [expression]}
        for key, value in data.items():

            temp_data = [[], [], [], [], []]

            for v in value:

                gene_exp = counts.get_gene_exp(v)

                if gene_exp is None:
                    continue

                for cluster, temp in gene_exp.items():
                    if temp:
                        temp1 = zscore(temp)

                        for z_idx, t in enumerate(temp):
                            temp_data[0].append(cluster)
                            temp_data[1].append(t)
                            temp_data[2].append(v)
                            temp_data[3].append(z_idx)
                            temp_data[4].append(temp1[z_idx])

            # data, out_dir, key, ncol, columns = args
            tasks.append([temp_data, out_dir, key, ["cluster", "value", "gene", "idx", "zscore"]])

        logger.info(len(tasks))
        with Pool(n_jobs) as p:
            p.map(_make_line_plot_, tasks)


def _make_heatmap_(args):
    u"""
    make heatmap using seurat
    :param args:
        - genes
        - title
        - output file
        - rds
        - group by
    :return:
    """

    genes, title, output, rds, group_by = args

    temp_genes = os.path.join(output + ".temp_genes")
    temp_r = os.path.join(output + ".r")

    with open(temp_genes, "w+") as w:
        w.write("\n".join(genes))

    rscript = """
    genes = "{0}"
    title = "{1}"
    output = "{2}"
    rds = "{3}"
    group_by = "{4}"
    """.format(temp_genes, title, output, rds, group_by)

    rscript += """
    load <- function(){
        library(ggplot2)
        library(Seurat)
        library(stringr)
    }
    suppressPackageStartupMessages(load())
      
    if(str_detect(packageVersion("Seurat"), "^3.\\\\d.\\\\d")) {
        detach("package:Seurat")
        devtools::install("Seurat", version="2.3.4")
        library("Seurat")
    }
    
    genes = read.table(genes)
    
    obj <- readRDS(rds)
    
    p <- DoHeatmap(
        obj, 
        group.by = group_by, 
        genes.use = genes[,1], 
        cex.col=0,
        slim.col.label = TRUE, 
        remove.key = TRUE,
        do.plot = F
    )
    
    height = length(genes[,1]) / 8
    
    if(height < 5) {
        height = 5
    } else if(height > 40){
        height = 40
    }
    
    ggsave(filename = output,
           plot = p,
           width = 6,
           height = height,
           limitsize = F,
           units = "in",
           dpi = 600)
    """
    with open(temp_r, "w+") as w:
        w.write(rscript)

    try:
        check_call("Rscript " + temp_r, shell=True)
    except CalledProcessError as err:
        logger.error(err)
    os.remove(temp_genes)
    os.remove(temp_r)


def make_heatmap(target_dir, output_dir, rds, n_jobs=10):
    data_labels = ["Imbalanced", "ENN", "RENN"]
    group_by = {
        "stage": "Stage",
        "cluster": "res.0.6"
    }

    tasks = []
    for i in data_labels:
        data = _read_module_(target_dir, i)

        out_dir = os.path.join(output_dir, i)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        for cluster, genes in data.items():
            if len(genes) >= 3:
                for key, group in group_by.items():
                    tasks.append([
                        genes,
                        cluster,
                        os.path.join(out_dir, "{0}_{1}_{2}.png".format(i, cluster, key)),
                        rds,
                        group
                    ])

    with Pool(n_jobs) as p:
        p.map(_make_heatmap_, tasks)


def _make_venn_(args):
    labels, genes, outfile = args

    fig, ax = plt.subplots()
    if len(labels) == 2:
        venn2(genes, set_labels=labels, ax=ax)
    elif len(labels) == 3:
        venn3(genes, set_labels=labels, ax=ax)
    else:
        logger.warning("length of different sets not in 2,3")

    plt.savefig(outfile, dpi=600)
    plt.close(fig)


def venn(target_dir: str, methods_num=4, n_jobs=10):
        u"""
        make upset plot
        :param target_dir:
        :param methods_num:
        :param group_spec:
        :param groups
        :return:
        """

        data_labels = ["Imbalanced", "ENN", "RENN"]

        # group_spec = pd.read_excel(group_spec, index_col=0)

        # select qualified markers
        # idx = []
        # for i, j in zip(group_spec["p_val_adj"] < pvalue, group_spec["avg_logFC"] > logfc):
        #     idx.append(i and j)
        #
        # if not isinstance(groups, tuple):
        #     groups = [groups]

        group_genes = {
            "Known genes": {
                "EGFR",
                "CTLA4",
                "CD28",
                "ALKAL1",
                "ALKAL2",
                "KRAS",
                "ROS1",
                "BRAF",
                "MET",
                "ERBB2",
                "PDCD1",
                "TP53",
                "CD274",
                "TNFRSF9"
            }
        }
        # for group in groups:
        #     group_genes["Cell cluster {0}".format(group)] = set(
        #         group_spec.loc[[str(x) == str(group) for x in group_spec["ident"]], "gene"]
        #     )

        tasks = []
        for i in data_labels:
            out_dir = os.path.join(target_dir, "venn/{0}".format(i))

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            logger.info("Make upset plot of {0}".format(i))
            data = {}
            out_files = glob(os.path.join(target_dir, "{0}_*_data.json".format(i)))

            for j in out_files:
                key = os.path.basename(j).split("_")[1]

                with open(j) as r:
                    temp = json.load(r)

                for k, v in temp.items():
                    if len(v) > 3:
                        data["{0}_{1}".format(key, k)] = v

            out_files = glob(os.path.join(target_dir, "{0}_*_data.txt".format(i)))

            for j in out_files:
                key = os.path.basename(j).split("_")[1]

                temp_data = {}
                with open(j) as r:
                    for line in r:
                        lines = line.split()
                        temp = temp_data.get(lines[1], [])
                        temp.append(lines[0])
                        temp_data[lines[1]] = temp

                for k, v in temp_data.items():
                    if len(v) > 3:
                        data["{0}_{1}".format(key, k)] = v

            # for_upset, for_upset_data, jaccard = [], [], []

            for k in range(1, methods_num + 1):
                for x in combinations(data.keys(), k):

                    if k == 1:
                        temp_genes = list(group_genes.values())
                        temp_labels = list(group_genes.keys())

                        temp_labels.append(x[0])
                        temp_genes.append(set(data[x[0]]))

                        tasks.append([temp_labels, temp_genes, os.path.join(out_dir, "{0}.png".format(x[0]))])

                        # fig, ax = plt.subplots()
                        # if len(temp_labels) == 2:
                        #     venn2(temp_genes, set_labels=temp_labels, ax=ax)
                        # elif len(temp_labels) == 3:
                        #     venn3(temp_genes, set_labels=temp_labels, ax=ax)
                        # else:
                        #     logger.warning("length of different sets not in 2,3")
                        #
                        # plt.savefig(os.path.join(out_dir, "{0}.png".format(x[0])), dpi=600)
                        # plt.close(fig)

            #         if len(set([i.split("_")[0] for i in x])) == k:
            #
            #             temp_data = None
            #             jaccard.append(set())
            #
            #             for y in x:
            #                 jaccard[-1] |= set(data[y])
            #                 if temp_data is None:
            #                     temp_data = set(data[y])
            #                 else:
            #                     temp_data &= set(data[y])
            #
            #             for_upset.append(x)
            #             for_upset_data.append(len(temp_data))
            #
            # with open(os.path.join(target_dir, "{0}_jaccard.txt".format(i)), "w+") as w:
            #     for idx in range(len(for_upset)):
            #         d, j, k = for_upset[idx], for_upset_data[idx], jaccard[idx]
            #
            #         w.write("{0}\t{1}\n".format("|".join(list(d)), j / len(k)))

        with Pool(n_jobs) as p:
            p.map(_make_venn_, tasks)


def sc_violin(target_dir, output_dir):
    u"""
    make stacked violin plot of scanpy style
    :param target_dir:
    :param output_dir:
    :return:
    """

    counts = pd.read_csv(os.path.join(target_dir, "normalized_counts.csv.gz"), index_col=0)
    meta = pd.read_csv(os.path.join(target_dir, "meta.csv.gz"), index_col=0)

    data = sc.AnnData(counts.transpose(), obs=meta)
    data.obs = meta

    genes = {
        "EGFR",
        "CTLA4",
        "CD28",
        "ALKAL1",
        "ALKAL2",
        "KRAS",
        "ROS1",
        "BRAF",
        "MET",
        "ERBB2",
        "PDCD1",
        "TP53",
        "CD274",
        "TNFRSF9"
    }

    data.obs["cluster"] = [str(x) for x in data.obs["res.0.6"]]
    used_genes = sorted(list(set(data.var_names) & set(genes)))

    fig = plt.figure(figsize=(6, 12))
    gs = gridspec.GridSpec(nrows=len(used_genes), ncols=1)

    axes = []
    for i in range(len(used_genes)):
        axes.append(plt.subplot(gs[i, :]))
        sc.pl.violin(
            data,
            used_genes[i],
            groupby="cluster",
            swap_axes=True,
            ax=axes[-1],
            show=False
        )
        axes[-1].set_ylabel(
            used_genes[i],
            rotation=0,
            labelpad=25,
            va="center"
        )
        axes[-1].set_xlabel("")

        if i < len(used_genes) - 1:
            axes[-1].set_xticks([])

    plt.savefig(os.path.join(output_dir, "known_genes_stacked_violin_cluster.png"))

    fig = plt.figure(figsize=(6, 12))
    gs = gridspec.GridSpec(nrows=len(used_genes), ncols=1)

    axes = []
    for i in range(len(used_genes)):
        axes.append(plt.subplot(gs[i, :]))
        sc.pl.violin(
            data,
            used_genes[i],
            groupby="Stage",
            swap_axes=True,
            ax=axes[-1],
            show=False
        )
        axes[-1].set_ylabel(
            used_genes[i],
            rotation=0,
            labelpad=25,
            va="center"
        )
        axes[-1].set_xlabel("")

        if i < len(used_genes) - 1:
            axes[-1].set_xticks([])

    plt.savefig(os.path.join(output_dir, "known_genes_stacked_violin_stage.png"))


if __name__ == '__main__':

    from fire import Fire
    Fire()

    pass
