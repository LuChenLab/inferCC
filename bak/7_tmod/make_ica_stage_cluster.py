#!/usr/bin/env python3
# -*- coding:utf- 8 -*-
u"""
@since 2019.06.08
"""
import json
import logging
import os
import warnings
from glob import glob
from multiprocessing import Pool
from subprocess import check_call, CalledProcessError


def warn(*args, **kwargs):
    pass

warnings.warn = warn


import matplotlib
import pandas as pd
from anndata import AnnData
from sklearn.decomposition import FastICA
from sklearn.linear_model import LassoLarsIC, ARDRegression
from tqdm import trange, tqdm

matplotlib.use("Agg")

from matplotlib import pyplot as plt
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.metrics import calinski_harabasz_score
from kneed import KneeLocator

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] - [%(name)s] - [%(levelname)s] - %(message)s')
logger = logging.getLogger(__name__)


__dir__ = os.path.abspath(os.path.dirname(__file__))
__RSCRIPT__ = os.path.join(__dir__, "make_tmod.R")


def load_from_csv(
        input_dir: str,
        counts_file: str="normalized_counts.csv.gz",
        low_expression=0
):
    u"""
    load data from csv files
    :param input_dir:
    :param counts_file:
    :param str
    :return:
    """
    logger.debug("Reading {0}".format(input_dir))

    input_file = os.path.join(input_dir, counts_file)

    # if not os.path.exists(input_file):
    #     input_file += ".gz"

    mtx = pd.read_csv(input_file, index_col=0)
    meta = pd.read_csv(os.path.join(input_dir, "meta.csv.gz"), index_col=0)
    meta = meta.loc[meta.index, :]

    # logger.info(mtx.shape)
    # filter low expressed genes
    genes_sum = [x / mtx.shape[1] > low_expression for x in mtx.sum(axis=1)]

    mtx = mtx.loc[genes_sum, :]

    # logger.info(mtx.shape)
    mtx = mtx.transpose()

    data = AnnData(mtx, obs=meta)
    data.obs = meta

    return data


def run_r(data):
    u"""
    execute r code
    :param data:
        - rds
        - prefix directory
        - cluster
        - gene
    :return:
    """
    rds, prefix, cluster, gene, group = data
    # print(cluster)

    rscript = """
        rds = "{0}"
        genes = strsplit("{1}", ",")[[1]]
        outfile = "{2}"
        title = "{3}"
        group = "{4}"
    """.format(
        rds,
        ",".join(gene),
        os.path.join(prefix, "{0}.png".format(cluster)),
        "cluster {0}".format(cluster),
        group
    )

    rscript += '''
    # scripts to make heatmap

    load <- function() {
        library(Seurat)
        library(ggplot2)
        library(stringr)
    }

    suppressMessages(load())

    if(str_detect(packageVersion("Seurat"), "^3.\\\\d.\\\\d")) {
        detach("package:Seurat")
        devtools::install("Seurat", version="2.3.4")
        library("Seurat")
    }


    obj <- readRDS(rds)


    p <- DoHeatmap(
        obj, 
        genes.use = genes, 
        disp.min = -2.5, 
        disp.max = 2.5, 
        group.by = group,
        group.order = NULL, 
        draw.line = TRUE, 
        col.low = "#FF00FF",
        col.mid = "#000000", 
        col.high = "#FFFF00", 
        slim.col.label = TRUE,
        remove.key = TRUE, 
        rotate.key = FALSE, 
        title = title, 
        do.plot = TRUE
    )



    height = length(genes) / 8

    if (height < 5) {
        height = 5
    } else if (height > 40) {
        height = 40
    }

    ggsave(filename = outfile, plot = p, width = 6, height = height, units = "in", dpi = 600)
    '''

    __dir__ = os.path.abspath(os.path.dirname(__file__))
    temp_r = os.path.join(__dir__, "tmp_rscript_{0}.R".format(cluster))
    with open(temp_r, "w+") as w:
        w.write(rscript)

    try:
        with open(os.devnull) as w:
            check_call("Rscript {0}".format(temp_r), shell=True, stdout=w, stderr=w)
    except CalledProcessError as err:
        print(err)
    finally:
        os.remove(temp_r)


def make_heatmap(data: str, group_by, rds, output, process):
    u"""
    make heatmap without scanpy, scanpy do work properly, and better make heatmap in same style
    just give up write one with matplotlib
    too much work for me

    just use Seurat instead
    :return:
    """
    with open(data) as r:
        data = json.load(r)

    tasks = []

    if not os.path.exists(output):
        os.makedirs(output)
    # rds, prefix, cluster, gene, group = data
    for k, v in data.items():
        tasks.append([rds, output, k, v, group_by])

    with Pool(processes=process) as p:
        p.map(run_r, tasks)


def estimate_best_k(
        data: pd.DataFrame,
        min_cluster=1,
        max_cluster=15,
        random_state=None,
        n_jobs=1,
        prefix=None
) -> (dict, int):
    u"""
    Using KMeans to cluster all the genes, and find best k through kneed
    :param data: data frame of gene expression matrix
    :param min_cluster: the min cluster
    :param max_cluster: the max cluster
    :param random_state: see sklearn.cluster.KMeans
    :param n_jobs: see sklearn.cluster.KMeans
    :param prefix: output prefix, if prefix is not None, make plots
    :return (dict, int): the km.labels_ of different cluster numbers, and estimated best k
    """

    sum_of_squared_distances = []
    max_cluster = min(max_cluster, data.shape[0])
    chs = []
    K = range(min_cluster, max_cluster + 1)
    res = {}
    for k in K:
        logger.debug("Estimate best k {0}/{1}".format(k, max_cluster))

        try:
            km = KMeans(n_clusters=k, random_state=random_state, n_jobs=n_jobs)
            km = km.fit(data)
            sum_of_squared_distances.append(km.inertia_)

            res[k] = list(km.labels_)
        except ValueError as err:
            logger.error(err)
            logger.error(k)
            logger.error(data.shape)
            break
        if k > 1:
            try:
                chs.append(calinski_harabasz_score(data, km.labels_))
            except ValueError as err:
                logger.error(err)
                chs.append(0)

        del km

    kn = KneeLocator(K, sum_of_squared_distances, curve='convex', direction='decreasing')

    if prefix is not None:
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(K, sum_of_squared_distances, 'bx-')
        ax.vlines(kn.knee, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
        plt.xlabel('k')
        plt.ylabel('Sum_of_squared_distances')
        plt.title('Elbow Method For Optimal k')
        plt.tight_layout()
        plt.xticks(K)
        plt.savefig("{0}_Elbow_for_k.png".format(prefix), dpi=600)

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(K[1:], chs, 'bx-')
        ax.vlines(kn.knee, plt.ylim()[0], plt.ylim()[1], linestyles='dashed')
        plt.xlabel('k')
        plt.ylabel('Calinski-Harabasz Index')
        plt.title('Calinski-Harabasz Index')
        plt.tight_layout()
        plt.xticks(K[1:])
        plt.savefig("{0}_calinski_harabasz.png".format(prefix), dpi=600)

    return res, kn.knee


def get_best_clusters(data: pd.DataFrame, km: dict, best_k: int):
    u"""
    construct clusters based on which k to use
    """
    clusters_of_best = {}
    for i, j in zip(data.index, km[best_k]):
        j = int(j)
        temp = clusters_of_best.get(j, [])
        temp.append(i)
        clusters_of_best[j] = temp

    return clusters_of_best


def estimate_ica(
        data: pd.DataFrame,
        n_comp=1,
        min_cluster=1,
        max_cluster=30,
        n_jobs=1,
        random_state=1,
        prefix=None
):
    u"""
    perform ICA and using KMeans to cluster, then using BIC to estimate the best components
    :param data:
    :param n_comp: the number of components
    :param min_cluster: see estimate_best_k
    :param max_cluster: see estimate_best_k
    :param prefix: see estimate_best_k
    :param n_jobs: see estimate_best_k
    :param random_state:
    :param doARD: whether to do ARDRegression
    :return:
    """
    logger.debug("ICA of {0}".format(n_comp))
    ica = FastICA(n_components=n_comp, random_state=random_state)
    S_ = ica.fit_transform(data)  # Reconstruct signals

    km, best_k = estimate_best_k(
        pd.DataFrame(S_),
        min_cluster=min_cluster,
        max_cluster=max_cluster,
        random_state=random_state,
        n_jobs=n_jobs,
        prefix=prefix
    )

    model_bic = LassoLarsIC(criterion='bic')
    try:
        model_bic.fit(data, km[best_k])
    except KeyError as err:
        print(err)

    return {"n_comp": n_comp, "ica": S_, "km": km, "best_k": best_k, "bic": model_bic.alpha_}


def perform_ica_on_data(
        data: AnnData,
        genes_use,
        n_pcs=50,
        min_cluster=1,
        max_cluster=30,
        prefix="",
        n_jobs=1,
        random_state=1,
        do_ard=False,
        group_by="Stage",
        rds=None
):
    u"""
    try to find the best number of components for ICA
    :param data:
    :param n_pcs: how many components to test in total
    :param group_spec: the markers by specific group
    :param min_cluster: see estimate_best_k
    :param max_cluster: see estimate_best_k
    :param prefix:
    :param n_jobs: see estimate_best_k
    :param random_state:
    :param do_ard: whether to do ARDRegression
    :param group_by:
    :param rds
    :return:
    """
    res = []
    mtx = data.to_df().transpose()

    genes_use = set(genes_use) & set(mtx.index)

    if genes_use:
        mtx = mtx.loc[genes_use, :]

    n_pcs = min(n_pcs, min(mtx.shape)) - 1
    for i in trange(1, n_pcs):
        try:
            res.append(estimate_ica(
                mtx,
                n_comp=i,
                min_cluster=min_cluster,
                max_cluster=min(max_cluster, min(mtx.shape) - 1),
                random_state=random_state,
                n_jobs=n_jobs,
                prefix=None
            ))
        except Exception as err:
            print(err)
            print("pass")

    bic_scores = [x["bic"] for x in res]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(range(1, len(bic_scores) + 1), bic_scores, "bx-")
    plt.xlabel("n_components")
    plt.ylabel("BIC")
    plt.tight_layout()
    plt.xticks(range(1, len(bic_scores) + 1))
    plt.savefig("{0}_bic.png".format(prefix), dpi=600)

    if do_ard:
        clf = ARDRegression(compute_score=True)
        for i in res:
            clf.fit(i["ica"], i["km"][i["best_k"]])
            i["ARD"] = clf.lambda_[0]

        scores = [x["ARD"] for x in res]

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(range(1, len(scores) + 1), scores, "bx-")
        plt.xlabel("n_components")
        plt.ylabel("ARD lambda_")
        plt.tight_layout()
        plt.xticks(range(1, len(scores) + 1))
        plt.savefig("{0}_ARD.png".format(prefix), dpi=600)

    # find best number of components
    # genes = {}
    # for idx, row in group_spec.iterrows():
    #     temp = genes.get(row["ident"], set())
    #     temp.add(idx)
    #     genes[row["ident"]] = temp

    idx = bic_scores.index(max(bic_scores))
    logger.debug("Estimate the best number of components is {0}".format(idx))

    clusters = get_best_clusters(mtx, res[idx]["km"], res[idx]["best_k"])
    # overlap_by_cluster, overlap_by_group = get_overlap_perc_of_most_variable_and_clusters(genes, clusters)

    output_json = "{0}_data.json".format(prefix)
    # print(output_json)
    with open(output_json, "w+") as w:
        json.dump(
            clusters,
            w, indent=4
        )

    if rds:
        make_heatmap(
            data=output_json,
            group_by=group_by,
            rds=rds,
            output=output_json.replace("data.json", "heatmap"),
            process=n_jobs
        )

    # overlap_by_cluster *= 100
    # overlap_by_group *= 100
    #
    # fig, ax = plt.subplots(figsize=(6, 6))
    # sns.heatmap(overlap_by_cluster, ax=ax, annot=True, fmt=".1f")
    # plt.tight_layout()
    # plt.savefig("{0}_heatmap_by_cluster.png".format(prefix), dpi=600)
    #
    # fig, ax = plt.subplots(figsize=(6, 6))
    # sns.heatmap(overlap_by_group, ax=ax, annot=True, fmt=".1f")
    # plt.tight_layout()
    # plt.savefig("{0}_heatmap_by_group.png".format(prefix), dpi=600)
    #
    # with open("{0}_data.pickle".format(prefix), "wb+") as w:
    #     pickle.dump(res, w)


def run(path, paga, output, rds, n_jobs=20):
    module = pd.read_excel(path)
    data = load_from_csv(paga)

    for i in tqdm(module["gene_module_id"].unique()):
        genes_use = set(module.loc[module["gene_module_id"] == i, "gene"])

        output_dir = os.path.join(output, "stage_{0}".format(i))

        try:
            perform_ica_on_data(
                data=data,
                genes_use=genes_use,
                prefix=output_dir,
                rds=rds,
                n_jobs=n_jobs
            )
        except Exception as err:
            print(err)
            print(output_dir)


def main(path, n_jobs):
    u"""
    main
    :return:
    """
    files = glob(os.path.join(os.path.abspath(path), "T_cells_SCC/mfuzz_gene_module/results.xlsx"))
    # tasks = []
    for f in sorted(files):
        print(f)
        dirname = os.path.dirname(os.path.dirname(f))
        output = os.path.join(dirname, "tmod_stage")

        if not os.path.exists(output):
            os.makedirs(output)

        rds = glob(os.path.join(dirname, "*.rds"))

        for r in rds:
            if "monocle" not in r and \
                "slingshot" not in r and \
                    "wgcna" not in r.lower():
                run(
                    f,
                    os.path.join(dirname, "paga"),
                    output=output,
                    rds=r,
                    n_jobs=n_jobs
                )


if __name__ == '__main__':
    from fire import Fire
    Fire(main)