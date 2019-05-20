#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Scripts to perform kmeans and ICA on marker genes,
trying to find out the best group of gene signatures
"""
import argparse as ap
import json
import logging
import os
import pickle
import sys
from glob import glob
from subprocess import check_call, CalledProcessError
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")

from sklearn.cluster import KMeans
from kneed import KneeLocator
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from anndata import AnnData
from imblearn.under_sampling import EditedNearestNeighbours, RepeatedEditedNearestNeighbours
from sklearn.metrics import calinski_harabasz_score
from tqdm import tqdm
from sklearn.decomposition import FastICA
from sklearn.linear_model import LassoLarsIC, ARDRegression

from itertools import combinations
from upsetplot import from_memberships, plot


logging.basicConfig(level=logging.INFO, format='[%(asctime)s] - [%(name)s] - [%(levelname)s] - %(message)s')
logger = logging.getLogger(__name__)


__date__ = "2019.05.16"
__author__ = "ygidtu"
__email__ = "ygidtu@gmail.com"
__dir__ = os.path.abspath(os.path.dirname(__file__))


def load_from_csv(input_dir: str, counts_file: str="normalized_counts.csv.gz") -> (AnnData, AnnData):
    u"""
    load data from csv files
    :param input_dir:
    :param counts_file:
    :return:
    """
    logger.info("Reading {0}".format(input_dir))
    mtx = pd.read_csv(os.path.join(input_dir, counts_file), index_col=0, engine="c")
    meta = pd.read_csv(os.path.join(input_dir, "meta.csv.gz"), index_col=0, engine="c")
    meta = meta.loc[meta.index, :]

    mtx = mtx.transpose()

    data = AnnData(mtx, obs=meta)
    data.obs = meta

    logger.info("Perform ENN")
    enn = EditedNearestNeighbours(n_jobs=10, return_indices=True)

    mtx_enn, group_enn, idx_enn = enn.fit_resample(mtx, meta["Stage"])

    data_enn = AnnData(mtx.iloc[list(idx_enn), :], meta.iloc[idx_enn, :])

    data_enn.obs = meta.iloc[idx_enn, :]

    logger.info("Perform RENN")
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
    print(cluster)

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
    chs = []
    K = range(min_cluster, max_cluster + 1)
    res = {}
    for k in tqdm(K):
        logger.debug("Estimate best k {0}/{1}".format(k, max_cluster))

        km = KMeans(n_clusters=k, random_state=random_state, n_jobs=n_jobs)
        km = km.fit(data)
        sum_of_squared_distances.append(km.inertia_)

        res[k] = list(km.labels_)

        if k > 1:
            chs.append(calinski_harabasz_score(data, km.labels_))

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


def get_overlap_perc_of_most_variable_and_clusters(most_variable_genes, clusters):
    u"""
    get overlap percentage between the new cluster and the grouped markers
    :param most_variable_genes: dict {group ident: [markers...]}
    :param clusters: dict {cluster id: [genes...]}
    :return:
    """
    overlap_perc_by_cluster = {}
    overlap_perc_by_group = {}
    for group, marker in most_variable_genes.items():
        for cluster, genes in clusters.items():
            temp = overlap_perc_by_cluster.get(group, {})
            temp[cluster] = len(set(marker) & set(genes)) / len(genes)
            overlap_perc_by_cluster[group] = temp

            temp = overlap_perc_by_group.get(group, {})
            temp[cluster] = len(set(marker) & set(genes)) / len(marker)
            overlap_perc_by_group[group] = temp

    return pd.DataFrame(overlap_perc_by_cluster), pd.DataFrame(overlap_perc_by_group)


def perform_kmeans_on_data(
        data: AnnData,
        group_spec: pd.DataFrame,
        prefix,
        random_state=1,
        n_jobs=1,
        min_cluster=2,
        max_cluster=20,
        group_by="Stage",
        rds=None
) -> None:
    u"""
    peform all the calculation and plotting on single AnnData
    :param data:
    :param group_spec:
    :param prefix:
    :param random_state:
    :param n_jobs:
    :param min_cluster:
    :param max_cluster:
    :param group_by:
    :param rds:
    :return:
    """

    genes = {}
    for idx, row in group_spec.iterrows():
        temp = genes.get(row["ident"], set())
        temp.add(idx)
        genes[row["ident"]] = temp

    mtx = data.to_df().transpose()
    mtx = mtx.loc[set(group_spec.index) & set(mtx.index), :]

    km, best_k = estimate_best_k(
        mtx,
        random_state=random_state, n_jobs=n_jobs,
        max_cluster=max_cluster, min_cluster=min_cluster,
        prefix=prefix
    )

    clusters = get_best_clusters(mtx.loc[set(group_spec.index) & set(mtx.index), :], km, best_k)

    overlap_by_cluster, overlap_by_group = get_overlap_perc_of_most_variable_and_clusters(genes, clusters)

    output_json = "{0}_data.json".format(prefix)
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

    overlap_by_cluster *= 100
    overlap_by_group *= 100

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(overlap_by_cluster, ax=ax, annot=True, fmt=".1f")
    plt.tight_layout()
    plt.savefig("{0}_heatmap_by_cluster.png".format(prefix), dpi=600)

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(overlap_by_group, ax=ax, annot=True, fmt=".1f")
    plt.tight_layout()
    plt.savefig("{0}_heatmap_by_group.png".format(prefix), dpi=600)


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
    logger.info("ICA of {0}".format(n_comp))
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
    model_bic.fit(data, km[best_k])

    return {"n_comp": n_comp, "ica": S_, "km": km, "best_k": best_k, "bic": model_bic.alpha_}


def perform_ica_on_data(
        data: AnnData,
        group_spec: pd.DataFrame,
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

    genes_use = set(group_spec.index) & set(mtx.index)

    if genes_use:
        mtx = mtx.loc[genes_use, :]

    for i in range(1, n_pcs + 1):
        res.append(estimate_ica(
            mtx,
            n_comp=i,
            min_cluster=min_cluster,
            max_cluster=max_cluster,
            random_state=random_state,
            n_jobs=n_jobs,
            prefix=None
        ))

    bic_scores = [x["bic"] for x in res]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(range(1, n_pcs + 1), bic_scores, "bx-")
    plt.xlabel("n_components")
    plt.ylabel("BIC")
    plt.tight_layout()
    plt.xticks(range(1, n_pcs + 1))
    plt.savefig("{0}_bic.png".format(prefix), dpi=600)

    if do_ard:
        clf = ARDRegression(compute_score=True)
        for i in res:
            clf.fit(i["ica"], i["km"][i["best_k"]])
            i["ARD"] = clf.lambda_[0]

        scores = [x["ARD"] for x in res]

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(range(1, n_pcs + 1), scores, "bx-")
        plt.xlabel("n_components")
        plt.ylabel("ARD lambda_")
        plt.tight_layout()
        plt.xticks(range(1, n_pcs + 1))
        plt.savefig("{0}_ARD.png".format(prefix), dpi=600)

    # find best number of components
    genes = {}
    for idx, row in group_spec.iterrows():
        temp = genes.get(row["ident"], set())
        temp.add(idx)
        genes[row["ident"]] = temp

    idx = bic_scores.index(max(bic_scores))
    logger.info("Estimate the best number of components is {0}".format(idx))

    clusters = get_best_clusters(mtx, res[idx]["km"], res[idx]["best_k"])
    overlap_by_cluster, overlap_by_group = get_overlap_perc_of_most_variable_and_clusters(genes, clusters)

    output_json = "{0}_data.json".format(prefix)
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

    overlap_by_cluster *= 100
    overlap_by_group *= 100

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(overlap_by_cluster, ax=ax, annot=True, fmt=".1f")
    plt.tight_layout()
    plt.savefig("{0}_heatmap_by_cluster.png".format(prefix), dpi=600)

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(overlap_by_group, ax=ax, annot=True, fmt=".1f")
    plt.tight_layout()
    plt.savefig("{0}_heatmap_by_group.png".format(prefix), dpi=600)

    with open("{0}_ica_temp.pickle".format(prefix), "wb+") as w:
        pickle.dump(res, w)


def perform_mfuzz(
        data: str,
        output: str,
        genes,
        rds:str,
        group_by="Stage",
        random_state=1,
):
    u"""
    perform mfuzz
    :param data: path to gzipped file
    :param output: output prefix
    :param genes: comman sperated gene list
    :param group_by:
    :param random_state:
    :param rds:
    :return:
    """
    rscript = """
        random_state = {0}
        data = "{1}"
        genes = "{2}"
        output = "{3}"
        group_by = "{4}"
        rds = "{5}"
        """.format(
            random_state,
            os.path.abspath(data),
            ",".join(genes),
            os.path.abspath(output),
            group_by,
            rds
        )

    rscript += """
    library(ggplot2)
    library(Mfuzz)
    library(Seurat)
    library(stringr)
      
    if(str_detect(packageVersion("Seurat"), "^3.\\\\d.\\\\d")) {
        detach("package:Seurat")
        devtools::install("Seurat", version="2.3.4")
        library("Seurat")
    }
    
    set.seed(random_state)
    obj <- readRDS(rds)
    r = gzfile(data)
    data = read.csv(r, row.names = 1)
    
    genes = strsplit(genes, ",")[[1]]
       
    expr = ExpressionSet(
        as.matrix(data[intersect(genes, rownames(data)),])
    )
    
    cl = mfuzz(expr, c=16, m=1.25)
    
    res = as.data.frame(cl$cluster)
    write.table(res, file = paste0(output, "_data.txt"), quote=F, col.name=F, sep = "\t")

    groups = unique(cl$cluster)
    
    image_dir = paste(output, "heatmap", sep = "_")
    dir.create(image_dir, showWarnings = F)
    
    for(i in groups) {
        temp = cl$cluster[cl$cluster == i]
        p <- DoHeatmap(
            obj, 
            group.by = group_by, 
            genes.use = names(temp), 
            cex.col=0,
            slim.col.label = TRUE, 
            remove.key = TRUE,
            do.plot = F
        )
        
        height = length(temp) / 8
        
        if(height < 5) {
            height = 5
        } else if(height > 40){
            height = 40
        }
        
        ggsave(filename = paste(image_dir, paste0(i, ".png"), sep = "/"),
               plot = p,
               width = 6,
               height = height,
               limitsize = F,
               units = "in",
               dpi = 600)
    }
    """

    temp = os.path.join(__dir__, "tmp_r.R")

    with open(temp, "w+") as w:
        w.write(rscript)

    try:
        check_call("Rscript {0}".format(temp), shell=True)
    except CalledProcessError as err:
        logger.error(err)
    finally:
        os.remove(temp)


def perform_wgcna(
        data: str,
        output: str,
        genes,
        rds,
        group_by = "Stage",
        random_state=1
):
    u"""

    :param data:
    :param output:
    :param genes:
    :param group_by:
    :param random_state:
    :param rds
    :return:
    """

    rscript = """
    random_state = {0}
    data = "{1}"
    genes = "{2}"
    output = "{3}"
    rds = "{4}"
    group_by = "{5}"
    """.format(
        random_state,
        os.path.abspath(data),
        ",".join(genes),
        os.path.abspath(output),
        os.path.abspath(rds),
        group_by
    )

    rscript += """
    library(WGCNA)
    library(ggplot2)
    library(Seurat)

    enableWGCNAThreads()
     
    set.seed(random_state)
        
    r = gzfile(data)
    data = read.csv(r, row.names = 1)
    
    genes = strsplit(genes, ",")[[1]]
    
    mtx = data[intersect(genes, rownames(data)),]
    
    
    sampleTree = hclust(dist(mtx), method = "average");
    png(file = paste0(output, "_sampleClustering.png"), width = 12, height = 9, res = 600, units = "in")
    par(cex = 0.6);
    par(mar = c(0,4,2,0))
    plot(sampleTree, 
         main = "Sample clustering to detect outliers", 
         sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
    dev.off()
    
    # 作图选择power
    # Call the network topology analysis function
    sft = pickSoftThreshold(t(mtx), verbose = 5)
    # Plot the results:
    
    png(file = paste0(output, "_softthreshold.png"), width = 12, height = 9, res = 600, units = "in")
    par(mfrow = c(1,2));
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=sft$powerEstimate,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.90,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=sft$powerEstimate, cex=cex1,col="red")
    dev.off()
    
    ## make network
    net = blockwiseModules(t(mtx), power = sft$powerEstimate,
                           TOMType = "unsigned", minModuleSize = 30,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = T,
                           saveTOMFileBase = paste0(output, "_wcgna"),
                           verbose = 3)
    
    
    res = as.data.frame(net$colors)
    write.table(res, file = paste0(output, "_data.txt"), quote=F, col.name=F, sep = "\t")
    
    
    # open a graphics window
    png(file = paste0(output, "_clustering.png"), width = 12, height = 9, res = 600, units = "in")
    # Convert labels to colors for plotting
    mergedColors = labels2colors(net$colors)
    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
                        
    dev.off()
    
    
    obj <- readRDS(rds)
    
    groups = unique(res[, 1])
    
    image_dir = paste(output, "heatmap", sep = "_")
    dir.create(image_dir, showWarnings = F)
    
    for(i in groups) {
        temp = res[res[, 1] == i, , drop=F]
        p <- DoHeatmap(
            obj, 
            group.by = group_by, 
            genes.use = rownames(temp), 
            cex.col=0,
            slim.col.label = TRUE, 
            remove.key = TRUE,
            do.plot = F
        )
        
        height = nrow(temp) / 8
        
        if(height < 5) {
            height = 5
        } else if(height > 40){
            height = 40
        }
        
        ggsave(filename = paste(image_dir, paste0(i, ".png"), sep = "/"),
               plot = p,
               width = 6,
               height = height,
               limitsize = F,
               units = "in",
               dpi = 600)
    }
    """
    temp = os.path.join(__dir__, "tmp_r.R")

    with open(temp, "w+") as w:
        w.write(rscript)

    try:
        check_call("Rscript {0}".format(temp), shell=True)
    except CalledProcessError as err:
        logger.error(err)
        exit(err)

    os.remove(temp)


def make_upset_plot(data_labels, target_dir):
    u"""
    make upset plot
    :param data_labels:
    :param target_dir
    :return:
    """

    for i in data_labels:
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

        for_upset, for_upset_data = [], []

        for i in range(1, len(data_labels) + 1):
            for j in combinations(data.keys(), i):
                if len(set([x.split("_")[0] for x in j])) == i:

                    temp_data = None

                    for k in j:
                        if temp_data is None:
                            temp_data = set(data[k])
                        else:
                            temp_data &= set(data[k])

                    for_upset.append(j)
                    for_upset_data.append(len(temp_data))

        for_upset = from_memberships(for_upset, data=for_upset_data)

        plot(for_upset)
        plt.savefig(os.path.join(target_dir, "{0}_upset.png".format(i)), dpi=600)


def command_line() -> ap.Namespace:
    u"""
    construct command line parameters
    :return:
    """
    parser = ap.ArgumentParser(description="Script to make gene signature")
    parser.add_argument(
        "-i",
        help="Path to input directory, the direcotry contains csv files from seurat",
        required=True,
        type=str
    )
    parser.add_argument(
        "-x",
        help="Path to xlsx file, contains the markers",
        required=True,
        type=str
    )
    parser.add_argument(
        "-o",
        help="Path to output directory",
        required=True,
        type=str
    )
    parser.add_argument(
        "-p",
        help="How many processes to use (default: %(default)s)",
        default=1,
        type=int
    )
    parser.add_argument(
        "--max-cluster",
        help="The max clusters for KMeans (default: %(default)s)",
        default=30,
        type=int
    )
    parser.add_argument(
        "--random-state",
        help="The random state to KMeans (default: %(default)s)",
        default=1,
        type=int
    )
    parser.add_argument(
        "--pvalue",
        help="The threshold to select marker genes (default: %(default)s)",
        default=0.05,
        type=float
    )
    parser.add_argument(
        "--logfc",
        help="The threshold to select marker genes (default: %(default)s)",
        default=0,
        type=float
    )
    parser.add_argument(
        "--n-comps",
        help="How many components to test in ICA (default: %(default)s)",
        default=50,
        type=int
    )
    parser.add_argument(
        "--ard",
        help="Whether to calculate the ARDRegression, but not used in estimate the best number of components",
        action="store_true"
    )
    parser.add_argument(
        "--file-name",
        default="normalized_counts.csv.gz",
        type=str
    )
    parser.add_argument(
        "--rds",
        help="Path to Seurat obj rds, if this is set, make heatmaps",
        type=str
    )
    parser.add_argument(
        "--specific-ident",
        help="Only using markers of specific group",
        type=str
    )

    if len(sys.argv) <= 1:
        parser.print_help()
        exit(0)

    try:
        return parser.parse_args(sys.argv[1:])
    except ap.ArgumentError as err:
        logger.error(err)
        parser.print_help()
        exit(err.__hash__())


def main(args: ap.Namespace):
    u"""
    main function
    :param args:
    :return:
    """
    if not os.path.exists(args.o):
        os.makedirs(args.o)

    logger.info("Reading {0}".format(args.x))
    group_spec = pd.read_excel(args.x, index_col=0)

    # select qualified markers
    idx = []
    for i, j in zip(group_spec["p_val_adj"] < args.pvalue, group_spec["avg_logFC"] > args.logfc):
        idx.append(i and j)

    if args.specific_ident:
        temp_idx, idx = idx, []

        for i, j in zip(temp_idx, group_spec["ident"]):
            idx.append(i and (str(j) == args.specific_ident))

    group_spec = group_spec.loc[idx, :]

    tsne = pd.read_csv(os.path.join(args.i, "tsne.csv.gz"), index_col=0, engine="c")

    data_labels = ["Imbalanced", "ENN", "RENN"]

    for i, data in zip(data_labels, load_from_csv(args.i, counts_file=args.file_name)):
        logger.info(i)

        temp_data = os.path.join(args.o, "{0}.csv.gz".format(i))
        data.to_df().transpose().to_csv(temp_data)

        make_dotplot(
            coord=tsne,
            data=data,
            filename="{0}.png".format(os.path.join(args.o, i))
        )
        perform_kmeans_on_data(
            data=data,
            group_spec=group_spec,
            prefix=os.path.join(args.o, "{0}_KMeans".format(i)),
            n_jobs=args.p,
            min_cluster=1,
            max_cluster=args.max_cluster,
            random_state=args.random_state,
            rds=args.rds,
            group_by="Stage" if args.x.endswith("stage.xlsx") else "res.0.6",
        )

        perform_ica_on_data(
            data=data,
            group_spec=group_spec,
            prefix=os.path.join(args.o, "{0}_ICA".format(i)),
            min_cluster=1,
            max_cluster=args.max_cluster,
            random_state=args.random_state,
            n_jobs=args.p,
            n_pcs=args.n_comps,
            do_ard=args.ard,
            rds=args.rds,
            group_by="Stage" if args.x.endswith("stage.xlsx") else "res.0.6",
        )

        if args.rds:
            perform_mfuzz(
                data=temp_data,
                group_by="Stage" if args.x.endswith("stage.xlsx") else "res.0.6",
                rds=args.rds,
                genes=group_spec.index,
                random_state=args.random_state,
                output=os.path.join(args.o, "{0}_Mfuzz".format(i))
            )

            perform_wgcna(
                data=temp_data,
                group_by="Stage" if args.x.endswith("stage.xlsx") else "res.0.6",
                rds=args.rds,
                genes=group_spec.index,
                random_state=args.random_state,
                output=os.path.join(args.o, "{0}_WGCNA".format(i))
            )

    logger.info("Make upset plot")
    make_upset_plot(data_labels, args.o)


if __name__ == '__main__':
    main(command_line())

