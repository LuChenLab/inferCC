#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.26

get and select part of good quality immu cells
"""

import logging
import os

import pandas as pd
from anndata import AnnData
from imblearn.under_sampling import RepeatedEditedNearestNeighbours
from tqdm import tqdm

CELLS = {
    # "Pulmonary alveolar type II cells": "Alveolar_II",
    "B cells": "B_cells",
    # "Basal cells": "Basal",
    "CD4+ T cells": "CD4",
    # "Ciliated cells": "Ciliated",
    # "Club cells": "Club",
    "Dendritic": "Dendritic",
    # "Endothelial cells": "Endothelial",
    "Erythroid-like and erythroid precursor cells": "Erythroid_precursor",
    "Exhaust T": "Exhaust_T",
    # "Fibroblasts": "",
    # "Airway goblet cells": "Goblet",
    "Granulocyte": "Granulocyte",
    "Mast": "Mast",
    "Monocytes": "Monocytes",
    # "Neuroendocrine cell": "",
    "NK": "NK",
    "T cells": "T_cells",
    "Treg": "Treg"
}

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] - [%(name)s] - [%(levelname)s] - %(message)s')
logger = logging.getLogger(__name__)


def load_from_csv(
        input_dir: str,
        counts_file: str="raw.csv.gz",
        n_jobs=1,
        low_expression=0.1
) -> (AnnData, AnnData):
    u"""
    load data from csv files
    :param input_dir:
    :param counts_file:
    :param n_jobs
    :param str
    :return:
    """
    
    logger.info("Reading {0}".format(input_dir))

    input_file = os.path.join(input_dir, counts_file)

    if not os.path.exists(input_file):
        return None

    mtx = pd.read_csv(input_file, index_col=0)
    meta = pd.read_csv(os.path.join(input_dir, "meta.csv.gz"), index_col=0)
    meta = meta.loc[meta.index, :]

    # logger.info(mtx.shape)
    # # filter low expressed genes
    # genes_sum = [x / mtx.shape[1] > low_expression for x in mtx.sum(axis=1)]

    # mtx = mtx.loc[genes_sum, :]

    logger.info(mtx.shape)
    mtx = mtx.transpose()
    
    data = AnnData(mtx, obs=meta)
    data.obs = meta

    # logger.info("Perform RENN")
    # renn = RepeatedEditedNearestNeighbours(n_jobs=n_jobs, return_indices=True)

    # mtx_renn, group_renn, idx_renn = renn.fit_resample(mtx, meta["Stage"])

    # data_renn = AnnData(mtx.iloc[list(idx_renn), :], meta.iloc[idx_renn, :])

    # data_renn.obs = meta.iloc[idx_renn, :]
    
    return data


def dump_to_mm(data, output: str, full_features=None):
    u"""
    dump dataframe to marketmatrix file
    :param output directory
    :param full_features a dict, key is gene symbol. convert gene symbol to ENSG00000001631 KRIT1   Gene Expression
    """
    data = data.to_dict()

    if not os.path.exists(output): 
        os.makedirs(output)
    
    features = set()
    barcodes = set()
    matrix = []
    
    for feature in tqdm(data.keys()):
        features.add(feature)
        temp = data[feature]
        
        for barcode in temp.keys():
            barcodes.add(barcode)

            if temp[barcode] > 0:
                matrix.append([feature, barcode, temp[barcode]])
                
    features = {y: x + 1 for x, y in enumerate(sorted(features))}
    barcodes = {y: x + 1 for x, y in enumerate(sorted(barcodes))}
    
    with open(os.path.join(output, "genes.tsv"), "w+") as w:
        if full_features is not None:
            for i in features:
                w.write(full_features[i] + "\n")
        else:
            w.write("\n".join(features) + "\n")
        
    with open(os.path.join(output, "barcodes.tsv"), "w+") as w:
        w.write("\n".join(barcodes) + "\n")
        
    with open(os.path.join(output, "matrix.mtx"), "w+") as w:
        w.write("%%MatrixMarket matrix coordinate integer general\n%metadata_json: {\"format_version\": 2, \"software_version\":\"3.0.1\"}\n")
        w.write("{0} {1} {2}\n".format(len(features), len(barcodes), len(matrix)))
        
        for i in matrix:
            w.write("{0} {1} {2}\n".format(features[i[0]], barcodes[i[1]], int(i[2])))


def get_id(id_, known_ids, num=0):
    u"""
    DUE TO R ALWAYS APPEND EXTRA NUMBER TO SAME ROWNAMES
    THEREFORE, USING THIS TO GET SAME GENE SYMBOL PATTERN
    """
    
    if num > 0:
        new_id_ = "{0}-{1}".format(id_, num)
    else:
        new_id_ = id_
    
    if new_id_ in known_ids:
        return get_id(id_, known_ids, num+1)
    else:
        return new_id_


def main(input_dir, output_dir, n_jobs = 20, raw=None):
    
    
    if raw is None:
        raw = None
    else:
        full_features = {}
        known_ids = set()
        
        with open(raw) as r:
            for line in r:
                lines = line.rstrip().split("\t")
                id_ = get_id(lines[1], known_ids)
                known_ids.add(id_)
                
                lines[1] = id_
                full_features[id_] = "\t".join(lines) 
    
    data = {}
    for cell in CELLS.values():
        for d in ["ADC", "SCC"]:
            key = "{0}_{1}".format(cell, d)
            path = os.path.join(input_dir, key, "paga")

            temp = load_from_csv(path, n_jobs=n_jobs)
            if temp:
                data[key] = temp
                
    res = []
    for i in data.values():
        res.append(i.to_df())
        
    res = pd.concat(res, axis=0)
    
    dump_to_mm(data = res, output="/mnt/raid62/Lung_cancer_10x/MetaCell/immu", full_features=full_features)
    
    output_dir = os.path.join("/mnt/raid62/Lung_cancer_10x/MetaCell/immu/anndata")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    for k, d in tqdm(data.items()):
        d.write(os.path.join(output_dir, k + ".h5ad"))
    
    