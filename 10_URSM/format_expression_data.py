#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.26

- single_cell_matrix: row is cell, col is gene
- single_cell_type: one col, 0..(K-1) -> indicate cell types
- bulk_rnaseq: row is sample, col is gene
"""
import os
from glob import glob
import random

import pandas as pd

from multiprocessing import Pool


def read_data(args):
    u"""
    read data into pd.DataFrame
    """
    path, seed, num = args
    print(path)
    data = pd.read_csv(os.path.join(path, "raw.csv.gz"), index_col=0)
    meta = pd.read_csv(os.path.join(path, "meta.csv.gz"), index_col=0)
    
    stage_groups = {}
    for idx, row in meta.iterrows():
        temp = stage_groups.get(row["Stage"], [])
        temp.append(idx)
        stage_groups[row["Stage"]] = temp
    
    cell = path.split("/")[-3]
    cells = []
    
    columns = []
    for stage, value in stage_groups.items():
        if num is not None and len(value) > num:
            value = random.Random(seed).choices(value, k=num)
            
        for i in value:
            cells.append("{0}_{1}".format(cell, stage))

        columns += value            
            
    return data[columns], cells


def main(input_dir, output_dir, bulk, n_jobs = 20, random_state=0, num_cells=50):
    u"""
    :param input_dir
    :param output_dir
    """
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    print(bulk)
    bulk = pd.read_csv(bulk, index_col=0)
    # bulk.index = bulk["gene_name"]
    print(bulk.shape)
    
            
    files = glob(os.path.join(input_dir, "*_SCC/paga/"))
    tasks = [[x, random_state, num_cells] for x in files]
    
    with Pool(n_jobs) as p:
        data = p.map(read_data, tasks)
    
    expr = []
    cells = []
    for i in data:
        expr.append(i[0])
        cells += i[1]
    
    print("concat")
    expr = pd.concat(expr, axis=1)
    print(expr.shape)
    
    print("common genes")
    genes = set(expr.index) & set(bulk.index)
    print(len(genes))
    
    expr = expr.loc[genes, :]
    bulk = bulk.loc[genes, :]
    
    print("dump bulk")
    bulk = bulk.transpose()
    bulk.to_csv(os.path.join(output_dir, "bulk_expr.csv"), header = False, index=False)
    
    with open(os.path.join(output_dir, "bulk_expr_sample.txt"), "w+") as w:
        for i in bulk.index:
            w.write("{}\n".format(i))
    
    print("dump single cell")
    expr.transpose().to_csv(os.path.join(output_dir, "single_cell_expr.csv"), header = False, index=False)

    with open(os.path.join(output_dir, "single_cell_expr_cells.txt"), "w+") as w:
        w.write("\n".join(expr.columns))
        
    with open(os.path.join(output_dir, "single_cell_expr_genes.txt"), "w+") as w:
        w.write("\n".join(expr.index))
        
    uniq_cells = sorted(set(cells))
    
    with open(os.path.join(output_dir, "single_cell_expr_cells_type.csv"), "w+") as w:
        for i in cells:
            w.write("{}\n".format(uniq_cells.index(i)))
    

if __name__ == '__main__':
    from fire import Fire
    Fire(main)
    
    

# from wget import download


# import os
# from subprocess import check_call, CalledProcessError


# def download(url, out):
#     try:
#         check_call("wget -c -O {1} {0} ".format(url, os.path.abspath(out)), shell = True)
#     except CalledProcessError as err:
#         download(url, out)


# for f in ["cgraph.guo2018_tpm_scaled_filt_Tumor_Normal_Blood.Rda",
#             "coclust.guo2018_tpm_scaled_filt_Tumor_Normal_Blood.Rda",
#             "gset.guo2018_lateral.Rda",
#             "gstat.guo2018_tpm_scaled_filt_Tumor_Normal_Blood.Rda",
#             "lfp_screenshot.png",
#             "mat.guo2018_tpm_scaled_filt_guo2018_lateral.Rda",
#             "mat.guo2018_tpm_scaled_filt_Tumor_Normal_Blood.Rda",
#             "mc.guo2018_tpm_scaled_filt_Tumor_Normal_Blood_outClean_nonAnn.Rda"]:
#     download(
#         "http://www.wisdom.weizmann.ac.il/~lubling/metac_data/metacell_annotation/" + f, out="lung_db/" + f)
    
