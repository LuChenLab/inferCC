#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.14
Make radar plots across all cell types

First to make a collect matrix of expression value and meta info
"""
import os
import re
from glob import glob
from multiprocessing import Pool

import pandas as pd

from tqdm import tqdm


def read(path):
    return pd.read_csv(path, index_col=0, engine='python')


def get_all_paga(path):
    u"""
    get paga path of each cells
    """
    pagas = glob(os.path.join(path, "*/paga"))
    
    res = []
    for i in pagas:
        if re.search(r"_(SCC|ADC|P\d+|Tumor)", i, re.I):
            continue
        
        res.append(i)
        
    return sorted(res)


def merge_counts(path):
    u"""
    merge scale.data not raw.data
    :param path: list of paga path
    """
    data = []
    for i in tqdm(path):
        data.append(pd.read_csv(os.path.join(i, "scale.csv.gz"), index_col=0))
    
    # path = [os.path.join(x, "scale.csv.gz") for x in path]
    
    # with Pool(n_jobs) as p:
    #     data = p.map(read, path)
        
    return pd.concat(data, axis=1)


def merge_metas(path):
    u"""
    merge meta.csv.gz
    """
    data = []
    for i in tqdm(path):
        temp = pd.read_csv(os.path.join(i, "meta.csv.gz"), index_col=0)
        temp = temp[["PatientID", "cell_name", "Stage"]]
        data.append(temp)
        
    return pd.concat(data, axis=0)



def merge_raw_counts(path):
    u"""
    merge scale.data not raw.data
    :param path: list of paga path
    """
    data = []
    for i in tqdm(path):
        data.append(pd.read_csv(os.path.join(i, "raw.csv.gz"), index_col=0))
    # path = [os.path.join(x, ) for x in path] 
    # with Pool(n_jobs) as p:
    #     data = p.map(read, path)
        
    return pd.concat(data, axis=1)    
    


if __name__ == '__main__':
    files = get_all_paga("/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/")
    
    scale = merge_counts(files)
    meta = merge_metas(files)
    raw = merge_raw_counts(files)
    
    scale.to_csv("/mnt/raid62/Lung_cancer_10x/tmod_stage/scale.csv")
    raw.to_csv("/mnt/raid62/Lung_cancer_10x/tmod_stage/raw.csv") 
    meta.to_csv("/mnt/raid62/Lung_cancer_10x/tmod_stage/meta.csv") 
