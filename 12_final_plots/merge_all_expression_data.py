#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Merge different expression data for next plot
"""
import os
from glob import glob

from multiprocessing import Pool

import pandas as pd


def read_data(path):
    u"""
    read data from paga
    """
    data = pd.read_csv(os.path.join(path, "scale.csv.gz"), index_col=0, engine="python")
    meta = pd.read_csv(os.path.join(path, "meta.csv.gz"), index_col=0, engine="python")
    try:
        meta = meta.loc[:, ["res.0.6", "Stage", "Disease"]]
    except KeyError as err:
        print(err)
        print(path)
        print(meta.head())
        
        raise KeyError(err)
    
    cell = os.path.basename(os.path.dirname(path))
    
    meta["Cell"] = [cell.split("_")[0] for x in range(meta.shape[0])]

    return data, meta


def main(input_dir, output_dir, n_jobs = 20):
    u"""
    main
    """
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    files = glob(os.path.join(input_dir, "*_ADC/paga"))
    files += glob(os.path.join(input_dir, "*_SCC/paga")) 
    
    with Pool(n_jobs) as p:
        res = p.map(read_data, files)
        
    data, meta = [], []
    
    for i in res:
        data.append(i[0])
        meta.append(i[1])
        
    pd.concat(data, axis=1).to_csv(os.path.join(output_dir, "scale.csv"))
    pd.concat(meta).to_csv(os.path.join(output_dir, "meta.csv")) 
    
    
if __name__ == '__main__':
    from fire import Fire
    
    Fire(main)
    

    