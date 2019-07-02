#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.07.01
I really do not like R

Still not working for R

But I gonna use this
"""

import os
from glob import glob
import pandas as pd
from multiprocessing import Pool




def call(args):
    u"""
    :param args
        path
        output
    """
    path, output_f = args
    # data = pd.read_csv(path, sep = "\t", index_col=0)
    
    # print(data.head())
    
    # data.index = data.iloc[:, 0]

    # data = data.iloc[:, range(1, data.shape[1])]
    
    # print(data.head()) 
    
    # data.columns = data.loc["mc_id", :]
    # data = data.drop(["mc_id", "n_cells", "group"])
    # data = data.drop("mc_id", axis = 1)

    # data.to_csv(output)
    
    row_idx = -1
    with open(output_f, "w+") as w:
        with open(path) as r:
            for line in r:
                line = line.strip()
                
                lines = line.split()
                
                for i in lines:
                    if i == "mc_id":
                        row_idx = lines.index("mc_id")
                        
                if row_idx >= 0:
                    w.write(",".join(lines[row_idx:]) + "\n")                    


def main(input_dir, output_dir, n_jobs = 10):
    u"""
    main
    """
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    tasks = []
    for i in glob(os.path.join(input_dir, "*/figs/*.log2_mc_fp.txt")):
        temp_out = os.path.join(output_dir, os.path.basename(i))
        
        if i != temp_out:
            tasks.append([i, temp_out])
    
    print(len(tasks))
    with Pool(n_jobs) as p:
        p.map(call, tasks)
    

if __name__ == '__main__':
    from fire import Fire
    
    Fire(main)
