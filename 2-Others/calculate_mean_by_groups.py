#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Calculate mean by groups
"""
import os
import gzip

import pandas as pd

from tqdm import tqdm


def main(meta, scaled_data, output):
    u"""
    Main function
    """
    
    meta = pd.read_csv(meta, index_col=0)
    meta = meta.loc[(meta["Batch"] != 3) & (meta["cell_name"] != ""), :]

    columns = ["cell_name"]  # "Disease", "Stage", "PatientID", 
    
    print("Format groups")
    groups = {}
    for cell, row in tqdm(meta.iterrows()):
        key = "{}|{}|{}".format(
            # row["PatientID"],
            row["Disease"],
            row["Stage"],
            row["cell_short"]
        )
           
        temp = groups.get(key, set())
        temp.add(cell)
        groups[key] = temp
        
    print("calculate mean")
    cell_idx = {}
    with open(output, "w+") as w:
        keys = sorted(groups.keys())
        
        w.write("gene\t{}\n".format("\t".join(keys)))
        
        with gzip.open(scaled_data, "rt") as r:
            for idx, line in enumerate(tqdm(r)):
                lines = line.replace('"', "").split(",")
                if idx == 0:
                    for i, j in enumerate(lines):
                        cell_idx[j] = i
                else:
                    res_lines = [lines[0]]
                    for key in keys:
                        temp = []
                        required_cells = groups[key]
                        
                        for cell in required_cells:
                            temp.append(float(lines[cell_idx[cell]]))
                        res_lines.append(str(sum(temp) / len(temp)))

                    w.write("\t".join(res_lines) + "\n")


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
                        
                
                