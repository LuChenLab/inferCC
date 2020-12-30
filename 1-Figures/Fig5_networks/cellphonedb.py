#!/usr/bin/env python3
import os
from subprocess import check_call
from fire import Fire
from shutil import rmtree



def main(input_dir, threads = 10, force = False):
    for parent, _, files in os.walk(input_dir):
        for f in files:
            meta = os.path.join(parent, "meta.txt")
            counts = os.path.join(parent, "counts.txt")
            if f == "meta.txt" and os.path.exists(counts):
                out = os.path.join(parent, "out")
                
                if force and os.path.exists(out):
                    rmtree(out)
                
                if not os.path.exists(out):
                    
                    os.makedirs(out)
                    
                    check_call(
                        f"cellphonedb method statistical_analysis {meta} {counts} --counts-data=gene_name --threads={threads} --output-path={out}",
                        shell=True
                    )
                    check_call(f"cellphonedb plot heatmap_plot --pvalues-path {out}/pvalues.txt --output-path {out} {meta}", shell=True)

if __name__ == '__main__':
    Fire(main)