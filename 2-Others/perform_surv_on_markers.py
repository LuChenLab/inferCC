#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2019.10.12

Perform surv on cluster markers
"""
import os

from glob import glob
from subprocess import check_call

import click
import pandas as pd

__dir__ = os.path.dirname(os.path.abspath(__file__))


@click.command(
    context_settings=dict(help_option_names=['-h', '--help']),
)
@click.version_option("0.0.1", message="Current version %(version)s")
@click.option(
    "-i",
    "--input-path",
    default=None,
    type=click.STRING,
    help="path to input directory",
    show_default=True
)
def main(input_path):
    u"""
    Main
    """
    
    exception = ["CD8LUAD", "CiliaLUAD", "CiliaLUSC", "ClubLUAD", "ECLUSC"]
    
    files = glob(os.path.join(input_path, "*/*/cluster_markers.csv"))
    
    for i in files:
        in_dir = os.path.dirname(i)
        cell = os.path.basename(os.path.dirname(in_dir))
        disease = os.path.basename(in_dir)
        
        print(cell, disease)
                
        try:
            data = pd.read_csv(i, index_col=0)
            data = data.loc[data["p_val_adj"] < 0.05, :]
            data = data.sort_values("avg_logFC", ascending=False)
            
            for ident in data["ident"].unique():
                genes = ",".join(data.loc[data["ident"] == ident, "gene"][:10])

                for c in [50, 75]:
                    out_file = os.path.join(in_dir, f"Cluster_surv_{ident}_{c}.pdf")
                    
                    if "{}{}".format(cell, disease) not in exception and os.path.exists(out_file):
                        continue
                    check_call(f"python {os.path.join(__dir__, 'gepia_survival.py')} -i {genes} -o {out_file} -d {disease} -c {c} --gene-set", shell=True)
        except KeyError as err:
            print(err)
        
if __name__ == '__main__':
    main()