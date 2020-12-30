#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@2019.02.07

Decompress cellranger compressed files
"""
import sys
import gzip
import os


def main(indir):
    u"""
    Main function
    :return:
    """
    files = {
        "barcodes.tsv.gz": ["barcodes.tsv"],
        "features.tsv.gz": ["features.tsv", "genes.tsv"],
        "matrix.mtx.gz": ["matrix.mtx"]
    }

    for i in os.listdir(indir):
        print(i)
        source_dir = os.path.join(indir, f"{i}/outs/filtered_feature_bc_matrix/")
        target_dir = os.path.join(indir, f"{i}/outs/filtered_feature_bc_matrix_for_seurat/")

        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        for ori, tar in files.items():
            open_files = [open(f"{target_dir}{x}", "w+") for x in tar]

            with gzip.open(f"{source_dir}{ori}", "rb") as r:
                content = r.read().decode("utf-8")

            [x.write(content) for x in open_files]
            [x.close() for x in open_files]


if __name__ == '__main__':

    if len(sys.argv) <= 1:
        print("Please provide the directory of cellranger results")
    else:
        main(sys.argv[1])
