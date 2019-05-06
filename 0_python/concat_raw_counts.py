#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu

This script is used to merge gene expression value into a raw count table (row -> genes, col -> cells)

"""
import gzip
import math
import os
import click
import pandas as pd
import matplotlib
matplotlib.use("svg")    # set the matplotlib backend, in case the matplotlib backend errors

import scanpy as sc
from tqdm import tqdm


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS
)
def main():
    u"""
    Main function
    :return:
    """
    pass


@click.command(
    context_settings=CONTEXT_SETTINGS,
    short_help="Merge cellranger output counts matrix"
)
@click.option(
    '-i',
    '--input-dir',
    help='Path to cellranger output directory',
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-o',
    '--output',
    help='Prefix of output',
    type=click.Path(),
    required=True
)
@click.option(
    "--rs",
    help="Threshold of rowSums, gene expression across all cells",
    default=0,
    type=click.IntRange(0)
)
@click.option(
    "--cs",
    help="Threshold of colSums, gene expression of single cells",
    default=0,
    type=click.IntRange(0)
)
@click.option(
    "--cells",
    help="Threshold of number of cells in single sample",
    default=200,
    type=click.IntRange(0)
)
@click.option(
    "--genes",
    help="Threshold of expressed genes in each cell",
    default=500,
    type=click.IntRange(0)

)
def merge(
        input_dir,
        output,
        rs,
        cs,
        cells,
        genes
):
    u"""
    Merge all the expression data from cellranger into raw count table

    \f
    :return:
    """

    if not os.path.exists(input_dir) or not os.path.isdir(input_dir):
        raise FileNotFoundError("{} not found".format(os.path.abspath(input_dir)))

    out_dir = os.path.dirname(input_dir)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    data = []
    for i in os.listdir(input_dir):
        mtx = os.path.join(input_dir, i)
        mtx = os.path.join(mtx, "outs/filtered_feature_bc_matrix")

        if not os.path.exists(mtx):
            print("{} not exists".format(mtx))
        else:
            print(mtx)

        try:
            temp = sc.read_10x_mtx(mtx, var_names='gene_symbols', cache=True)
            temp = temp.to_df().transpose()
            temp.columns = ["{0}_{1}".format(i, x) for x in temp.columns]

            drop = []
            for i in tqdm(temp.columns, desc="Finding low expression cells"):
                if temp[i].sum() < cs:
                    drop.append(i)
                    continue

                if sum(temp[i] > 0) < genes:
                    drop.append(i)

            temp = temp.drop(drop, axis=1)

            if temp.shape[1] < cells:
                continue

            data.append(temp)
        except Exception as err:
            # there 2 or 3 different exceptions that scanpy could throw, can't remember all
            # ZeroDivisionError or sth else
            print("read {0} throws {1}".format(mtx, err))
            continue

    print("concat")
    data = pd.concat(data, axis=1, sort=True)

    print("float to int")
    data = data.astype("int")

    print("drop genes by expression")
    drop = []
    for i, j in tqdm(data.iterrows(), desc="Finding low expression genes"):
        if j.sum() < rs:
            drop.append(i)
    data = data.drop(drop)

    print("Write to {0}".format(output + ".csv.gz"))

    data.to_csv(output + ".csv", na_rep=0, chunksize=10)


@click.command(
    context_settings=CONTEXT_SETTINGS,
    short_help="Split raw counts table by tissue"
)
@click.option(
    '-i',
    '--input-file',
    help='Path to raw counts table',
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '-m',
    "--meta",
    help="Path to meta info by cell csv",
    type=click.Path(exists=True),
    required=True
)
@click.option(
    '--column',
    help="Which columns to split",
    type=click.IntRange(1),
    default=1
)
@click.option(
    '-o',
    '--output',
    help='Path to output directory',
    type=click.Path(),
    required=True
)
def split(
        input_file,
        column,
        meta,
        output
):
    u"""
    Merge all the expression data from cellranger into raw count table

    \f
    :return:
    """
    data = {}

    meta = pd.read_csv(meta, index_col=0)
    for key, value in tqdm(meta.iloc[:, column - 1].to_dict().items()):

        if isinstance(key, float) and math.isnan(key):
            continue

        temp = data.get(value, [])
        temp.append(key)
        data[value] = temp

    infile = gzip.open(input_file, "rb")

    file_handlers = {}
    for k, v in data.items():
        if isinstance(k, float):
            k = "Unkown"
        file_handlers[k] = gzip.open(os.path.join(output, k + ".csv.gz"), "wt+")

        try:
            file_handlers[k].write(",{0}\n".format(",".join(v)))
        except TypeError:
            print(v)
            exit(0)

    file_index = {}
    for idx, line in tqdm(enumerate(infile)):
        lines = line.decode("utf-8").rstrip().split(",")
        if idx == 0:

            temp_indexes = {j: i for i, j in enumerate(lines)}

            for k, values in tqdm(data.items()):
                temp = []
                for v in values:
                    temp.append(temp_indexes[v])
                file_index[k] = temp

        else:
            for k, indexes in file_index.items():

                file_handlers[k].write("{0},{1}\n".format(lines[0], ",".join([lines[x] for x in indexes])))

    for v in file_handlers.values():
        v.close()


if __name__ == '__main__':
    main.add_command(merge)
    main.add_command(split)
    main()




