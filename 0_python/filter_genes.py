#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Filter out low expression genes
"""
import sys
import os
import gzip
import numpy as np
from tqdm import tqdm


def filter_by_expression(infile, outfile, threshold=1000):
    u"""
    filter by rowSums
    :param infile: path to input file
    :param outfile: path to output file
    :param threshold: int
    :return:
    """
    qualified = set()
    pbar = tqdm(total=os.path.getsize(infile))
    read = 0
    with open(outfile, "w+") as w:
        with open(infile) as r:
            line = r.readline()
            w.write(line)

            while True:
                line = r.readline()

                if not line:
                    break

                pbar.update(r.tell() - read)
                read = r.tell()

                lines = line.split()

                if sum([int(x) for x in lines[1:]]) >= 1000:
                    qualified.add(lines[0])

                w.write(line)

    with open(outfile + ".qualified", "w+") as w:
        for i in qualified:
            w.write(f"{i}\n")


def calculate_rs_and_cs(infile, outfile):
    rowSums = {}
    columns = None
    colSums = None

    with open(infile) as r:
        for idx, line in tqdm(enumerate(r)):
            lines = line.split()[1:]
            if idx == 0:
                columns = lines[1:]
                colSums = np.zeros(len(columns))
            else:

                value = np.array(lines[1:], dtype=np.int)
                rowSums[lines[0]] = value.sum()
                colSums += value

    with open(outfile + ".rowSums", "w+") as w:
        for k, v in tqdm(rowSums.items()):
            w.write(f"{k}\t{v}\n")

    with open(outfile + ".colSums", "w+") as w:
        for k, v in tqdm(zip(columns, colSums)):
            w.write(f"{k}\t{v}\n")


def merge_tow_counts(rowSums, colSums, int_gzip, counts, outfile):
    u"""

    :param rowSums:
    :param colSums:
    :param int_gzip:
    :param counts:
    :param outfile:
    :return:
    """

    rows = set()
    with open(rowSums) as r:
        for idx, line in tqdm(enumerate(r), desc="Reading RowSums"):
            if idx == 0:
                continue

            rows.add(line.split()[0])

    cols = set()
    with open(colSums) as r:
        for idx, line in tqdm(enumerate(r), desc="Reading ColSums"):
            if idx == 0:
                continue

            cols.add(line.split()[0])

    data = {}
    col_idx = {}
    with gzip.open(int_gzip) as r:
        for idx, line in tqdm(enumerate(r), desc="Reading int"):
            line = line.decode("utf-8")
            if idx == 0:
                for index, col in enumerate(line.split()):
                    if col in cols:
                        col_idx[index] = col
            else:
                lines = line.split()

                if lines[0] not in rows:
                    continue

                temp = {}
                for index, col in enumerate(lines[1:]):

                    try:
                        temp[col_idx[index]] = col
                    except KeyError:
                        continue
                data[lines[0]] = temp

    col_idx = {}
    with open(counts) as r:
        for idx, line in tqdm(enumerate(r), desc="Reading zemin"):
            lines = line.split()[1:]
            if idx == 0:
                for index, col in enumerate(lines[1:]):
                    if col in cols:
                        col_idx[index] = col
            else:
                if lines[0] not in rows:
                    continue

                try:
                    temp = data[lines[0]]
                except KeyError:
                    continue
                for index, col in enumerate(lines[1:]):
                    try:
                        temp[col_idx[index]] = col
                    except KeyError:
                        continue
                data[lines[0]] = temp

    cols = sorted(list(cols))
    with open(outfile, "w+") as w:
        cols_str = '\t'.join(cols)
        w.write(f"{cols_str}\n")
        for row in rows:
            try:
                new_line = [data[row][col] for col in cols]
                new_line = "\t".join(new_line)
                w.write(f"{row}\t{new_line}\n")
            except KeyError:
                continue


if __name__ == '__main__':
    # filter_by_expression(
    #     sys.argv[1],
    #     sys.argv[2],
    #     1000
    # )

    # calculate_rs_and_cs(
    #     sys.argv[1],
    #     sys.argv[2]
    # )

    merge_tow_counts(
        rowSums="rowSums",
        colSums="colSums",
        int_gzip="raw_counts_int.txt.gz",
        counts="count.txt",
        outfile=".txt"
    )
