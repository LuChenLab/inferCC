#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@Created at 2019.05.06+

"""

import logging
import os
import sys
from argparse import ArgumentParser, ArgumentError, Namespace
from glob import glob
from multiprocessing import Pool, Manager
from typing import List

import anndata
import pandas as pd
import peewee as pw
from tqdm import trange, tqdm

__dir__ = os.path.abspath(os.path.dirname(__file__))

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

STR_LIST = List[str]


db = pw.SqliteDatabase(os.path.join("/mnt/data1/zhangyiming/CellRanger", "velocyto.db"))


class BaseModel(pw.Model):
    u"""
    Base model of peewee database
    """

    class Meta:
        database = db


class Feature(BaseModel):

    feature_id = pw.CharField(index=True, unique=True)
    feature_name = pw.CharField(index=True)
    chromosome = pw.CharField()
    start = pw.IntegerField()
    end = pw.IntegerField()
    strand = pw.CharField()


class Barcode(BaseModel):
    barcode = pw.CharField()
    clusters = pw.IntegerField()
    x = pw.IntegerField()
    y = pw.IntegerField()


class Layers(BaseModel):

    feature = pw.ForeignKeyField(Feature)
    barcode = pw.ForeignKeyField(Barcode)

    value = pw.FloatField()
    type = pw.CharField()


class Matrix(BaseModel):

    feature = pw.ForeignKeyField(Feature)
    barcode = pw.ForeignKeyField(Barcode)

    value = pw.FloatField()


class Ambiguous(BaseModel):

    feature = pw.ForeignKeyField(Feature)
    barcode = pw.ForeignKeyField(Barcode)

    value = pw.FloatField()



class Spliced(BaseModel):

    feature = pw.ForeignKeyField(Feature)
    barcode = pw.ForeignKeyField(Barcode)

    value = pw.FloatField()


class Unspliced(BaseModel):

    feature = pw.ForeignKeyField(Feature)
    barcode = pw.ForeignKeyField(Barcode)

    value = pw.FloatField()


def loom_to_csv(args) -> None:
    u"""
    Convert velocyte loom file to csv
    :param loom_path:
    :param output:
    :return:
    """
    loom_path, output = args
    if not os.path.exists(output):
        os.makedirs(output)

    logger.info("Loading from {0}".format(loom_path))
    data = anndata.read_loom(os.path.abspath(loom_path))

    logger.info("Feature of {0}".format(os.path.basename(loom_path)))
    data.var.to_csv(os.path.join(output, "var.csv.gz"))

    logger.info("Barcode of {0}".format(os.path.basename(loom_path)))
    data.obs.to_csv(os.path.join(output, "obs.csv.gz"))

    for i in ["matrix", "ambiguous", "spliced", "unspliced"]:
        logger.info("{1} of {0}".format(os.path.basename(loom_path), i))

        temp = pd.DataFrame(data.layers[i].todense())
        temp.columns = data.var.index
        temp["index"] = data.obs.index
        temp = temp.melt(id_vars="index")
        temp = temp.loc[temp["value"] > 0, :]
        temp.to_csv(os.path.join(output, "{0}.csv.gz".format(i)))


def insert_loom(path: str):
    u"""
    Load loom
    :param path: path to loom file
    :return:
    """
    path = os.path.abspath(path)

    logger.info("Loading from {0}".format(path))
    data = anndata.read_loom(os.path.abspath(path))

    temp_feature_idx = {}
    temp_barcode_idx = {}

    # add features
    features = []
    query = Feature.select()
    db_idx = query.count() + 1
    curr_idx = 0

    exists_features = {}
    for i in query:
        exists_features[i.feature_id] = i.id

    for idx, row in data.var.iterrows():

        if row["Accession"] in exists_features.keys():
            temp_feature_idx[curr_idx] = exists_features[row["Accession"]]
        else:
            features.append({
                "id": db_idx,
                "feature_id": row["Accession"],
                "feature_name": idx,
                "chromosome": row["Chromosome"],
                "start": row["Start"],
                "end": row["End"],
                "strand": row["Strand"]
            })

            temp_feature_idx[curr_idx] = db_idx
            db_idx += 1

        curr_idx += 1

    with db.atomic():
        for i in trange(0, len(features), 100):
            Feature.insert_many(features[i:i + 100]).execute()

    del features

    # add barcode
    barcodes = []
    db_idx = Barcode.select().count() + 1
    curr_idx = 0
    for idx, row in data.obs.iterrows():
        barcodes.append({
            "barcode": idx,
            "clusters": row["Clusters"],
            "x": row["_X"],
            "y": row["_Y"]
        })

        temp_barcode_idx[curr_idx] = db_idx
        curr_idx += 1
        db_idx += 1

    with db.atomic():
        for i in trange(0, len(barcodes), 100):
            Barcode.insert_many(barcodes[i:i + 100]).execute()

    del barcodes

    # add each layer
    for i in ["matrix", "ambiguous", "spliced", "unspliced"]:

        temp_data = []
        temp_layer = data.layers[i]

        for row in range(temp_layer.shape[0]):
            cols = temp_layer[row].nonzero()[1]
            temp_layer_data_single_barcode = {}
            barcode = temp_barcode_idx[row]
            for col in cols:
                feature = temp_feature_idx[col]
                value = temp_layer[row, col]

                temp_data.append({
                    "feature": feature,
                    "barcode": barcode,
                    "value": value,
                    "type": i
                })

        with db.atomic():
            for i in trange(0, len(temp_data), 100):
                Layers.insert_many(temp_data[i: i + 100]).execute()


def insert_barcodes(barcods: STR_LIST):
    u"""
    insert barcodes into database
    :param barcods: list of path of barcodes
    :return:
    """
    if Barcode.table_exists():
        Barcode.drop_table()
    Barcode.create_table()

    data = []
    for i in tqdm(barcods):
        temp = pd.read_csv(i, index_col=0)

        for idx, row in temp.iterrows():
            data.append({
                "barcode": idx,
                "clusters": row["Clusters"],
                "x": row["_X"],
                "y": row["_Y"]
            })

    with db.atomic():
        for i in trange(0, len(data), 100):
            Barcode.insert_many(data[i: i + 100]).execute()


def insert_features(features: STR_LIST):
    u"""
    insert features into database
    :param features: list of path of features
    :return:
    """
    if Feature.table_exists():
        Feature.drop_table()
    Feature.create_table()

    data = {}
    for i in tqdm(features):
        temp = pd.read_csv(i, index_col=0)

        for idx, row in temp.iterrows():
            data[row["Accession"]] = {
                "feature_id": row["Accession"],
                "feature_name": idx,
                "chromosome": row["Chromosome"],
                "start": row["Start"],
                "end": row["End"],
                "strand": row["Strand"]
            }

    data = list(data.values())
    with db.atomic():
        for i in trange(0, len(data), 100):
            Feature.insert_many(data[i: i + 100]).execute()


def insert(args):
    f, barcodes, features, table, lock = args
    logger.info(f)
    temp = pd.read_csv(f)

    data = []
    for _, row in temp.iterrows():
        data.append({
            "barcode": barcodes[row["index"]],
            "feature": features[row["variable"]],
            "value": row["value"]
        })

    lock.acquire()
    with db.atomic():
        for j in trange(0, len(data), 100):
            table.insert_many(data[j: j + 100]).execute()
    lock.release()


def insert_value(table: BaseModel, files: STR_LIST, processes: int):
    u"""
    insert layer value
    :param table: target table
    :param files: list of file path
    :param processes: use how many processes
    :return:
    """
    if table.table_exists():
        table.drop_table()
    table.create_table()

    querys = Barcode.select()
    barcodes = {}
    for i in tqdm(querys):
        barcodes[i.barcode] = i.id

    querys = Feature.select()
    features = {}
    for i in tqdm(querys):
        features[i.feature_name] = i.id

    m = Manager()
    lock = m.Lock()

    args = [[x, barcodes, features, table, lock] for x in files]

    with Pool(processes) as p:
        p.map(insert, args)


def main(args: Namespace):
    u"""

    :param args:
    :return:
    """

    logger.info("Convert loom files to csv")

    looms = glob(os.path.join(args.input, "*/velocyto/*.loom"))

    assert len(looms) > 0, "cannot find any loom file"

    csv_dir = os.path.join(args.output, "csv")

    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)

    cmds = []
    for i in looms:
        cmds.append([
            i,
            os.path.join(csv_dir, os.path.basename(i).split(".")[0])
        ])

    with Pool(args.processes) as p:
        p.map(loom_to_csv, cmds)

    logger.info("Insert csv to database")

    barcodes = glob(os.path.join(csv_dir, "*/obs.csv.gz"))
    insert_barcodes(barcodes)

    features = glob(os.path.join(csv_dir, "*/var.csv.gz"))
    insert_features(features)

    matrix = glob(os.path.join(csv_dir, "*/matrix.csv.gz"))
    ambiguous = glob(os.path.join(csv_dir, "*/ambiguous.csv.gz"))
    spliced = glob(os.path.join(csv_dir, "*/spliced.csv.gz"))
    unspliced = glob(os.path.join(csv_dir, "*/unspliced.csv.gz"))

    insert_value(Ambiguous, ambiguous, args.processes)
    insert_value(Spliced, spliced, args.processes)
    insert_value(Unspliced, unspliced, args.processes)
    insert_value(Matrix, matrix, args.processes)


if __name__ == '__main__':
    parser = ArgumentParser(description="Combine loom from different sample by cell type")

    parser.add_argument(
        "-i",
        "--input",
        help="Path to CellRanger output directory",
        type=str,
        required=True
    )

    parser.add_argument(
        "-m",
        "--meta",
        help="Path to meta data, csv file",
        type=str,
        required=True
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Path to output directory",
        required=True,
        type=str
    )

    parser.add_argument(
        "-p",
        "--processes",
        help="How many processes to use",
        default=1,
        type=int
    )

    if len(sys.argv) <= 1:
        parser.print_help()
        exit(0)

    try:
        main(parser.parse_args(sys.argv[1:]))
    except ArgumentError as err:
        logger.error(err)
        parser.print_help()


