#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.25
This scripts is used to merge different cellranger output matrix toghether
"""

import os
import pandas as pd
from glob import glob
from tqdm import tqdm
from multiprocessing import Pool


class Matrix(object):
    u"""
    object to load features
    """
    
    def __init__(self, path: str):
        u"""
        init this class
        """
        print(path)
        self.path = os.path.abspath(path)
        self.cell = os.path.abspath(path).split("/")[-3]

        self.matrix = []
        
        self.__load__()
        
    def __load__(self):
        u"""
        load data
        """
        barcodes = []
        features = []
        with open(os.path.join(self.path, "barcodes.tsv")) as r:
            for line in r:
                line = line.strip()
                
                barcodes.append("{0}-{1}".format(self.cell, line))

        with open(os.path.join(self.path, "features.tsv")) as r:
            for line in r:
                line = line.strip()
                features.append(line)

        header = True
        with open(os.path.join(self.path, "matrix.mtx")) as r:
            for line in r:
                if line.startswith("%"):
                    continue
                
                if header:
                    header = False
                    continue
                
                lines = line.split()
                
                self.matrix.append({
                    "barcode": barcodes[int(lines[1]) - 1],
                    "feature": features[int(lines[0]) - 1],
                    "value": int(lines[2])
                })
                
    def __add__(self, other):
        u"""
        merge two object together
        """
        
        self.matrix += other.matrix
        
        return self
        
    def dump(self, output):
        u"""
        dump matrix to marketmatrix
        :param output path to output file
        """
        
        if not os.path.exists(output):
            os.makedirs(output)
        
        barcodes = set()
        features = set()
        total = 0
        for i in self.matrix:
            barcodes.add(i["barcode"])
            features.add(i["feature"])
            total += i["value"]
        
        barcodes = {y: x + 1 for x, y in enumerate(sorted(barcodes))}
        
        
        features = {y: x + 1 for x, y in enumerate(sorted(features))}
        
        with open(os.path.join(output, "barcodes.tsv"), "w+") as w:
            w.write("\n".join(barcodes.keys()) + "\n")
            
        with open(os.path.join(output, "features.tsv"), "w+") as w:
            w.write("\n".join(features.keys()) + "\n")
            
        with open(os.path.join(output, "matrix.mtx"), "w+") as w:
            w.write("%%MatrixMarket matrix coordinate integer general\n%metadata_json: {\"format_version\": 2, \"software_version\":\"3.0.1\"}\n")
            w.write("{0} {1} {2}\n".format(len(features), len(barcodes), len(self.matrix)))
            for i in tqdm(self.matrix):
                w.write("{0} {1} {2}\n".format(
                    features[i["feature"]], 
                    barcodes[i["barcode"]], 
                    i["value"]
                ))


class Raw(object):
    u"""
    dump pandas frame to marketmatrix file
    """
    
    def __init__(self, path):
        u"""
        init
        """
        data = pd.read_csv(path, index_col=0)
        # print("convert to dict")
        self.data = data.to_dict()
        
                
    def dump_to_mm(self, output: str, full_features=None):
        u"""
        dump dataframe to marketmatrix file
        :param output directory
        :param full_features a dict, key is gene symbol. convert gene symbol to ENSG00000001631 KRIT1   Gene Expression
        """

        if not os.path.exists(output): 
            os.makedirs(output)
        
        features = set()
        barcodes = set()
        matrix = []
        
        for barcode in tqdm(self.data.keys()):
            barcodes.add(barcode)
            temp = self.data[barcode]
            
            for feature in temp.keys():
                features.add(feature)

                if temp[feature] > 0:
                    matrix.append([feature, barcode, temp[feature]])
                    
        features = {y: x + 1 for x, y in enumerate(sorted(features))}
        barcodes = {y: x + 1 for x, y in enumerate(sorted(barcodes))}
        
        with open(os.path.join(output, "genes.tsv"), "w+") as w:
            if full_features is not None:
                for i in features:
                    w.write(full_features[i] + "\n")
            else:
                w.write("\n".join(features) + "\n")
            
        with open(os.path.join(output, "barcodes.tsv"), "w+") as w:
            w.write("\n".join(barcodes) + "\n")
            
        with open(os.path.join(output, "matrix.mtx"), "w+") as w:
            w.write("%%MatrixMarket matrix coordinate integer general\n%metadata_json: {\"format_version\": 2, \"software_version\":\"3.0.1\"}\n")
            w.write("{0} {1} {2}\n".format(len(features), len(barcodes), len(matrix)))
            
            for i in matrix:
                w.write("{0} {1} {2}\n".format(features[i[0]], barcodes[i[1]], int(i[2])))


def wrapper(args):
    path, full_features = args
    
    print(path)
    
    data = Raw(path)
    
    output = os.path.dirname(path)

    data.dump_to_mm(output, full_features=full_features)


def get_id(id_, known_ids, num=0):
    u"""
    DUE TO R ALWAYS APPEND EXTRA NUMBER TO SAME ROWNAMES
    THEREFORE, USING THIS TO GET SAME GENE SYMBOL PATTERN
    """
    
    if num > 0:
        new_id_ = "{0}-{1}".format(id_, num)
    else:
        new_id_ = id_
    
    if new_id_ in known_ids:
        return get_id(id_, known_ids, num+1)
    else:
        return new_id_



def main(input_dir, output_dir=None, n_jobs=20, raw=None):
    u"""
    main function
    :param input_dir: 
    :param output_dir: output file path
    :param n_jobs:
    """
    if raw:
        full_features = {}
        known_ids = set()
        
        with open(raw) as r:
            for line in r:
                lines = line.rstrip().split("\t")
                id_ = get_id(lines[1], known_ids)
                known_ids.add(id_)
                
                lines[1] = id_
                full_features[id_] = "\t".join(lines)
                
        files = sorted(glob(os.path.join(input_dir, "*/paga/raw.csv.gz")))
        
        tasks = [[x, full_features] for x in files][:1]
        
        with Pool(n_jobs) as p:
            p.map(wrapper, tasks)
        
    else:
        files = glob(os.path.join(input_dir, "*/outs/filtered_feature_bc_matrix_for_seurat"))
        if files:
            data = Matrix(files[0])
            
            if len(files) > 1:
                for i in tqdm(files[1:]):
                    data += Matrix(i)
            
            # print(len(files))
            
            # with Pool(n_jobs) as p:
            #     data = p.map(wrapper, list(chunks(files, n_jobs)))
                
            # curr = data[0]
            # for i in data[1:]:
            #     curr += i
            
            output_dir = os.path.abspath(output_dir)
            data.dump(output_dir)


if __name__ == '__main__':
    from fire import Fire
    
    Fire(main)


