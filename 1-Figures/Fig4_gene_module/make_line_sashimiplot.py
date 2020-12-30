#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2019.11.28

chr2    HAVANA  transcript      218781749       218815288

python pysashimi/main.py line \
    -b LungCancer10x/tests/all_bam.txt \
    -e chr2:218780000-218786749:+ \
    -g /mnt/raid63/LungCancerXieDanData/genome/gencode.v32.annotation.gtf  \
    --color-factor 3 \
    -o LungCancer10x/tests/ENST00000494263_ATAC_not_y.pdf \
    -p 20 \
    --share-y \
    --plot-by 4 \
    --sep-by-color
"""
import os
from multiprocessing import Pool
from subprocess import check_call, CalledProcessError


SCRIPT = "pysashimi/main.py"


class Gene(object):

    def __init__(self, chrom, start, end, strand, gene):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.gene = gene
        self.strand = strand

    def __hash__(self):
        return hash(self.gene)

    def promoter(self, up = 2000, down = 500):
        u"""
        return a promoter region
        """

        return "{}:{}-{}:{}".format(
            self.chrom,
            self.start - up if self.strand == "+" else self.end + up,
            self.start + down if self.strand == "+" else self.end - down,
            self.strand
        )


def load_reference(bed):
    u"""
    Load reference
    :param bed: path to a bed file
    """
    print("Loading %s" % bed)
    data = {}
    with open(bed) as r:
        for line in r:
            if line.startswith("M"):
                continue

            lines = line.split()

            data[lines[3]] = Gene(lines[0], lines[1], lines[2], lines[5], lines[3])

    return data


def call(cmd):
    try:
        with open(os.devnull, "w+") as w:
            check_call(cmd, shell = True, stdout = w, stderr = w)
    except CalledProcessError as err:
        print(err)


def main(bed, sample, bam_list, gtf, output, n_jobs = 10):
    u"""
    Main
    """
    if not os.path.exists(output):
        os.makedirs(output)
    print("Reading sample")
    # with open(sample) as r:
    #     genes = set([x.strip() for x in r])
    
    genes = sample.split(",")

    # genes = [
    #    "CYP2B7P", "PIGR", "MUC1",
    #    "SPINK1", "CKS1B", "SFTPA1",
    #    "SFTPA2", "IGHG1", "CCND1",
    #    "UBB", "H2AFZ", "IGFBP7",
    #    "VIM"
    #]

    #print(genes)
    reference = load_reference(bed)

    tasks = []
    for i in genes:
        if i not in reference.keys():
            continue
        print(i)
        tasks.append("python /mnt/raid61/Personal_data/zhangyiming/code/pysashimi/main.py normal -e {} -b {} -g {} --color-factor 3 -o {} -p {} --remove-empty-gene --no-gene --title {} --share-y".format(  # --plot-by 4  --sep-by-color --distance-ratio 0.35
            reference[i].promoter(),
            bam_list, gtf,
            os.path.join(output, "{}.txt".format(i)), n_jobs,
            i
        ))


    with Pool(min(n_jobs, 6)) as p:
        p.map(call, tasks)


if __name__ == '__main__':
    import pandas as pd
    
    # data = pd.read_excel("LungCancer10x/09_bulk/RNA_seq/DEGs/Res/tf.xlsx")
    # data = "LungCancer10x/11_CNV/each_cells/ATII/cluster_markers_LUAD.csv"
    # data = pd.read_csv(data)
    # data = data.loc[data["ident"].isin([1, 2, 3, 4]), :]
    
    # data.sort_values(['ident,avg_logFC'], ascending=[1,0], inplace=True)
    
    # data = data.groupby("ident").head(50)
    
    # print(data.head())
    
    main(
        bed = "LungCancerData/genome/gencode.v32.annotation.bed",
        sample = "KLF5", # ",".join(set(data["gene"])) ADAR,DSP,ELF3,ETS2,H1F0,HES1,HMGA1,HMGB2,HMGB3,HOPX,JUND,KLF6,LGR4,MYC,NKX2-1,PLXNB2,RNF213,SOX4,TAX1BP1,ID1,ID2,ID3,IFI16,HLA-DRB1,ZFP36,ZFP36L1,ZFP36L2,TOX4",
        bam_list = "LungCancer10x/tests/all_bam.txt",
        gtf = "/mnt/raid63/LungCancerData/genome/gencode.v32.annotation.gtf",
        output = "LungCancer10x/11_CNV/gnet"
    )