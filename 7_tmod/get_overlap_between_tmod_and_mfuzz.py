#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.12
This script is used to get overlap between tmod msig module and the mfuzz results
"""

import os
from glob import glob
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import pandas as pd
import xlrd


def read_mfuzz(path):
    u"""
    read mfuzz to dict
    :param path: path to mfuzz results xlsx file
    :return: {cell_name: {id: [genes]}}
    """
    data = pd.read_excel(path)

    cell_name = os.path.basename(os.path.dirname(os.path.dirname(path)))

    res = {}
    for _, row in data.iterrows():
        temp = res.get(row[0], [])
        temp.append(row[1])
        res[row[0]] = temp

    return {cell_name: res}


def read_msig(path: str):
    u"""
    read extracted msig file, first column is MODULE ID, second column is gene symbol

    no header
    :param path: path to file
    :return: dict {msig module id: [ gene symbols ]}
    """

    res = {}
    with open(path) as r:
        for line in r:
            lines = line.split()

            temp = res.get(lines[0], [])
            temp.append(lines[1])
            res[lines[0]] = temp

    return res


def read_tmod(path: str, auc=0.5, p_val=0.05):
    u"""
    read tmod results
    :param path:
    :param auc: low threshold of AUC
    :param p_val: upper threshold of adj.P.Value
    :return: {cell: {stages: {msig module id: msig module title}}
    """
    sheets = xlrd.open_workbook(path, on_demand=True).sheet_names()

    res = {}
    for sheet in tqdm(sheets):
        data = pd.read_excel(path, sheet_name=sheet)

        # ID	Title	U	N1	AUC	P.Value	adj.P.Val	ident	cell
        for _, row in data.iterrows():
            
            if row["AUC"] > auc and row["adj.P.Val"] < p_val:
                temp_cell = res.get(row["cell"], {})
                temp_modules = temp_cell.get(row["ident"], {})
                temp_modules[row["ID"]] = row["Title"]
                temp_cell[row["ident"]] = temp_modules
                res[row["cell"]] = temp_cell

    return res


class GeneModule(object):
    u"""
    Class to record gene module and prepare for further merge
    """
    
    def __init__(self, cell: str, msig_id, stage: str, mfuzz_id: str, genes):
        u"""
        init this class
        :param cell: cell name
        :param msig_id: msig module id and title in tuple()
        :param stage: stage
        :param mfuzz_id:
        :param genes: list or set of gene symbols
        """
        self.cell = cell
        self.msig_id = set([msig_id])
        self.stage = stage
        self.mfuzz_id = set([mfuzz_id])
        self.genes = set(genes)
        
    def __str__(self):
        u"""
        Convert this class to string
        """
        msig_id = sorted(self.msig_id, key=lambda x: x[0])
        
        return "\"{0}\",\"{1}\",\"{2}\",\"{3}\",\"{4}\",\"{5}\",\"{6}\"".format(
            self.cell, 
            "|".join(sorted(set([x[0] for x in msig_id]))),
            "|".join(sorted(set([x[1] for x in msig_id]))),    
            self.stage, 
            "|".join([str(x) for x in sorted(set(self.mfuzz_id))]),
            "|".join(sorted(self.genes)),
            len(self.genes)
        )
        
    def to_list(self):
        
        msig_id = sorted(self.msig_id, key=lambda x: x[0])
        
        return [
            self.cell, 
            "|".join(sorted(set([x[0] for x in msig_id]))),
            "|".join(sorted(set([x[1] for x in msig_id]))),    
            self.stage, 
            "|".join([str(x) for x in sorted(set(self.mfuzz_id))]),
            "|".join(sorted(self.genes)),
            len(self.genes) 
        ]

    def __gt__(self, other):
        u"""
        greater than
        By cell, stage and number of genes
        """
        if self.cell != other.cell:
            return self.cell > other.cell
        
        if self.stage != other.stage:
            return self.stage > other.stage
        
        if len(self.genes) == len(other.genes):
            return "|".join(sorted(self.genes)) > "|".join(sorted(other.genes))
        
        return len(self.genes) > len(other.genes)
    
    def __lt__(self, other):
        if self.cell != other.cell:
           return self.cell < other.cell
        
        if self.stage != other.stage:
           return self.stage < other.stage
        
        if len(self.genes) == len(other.genes):
           return "|".join(sorted(self.genes)) < "|".join(sorted(other.genes))
        
        return len(self.genes) < len(other.genes)
    
    def merge(self, other):
        u"""
        Merge other gene module to this one
        """
        if self.cell != other.cell or self.stage != other.stage:
            raise ValueError("cell or stage not match")
        
        self.msig_id |= other.msig_id
        self.mfuzz_id |= other.mfuzz_id
        self.genes |= other.genes
    
    def is_in_this(self, other):
        u"""
        Check if genes from other is all in this
        """
        return self.is_same_cell_stage(other) and len(self.genes & other.genes) == len(other.genes)
    
    def is_same_cell_stage(self, other):
        return self.cell == other.cell and self.stage == other.stage
    
def mergeGeneModule(module_list):
        u"""
        merge gene module by the genes
        
        if some gene modules genes all locate in another, then merge these two together
        :param module_list list of GeneModule 
        """
        res = []
        
        module_list.sort()
        
        for i in range(len(module_list)):
            curr_i = module_list[i]
            res.append(curr_i)
            for j in range(i + 1, len(module_list)):
                curr_j = module_list[j]
                if curr_j.is_same_cell_stage(curr_i):
                    if curr_j.is_in_this(curr_i):
                        curr_j.merge(curr_i)
                        res = res[:-1]
                        break
                else:
                    break
                
        return res
    
def mergeSecondRound(module_list):
    u"""
    Merge Second round, merge genes from same cell and same stage togheter
    """
    res = []
    
    res = {}
    for i in module_list:
        key = "{0}{1}".format(i.cell, i.stage)
        
        temp = res.get(key)
        if temp is None:
            temp = i
        else:
            temp.merge(i)
        res[key] = temp

    return sorted(list(res.values()))
    

def main(input_dir, msig, tmod, output, auc=0.5, p_val=0.05, n_jobs=10):
    u"""

    :param input_dir:
    :param msig:
    :param tmod:
    :param auc:
    :param p_val:
    :param output:
    :return:
    """
    files = glob(os.path.join(input_dir, "*/mfuzz_gene_module/results.xlsx"))
    mfuzz = {}

    with Pool(n_jobs) as p:
        temp = p.map(read_mfuzz, files)

    p.close()
    p.join()
    
    for i in temp:
        mfuzz.update(i)
    
    print("Read {0}".format(tmod))
    tmod = read_tmod(path=tmod, auc=auc, p_val=p_val)

    print("Read {0}".format(msig))
    msig = read_msig(msig)


    print("Start to combine")
    res = []

    for cell, data in tqdm(tmod.items()):
        for stage, module in tqdm(data.items()):
            for m, m_title in tqdm(module.items()):
                module_genes = msig.get(m, set())
                
                mfuzz_data = mfuzz.get(cell, {})

                for mfuzz_id, mfuzz_gene in mfuzz_data.items():
                    temp_genes = set(module_genes) & set(mfuzz_gene)

                    if len(temp_genes) > 0:
                        res.append(GeneModule(
                            cell=cell,
                            msig_id=(m, m_title), 
                            stage=stage,
                            mfuzz_id=mfuzz_id, 
                            genes=temp_genes
                        ))
                        
                        # res.append({
                        #     "cell": cell,
                        #     "msig": m,
                        #     "stage": stage,
                        #     "mfuzz": mfuzz_id,
                        #     "genes": "|".join(sorted(temp_genes))
                        # })

    col_names = [
        "Cell_name", "Msig_ID", 
        "Msig Title", "Stage", 
        "Mfuzz_ID", "Genes", 
        "Gene_num"
    ]

    first_round = pd.DataFrame([x.to_list() for x in mergeGeneModule(res)])
    second_round = pd.DataFrame([x.to_list() for x in mergeSecondRound(res)])
    
    first_round.columns = col_names
    second_round.columns = col_names
    
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    first_round.to_excel(writer, sheet_name='First round', index=False)
    second_round.to_excel(writer, sheet_name='Second round', index=False)
    
    writer.save()

    # print("Write {0}".format(output))
    # with open(output, "w+") as w:
    #     w.write()
        
    #     for i in tqdm(mergeGeneModule(res)):
    #         # w.write("{cell}\t{msig}\t{stage}\t{mfuzz}\t{genes}\n".format(**i))
    #         w.write("{0}\n".format(i))


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
