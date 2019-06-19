#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.19
make tmod per stage. multiple results required

- Utest
- CERNO
- Ztest
"""

import os
from glob import glob
import sys
from multiprocessing import Pool
from subprocess import check_call

import pandas as pd

import xlrd
from tqdm import tqdm


def call(args):
    u"""
    Run rscript
    :param args: input xlsx, output dir, path to msig.xml
    """
    
    xlsx, output, msig = args

    rscript = """
    xlsx = "{0}"
    output = "{1}"
    msig = "{2}"
    """.format(
        xlsx, output, msig
    )

    rscript += """
    library(tmod)
    library(openxlsx)
    
    msig = msig <- tmodImportMSigDB(msig)
    sel <- msig$MODULES$Category %in% c("H","C5","C2")
    
    
    clt = read.xlsx(xlsx)
    
    
    resU = NULL
    resZ = NULL
    resC = NULL
    
    for (i in unique(clt$ident)) {
        temp = clt[clt$ident == i, ]
        
        temp_genes = temp[order(temp$avg_logFC, decreasing = T), "gene"]
                
        temp = tmodUtest(temp_genes, mset = msig[sel], qval = Inf)    
        
        if(nrow(temp) > 0) {
            temp$ident = i
            resU = rbind(resU, temp)
        }
        
        
        temp = tmodCERNOtest(temp_genes, mset = msig[sel], qval = Inf)    
    
        if(nrow(temp) > 0) {
            temp$ident = i
            resC = rbind(resC, temp)
        }
        
        temp = tmodZtest(temp_genes, mset = msig[sel], qval = Inf)    
    
        if(nrow(temp) > 0) {
            temp$ident = i
            resZ = rbind(resZ, temp)
        }
    }
    
    wb = createWorkbook()
    addWorksheet(wb, "Utest")
    writeData(wb, 1, resU)
    
    addWorksheet(wb, "CERNO")
    writeData(wb, 2, resC)
    
    addWorksheet(wb, "Ztest")
    writeData(wb, 3, resZ)
    
    saveWorkbook(wb, file = output, overwrite = T)
    """    
    
    temp = output + ".R"
    
    with open(temp, "w+") as w:
        w.write(rscript)
        
    check_call("Rscript {0}".format(temp), shell = True)
    
    os.remove(temp)
    
    
def read_for_merge(path, output):
    u"""
    merge results from different cells into one
    """
    res = {}
    
    cell = os.path.basename(i).split(".")[0]    
    sheets = xlrd.open_workbook(i, on_demand=True).sheet_names()

    for sheet in sheets:
        data = pd.read_excel(i, sheet_name=sheet, engine="python")
        
        data["cell"] = [cell for x in range(data.shape[0])]
        
        res[sheet] = data
        
    return res
        
    
def main(input_dir, output_file, msig, n_jobs=10):
    u"""
    Main function
    :param input_dir: path to input directory
    :param output_file: path to output file
    :param n_jobs: use how many cpus
    """
    
    # files = glob(os.path.join(os.path.abspath(input_dir), "*_ADC/annotation_results_by_stage.xlsx"))
    # files += glob(os.path.join(os.path.abspath(input_dir), "*_SCC/annotation_results_by_stage.xlsx")) 
    
    output_file = os.path.abspath(output_file)
    temp_output_dir = output_file + ".temp"
    
    if not os.path.exists(temp_output_dir):
        os.makedirs(temp_output_dir)
    
    # tasks = []
    # for i in files:
        
    #     cell_name = os.path.basename(os.path.dirname(i))
        
    #     tasks.append([
    #         os.path.abspath(i),
    #         os.path.join(temp_output_dir, cell_name + ".xlsx"),
    #         os.path.abspath(msig)
    #     ])
    
    # with Pool(n_jobs) as p:
    #     p.map(call, tasks)
            
    files = [os.path.join(path, x) for x in os.listdir(path)]
    
    with Pool(n_jobs) as p:
        temp = p.map(read_for_merge, files)
        
    res = {}
    
    for i in temp:
        for k, v in i.items():
            temp_data = res.get(k, [])
            temp_data.append(v)
            res[k] = temp_data
            
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    
    for key, data in res.items():
        pd.concat(data).to_excel(writer, sheet_name = key)
    
    writer.save()
    
    
    # os.removedirs(temp_output_dir)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)
