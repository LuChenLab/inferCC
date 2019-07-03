#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Get all different expressed genes from GEPIA


http://gepia2.cancer-pku.cn/assets/PHP2/differential_genes.php?methodoption=ANOVA&dataset=ACC&fccutoff=0&type=gene&qcutoff=1
http://gepia2.cancer-pku.cn/assets/PHP2/differential_genes.php?methodoption=LIMMA&dataset=ACC&fccutoff=0&type=gene&qcutoff=1
http://gepia2.cancer-pku.cn/assets/PHP2/differential_genes.php?methodoption=TOP10&dataset=ACC&fccutoff=-5&type=gene&qcutoff=-5
"""
import os
import json
from bs4 import BeautifulSoup
from tqdm import tqdm
import pandas as pd

import requests
import gepia


THEAD = [
    "Symbol",
    "ID",
    "Median_normal",
    "Median_cancer",
    "logFC",
    "adjp",
    "disease"
]


def decode_data(methods, disease):
    u"""
    decode url
    """
    res = []
    
    if methods == "TOP10":
        url = "http://gepia2.cancer-pku.cn/assets/PHP2/differential_genes.php?methodoption=TOP10&dataset={0}&fccutoff=-5&type=gene&qcutoff=-5".format(disease)
    else:
        url = "http://gepia2.cancer-pku.cn/assets/PHP2/differential_genes.php?methodoption={0}&dataset={1}&fccutoff=-5&type=gene&qcutoff=1".format(methods, disease)
    
    data = requests.get(url)
    try:
        
        soup = BeautifulSoup(data.json()["output"], "lxml")
        tbody = soup.find("tbody")
        
        for i in tbody.find_all("tr"):
            temp = {k: v.get_text() for k, v in zip(THEAD, i.find_all("td"))}
            temp["method"] = methods
            temp["disease"] = disease
            
            res.append(temp)
    except json.decoder.JSONDecodeError as err:
        tqdm.write(str(err))
        
    return res
            
def main(output_data):
    u"""
    
    """
    res = []
    for i in tqdm(gepia.CANCERType):
        tqdm.write(i)
        for j in ["ANOVA", "LIMMA", "TOP10"]:
            res += decode_data(j, i)
    
    with open(output_data, "w+") as w:
        json.dump(res, w, indent=4)
    
    # with open(output_data) as r:
    #     res = json.load(r)
            
    anova = [x for x in res if x["method"] == "ANOVA"]
    anova = pd.DataFrame(anova)
    
    limma = [x for x in res if x["method"] == "LIMMA"]
    limma = pd.DataFrame(limma)
    
    top10 = [x for x in res if x["method"] == "TOP10"]
    top10 = pd.DataFrame(top10) 
    
    genes = [
        "ALKBH3",
        "ALOX5AP",
        "ALKBH3",
        "ALOX5AP",
        "BRAF",
        "MET",
        "ESRP1",
        "KRAS",
        "LSP1",
        "PSMB4",
        "RBPJ",
        "SSNA1",
        "TIGIT",
        "TNFRSF18",
        "TNFRSF4",
    ]
    
    temp_head = [x for x in THEAD]
    temp_head[-2] = "Percentage"
        
    anova_genes = [x for x in res if x["method"] == "ANOVA" and x["disease"] in ["LUAD", "LUSC"] and x["Symbol"] in genes]
    limma_genes = [x for x in res if x["method"] == "LIMMA" and x["disease"] in ["LUAD", "LUSC"] and x["Symbol"] in genes]
    top10_genes = [x for x in res if x["method"] == "TOP10" and x["disease"] in ["LUAD", "LUSC"] and x["Symbol"] in genes]
    
    writer = pd.ExcelWriter(output_data.replace("json", "xlsx"))
    anova[THEAD].to_excel(writer,'ANOVA', index=False)
    limma[THEAD].to_excel(writer,'LIMMA', index=False)
    top10 = top10[THEAD]
    top10.columns = temp_head
    top10.to_excel(writer, "TOP10", index=False)
    
    pd.DataFrame(anova_genes)[THEAD].to_excel(writer, "ANOVA_required", index=False)
    pd.DataFrame(limma_genes)[THEAD].to_excel(writer, "LIMMA_required", index=False)
    top10_genes = pd.DataFrame(top10_genes)[THEAD]
    top10_genes.columns = temp_head
    top10_genes.to_excel(writer, "TOP10_required", index=False)
    
    writer.save()
        

if __name__ == '__main__':
    from fire import Fire
    Fire(main)

