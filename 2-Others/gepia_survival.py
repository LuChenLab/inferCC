#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import json
import os
import time

import click
import requests

from homura import download

VALID = "http://gepia2.cancer-pku.cn/assets/PHP4/valid.php"
SURVIVAL = "http://gepia2.cancer-pku.cn/assets/PHP4/survival_zf.php"
OUTPUT = "http://gepia2.cancer-pku.cn/tmp/"


def valid(genes, disease):
    
    params = [
        ("q[0][]", "gene"),
        ("q[0][]", "\n".join(genes)),
        ("q[1][]", "gene"),
        ("q[1][]", ""),
        ("q[2][]", "ds1"),
        ("q[2][]", disease)
    ]
    req = requests.post(VALID, data = params)
    
    data = req.json()
    
    for i in data:
        if not i['success']:
            return i["failure"]
    return None

def post(genes, disease, groupcutoff1 = 75, groupcutoff2 = 25):
    u"""

    """
    
    params = {
            "methodoption": "os",
            "dataset": "\n".join(disease),
            "signature": "\n".join(genes),
            "highcol": "#ff0000",
            "lowcol": "#0000ff",
            "groupcutoff1": groupcutoff1,
            "groupcutoff2": groupcutoff2,
            "axisunit": "month",
            "ifhr": "hr",
            "ifconf": "conf",
            "signature_norm": "",	
            "is_sub": False,
            "subtype": ""
    }
    
    data = requests.post(SURVIVAL, params).json()
        
    return data


@click.command(
    context_settings=dict(help_option_names=['-h', '--help']),
)
@click.version_option("0.0.1", message="Current version %(version)s")
@click.option(
    "-i",
    "--input-genes",
    default=None,
    type=click.STRING,
    help="path to gene list or comma seperated gene list",
    show_default=True
)
@click.option(
    "-o",
    "--output",
    default=os.path.dirname(__file__),
    type=click.STRING,
    help="Output directory",
    show_default=True
)
@click.option(
    "-d",
    "--disease",
    default="LUAD",
    type=click.STRING,
    help="comma seperated disease",
    show_default=True
)
@click.option(
    "-c",
    "--cutoff",
    default=75,
    type=click.IntRange(min=0, max=100),
    help="expression cutoff",
    show_default=True
)
@click.option(
    '--gene-set', 
    is_flag=True,
    default=False,
    help="By gene set",
    show_default=True
)
def main(input_genes, output, disease, cutoff, gene_set):
    u"""
    main function
    """    
    is_input_file = None
    if os.path.exists(input_genes):
        is_input_file = input_genes
        with open(input_genes) as r:
            input_genes = [x.strip() for x in r.read().split("\n")]
            input_genes = [x for x in input_genes if x]
    else:
        input_genes = input_genes.split(",")
    
    print("Input genes: {}".format(len(input_genes)))
    check = valid(input_genes, disease)
    while check is not None:
        input_genes.remove(check)
        print("Input genes: {}; {} is invalid".format(len(input_genes), check))
        time.sleep(0.5)
        check = valid(input_genes, disease)
        
    if is_input_file is not None:
        with open(is_input_file, "w+") as w:
            for i in sorted(input_genes):
                w.write("{}\n".format(i))
    
    if not input_genes:
        print("No valid genes")
        exit(0)
    
    if gene_set:
        data = post(genes = input_genes, disease = disease.split(","), groupcutoff1=cutoff, groupcutoff2 = 100 - cutoff)
        
        if data["outdir"] == "fail":
            print("Failed {}".format(data["outdir"]))
            exit(data)
        
        outdir = os.path.dirname(output)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        if os.path.exists(output):
            os.remove(output)
    
        download(url = os.path.join(OUTPUT, data["outdir"]), path = output)

    else:
        if not os.path.exists(output):
            os.makedirs(output)
            
        for i in input_genes:
            print(i)
            data = post(genes = [i], disease = disease.split(","), groupcutoff1=cutoff, groupcutoff2 = 100 - cutoff)
            
            if data["outdir"] == "fail":
                print("Failed {}".format(data["outdir"]))
                continue
            
            outfile = os.path.join(output, data["outdir"])
            if os.path.exists(outfile):
                os.remove(outfile)

            download(url = os.path.join(OUTPUT, data["outdir"]), path = outfile)
    

if __name__ == '__main__':
    main()
    