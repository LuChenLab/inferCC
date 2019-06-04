#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.03
@author Zhang yiming
A simple crawler of MedGen to collect phenotype or a disorder associate with gene
"""
import json
import os
from glob import glob
from multiprocessing import Pool

import bs4
import pandas as pd
import requests
from fire import Fire
from tqdm import tqdm


# from lxml import etree


def getPage(gene) -> list:
    u"""
    get data from NCBI
    :param gene:
    :return:
    """
    # print(gene)
    content = requests.get("https://www.ncbi.nlm.nih.gov/medgen/?term={0}".format(gene))

    return [{"gene": gene, "data": pageToJson(content.content)}]


def pageToJson(content) -> list:
    u"""
    Collect a page data into json
    :return: dict
    """
    soup = bs4.BeautifulSoup(str(content), "lxml")

    rprts = soup.find_all("div", {"class": "rslt"})

    # tree = etree.fromstring(soup.prettify())
    #
    # rprts = tree.xpath("//div[class='rprt']/div[class='rslt']")

    data = []
    for i in rprts:
        # title = i.xpath("/p[class='title']/a/text()")

        title = i.find("p", {"class": "title"}).find("a")

        ids = i.find("dl", {"class": "rprtid"}).findChildren(recursive=False)

        temp = {
            "title": title.text,
            "title_href": title.href,
            ids[0].text.replace(":", ""): ids[1].text,
            ids[2].find(text = True, recursive=False)[1].replace(":", "").strip(): ids[3].text
        }

        data.append(temp)

    return data


def main(path: str, output: str, n_jobs=3):
    u"""
    :param path: input directory
    :param output:
    :param n_jobs:
    :return:
    """
    files = glob(os.path.join(os.path.abspath(path), "*/*_gene_module/results.xlsx"))

    exists_data = []
    exists_genes = set()
    if os.path.exists(output):
        with open(output) as r:
            exists_data += json.load(r)
            exists_genes |= set([x["gene"] for x in exists_data])

    genes = set()
    for i in files:
        print(i)
        data = pd.read_excel(i)
        genes = genes | set(data["gene"])

    genes = genes - exists_genes

    with Pool(n_jobs) as p:
        temp_res = list(tqdm(p.imap_unordered(getPage, list(genes)), total=len(genes)))

    p.close()
    p.join()

    res = []
    for i in temp_res:
        res += i

    res += exists_data
    with open(output, "w+") as w:
        json.dump(res, w, indent=4)


if __name__ == '__main__':
    Fire(main)


