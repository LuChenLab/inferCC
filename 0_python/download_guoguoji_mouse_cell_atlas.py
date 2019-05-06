#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Create by ygidtu@gmail.com at 2019.01.18

"""
import argparse as ap
import os
import sys
import json
import requests
from bs4 import BeautifulSoup
import wget
from urllib.error import HTTPError
import pandas as pd


base_url = "http://bis.zju.edu.cn/MCA/search.html"
json_url = "http://bis.zju.edu.cn/MCA/data/tissues/%(tissue)s/mca_top_markers_%(tissue)s.json"

def main(args):
    u"""
    Main function
    :param args:
    :return:
    """

    output = os.path.abspath(args.output)

    out_dir = os.path.dirname(output)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_json = os.path.join(out_dir, "json")
    if not os.path.exists(out_json):
        os.makedirs(out_json)

    page_content = requests.get(base_url)
    soup = BeautifulSoup(page_content.content, "lxml")

    for tissue in soup.find("select", {"id": "tissue_input"}).find_all("option"):
        print(tissue)
        tissue = tissue["value"]

        tmp = json_url % {"tissue": tissue}

        out_path = os.path.join(out_json, os.path.basename(tmp))

        if os.path.exists(out_path):
            os.remove(out_path)

        try:
            wget.download(tmp, out_path)
        except HTTPError:
            continue

    data = []
    for i in os.listdir(out_json):
        input_json = os.path.join(out_json, i)

        print(input_json)

        try:
            with open(input_json) as w:
                tmp_data = json.load(w)
        except json.decoder.JSONDecodeError:
            continue

        for j in tmp_data["data"]:
            j["Tissue"] = i.split("_")[-1].replace(".json", "")
            data.append(j)

    data = pd.DataFrame(data)

    data.to_excel(output)


if __name__ == '__main__':

    parser = ap.ArgumentParser("Download Mouse Cell Atlas")
    parser.add_argument(
        "-o",
        "--output"
    )

    if len(sys.argv) < 2:
        parser.print_help()

    else:
        main(parser.parse_args(sys.argv[1:]))

