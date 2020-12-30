#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys
import os
import pandas as pd
from glob import glob
from tqdm import tqdm


def main(param=sys.argv[1:]):
    xlsx = glob(os.path.join(param[0], "*.xlsx"))

    msd, logfc, pval = [], [], []
    for i in tqdm(xlsx):
        label = os.path.basename(i).replace(".xlsx", "")
        temp = pd.read_excel(i, sheet_name=1)
        temp["module"] = [label for _ in range(temp.shape[0])]
        msd.append(temp)

        temp = pd.read_excel(i, sheet_name=2)
        temp["module"] = [label for _ in range(temp.shape[0])]
        logfc.append(temp)

        temp = pd.read_excel(i, sheet_name=3)
        temp["module"] = [label for _ in range(temp.shape[0])]
        pval.append(temp)

    writer = pd.ExcelWriter(param[1], engine='xlsxwriter')

    pd.concat(msd).to_excel(writer, sheet_name="MSD")
    pd.concat(logfc).to_excel(writer, sheet_name="logFC")
    pd.concat(pval).to_excel(writer, sheet_name="p_adj_val")

    writer.save()


if __name__ == '__main__':
    main()

