#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by Zhang yiming at 2019.01.15

Run Seurat in batch
"""
import argparse as ap
import os
from glob import glob
import sys
import subprocess as sp
import multiprocessing as mp
import gzip
from tqdm import tqdm
from shutil import copytree, rmtree, copyfile
from openpyxl import load_workbook
from pypinyin import lazy_pinyin


__dir__ = os.path.dirname(os.path.abspath(__file__))
__script__ = os.path.join(__dir__, "run_seurat.R")
__script2__ = os.path.join(__dir__, "run_seurat2.R")
__static__ = os.path.join(__dir__, "static")


class RunSeurat:
    u"""
    Class to handle all issues
    """

    def __init__(self, args):
        u"""
        init this class
        :param args: argparse namespace
        """

        self.input = args.input

        if not os.path.exists(self.input):
            raise ValueError("%s not exists" % self.input)

        self.output = args.output

        if not os.path.exists(self.output):
            os.makedirs(self.output)

        self.xlsx = args.xlsx
        self.processes = args.processes

    @staticmethod
    def run_seurat_clustering_single(args):
        u"""
        run seurat
        :param args: {"input": , "output":}
        :return:
        """
        if not os.path.exists(args["output"]):
            os.makedirs(args["output"])

        cmds = [
            "Rscript",
            __script__,
            "-i {input} -o {output} -l {label}".format(**args)
        ]

        with open(os.devnull, "w") as w:
            try:
                sp.check_call(" ".join(cmds), shell=True, stdout=w, stderr=w)
            except sp.CalledProcessError:
                print(" ".join(cmds))

    def run_seurat_clustering(self):
        u"""
        Run seurat in batch
        :return:
        """
        args = []
        for i in os.listdir(self.input):
            args.append({
                "input": os.path.abspath(os.path.join(self.input, i)),
                "output": os.path.abspath(os.path.join(self.output, i)),
                "label": i
            })

        with mp.Pool(processes=self.processes) as pool:
            list(tqdm(pool.imap_unordered(self.run_seurat_clustering_single, args), total=len(args)))

    @staticmethod
    def run_seurat_comparison_single(args):
        cmds = [
            "Rscript",
            __script2__,
            "-c {cancer} -n {normal} -o {output} -l {label}".format(**args)
        ]

        with open(os.devnull, "w") as w:
            try:
                sp.check_call(" ".join(cmds), shell=True, stdout=w, stderr=w)
            except sp.CalledProcessError:
                print(" ".join(cmds))

    def run_seurat_comparison(self):
        u"""
        Run seurat comparison
        :return:
        """
        wb = load_workbook(self.xlsx)
        ws = wb[wb.sheetnames[0]]

        data = {}

        header = {
            "Name": None,
            "tissue_type": None,
            "status": None,
            "patient": None
        }
        for row_idx, row in enumerate(ws.rows):
            if row_idx == 0:
                for col_idx, col in enumerate(row):
                    if col.value in header.keys():
                        header[col.value] = col_idx
            elif not row[0].value:
                continue
            else:
                if row[header["status"]].value == 0:
                    continue

                if row[header["patient"]].value not in data.keys():
                    tmp = {}
                else:
                    tmp = data[row[header["patient"]].value]

                tmp[row[header["tissue_type"]].value.strip()] = row[header["Name"]].value

                data[str(row[header["patient"]].value.strip())] = tmp

        inputs = {}
        for i in os.listdir(self.input):
            inputs[i] = os.path.abspath(os.path.join(self.input, i))

        args = []
        for key, value in data.items():
            if len(value) < 2:
                continue

            key = "_".join(lazy_pinyin(key))

            out_dir = os.path.join(self.output, key)

            try:
                normal = value["normal"]

                for k, v in value.items():
                    if k == "normal":
                        continue

                    k = k.replace(" ", "_")

                    real_out_dir = os.path.join(out_dir, k)
                    if not os.path.exists(real_out_dir):
                        os.makedirs(real_out_dir)

                    args.append({
                        "normal": inputs[normal],
                        "cancer": inputs[v],
                        "label": k,
                        "output": real_out_dir
                    })

            except KeyError:
                continue

        with mp.Pool(processes=self.processes) as pool:
            list(tqdm(pool.imap_unordered(self.run_seurat_comparison_single, args), total=len(args)))


class Report:
    u"""
    Make html
    """

    def __init__(self, args):

        self.input = args.input

        if not os.path.exists(self.input):
            raise ValueError("%s not exists" % self.input)

        self.output = args.output

        if not os.path.exists(self.output):
            os.makedirs(self.output)
        self.xlsx = args.xlsx

    def generate_report(
            self,
            plots,
            columns
    ):
        u"""
        Generate report
        :return:
        """

        target = os.path.join(self.output, "static")

        if os.path.exists(target):
            rmtree(target)

        copytree(__static__, target)

        home = self.__generate_home_page__(columns)

        with open(os.path.join(self.output, "index.html"), "w+") as w:
            w.write(home)

        out_html = os.path.join(self.output, "html")
        if not os.path.exists(out_html):
            os.makedirs(out_html)

        out_img = os.path.join(self.output, "img")
        if not os.path.exists(out_img):
            os.makedirs(out_img)

        for i in os.listdir(self.input):
            full_path = os.path.join(self.input, i)

            for j in os.listdir(full_path):
                if j.endswith("png"):

                    target = os.path.join(out_img, i)

                    if not os.path.exists(target):
                        os.makedirs(target)

                    copyfile(
                        os.path.join(full_path, j),
                        os.path.join(target, j)
                    )

            with open(os.path.join(out_html, i + ".html"), "w+") as w:
                w.write(self.__generate_single_page__(i, plots))

    def __generate_home_page__(self, columns):
        u"""
        :param columns: list of ints, to which columns are kept of xlsx
        :return:
        """
        head = []
        for i in sorted([x for x in os.listdir(__static__) if os.path.isfile(os.path.join(__static__, x))]):
            if i.endswith("js"):
                head.append(
                    '<script type="text/javascript" src="./static/%s"></script>' % i
                )
            else:
                head.append(
                    '<link rel="stylesheet" href="./static/%s">' % i
                )
        head = "\n".join(head)

        html = """
        <html>
            <head>
                <meta charset="utf-8">
                <title>Seurat report</title>
                %(head)s
            </head>
    
            <body>
                <div class="navbar navbar-expand-lg fixed-top navbar-dark bg-primary">
                    <div class="container">
                        <a href="." class=".">Home</a> 
                    </div>
                </div>
                <div class="container">
                    <hr>
                    <div class="bs-docs-section">
                        <div class="page-header">
                            <div class="row">
                                <div class="col-lg-12">
                                    <h3>统计信息</h3>
                                </div>
                            </div>
                        </div>
                    
                        <div class="row">
                            <div class="col-lg-12">
                                <div class="page-header">
                                    <h4 id="tables">测序质量</h4>
                                </div>
                                <div class="bs-component">
                                    <table id="table" class="table table-hover">
                                    <thead>
                                        <tr>
                                        %(thead)s
                                        </tr>
                                    </thead>
                                    <tbody>
                                        %(tbody)s
                                    </tbody>
                                    </table>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <script type="text/javascript">
                    $(document).ready(function() {
                        $("#table").DataTable();
                    } );
                </script>
            </body>
        </html>
        """

        wb = load_workbook(self.xlsx)
        ws = wb[wb.sheetnames[0]]

        data = {"thead": "", "tbody": "", "head": head}

        for row in range(1, ws.max_row + 1):

            if not ws[row][0].value:
                continue

            tmp = "<tr>\n%s</tr>\n"
            tbody_row = []
            for col in columns:
                col -= 1

                if row == 1:
                    data["thead"] += "<th>%s</th>" % ws[row][col].value
                else:
                    if col == 0:
                        tbody_row.append(
                            '<td><a href="./html/{0}.html">{0}</a></td>'.format(ws[row][col].value)
                        )
                    else:
                        tbody_row.append(
                            "<td>%s</td>" % ws[row][col].value
                        )
            if tbody_row:
                data["tbody"] += tmp % "\n".join(tbody_row)
        return html % (data)

    @staticmethod
    def __generate_single_page__(indir, plots):
        u"""
        :param plots: {fig name: fig alias}
        :return:
        """
        head = []
        for i in sorted([x for x in os.listdir(__static__) if os.path.isfile(os.path.join(__static__, x))]):
            if i.endswith("js"):
                head.append(
                    '<script type="text/javascript" src="../static/%s"></script>' % i
                )
            else:
                head.append(
                    '<link rel="stylesheet" href="../static/%s" >' % i
                )

        head = "\n".join(head)

        label = os.path.basename(indir)

        html = """
            <html>
                <head>
                    <meta charset="utf-8">
                    <title>%(label)s</title>
                    %(head)s
                </head>
    
                <body>
                    <div class="navbar navbar-expand-lg fixed-top navbar-dark bg-primary">
                      <div class="container">
                         <a href="../index.html" class=".">Home</a> 
                      </div>
                    </div>
                    <hr>
                    <div class="container">
                        <div class="page-header">
                            <div class="row">
                                <div class="col-lg-12">
                                    <h3>%(label)s</h3>
                                </div>
                            </div>
                        </div>
                        <div class="bs-docs-section">
                            <div class="row">
                                <div class="col-lg-12">
                                    <div class="page-header">
                                        <h4>%(label)s</h4>
                                    </div>
                                    <div class="bs-component">
                                        %(content)s
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </body>
            </html>
        """
        content = []
        for i in plots.keys():
            content.append(
                """
                <div class="bs-docs-section">
                    <div class="row">
                        <div class="col-lg-12">
                            <div class="page-header">
                                <h5>%(label)s</h5>
                            </div>
                            <div class="bs-component">
                                %(img)s
                            </div>
                        </div>
                    </div>
                </div>
                """ % ({
                    "img": '<img width="90%" height="800" src="../img/{0}/{1}" />'.format(label, i),
                    "label": plots[i]
                })
            )

            content.append("<hr />")

        return html % ({
            "head": head,
            "label": label,
            "content": "\n".join(content)
        })


def __decompress__(args):
    f, outdir, lock = args

    sample_id = os.path.basename(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(f)
            )
        )
    )

    out_dir = os.path.join(outdir, sample_id)

    lock.acquire()
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    lock.release()

    basename = os.path.basename(f).replace(".gz", "")
    if basename == "features.tsv":
        basename = "genes.tsv"
    outfile = os.path.join(out_dir, basename)

    with open(outfile, "w+") as w:
        with gzip.open(f, "rb") as r:
            w.write(str(r.read().decode("utf-8")))


def decompress(indir, outdir, processes):
    u"""
    Compress data
    :param indir:
    :param outdir:
    :return:
    """
    input_files = glob(os.path.join(indir, "*/outs/filtered_feature_bc_matrix/*"))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    args = []
    m = mp.Manager()
    lock = m.Lock()
    for f in input_files:
        args.append([f, outdir, lock])

    with mp.Pool(processes) as pool:
        list(tqdm(pool.imap_unordered(__decompress__, args), total=len(args)))


if __name__ == '__main__':
    parser = ap.ArgumentParser("Run Seurat")
    parser.add_argument(
        "-i", "--input",
        help="Path to input directory",
        type=str,
        required=True
    )
    parser.add_argument(
        "-o", "--output",
        help="Path to output directory",
        type=str,
        required=True
    )
    parser.add_argument(
        "-p", "--processes",
        default=1,
        help="CPU usage",
        type=int,
    )
    parser.add_argument(
        "-x",
        "--xlsx",
        help="Path to meta info xlsx",
        required=True,
        type=str
    )
    parser.add_argument(
        "--decompress",
        action="store_true",
        default=False,
        help="Compress files"
    )
    parser.add_argument(
        "--clustering",
        action="store_true",
        default=False,
        help="Do clustering rather than comparision"
    )
    parser.add_argument(
        "--html",
        action="store_true",
        default=False,
        help="Collect results into output html"
    )

    if len(sys.argv) <= 1:
        parser.print_help()
    else:
        args = parser.parse_args(sys.argv[1:])

        default_columns = [1, 2, 3, 4, 9, 10, 11, 12, 13, 14]

        if args.clustering:
            default_plots = {
                "VlnPlot.png": "Vln Plot",
                "GenePlot.png": "Gene Plot",
                "Variable.png": "Variable Genes",
                "VizPCA.png": "Viz PCA",
                "PCAPlot.png": "PCA Plot",
                "PCHeatmap.png": "PCA Heatmap",
                "JackStrawPlot.png": "Jack Straw Plot",
                "PCElbowPlot.png": "PCElBowPlot",
                "tSNE.png": "tSNE"
            }
        else:
            default_plots = {
                "CCA.png": "CCA",
                "MetageneBiorPlot.png": "Meta gene Bior Plot",
                "DimHeatmap.png": "Dim Heatmap",
                "AlignedCCAVlnPlot.png": "Aligned CCA Vln Plot",
                "tSNE.png": "tSNE",
                "feature.png": "feature",
                "tSNE_with_cell_name.png": "tSNE with cell name",
                "SpliceDotPlotGG.png": "Splice Dot Plot GG"
            }

        if args.html:
            Report(args).generate_report(
                columns= default_columns,
                plots=default_plots
            )
        elif args.decompress:
            decompress(
                args.input,
                args.output,
                args.processes
            )
        else:
            tmp = RunSeurat(args)

            if args.clustering:
                tmp.run_seurat_clustering()
            else:
                tmp.run_seurat_comparison()


