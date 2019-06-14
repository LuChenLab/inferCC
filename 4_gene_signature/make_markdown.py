#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Collect image to different markdown and to pdf file
"""
import os
import glob
from fire import Fire
from subprocess import check_call
from multiprocessing import Process, Pool


class Directory2Md(object):
    u"""
    directory file to markdown
    """

    def __init__(self, path:str, outfile:str):
        u"""
        path to file
        :param path:
        """

        self.path = path

        # "Imbalanced", "ENN",
        self.methods = ["ICA", "Mfuzz", "WGCNA"]
        self.clean = ["Imbalanced", "ENN", "RENN"]
        self.outfile = outfile

    def image_to_table(self, clean, method, cluster):
        u"""
        convert image to markdown
        :param clean:
        :param method:
        :param cluster:
        :return:
        """
        content = ["## {0} {1}".format(method, cluster)]

        heatmap1 = os.path.join(
            self.path,
            "heatmaps/{0}/{0}_{1}_{2}_cluster.png".format(clean, method, cluster)
        )
        heatmap2 = os.path.join(
            self.path,
            "heatmaps/{0}/{0}_{1}_{2}_stage.png".format(clean, method, cluster)
        )
        violin = os.path.join(
            self.path,
            "violin/{0}/{1}_{2}_violin.png".format(clean, method, cluster)
        )
        line = os.path.join(
            self.path,
            "violin/{0}/{1}_{2}_lineplot.png".format(clean, method, cluster)
        )
        venn = os.path.join(
            self.path,
            "venn/{0}/{1}_{2}.png".format(clean, method, cluster)
        )
        violin1 = os.path.join(
            self.path,
            "violin/{0}/{1}_{2}_violin_value.png".format(clean, method, cluster)
        )

        if not os.path.exists(line) or not os.path.exists(violin):
            return ""

        content.append("### Heatmaps")
        content.append("| By cluster | By Stage |")
        content.append("| :--- | ---: |")
        content.append("| ![]({0}) | ![]({1}) |".format(heatmap1, heatmap2))
        content.append("")

        content.append("### Venn")
        content.append("Overlap with known gene set")
        content.append("![]({0})".format(venn))

        content.append("### Violin plots")
        content.append("| zscore | Normalized_counts |")
        content.append("| :--- | ---: |")
        content.append("| ![]({0}) | ![]({1}) |".format(violin, violin1))
        content.append("")

        content.append("### line plots")
        content.append("![]({0})".format(line))

        return "\n".join(content)

    def to_md(self, outfile=None):

        if outfile is None:
            outfile = self.outfile

        temp_mds = []
        temp_htmls = []

        for i in self.clean:
            content = []
            content.append("# {0}".format(i))

            for j in self.methods:
                heatmap = os.path.join(self.path, "{0}_{1}_heatmap".format(i, j))

                if not os.path.exists(heatmap):
                    continue

                for k in sorted([int(x.strip().split(".")[0]) for x in os.listdir(heatmap) if not x.startswith(".")]):

                    content.append(self.image_to_table(clean=i, method=j, cluster=k))

                content.append("---")
                content.append("")

                content.append("---")
                content.append("")

            temp_mds.append("{0}_{1}.md".format(outfile, i))
            temp_htmls.append("{0}_{1}.html".format(outfile, i))
            with open(temp_mds[-1], "w+") as w:
                w.write("\n".join(content))

            check_call(
                "python /mnt/raid61/Personal_data/zhangyiming/mdconverter/converter.py -md {0} -html {1} -pdf {2}".format(
                    temp_mds[-1],
                    temp_htmls[-1],
                    "{0}_{1}.pdf".format(outfile, i)
                ),
                shell=True
            )

            os.remove(temp_mds[-1])
            os.remove(temp_htmls[-1])
        # os.rename(os.path.join(self.path, "images.pdf"), outfile)


def mfuzz_gene_module_to_md(args):
    u"""

    :param input_dir:
    :param output_path:
    :return:
    """
    input_dir, output_path = args
    cell_name = os.path.basename(os.path.dirname(input_dir))

    heatmaps = glob.glob(os.path.join(input_dir, "heatmap_*.png"))

    dotplots = glob.glob(os.path.join(input_dir, "dotplot_*.png"))
    # print(heatmaps)
    # heatmaps = sorted(heatmaps, key=lambda x: int(os.path.basename(x).replace("heatmap_", "").split(".")[0]))
    # dotplots = sorted(dotplots, key=lambda x: int(os.path.basename(x).replace("dotplot_", "").split(".")[0]))

    md = os.path.join(input_dir, "temp.md")
    html = os.path.join(input_dir, "temp.html")

    with open(md, "w+") as w:
        w.write(f"# {cell_name}\n\n")

        for i in sorted(heatmaps):
            w.write(f"![]({i})\n\n")

        for i in sorted(dotplots):

            # if "stage" in i or "cluster" in i:
            #     continue

            w.write(f"![]({i})\n\n")

        w.write(f"![]({os.path.join(input_dir, 'kegg.png')})\n\n")

        w.write(f"![]({os.path.join(input_dir, 'BP.png')})\n\n")
        w.write(f"![]({os.path.join(input_dir, 'CC.png')})\n\n")
        w.write(f"![]({os.path.join(input_dir, 'MF.png')})\n\n")

        w.write(f"![]({os.path.join(input_dir, 'do.png')})\n\n")

    check_call(
        "python /mnt/raid61/Personal_data/zhangyiming/mdconverter/converter.py -md {0} -html {1} -pdf {2}".format(
            md,
            html,
            "{0}.pdf".format(os.path.join(output_path, cell_name))
        ),
        shell=True
    )

    os.remove(md)
    os.remove(html)


def call(path, output):
    Directory2Md(path, output).to_md()


def main(path, output):
    files = glob.glob(os.path.join(os.path.abspath(path), "*/cluster_gene_module"))

    output = os.path.abspath(output)
    if not os.path.exists(output):
        os.makedirs(output)

    tasks = []
    for i in files:
        print(i)
        cell_name = os.path.basename(os.path.dirname(i))

        # if cell_name != "CD4_SCC":
        #     continue

        print(i)
        # if cell_name == "Alveolar_II":
        # Directory2Md(i).to_md(os.path.join(output, cell_name) + ".pdf")
        # tasks.append(Process(target=mfuzz_gene_module_to_md, args=[i, output]))
        # tasks.append(Directory2Md()

        tasks.append([i, output])

        # mfuzz_gene_module_to_md(i, output)

    # tasks = [Process(x.to_md) for x in tasks]
    # [x.start() for x in tasks]
    # [x.join() for x in tasks]

    with Pool(10) as p:
        p.map(mfuzz_gene_module_to_md, tasks)


if __name__ == '__main__':
    Fire(main)

