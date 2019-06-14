#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2019.04.12
"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
from collections.abc import Iterable
from itertools import chain
import numpy as np
from tqdm import tqdm


class BasicBar(object):
    u"""
    负责做单个的百分比柱子
    """

    def __init__(
            self,
            ax,
            values,
            label="",
            ylabel="",
            colors=None,
            no_x_axis=True,
            boxplot=False,
            height=1,
            font_size=10
    ):
        u"""
        初始化该类
        :param ax: matplotlib axis
        :param values: a list of int/float or a single float/int
        :param label: the label show in legend
        :param ylabel: the lable shoed in y axis
        :param colors: a list of RGB color or a string represented RGB value
        :param no_x_axis: whether to plot x axis
        :param boxplot: whether to make boxplot rather than bar plot
        :param height: height of single bar, default 1
        :param xlim: single value, set the x axis limitation
        """
        self.values = values

        if colors is not None:
            if isinstance(colors, str):
                colors = [colors]
            elif len(colors) < len(values):
                raise ValueError("Input colors less than values")

        self.colors = colors

        self.no_x_axis = no_x_axis
        self.boxplot = boxplot
        self.ax = ax
        self.label = label
        self.height = height
        self.font_size = font_size

        self.plot()

    def legend(self, *args, **wargs):
        u"""
        wrapper of matplotlib.axis
        """
        self.ax.legend(*args, **wargs)

    def show_x_axis(self):
        u"""
        show x axis
        """
        self.no_x_axis = False
        self.set_border()

    def set_xlim(self, xmax):
        u"""
        Set x limitation
        :param xmax: int/float
        """
        self.ax.spines['bottom'].set_bounds(0, xmax)

        self.ax.set_xlim(0, xmax)

    def set_ylabel(self, distance=None):
        u"""
        set ylabel
        :param distance: distance between plot and y label
        """
        self.ax.set_ylabel(
            self.label,
            rotation=0,
            va="center",
            fontsize=self.font_size,
            labelpad=distance
        )

    def set_border(self):
        u"""
        as function name said
        """
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['bottom'].set_visible(not self.no_x_axis)
        self.ax.spines['left'].set_visible(False)

        self.ax.tick_params(
            top=False,
            bottom=not self.no_x_axis,
            left=False,
            right=False,
            labelleft=False,
            labelbottom=not self.no_x_axis
        )

    def make_percentage_barplot(self):
        u"""
        make percentage stacked barplot
        """

        idx, curr = 0, 0
        for k, v in self.values.items():
            if self.colors:
                self.ax.fill_between([curr, curr + v], y1=0, y2=self.height, color=self.colors[idx], label=k)
            else:
                self.ax.fill_between([curr, curr + v], y1=0, y2=self.height, label=k)
            idx += 1
            curr += v

    def make_barplot(self):
        u"""
        make barplot
        :return:
        """
        if self.colors:
            self.ax.fill_between([0, self.values], y1=0, y2=self.height, color=self.colors[0])
        else:
            self.ax.fill_between([0, self.values], y1=0, y2=self.height)

    def make_boxplot(self):
        u"""
        """
        sns.boxplot(x=self.values, ax=self.ax, color="yellow")

    def plot(self):
        u"""
        Make plot base on input values
        """
        if self.boxplot:
            self.make_boxplot()
        elif isinstance(self.values, dict):
            self.make_percentage_barplot()
        elif isinstance(self.values, int) or isinstance(self.values, float):
            self.make_barplot()

        self.set_border()


def format_data(data, first="res.0.6", second="cell_name", third="Tissue", boxplot=False):
    u"""
    pandas is kind of tricky
    dict to collect all the values I needs
    first layer is cell_names
    second layer is {cluster id: counts}

    :param data: pandas data
    :param first: str
    :param second: str or None
    :param third: str or None
    :param boxplot: whether to collect data for boxplot
    :return:
    """
    res = {}

    first = data[first].tolist()
    second = data[second].tolist()

    if third is not None:
        third = data[third].tolist()

        for i, j, k in tqdm(zip(first, second, third)):
            temp = res.get(j, {})

            temp2 = temp.get(i, {})

            temp2[k] = temp2.get(k, 0) + 1

            temp[i] = temp2

            res[j] = temp
    else:
        for i, j in tqdm(zip(first, second)):
            temp = res.get(j, {})
            temp[i] = temp.get(i, 0) + 1
            res[j] = temp

    # if third is not None, convert int to percentage
    if third is not None:
        if not boxplot:
            for key, values in res.items():
                temp = {}
                for cluster, val in values.items():
                    temp1 = {}
                    for k, v in val.items():
                        temp1[k] = v / sum(val.values())
                    temp[cluster] = temp1

                res[key] = temp
        else:
            for key, values in res.items():
                for cluster, val in values.items():
                    values[cluster] = [np.log10(x + 1) for x in val.keys()]

    return res


def make_plots(gs, data, column=0, boxplot=False, ylabel_dist=None):
    u"""
    make plots
    :param gs: grid spec
    :param data: list of values
    :param column: which column to put that figures
    :param boxplot: whether plot boxplot
    :param percentage: whether to make percentage plot
    :param ylabel_dist: int -> then show ylabel
    :return: None
    """
    xmax = 0
    axes = []
    curr = 0
    for key, values in data.items():

        for k, v in values.items():

            if isinstance(v, dict):
                xmax = max(xmax, max(v.values()))
            elif isinstance(v, list):
                xmax = max(xmax, max(v))
            else:
                xmax = max(xmax, v)

            ax = plt.subplot(gs[curr, column])

            ax = BasicBar(
                ax=ax,
                values=v,
                label=k,
                no_x_axis=True,
                boxplot=boxplot
            )

            axes.append(ax)

            curr += 1

        curr += 1

    axes[-1].show_x_axis()

    for i in axes:
        i.set_xlim(xmax)

        if ylabel_dist is not None:
            i.set_ylabel(xmax * ylabel_dist)

    return axes


def test():
    import random

    plt.figure(1, figsize=(9, 3))

    gs = gridspec.GridSpec(3, 3)

    # percentage
    ax = plt.subplot(gs[0, 0])

    BasicBar(ax=ax, values=[0.5, 0.3, 0.2], label="1")

    ax = plt.subplot(gs[2, 0])
    BasicBar(ax=ax, values=[0.4, 0.3, 0.3], label="2", no_x_axis=False)

    # barplot
    ax = plt.subplot(gs[0, 1])
    BasicBar(ax=ax, values=5, xlim=5)

    ax = plt.subplot(gs[2, 1])
    BasicBar(ax=ax, values=4, no_x_axis=False, xlim=5)

    # boxplot
    ax = plt.subplot(gs[0, 2])
    BasicBar(ax=ax, values=[random.randint(0, 100) for x in range(200)], boxplot=True, xlim=150)

    ax = plt.subplot(gs[2, 2])
    BasicBar(ax=ax, values=[random.randint(50, 150) for x in range(200)], boxplot=True, xlim=150, no_x_axis=False)

    plt.subplots_adjust(hspace=.1, wspace=.7)

    plt.show()


if __name__ == '__main__':
    '''
    Index(['nGene', 'nUMI', 'SampleID', 'Batch', 'PatientID', 'Age', 'Sex',
       'Stage', 'Tissue', 'Disease', 'Phage', 'Additional', 'Cells',
       'orig.ident', 'batch', 'percent.mito', 'res.0.6', 'res.0.8', 'tSNE_1',
       'tSNE_2', 'UMAP1', 'UMAP2'],
      dtype='object')
    '''

    data = pd.read_csv(
        "/Users/zhangyiming/Downloads/lung_cancer_10X/讨论/2019.04.12/20190412_seurat_meta.csv",
        index_col=0,
        dtype={
            "nGene": np.int32,
            "nUMI": np.int32,
            "SampleID": str,
            "Batch": str,
            "PatientID": str,
            "Age": np.float,
            "Sex": str,
            "orig.ident": str,
            "batch": str,
            "percent.mito": np.float,
            "res.0.6": np.int,
            "res.0.8": np.int,
            "tSNE_1": np.float,
            "tSNE_2": np.float,
            "UMAP1": np.float,
            "UMAP2": np.float
        }
    )

    len(data["res.0.6"].unique()) + len(data["cell_name"].unique())

    temp = format_data(data)

    #%% prepare data

    plt.figure(1, figsize=(9, 3))

    gs = gridspec.GridSpec(len(data["res.0.6"].unique()) + len(data["cell_name"].unique()), 3)

    #%% generate grid spec

    make_plots(gs, data=temp)

    temp = format_data(data, third="nUMI", boxplot=True)

