#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
@since 2019.06.24

A randomForest classifier on stage module
"""
import os
import math

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import seaborn as sns
from anndata import AnnData
from matplotlib import gridspec as gridspec
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.neighbors.classification import KNeighborsClassifier
from sklearn.model_selection import cross_validate



STAGE = {
    "I": 1,
    "II": 2,
    "III": 3,
    "IV": 4
}


def load_data(path:str):
    u"""
    load data from specific path
    
    :param path: path to directory
    """
    
    # raw = pd.read_csv(os.path.join(path, "raw.csv.gz"), index_col=0, engine="python")
    scale = pd.read_csv(os.path.join(path, "scale.csv.gz"), index_col=0, engine="python")
    meta = pd.read_csv(os.path.join(path, "meta.csv.gz"), index_col=0, engine="python")
    
    data = AnnData(X=scale.transpose(), obs = meta)
    data.obs = meta
    # data.raw = raw
    
    return data
    
    
def perform_pca(data: AnnData, genes_use, n_comp=10):
    u"""
    Perform PCA
    :param data: AnnData
    :param genes_use: list of genes
    """ 
    pca = PCA(n_components=n_comp)
    
    mat = data.to_df()[genes_use].transpose()

    pca_res = pca.fit(mat)
    
    res = []
    
    for idx, comp in enumerate(pca_res.components_):
        comp = pd.DataFrame(comp)
        comp.index = mat.columns
        comp.columns = ["PC{}".format(idx + 1)]
        
        res.append(comp)
        
    return pd.concat(res, sort=False, axis=1)
    

def make_plots(data: AnnData, pca: pd.DataFrame, n_comp=10, ncol=3, dpi=300, output=None):
    u"""
    make combined scatter plots
    :param data: AnnData
    :param pca: DataFrame, columns is different PCs
    :param ncol: how many subplots in single row
    :param output: if None, plt.show
    """
    n_comp = min(n_comp, pca.shape[1])
    nrow = math.ceil((n_comp - 1) / ncol)
    
    fig = plt.figure(figsize=(nrow * 6, ncol * 6))
    gs = gridspec.GridSpec(nrow, ncol, figure=fig)
                           
    for i in range(n_comp - 1):
        curr_row = i // ncol
        curr_col = i % ncol
        
        curr_ax = plt.subplot(gs[curr_row, curr_col])
        
        temp = pca[["PC{}".format(i + 1), "PC{}".format(i + 2)]]
        temp = pd.concat([temp, data.obs["Stage"]], axis=1)

        sns.scatterplot(
            x="PC{}".format(i + 1), 
            y="PC{}".format(i + 2), 
            hue="Stage", 
            ax=curr_ax, 
            data=temp,
        )
        
    plt.tight_layout()
    if output:
        plt.savefig(output, dpi=dpi)
    else:
        plt.show()


def  train_random_forest(data: AnnData, genes_use=None, n_estimitors=100, max_depth=2, random_state=0, n_jobs=10):
    u"""
    train random forest by mfuzz
    
    :param data
    """
    
    if genes_use is None:
        genes_use = data.to_df().column

    clf = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0, n_jobs=n_jobs)
    
    mat = data.to_df()[genes_use]
    stage = data.obs.loc[mat.index, ["Stage"]]

    stage = [STAGE[x] for x in stage["Stage"]]
    
    clf.fit(mat, stage)
    
    features = pd.DataFrame(
        data = clf.feature_importances_,
        index = mat.columns,
        columns=['importance']
    ).sort_values('importance', ascending=False)
    
    return clf, features


def find_best_features(
    data: AnnData, 
    test_set: AnnData, 
    disease: str, 
    init_features=None, 
    n_iter=10, 
    random_state=0, 
    n_estimitors=100, 
    max_depth=2, 
    n_jobs=10
):
    u"""
    find best features
    """

    clf, features = train_random_forest(
        data, 
        init_features, 
        random_state=random_state, 
        n_jobs=n_jobs,
        n_estimitors=100, 
        max_depth=2
    )

    res = []
    for i in range(n_iter):
        test_Y = test_set.obs.loc[test_set.obs["Disease"] == disease, ["Stage"]]
        test_X = test_set.to_df().loc[test_Y.index, features.index]

        temp = clf.predict(test_X)

        fpr, tpr, thresholds = metrics.roc_curve([STAGE[x] for x in test_Y["Stage"]], temp.tolist(), pos_label=2)

        s_ = [
            features,
            metrics.auc(fpr, tpr),
            metrics.accuracy_score([STAGE[x] for x in test_Y["Stage"]], temp.tolist())
        ]
        
        res.append(s_)
        clf, features = train_random_forest(
            data, 
            features.loc[features["importance"] > 0, :].index, 
            random_state=random_state, 
            n_jobs=n_jobs,
            n_estimitors=100, 
            max_depth=2
        )
    
    for i in res:
        print([i[0].shape, i[1], i[2]])
        
    return res



# data = load_data("paga")


mfuzz = pd.read_excel(path)

# pca_res = perform_pca(data, mfuzz["gene"].unique())

# make_plots(data, pca=pca_res, output="/mnt/raid62/Lung_cancer_10x/project_scripts/pca.png")

# test_set = load_data("/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Alveolar_II/paga")


best_features = find_best_features(
    train, test, disease="ADC",
    init_features=mfuzz["gene"].unique()
)
        

pca_best_features = perform_pca(data, best_features[3][0].index)
# make_plots(data, pca=pca_best_features, output="/mnt/raid62/Lung_cancer_10x/project_scripts/pca_best_features.png")


def make_violin_of_test_score(data, figsize=None, output=None):
    u"""
    make violin plot of score
    """
    
    res = []
    for i in data:
        for j in i["test_score"]:
            res.append([i["iter"], j])
    
    res = pd.DataFrame(res, columns=["ident", "score"])
    
    if figsize:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig, ax = plt.subplots()
    sns.violinplot(x="ident", y="score", ax=ax, data=res)
    
    plt.tight_layout()
    if output:
        plt.savefig(output, dpi = 600)
    else:
        plt.show()
    


def find_best_features_by_cv(
    data: AnnData, 
    init_features=None, 
    n_iter=10, 
    random_state=0, 
    n_estimitors=100, 
    max_depth=2, 
    n_jobs=10,
    cv=10
):
    u"""
    Find best features by crossvalidation
    :param data: 
    """

    rf = RandomForestClassifier(n_estimators=n_estimitors, max_depth=max_depth, random_state=random_state, n_jobs=n_jobs)
    
    if init_features is None:
        init_features = data.to_df().columns
        
    res = []
    for i in range(n_iter):
        mat = data.to_df().loc[:, init_features]
        cv_results = cross_validate(rf, X=mat, y=data.obs["Stage"], cv=cv, n_jobs=n_jobs, return_estimator=True)
    
        res.append(
            {
                "iter": i,
                "test_score": cv_results["test_score"],
                "features": init_features
            }
        )
        
        importance = [] 
        for idx, estimator in enumerate(cv_results["estimator"]):
            importance.append(
                pd.DataFrame(
                    estimator.feature_importances_,
                    index=init_features,
                )
            )


        importance = pd.concat(importance, axis=1)
        
        importance = importance > 0
        importance = importance.sum(axis=1)
        
        init_features = importance[importance > (cv / 2)].index
        
    return res


def find_best_features_by_cv_knn(
    data: AnnData, 
    init_features=None, 
    max_neighbors=30,
    n_iter=10, 
    random_state=0, 
    n_jobs=10,
    cv=10
):
    u"""
    Find best features by crossvalidation
    :param data: 
    """

    
    
    if init_features is None:
        init_features = data.to_df().columns
        
    res = []
    mat = data.to_df().loc[:, init_features]
    
    for i in range(3, max_neighbors + 1):
        knn = KNeighborsClassifier(n_neighbors=i)
        cv_results = cross_validate(knn, X=mat, y=data.obs["Stage"], cv=cv, n_jobs=n_jobs, return_estimator=False)

        res.append(
            {
                "iter": i,
                "test_score": cv_results["test_score"],
                # "features": init_features
            }
        )
                
    return res


features = [
    'GPC4', 'CYBA', 'ELF3', 'AZGP1', 'FKBP5', 'MUC1', 'TSC22D3', 'C4BPA',
       'MT2A', 'XIST', 'SCGB3A1', 'FAM84B', 'MBNL3', 'SPINK1', 'HLA-DQA2',
       'MT-RNR1', 'CHCHD2', 'HLA-DRB5', 'TAGLN2', 'POLR2J3-1', 'ERBB3',
       'CEACAM6', 'TFPI', 'SFTPC', 'QPRT', 'RNASE1', 'MT1E', 'SEC61G', 'EPS8',
       'CEACAM5', 'RPS4Y1', 'FABP5', 'MT-ND5', 'RPS19', 'MUC21', 'MT-ND4L',
       'AREG', 'SOX4', 'DRAM1', 'CTGF', 'DDIT4', 'EGR1', 'NET1', 'B3GNT7',
       'SCGB3A2', 'CD74', 'TRIB1', 'CP', 'RPS26', 'HMGN3', 'APLP2', 'PRR15L',
       'MYO6', 'MYLIP', 'PRSS8', 'C3', 'EFNA1', 'FOSB', 'BCAP31', 'JUNB',
       'S100A6', 'WNK1', 'PMAIP1'
]


cv_results_all = find_best_features_by_cv(data, init_features=None, n_iter=50, random_state=0, n_estimitors=100, max_depth=None, n_jobs=10, cv=3)
for i in cv_results_all:
    print(sum(i["test_score"]) / len(i["test_score"]), max(i["test_score"]))
    
make_violin_of_test_score(cv_results_all, figsize=(10, 6) ,output="/mnt/raid62/Lung_cancer_10x/project_scripts/RF_all.png")

cv_results_mfuzz = find_best_features_by_cv(data, init_features=mfuzz["gene"].unique(), n_iter=1, random_state=0, n_estimitors=100, max_depth=None, n_jobs=10, cv=3)
for i in cv_results_mfuzz:
    print(sum(i["test_score"]) / len(i["test_score"]), max(i["test_score"]))

cv_results = find_best_features_by_cv(data, init_features=features, n_iter=1, random_state=0, n_estimitors=100, max_depth=None, n_jobs=10, cv=3)
for i in cv_results:
    print(sum(i["test_score"]) / len(i["test_score"]), max(i["test_score"]))


knn_res_all = find_best_features_by_cv_knn(data, max_neighbors=40)
for i in knn_res_all:
    print(sum(i["test_score"]) / len(i["test_score"]), max(i["test_score"]))

make_violin_of_test_score(knn_res_all, figsize=(10, 6) ,output="/mnt/raid62/Lung_cancer_10x/project_scripts/KNN_all.png")

knn_res_mfuzz = find_best_features_by_cv_knn(data, init_features=mfuzz["gene"].unique(), max_neighbors=40)
for i in knn_res_mfuzz:
    print(sum(i["test_score"]) / len(i["test_score"]), max(i["test_score"]))

make_violin_of_test_score(knn_res_mfuzz, figsize=(10, 6) ,output="/mnt/raid62/Lung_cancer_10x/project_scripts/KNN_mfuzz.png")

knn_res = find_best_features_by_cv_knn(data, init_features=features, max_neighbors=40)
for i in knn_res:
    print(sum(i["test_score"]) / len(i["test_score"]), max(i["test_score"]))
    
make_violin_of_test_score(knn_res, figsize=(10, 6) ,output="/mnt/raid62/Lung_cancer_10x/project_scripts/KNN.png")
