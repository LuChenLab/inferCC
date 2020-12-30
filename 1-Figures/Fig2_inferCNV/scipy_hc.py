from scipy.cluster.hierarchy import linkage, dendrogram

def hclust_scipy(X, method="complete"):
    return [int(i) for i in dendrogram(linkage(X, method=method))['ivl']]