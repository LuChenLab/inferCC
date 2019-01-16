###
### Created by ygidtu@gmail.com at 2019.01.15
### This scripts is used to run Seurat on Lung Cancer 10X Genome sequencing data
### 
load_pacages <- function() {
    if(!require("argparse")) {
        install.packages("argparse")
        
        library("argparse")
    }
    
    if(!require("Seurat")) {
        install.packages("Seurat")
        
        library("Seurat")
    }
}

load_pacages()
# suppressPackageStartupMessages(load_pacages)
# suppressPackageStartupMessages(library(Seurat))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument(
    "-i", 
    "--input", 
    type="character",
    help="Path to input directory [default %(default)s]"
)
parser$add_argument(
    "-o", 
    "--output", 
    help="Path to output directoy"
)
parser$add_argument(
    "-l",
    "--label",
    help=""
)


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

if (is.na(args$input)) {
    stop("no input")
} else if (is.na(args$input)) {
    stop("no output")
}

input = args$input
output = args$output
label = args$label


# input = "/Volumes/WD/Seurat/filtered_gene_bc_matrices/tests"
# output = "/Volumes/WD/Seurat/filtered_gene_bc_matrices/output"
# label = "test"

setwd(output)

##### Load Data

# Load the dataset
# The cell ranger outs filtered features
# must decompressed and the features.tsv renamed to genes.tsv

print("Reading")
data <- Read10X(data.dir = input)

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = data))

sparse.size <- object.size(x = data)


# Initialize the Seurat object with the raw (non-normalized data).  
# Keep all genes expressed in >= 3 cells (~0.1% of the data). 
# Keep all cells with at least 200 detected genes
# 保证每个基因至少在三个细胞系中有表达
# 保证每个细胞至少能检测到两百个基因
print("Create SeuratObject")
pbmc <- CreateSeuratObject(
    raw.data = data, 
    min.cells = 3, 
    min.genes = 200, 
    project = label
)


##### QC

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat. 
# For non-UMI data, nUMI represents the sum of the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and non-log-normalized counts The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.

# 计算每个细胞中的基因数和UMI数量，并且计算其中线粒体基因的比例
print("VlnPlot")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
# nGenes: 每个细胞中的基因数
# nUMI: 每个细胞中的UMI数
# percent.mito: 每个细胞中线粒体基因的比例
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

p <- VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
ggsave(filename="VlnPlot.png", plot=p)


# GenePlot is typically used to visualize gene-gene relationships,
# but can be used for anything calculated by the object, i.e. columns in object@meta.data, PC scores etc.  
# Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well

# 计算UMI数量分别与线粒体基因数量和基因数量的相关性
print("GenePlot")
png(filename="GenePlot.png", res=300, width=1800, height = 1200)
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
dev.off()


# We filter out cells that have unique gene counts over 2,500 or less than 200 
# Note that low.thresholds and high.thresholds are used to define a 'gate'.
# -Inf and Inf should be used if you don't want a lower or upper threshold.
# 按照每个细胞中超过2500个基因或者少于200个基因的标准进行筛选
print("Fitler")
pbmc <- FilterCells(
    object = pbmc, 
    subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf),
    high.thresholds = c(2500, 0.05)
)


#### Normalization
print("Normalize")
pbmc <- NormalizeData(
    object = pbmc,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)


#### Detect variable genes
print("VariableGenes")
png("Variable.png", res = 300, width = 1800, height = 1800)
pbmc <- FindVariableGenes(
    object = pbmc, 
    mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, 
    x.high.cutoff = 3, 
    y.cutoff = 0.5
)
dev.off()

#### Scale data and removing unwanted sources of variation
print("Scale")
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))


#### PCA
print("PCA")
pbmc <- RunPCA(
    object = pbmc, 
    pc.genes = pbmc@var.genes, 
    do.print = FALSE, 
    pcs.print = 1:5, 
    genes.print = 5
)
# PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

print("VizPCA")
png(filename="VizPCA.png", res=300, width=1800, height = 1800)
p <- VizPCA(object = pbmc, pcs.use = 1:2)
# ggsave(filename="VizPCA.png", plot = p)
p
dev.off()


print("PCAPlot")
p <- PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
ggsave(filename="PCAPlot.png", plot=p)

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)

print("PCheatmap")
png(filename="PCHeatmap.png", res=300, width=2400, height = 2400)

if (length(pbmc@meta.data$nGene) > 500) {
    PCHeatmap(
        object = pbmc,
        pc.use = 1:12, 
        do.balanced = TRUE, 
        label.columns = FALSE,
        use.full = TRUE,
        cells.use = 500
    )
} else {
    PCHeatmap(
        object = pbmc,
        pc.use = 1:12, 
        do.balanced = TRUE, 
        label.columns = FALSE,
        use.full = TRUE
    )
}

dev.off()


#### 

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)

print("JackStrawPlot")
png(filename="JackStrawPlot.png", res=300, width = 1800, height = 1800)
JackStrawPlot(object = pbmc, PCs = 1:12)
dev.off()

print("PCElbowPlot")
png(filename="PCElbowPlot.png", res=300, width = 1800, height = 1800)
PCElbowPlot(object = pbmc)
dev.off()

saveRDS(pbmc, file = "1_before_cluster.rds")

#### cluster

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
print("Cluster")
pbmc <- FindClusters(
    object = pbmc, 
    reduction.type = "pca", 
    dims.use = 1:10, 
    resolution = 0.6, 
    print.output = 0, 
    save.SNN = TRUE
    )

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10)

print("tSNE")
ggsave(
    filename="tSNE.png", 
    plot = TSNEPlot(object = pbmc)
    )

saveRDS(pbmc, file = "2_clustered.rds")


