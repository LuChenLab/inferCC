###
### Created by ygidtu@gmail.com at 2019.01.16
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
    
    if (!require("openxlsx")) {
        install.packages("openxlsx")
        
        library("openxlsx")
    }
}

load_pacages()


### set command line arguments

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
    help="label the output figures"
)
parser$add_argument(
    "-x",
    "--xlsx",
    type="character",
    help="Path to xlsx meta info"
)


input=parser$input
output=parser$output
label=parser$label
xlsx=parser$xlsx


#### testing

xlsx = "/Volumes/WD/Seurat/results/测序质量.xlsx"
label = "test"


if (is.na(input)) {
    stop("-i/--input required")
} else if (is.na(output)) {
    stop("-o/--output required")
} else if (is.na(label)) {
    stop("-l/--label required")
} else if (is.na(xlsx)) {
    stop("-x/--xlsx required")
}


#### load meta info
meta <- read.xlsx(xlsx)
meta <- meta[meta$status == 1, c("Name", "tissue_type")]

####

ctrl.data <- read.table("~/Downloads/immune_control_expression_matrix.txt.gz", 
                        sep = "\t")
stim.data <- read.table("~/Downloads/immune_stimulated_expression_matrix.txt.gz", 
                        sep = "\t")

# Set up control object
ctrl <- CreateSeuratObject(raw.data = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl@meta.data$stim <- "CTRL"
ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl, display.progress = F)
# Set up stimulated object
stim <- CreateSeuratObject(raw.data = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim@meta.data$stim <- "STIM"
stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
stim <- NormalizeData(stim)
stim <- ScaleData(stim, display.progress = F)

# Gene selection for input to CCA
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))


#### Run CCA
print("Run CCA")
immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = immune.combined, reduction.use = "cca", group.by = "stim", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "stim", 
              do.return = TRUE)
p <- plot_grid(p1, p2)

ggsave(filename="CCA.png", plot=p)


#### MetageneBicorPlot
print("MetageneBicorPlot")
p <- MetageneBicorPlot(
    immune.combined, 
    grouping.var = "stim", 
    dims.eval = 1:30, 
    display.progress = FALSE
)
ggsave(filename="MetageneBiorPlot.png", plot = p)


#### DimHeatmap

p <- DimHeatmap(
    object = immune.combined, 
    reduction.type = "cca", 
    cells.use = 500, 
    dim.use = 1:9, 
    do.balanced = TRUE
)

ggsave(filename = "DimHeatmap.png", plot = p)


#### align CCA subspaces
immune.combined <- AlignSubspace(
    immune.combined, 
    reduction.type = "cca",
    grouping.var = "stim", 
    dims.align = 1:20
)

print("Aligned CCA VlnPlot")
p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "stim", 
              do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "stim", 
              do.return = TRUE)
p <- plot_grid(p1, p2)

ggsave(filename = "AlignedCCAVlnPlot.png", plot = p)

#### tSNE and clustering
print("tSNE and clustering")
immune.combined <- RunTSNE(
    immune.combined, 
    reduction.use = "cca.aligned",
    dims.use = 1:20, 
    do.fast = T
)

immune.combined <- FindClusters(
    immune.combined, 
    reduction.type = "cca.aligned", 
    resolution = 0.6, 
    dims.use = 1:20
)


print("tSNE plot")
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
p <- plot_grid(p1, p2)

ggsave(filename = "tSNE.png", plot = p)

#### conserved marker

nk.markers <- FindConservedMarkers(
    immune.combined, 
    ident.1 = 7, 
    grouping.var = "stim", 
    print.bar = FALSE
)

####
#### 这里需要不同的基因list来确定具体细胞群的类别
####

p <- FeaturePlot(
    object = immune.combined, 
    features.plot = c(
        "CD3D", "SELL", "CREM", 
         "CD8A", "GNLY", "CD79A", 
        "FCGR3A", "CCL2", "PPBP"
    ), 
    min.cutoff = "q9", 
    cols.use = c("lightgrey", "blue"), 
    pt.size = 0.5
)

ggsave(filename = "feature.png", plot = p)



