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
    
    if(!require("logging")) {
        install.packages("logging")
        
        library("logging")
    }
}

load_pacages()

basicConfig()
addHandler(writeToConsole)

### set command line arguments

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
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
    "-n",
    "--normal",
    type="character",
    help="Path to normal rds"
)
parser$add_argument(
    "-c",
    "--cancer",
    type="character",
    help="Path to cancer rds"
)


args <- parser$parse_args()
output=args$output
label=args$label
normal = args$normal
cancer = args$cancer


#### testing

# setwd("/Volumes/WD/Seurat/results/tests")
# 
# normal <- Read10X(data.dir = "/Volumes/WD/Seurat/results/filtered_feature_bc_matrix/2018jz39")
# cancer <- Read10X(data.dir = "/Volumes/WD/Seurat/results/filtered_feature_bc_matrix/2018jz38")

if (is.na(output)) {
    stop("-o/--output required")
} else if (is.na(label)) {
    stop("-l/--label required")
} else if (is.na(normal)) {
    stop("-n/--normal required")
} else if (is.na(cancer)) {
    stop("-c/--cancer required")
}

setwd(output)

normal <- Read10X(data.dir = normal)
cancer <- Read10X(data.dir = cancer)


####

normal <- CreateSeuratObject(raw.data = normal, min.cells = 5, project = "10X")
cancer <- CreateSeuratObject(raw.data = cancer, min.cells = 5, project = "10X")


normal@meta.data$stim <- "Normal"
normal <- FilterCells(normal, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
normal <- NormalizeData(normal)
normal <- ScaleData(normal, display.progress = F)
# Set up stimulated object
cancer@meta.data$stim <- "Cancer"
cancer <- FilterCells(cancer, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
cancer <- NormalizeData(cancer)
cancer <- ScaleData(cancer, display.progress = F) 


normal <- FindVariableGenes(normal, do.plot = F)
cancer <- FindVariableGenes(cancer, do.plot = F)

g.1 <- head(rownames(normal@hvg.info), 1000)
g.2 <- head(rownames(cancer@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(normal@scale.data))
genes.use <- intersect(genes.use, rownames(cancer@scale.data))


#### Run CCA
loginfo("Run CCA")
immune.combined <- RunCCA(
    normal, 
    cancer, 
    genes.use = unique(genes.use), 
    num.cc = 30,
    add.cell.id1 = "Normal",
    add.cell.id2 = "Cancer"
    )

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(
    object = immune.combined, 
    reduction.use = "cca", 
    group.by = "stim", 
    pt.size = 0.5, 
    do.return = TRUE
)

p2 <- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "stim", 
              do.return = TRUE)
p <- plot_grid(p1, p2)

ggsave(filename="CCA.png", plot=p)


#### MetageneBicorPlot
loginfo("MetageneBicorPlot")
p <- MetageneBicorPlot(
    immune.combined, 
    grouping.var = "stim", 
    dims.eval = 1:30, 
    display.progress = FALSE
)
ggsave(filename="MetageneBiorPlot.png", plot = p)


#### DimHeatmap
png("DimHeatmap.png", res = 300, width = 2700, height = 2700)
DimHeatmap(
    object = immune.combined, 
    reduction.type = "cca", 
    cells.use = 500, 
    dim.use = 1:9, 
    do.balanced = TRUE
)
dev.off()


#### align CCA subspaces
immune.combined <- AlignSubspace(
    immune.combined, 
    reduction.type = "cca",
    grouping.var = "stim", 
    dims.align = 1:20
)

loginfo("Aligned CCA VlnPlot")
p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "stim", 
              do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "stim", 
              do.return = TRUE)
p <- plot_grid(p1, p2)

ggsave(filename = "AlignedCCAVlnPlot.png", plot = p)

#### tSNE and clustering
loginfo("tSNE and clustering")
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


loginfo("tSNE plot")
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
p <- plot_grid(p1, p2)

ggsave(filename = "tSNE.png", plot = p)

#### conserved marker

nk.markers <- FindConservedMarkers(
    immune.combined, 
    ident.1 = 7, 
    grouping.var = "stim", 
    loginfo.bar = FALSE
)

####
#### 这里需要不同的基因list来确定具体细胞群的类别
#### 此处先采用血细胞的buf
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


#### Draw tSNE based on clustering

# new.ident <- c(
#     "CD14 Mono",
#     "CD4 Naive T",
#     "CD4 Memory T",
#     "B",
#     "CD16 Mono",
#     "T activated",
#     "CD8 T",
#     "NK",
#     "DC",
#     "B activated",
#     "Mk",
#     "pDC",
#     "Eryth"
# )

## assign cell type to cluster
# for (i in 0:11) {
#     immune.combined <- RenameIdent(
#         object = immune.combined,
#         old.ident.name = i,
#         new.ident.name = i # new.ident[i + 1]
#     )
# }

p <- TSNEPlot(immune.combined, do.label = T, pt.size = 0.5)
ggsave(filename = "tSNE_with_cell_name.png", plot = p)


# immune.combined@ident <- factor(
#     immune.combined@ident,
#     levels = (
#         c(
#             "pDC", "Eryth",
#             "Mk", "DC",
#             "CD14 Mono",
#             "CD16 Mono",
#             "B activated",
#             "B", "CD8 T",
#             "NK", "T activated",
#             "CD4 Naive T",
#             "CD4 Memory T"
#         )
#     )
# )
# 
# markers.to.plot <- c(
#     "CD3D", "CREM",
#     "HSPH1", "SELL",
#     "GIMAP5", "CACYBP",
#     "GNLY", "NKG7",
#     "CCL5", "CD8A",
#     "MS4A1", "CD79A",
#     "MIR155HG", "NME1",
#     "FCGR3A", "VMO1",
#     "CCL2", "S100A9",
#     "HLA-DQA1", "GPR183",
#     "PPBP", "GNG11",
#     "HBA2", "HBB",
#     "TSPAN13", "IL3RA"
# )


sdp <- SplitDotPlotGG(
    immune.combined,
    genes.plot = genes.use[1:10],
    cols.use = c("blue", "red"),
    x.lab.rot = T,
    plot.legend = T,
    dot.scale = 8,
    do.return = T,
    grouping.var = "stim"
)

ggsave(filename="SpliceDotPlotGG.png", plot = sdp) 


saveRDS(immune.combined, file = "immune_combined.rds")

