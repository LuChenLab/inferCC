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
    
    if(!require("dplyr")) {
        install.packages("dplyr")
        
        library("dplyr")
    }
}

load_pacages()

basicConfig()
# addHandler(writeToConsole)

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


normal@meta.data$stim <- basename(args$normal)
normal <- FilterCells(normal, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
normal <- NormalizeData(normal)
normal <- ScaleData(normal, display.progress = F)
# Set up stimulated object
cancer@meta.data$stim <- basename(args$cancer)
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


# customize_Seurat_FeaturePlot <- function(p, alpha.use = 1, gradient.use = c("yellow", "red"), expression.threshold = 0, is.log1p.transformed = F) {
#     
#     #### Main function ####
#     main_function <- function(p = p, alpha.use = alpha.use, gradient.use = gradient.use, expression.threshold = expression.threshold, is.log1p.transformed = is.log1p.transformed) {
#         
#         # Order data by gene expresion level
#         p$data <- p$data[order(p$data$gene),]
#         
#         # Define lower limit of gene expression level
#         if (isTRUE(is.log1p.transformed)) {
#             expression.threshold <- expression.threshold
#         } else {
#             expression.threshold <- log1p(expression.threshold)
#         }
#         
#         # Compute maximum value in gene expression
#         max.exp <- max(p$data$gene)
#         
#         # Fill points using the gene expression levels
#         p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
#         
#         # Define transparency of points
#         p$layers[[1]]$mapping$alpha <- alpha.use
#         
#         # Change fill and colour gradient values
#         p <- p + scale_colour_gradientn(colours = gradient.use, guide = F, limits = c(expression.threshold, max.exp), na.value = "grey") +
#             scale_fill_gradientn(colours = gradient.use, name = expression(atop(Expression, (log))), limits = c(expression.threshold, max.exp), na.value = "grey") +
#             scale_alpha_continuous(range = alpha.use, guide = F)
#     }
#     
#     #### Execution of main function ####
#     # Apply main function on all features
#     p <- lapply(X = p, alpha.use = alpha.use, gradient.use = gradient.use, 
#                 expression.threshold = expression.threshold, is.log1p.transformed = is.log1p.transformed,
#                 FUN = main_function)
#     
#     # Arrange all plots using cowplot
#     # Adapted from Seurat
#     # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
#     # ncol argument adapted from Josh O'Brien
#     # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
#     # cowplot::plot_grid(plotlist = p, ncol = ceiling(sqrt(length(p))))
#     marrangeGrob(p, nrow=3, ncol=3)
# }


nk.markers <- FindAllMarkers(
    immune.combined,
    only.pos = TRUE, 
    min.pct = 0.25, 
    thresh.use = 0.25
)

nk.markers <- nk.markers[nk.markers$p_val_adj<0.05,]
nk.markers %>% group_by(cluster) %>% top_n(10) -> top10
nk.markers %>% group_by(cluster) %>% top_n(2) -> top2

saveRDS(immune.combined, file = "1_immune_combined.rds")


####
#### 这里需要不同的基因list来确定具体细胞群的类别
#### 此处先采用血细胞的buf
####

png("feature.png", res = 300, width = 6400, height = 6400)
FeaturePlot(
    object = immune.combined, 
    features.plot = sort(unique(top10$gene))[1:10], 
    min.cutoff = "q9", 
    cols.use = c("lightgrey", "blue"), 
    pt.size = 0.5
)
dev.off()



#### Draw tSNE based on clustering

p <- TSNEPlot(immune.combined, do.label = T, pt.size = 0.5)
ggsave(filename = "tSNE_with_cell_name.png", plot = p)


sdp <- SplitDotPlotGG(
    immune.combined,
    genes.plot = sort(unique(top10$gene)),
    cols.use = c("blue", "red"),
    x.lab.rot = T,
    plot.legend = T,
    dot.scale = 8,
    do.return = T,
    grouping.var = "stim"
)

ggsave(filename="SpliceDotPlotGG.png", plot = sdp, width = 15) 



#### 

LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
    for (i in genes) {
        x1 <- exp.mat[i, 1]
        y1 <- exp.mat[i, 2]
        plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                                label = i, size = text.size)
        plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                                    adj.y.s, yend = y1, size = segment.size)
    }
    return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                      adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                      adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}


t.cells <- SubsetData(immune.combined, ident.use = "Normal", subset.raw = T)
t.cells <- SetAllIdent(t.cells, id = "stim")
avg.t.cells <- log1p(AverageExpression(t.cells, show.progress = FALSE))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- SubsetData(immune.combined, ident.use = "Cancer", subset.raw = T)
cd14.mono <- SetAllIdent(cd14.mono, id = "stim")
avg.cd14.mono <- log1p(AverageExpression(cd14.mono, show.progress = FALSE))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

p1 <- ggplot(avg.t.cells, aes(normal, cancer)) + 
    geom_point() + 
    ggtitle("Cancer")

p1 <- LabelUR(
    p1, 
    genes = c(genes.to.label1, genes.to.label2), 
    avg.t.cells, 
    adj.u.t = 0.3, 
    adj.u.s = 0.23
)

p1 <- LabelUL(
    p1, 
    genes = top10$gene, 
    avg.t.cells, 
    adj.u.t = 0.5, 
    adj.u.s = 0.4, 
    adj.l.t = 0.25, 
    adj.l.s = 0.25
)

p2 <- ggplot(
    avg.cd14.mono, 
    aes(normal, cancer)) + 
    geom_point() + 
    ggtitle("Cancer")

p2 <- LabelUR(
    p2, 
    genes = top10$gene, 
    avg.cd14.mono, 
    adj.u.t = 0.3, 
    adj.u.s = 0.23
)

p2 <- LabelUL(
    p2, 
    genes = genes.to.label2, 
    avg.cd14.mono, 
    adj.u.t = 0.5, 
    adj.u.s = 0.4, 
    adj.l.t = 0.25, 
    adj.l.s = 0.25
)

p <- plot_grid(p1, p2)

ggsave("differential_expressed_genes.png", plot = p)


immune.combined@meta.data$celltype.stim <- paste0(
    immune.combined@ident, 
    "_", 
    immune.combined@meta.data$stim
)

immune.combined <- StashIdent(
    immune.combined, 
    save.name = "celltype"
)
immune.combined <- SetAllIdent(
    immune.combined, 
    id = "celltype.stim"
)
b.interferon.response <- FindMarkers(
    immune.combined, 
    ident.1 = "Normal", 
    ident.2 = "Cancer", 
    print.bar = FALSE
)

png("feature_heatmap.png", width = 2400, height = 2400, res= 300)
FeatureHeatmap(
    immune.combined, 
    features.plot = top10$gene, 
    group.by = "stim", 
    pt.size = 0.25, 
    key.position = "top", 
    max.exp = 3
)
dev.off()


saveRDS(immune.combined, file = "immune_combined.rds")

