#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
u"""
@since 2019.06.24
"""

library("bigSCale")
library("Seurat")
library("openxlsx")
library(pheatmap)

setwd("/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Alveolar_II_ADC/")

obj <- readRDS("Alveolar_II.rds")

# de = read.xlsx("annotation_results_by_stage.xlsx")
# de = de[de$p_val_adj < 0.05, ]


# network = list()
# label = c()
# for(i in unique(de$ident)) {
#     tryCatch(
#         {
#             network[[length(network) + 1]] = compute.network(expr.data = obj@raw.data, gene.names = de[de$ident == i, "gene"], clustering = 'direct')
#             label = c(label, i)
#         }, error=function(err) {

#         }
#     )
# }

# network.centrality = list()
# for(i in network) {
#     network.centrality[[length(network.centrality) + 1]] = i$centrality
# }


res = compute.network(expr.data = obj@raw.data, gene.names = rownames(obj@raw.data), clustering = 'direct')

png(filename = "bigscale_corr_heatmap.png", type="cairo", width = 36, height = 36, units = "in", res = 300)
pheatmap(as.matrix(res$correlations@Data), show_rownames = F)
dev.off()


png(filename = "newtwork.png", type="cairo", width = 8, height = 8, units = "in", res = 300)
plot(res$graph)
dev.off()


for(i in 1:length(res$graph)) {
    output = names(res$graph[[i]])

    print(output)

    png(filename = paste0(output, ".png"), type="cairo", width = 8, height = 8, units = "in", res = 300)
    rglplot(res$graph[[i]][[output]])
    dev.off() 
}



comparison = compare.centrality(network.centrality, label)

# save results
wb = createWorkbook()
addWorksheet(wb, "Degree")
writeData(wb, 1, comparison$Degree, rowNames = T)

addWorksheet(wb, "Betweenness")
writeData(wb, 2, comparison$Betweenness, rowNames = T)

addWorksheet(wb, "Closeness")
writeData(wb, 3, comparison$Closeness, rowNames = T)

addWorksheet(wb, "PAGErank")
writeData(wb, 4, comparison$PAGErank, rowNames = T)

saveWorkbook(wb, "bigscale.xlsx", overwrite=T)

saveRDS(comparison, "bigscale_DE.rds")

# function to get top n genes, which intersect of 4 different methods
# :param x comparison data
# :param top: top n
# :param methods: get results from which methods
find_top_n <- function(x, top=1000, methods=NULL) {

    if (is.null(methods)) {
        methods = names(x)
    }

    genes = c()
    for (i in methods) {
        temp = x[[i]]
        temp = temp[order(temp$Ranking), ]
        genes = c(genes, rownames(temp[1:top, ]))
    }

    genes = table(genes)
    genes = names(genes[genes > 1])

    return(genes)
}

genes = find_top_n(comparison)

p <- DoHeatmap(obj, genes.use = genes, group.by = "Stage", slim.col.label=TRUE, do.plot = TRUE)
ggsave(
    filename = "bigscale_heatmap_DE.png",
    plot = p,
    width = 6,
    height = 8,
    units = "in",
    dpi = 600,
    type = "cairo-png"
)


p <- DoHeatmap(obj, genes.use = find_top_n(comparison, top=500), group.by = "Stage", slim.col.label=TRUE, do.plot = TRUE)
ggsave(
    filename = "bigscale_heatmap_DE_top500.png",
    plot = p,
    width = 6,
    height = 8,
    units = "in",
    dpi = 600,
    type = "cairo-png"
)

mfuzz <- read.xlsx("mfuzz_gene_module/results.xlsx")


library(VennDiagram)
upset_list = list()
upset_list[["bigscale"]] = genes
upset_list[["bigscale_top500"]] = find_top_n(comparison, top=500)

# upset_list[["mfuzz"]] = unique(mfuzz$gene)
for(i in unique(mfuzz$gene_module_id)) {
    upset_list[[paste0("mfuzz_", i)]] = mfuzz[mfuzz$gene_module_id == i, "gene"]
}

png(filename = "bigscale_upset.png", type="cairo", width = 12, height = 6, units = "in", res = 600)
upset(
    fromList(upset_list), 
    intersections = list(
        c("bigscale", "mfuzz_1"),
        c("bigscale", "mfuzz_2"), 
        c("bigscale", "mfuzz_3"),
        c("bigscale_top500", "mfuzz_1"),
        c("bigscale_top500", "mfuzz_2"),
        c("bigscale_top500", "mfuzz_3"),
        c("bigscale", "mfuzz_1", "mfuzz_2", "mfuzz_3"),
        c("bigscale_top500", "mfuzz_1", "mfuzz_2", "mfuzz_3") 
    )
)
dev.off()
