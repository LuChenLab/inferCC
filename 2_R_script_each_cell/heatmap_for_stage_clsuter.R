# scripts to make plots of different percentage of different stage and clusters
# using heatmap

args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(pheatmap)
library(reshape2)

rds = args[1]
name = args[2]

setwd(dirname(rds))

obj <- readRDS(basename(rds))

temp <- obj@meta.data %>% 
    dplyr::select(res.0.6, Stage) %>% 
    dplyr::group_by(res.0.6, Stage) %>% 
    dplyr::add_tally() %>%
    unique() %>%
    dplyr::group_by(res.0.6) %>%
    dplyr::mutate(freq = n / sum(n) * 100)

stages = sort(unique(obj@meta.data$Stage))
cluster = sort(unique(obj@meta.data$res.0.6))

temp <- dcast(temp, res.0.6~Stage, value.var = "freq")
rownames(temp) <- temp[,1]
temp <- temp[, 2:ncol(temp)]


png("heatmap_stage_cluster.png", width = length(stages), height = length(cluster), res = 600, units = "in")

tryCatch({
    pheatmap(
        temp, 
        cluster_rows = F, 
        display_numbers = T, 
        main=gsub(pattern = "_", replacement = " ", x = name, perl = T)
    )
}, error = function(e) {
    pheatmap(
        temp, 
        cluster_rows = F, 
        cluster_cols = F,
        display_numbers = T, 
        main=gsub(pattern = "_", replacement = " ", x = name, perl = T)
    )
}
)

dev.off()
