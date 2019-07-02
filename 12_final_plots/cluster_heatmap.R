#!/usr/bin/env Rscript
# @since 2019.06.30
# This scripts is used to make complex heatmap of cluster module


# prepare variables
library(argparse)


parser <- ArgumentParser(
    description='Process some integers'
)
parser$add_argument(
    '-i', 
    type="character", 
    nargs='+',
    help='path to rds'
)
parser$add_argument(
    '-o', 
    type="character",
    help='path to output'
)

parser$add_argument(
    "-d", 
    action="store_true", 
    default=TRUE,
    help="Print extra output [default]"
)

args <- parser$parse_args()

rds = args$i
output_dir = args$o

is_disease = args$d


# rds = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Alveolar_II_ADC/Alveolar_II.rds"

input_dir = dirname(rds)
cluster_module = paste(input_dir, "cluster_gene_module/results.xlsx", sep = "/")

if (!file.exists(cluster_module)) {
    stop("no")
}


# load packages
library(openxlsx)
library(dplyr)
library(ComplexHeatmap)


custom_colors <- c(
    "#1B62A3", "#9CB8DD", "#EE6719", "#F8A966", "#278F36", 
    "#8CC873", "#C70F1F", "#EF7E82", "#7D509A", "#B59BC7",
    "#76423A", "#B48880", "#CB60A1", "#F1A3C5", "#7E7F83", "#BDC3C4", 
    "#ABB126", "#D1D57A", "#1DAFC3", "#8BCFDA"
)


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


# extract temp information
stage_colors = c("Benign"="#3E6596", "I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E")
disease_colors = c("ADC" = "#063373", "LCC"="#FECC1B", "Normal"="#73BDFF", "Other"="#DA6906", "SCC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B", "LBT"="#0084D1")
tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")
patient_colors = gg_color_hue(50)
names(patient_colors) <- sapply(1:50, function(x){sprintf("P%02d", x)})


# prepare data
obj <- readRDS(rds)

cluster_module = read.xlsx(cluster_module)

cluster_module = as.data.frame(cluster_module %>% group_by(ident) %>% top_n(10, wt=avg_logFC))


# make_heatmap

## first get the expression matrix

meta = obj@meta.data

meta = meta[order(as.numeric(meta$res.0.6)), ]


if (is_disease) {
    ha = HeatmapAnnotation(
        cluster = as.factor(meta$res.0.6), 
        stage = factor(meta$Stage, levels = sort(unique(meta$Stage))),
        patient = factor(meta$PatientID, levels = names(patient_colors)),
        col = list(
            stage = sapply(meta$Stage, function(x) {stage_colors[[x]]}),
            patient = sapply(meta$PatientID, function(x) {patient_colors[[x]]})
        )
    )
} else {
    ha = HeatmapAnnotation(
        cluster = as.factor(meta$res.0.6), 
        stage = factor(meta$Stage, levels = sort(unique(meta$Stage))),
        patient = factor(meta$PatientID, levels = names(patient_colors)),
        disease = meta$Disease,
        col = list(
            stage = sapply(meta$Stage, function(x) {stage_colors[[x]]}),
            patient = sapply(meta$PatientID, function(x) {patient_colors[[x]]}),
            disease = sapply(meta$Disease, function(x) {disease_colors[[x]]})
        )
    )
}


png(filename = output_dir, width = ncol(obj@raw.data) / 100 + 10, height = nrow(cluster_module) / 10 + 4, res = 600, units = "in")
Heatmap(
    obj@scale.data[cluster_module$gene[order(cluster_module$ident)], rownames(meta)],
    top_annotation = ha,
    cluster_rows = F,
    cluster_columns = F,
    show_column_names = F,
    name = ""
)

dev.off()
