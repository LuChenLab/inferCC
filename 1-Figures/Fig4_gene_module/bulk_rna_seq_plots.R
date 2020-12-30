# Create at 2019.09.26
# make volcano plots and compare with scRNA



library(ggplot2)
library(dplyr)
library(ggrepel)
library(openxlsx)
library(stringr)
library(ggpubr)
library(cowplot)
library(doMC)
library(here)

set_here("LungCancer10x/")


################################
## load expression data to make sure of logFC

expr = read.xlsx(here("09_bulk/RNA_seq/RSEM.xlsx"), rowNames = T)


luad_samples = colnames(expr)[str_detect(colnames(expr), "LUAD")]
lusc_samples = colnames(expr)[str_detect(colnames(expr), "LUSC")]
bspn_samples = colnames(expr)[str_detect(colnames(expr), "BSPN")]


################################
## Make valcano plots


wb = createWorkbook()

# 1. LUAD vs LUSC

data = read.xlsx(here("09_bulk/RNA_seq/Disease_DEGs.xlsx"), sheet = 1)

temp_data = data[abs(data$log2FoldChange) > 1 & data$padj < 0.05, ]
temp_data = na.omit(temp_data[order(temp_data$log2FoldChange), ])
gene.label = c(temp_data[1:10, "symbol"], temp_data[(nrow(temp_data) - 9):nrow(temp_data), "symbol"])


## Set label
# print(data[1, "log2FoldChange"])
# print(log2(mean(as.numeric(expr[data[1, "gene_id"], lusc_samples])) / mean(as.numeric(expr[data[1, "gene_id"], luad_samples]))))
# print(mean(as.numeric(expr[data[1, "gene_id"], lusc_samples])) < mean(as.numeric(expr[data[1, "gene_id"], luad_samples])))

# log2FoldChange > 1 is LUSC
# log2FoldChange < -1 is LUAD


data$Type = apply(data, 1, function(row) {
    if (is.na(row["padj"])) {
        row["padj"] = 1
    }
    
    if (as.numeric(row["log2FoldChange"]) > 1 && as.numeric(row["padj"]) < 0.05) {
        return("LUSC")
    } else if (as.numeric(row["log2FoldChange"]) < -1 && as.numeric(row["padj"]) < 0.05) {
        return("LUAD")
    } 
    return("Not")
})


## save LUAD vs LUSC data
addWorksheet(wb, "LUAD_vs_LUSC")
writeData(wb, 1, data)

gene.label = c()
p <- ggplot(data, aes(x=log2FoldChange, y=-log10(padj), color=Type)) + 
    geom_point() +
    scale_color_manual(values=c("LUAD"="#0084D1", "BSPN"="#9CD7D7", "LUSC"="#A0C807", "Not"="grey")) +
    geom_label_repel(
        data = data[data$symbol %in% c(gene.label, "DNAJB1", "GNLY"), ], 
        aes(x=log2FoldChange, y=-log10(padj), label=symbol),
        size = 5,
        color="black"
    ) +
    labs(color = "", title = "RNA-seq LUAD vs. LUSC") +
    theme_bw() + 
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = c(0.85, 0.15),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill=alpha('white', 0))
    ) +
    guides(color = guide_legend(override.aes = list(size = 8)))


p

ggsave(
    filename = here("09_bulk/RNA_seq/volcano_plot_LUAD_vs_LUSC.pdf"),
    plot = p,
    width = 6,
    height = 6
)


###################################
# 2. LUAD vs BSPN

data = read.xlsx(here("09_bulk/RNA_seq/Disease_DEGs.xlsx"), sheet = 2)

temp_data = data[abs(data$log2FoldChange) > 1 & data$padj < 0.05, ]
temp_data = na.omit(temp_data[order(temp_data$log2FoldChange), ])
gene.label = c(temp_data[1:10, "symbol"], temp_data[(nrow(temp_data) - 9):nrow(temp_data), "symbol"])


## Set label
print(data[1, "log2FoldChange"])
print(log2(mean(as.numeric(expr[data[1, "gene_id"], luad_samples])) / mean(as.numeric(expr[data[1, "gene_id"], bspn_samples]))))
print(mean(as.numeric(expr[data[1, "gene_id"], luad_samples])) < mean(as.numeric(expr[data[1, "gene_id"], bspn_samples])))

# log2FoldChange > 1 is LUAD
# log2FoldChange < -1 is BSPN


data$Type = apply(data, 1, function(row) {
    if (is.na(row["padj"])) {
        row["padj"] = 1
    }
    
    if (as.numeric(row["log2FoldChange"]) > 1 && as.numeric(row["padj"]) < 0.05) {
        return("LUAD")
    } else if (as.numeric(row["log2FoldChange"]) < -1 && as.numeric(row["padj"]) < 0.05) {
        return("BSPN")
    } 
    return("Not")
})


## save LUAD vs LUSC data
addWorksheet(wb, "LUAD_vs_BSPN")
writeData(wb, 2, data)


p <- ggplot(data, aes(x=log2FoldChange, y=-log10(padj), color=Type)) + 
    geom_point() +
    scale_color_manual(values=c("LUAD"="#0084D1", "BSPN"="#9CD7D7", "LUSC"="#A0C807", "Not"="grey")) +
    geom_text_repel(
        data = data[data$symbol %in% gene.label, ], 
        aes(x=log2FoldChange, y=-log10(padj), color=Type, label=symbol),
        size = 5
    ) +
    labs(color = "", title = "RNA-seq LUAD vs. BSPN") +
    theme_bw() + 
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = c(0.85, 0.15),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill=alpha('white', 0))
    ) +
    guides(color = guide_legend(override.aes = list(size = 8)))

p

ggsave(
    filename = here("09_bulk/RNA_seq/volcano_plot_LUAD_vs_BSPN.pdf"),
    plot = p,
    width = 6,
    height = 6
)



###################################
# 3. LUSC vs BSPN

data = read.xlsx(here("09_bulk/RNA_seq/Disease_DEGs.xlsx"), sheet = 3)

temp_data = data[abs(data$log2FoldChange) > 1 & data$padj < 0.05, ]
temp_data = na.omit(temp_data[order(temp_data$log2FoldChange), ])
gene.label = c(temp_data[1:10, "symbol"], temp_data[(nrow(temp_data) - 9):nrow(temp_data), "symbol"])


## Set label
print(data[1, "log2FoldChange"])
print(log2(mean(as.numeric(expr[data[1, "gene_id"], lusc_samples])) / mean(as.numeric(expr[data[1, "gene_id"], bspn_samples]))))
print(mean(as.numeric(expr[data[1, "gene_id"], lusc_samples])) < mean(as.numeric(expr[data[1, "gene_id"], bspn_samples])))

# log2FoldChange > 1 is LUAD
# log2FoldChange < -1 is BSPN


data$Type = apply(data, 1, function(row) {
    if (is.na(row["padj"])) {
        row["padj"] = 1
    }
    
    if (as.numeric(row["log2FoldChange"]) > 1 && as.numeric(row["padj"]) < 0.05) {
        return("LUSC")
    } else if (as.numeric(row["log2FoldChange"]) < -1 && as.numeric(row["padj"]) < 0.05) {
        return("BSPN")
    } 
    return("Not")
})


## save LUAD vs LUSC data
addWorksheet(wb, "LUSC_vs_BSPN")
writeData(wb, 3, data)


p <- ggplot(data, aes(x=log2FoldChange, y=-log10(padj), color=Type)) + 
    geom_point() +
    scale_color_manual(values=c("LUAD"="#0084D1", "BSPN"="#9CD7D7", "LUSC"="#A0C807", "Not"="grey")) +
    geom_text_repel(
        data = data[data$symbol %in% gene.label, ], 
        aes(x=log2FoldChange, y=-log10(padj), color=Type, label=symbol),
        size = 5
    ) +
    labs(color = "", title = "RNA-seq LUAD vs. BSPN") +
    theme_bw() + 
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = c(0.85, 0.15),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill=alpha('white', 0))
    ) +
    guides(color = guide_legend(override.aes = list(size = 8)))

p

ggsave(
    filename = here("09_bulk/RNA_seq/volcano_plot_LUSC_vs_BSPN.pdf"),
    plot = p,
    width = 6,
    height = 6
)


saveWorkbook(wb, here("09_bulk/RNA_seq/Disease_DEGs.xlsx"), overwrite = T)


###############################################
## Make volcano plots of Stage specific
stages = c("I", "II", "III")

# LUAD
data = read.xlsx(here("09_bulk/RNA_seq/Stage_DEGs.xlsx"), sheet = 1)

for (i in sort(unique(data$ident))) {
    data_stage = data[data$ident == i, ]
    
    temp_data = data[abs(data$log2FoldChange) > 1 & data$padj < 0.05, ]
    temp_data = na.omit(temp_data[order(temp_data$log2FoldChange), ])
    gene.label = c(temp_data[1:10, "symbol"], temp_data[(nrow(temp_data) - 9):nrow(temp_data), "symbol"])
    
    
    data_stage$Type = apply(data_stage, 1, function(row) {
        if (is.na(row["padj"])) {
            row["padj"] = 1
        }
        
        if (as.numeric(row["log2FoldChange"]) > 1 && as.numeric(row["padj"]) < 0.05) {
            return("Up")
        } else if (as.numeric(row["log2FoldChange"]) < -1 && as.numeric(row["padj"]) < 0.05) {
            return("Down")
        } 
        return("Not")
    })
    
    p <- ggplot(data_stage, aes(x=log2FoldChange, y=-log10(padj), color=Type)) + 
        geom_point() +
        scale_color_manual(values=c("Up"="red", "Down"="green", "Not"="grey")) +
        geom_text_repel(
            data = data_stage[data_stage$symbol %in% gene.label, ], 
            aes(x=log2FoldChange, y=-log10(padj), color=Type, label=symbol),
            size = 5
        ) +
        labs(color = "", title = paste("LUAD", i)) +
        theme_bw() + 
        theme(
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.position = c(0.85, 0.15),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.background = element_rect(fill=alpha('white', 0))
        ) +
        guides(color = guide_legend(override.aes = list(size = 8)))
    
    ggsave(
        filename = here(paste0("09_bulk/RNA_seq/volcano_plot_LUAD_", i, ".pdf")),
        plot = p,
        width = 6,
        height = 6
    )
    
}



# LUSC
data = read.xlsx(here("09_bulk/RNA_seq/Stage_DEGs.xlsx"), sheet = 2)

for (i in sort(unique(data$ident))) {
    data_stage = data[data$ident == i, ]
    
    temp_data = data[abs(data$log2FoldChange) > 1 & data$padj < 0.05, ]
    temp_data = na.omit(temp_data[order(temp_data$log2FoldChange), ])
    gene.label = c(temp_data[1:10, "symbol"], temp_data[(nrow(temp_data) - 9):nrow(temp_data), "symbol"])
    
    
    data_stage$Type = apply(data_stage, 1, function(row) {
        if (is.na(row["padj"])) {
            row["padj"] = 1
        }
        
        if (as.numeric(row["log2FoldChange"]) > 1 && as.numeric(row["padj"]) < 0.05) {
            return("Up")
        } else if (as.numeric(row["log2FoldChange"]) < -1 && as.numeric(row["padj"]) < 0.05) {
            return("Down")
        } 
        return("Not")
    })
    
    p <- ggplot(data_stage, aes(x=log2FoldChange, y=-log10(padj), color=Type)) + 
        geom_point() +
        scale_color_manual(values=c("Up"="red", "Down"="green", "Not"="grey")) +
        geom_text_repel(
            data = data_stage[data_stage$symbol %in% gene.label, ], 
            aes(x=log2FoldChange, y=-log10(padj), color=Type, label=symbol),
            size = 5
        ) +
        labs(color = "", title = paste("LUAD", i)) +
        theme_bw() + 
        theme(
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.position = c(0.85, 0.15),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.background = element_rect(fill=alpha('white', 0))
        ) +
        guides(color = guide_legend(override.aes = list(size = 8)))
    
    ggsave(
        filename = here(paste0("09_bulk/RNA_seq/volcano_plot_LUSC_", i, ".pdf")),
        plot = p,
        width = 6,
        height = 6
    )
    
}



##############################################
## make overlap with scRNA-seq

overlap_degs <- c()

scrna <- foreach(f = list.files(here("03_each_cells/total/"), pattern="markers.xlsx", full.names = T, recursive = T), .combine="rbind") %dopar% {
    
    if (basename(dirname(f)) == "disease") {
        temp = read.xlsx(f, rowNames = T)
        temp$cell = basename(dirname(dirname(f)))
        temp$gene = as.character(rownames(temp))
        return(temp)
    }
    
    return(NULL)
}

scrna$avg_logFC <- as.numeric(scrna$avg_logFC)
scrna$p_val_adj <- as.numeric(scrna$p_val_adj)

scrna$Type = "Not"
scrna$Type[scrna$avg_logFC > 0.25 & scrna$p_val_adj < 0.05] <- "LUAD"
scrna$Type[scrna$avg_logFC < -0.25 & scrna$p_val_adj < 0.05] <- "LUSC"

table(scrna$Type)


bulk_degs <- read.xlsx(here("09_bulk/RNA_seq/Disease_DEGs.xlsx"), sheet = 1)
bulk_degs <- bulk_degs[bulk_degs$Type != "Not", "symbol"]
overlap_degs <- c(overlap_degs, bulk_degs)

scrna$Bulk = NA
scrna[scrna$Type != "Not", "Bulk"] = "scRNA-seq only"
scrna[is.na(scrna$Bulk), "Bulk"] = "Bulk Only"
scrna[scrna$gene %in% bulk_degs, "Bulk"] = "Both"

degs <- c(bulk_degs, scrna$gene[scrna$Type != "Not"])

scrna = scrna[scrna$gene %in% degs, ]


res = scrna %>% group_by(cell, Bulk) %>% add_tally()
res = res %>% dplyr::select(cell, Bulk, n) %>% unique()

level = c("Bulk Only", "Both", "scRNA-seq only")
res$Bulk = factor(res$Bulk, levels = level)

res = res %>% group_by(cell) %>% arrange(cell, desc(Bulk)) %>% mutate(ypos = cumsum(n))
res = as.data.frame(res)
res$cell = str_replace_all(res$cell, "_", " ")

non_immu = c(
    "Alveolar II",
    "Basal",
    "Ciliated",
    "Club",
    "Epithelial",
    "Fibroblasts",
    "Neuroendocrine"
)

res$Type = "Immu"
res$Type[res$cell %in% non_immu] = "Non-Immu"

res = res[order(res$n, decreasing = T), ]
cell_order <- res$cell[res$Bulk == "Both"]


short_names <- c(
    "APII"="ATII",
    "Alveolar II"="ATII",
    "Basal"="Basal",
    "B cells"="B",
    "CD4+"="CD4",
    "CD4"="CD4",
    "CD8+"="CD8",
    "CD8"="CD8",
    "Ciliated"="Cilia",
    "Club"="Club",
    "DC"="DC",
    "Dendritic"="DC",
    "Endothelial"="EC",
    "Epithelial"="Epi",
    "Exhaust T"="Exh T",
    "Fibroblasts"="Fib",
    "Granulocyte"="Gran",
    "Mast"="Mast",
    "Mascrophages"="Mφ",
    "Macrophages"="Mφ",
    "Monocytes"="Mo",
    "Neuroendocrine"="NE",
    "NK"="NK",
    "Tregs"="Tregs"
)

res$cell = sapply(res$cell, function(x) {
    if (!x %in% names(short_names)) {
        x
    } else {
        short_names[x]
    }
})

res <- res[res$cell != "Mφ", ]

svg(here("09_bulk/RNA_seq/DEGs_overlap_counts_bar_DESeq2.svg"), width = 8, height = 6)
ggplot(res, aes(x=cell, y=n, fill = Bulk, label = n)) +
    geom_bar(stat = "identity") +
    geom_text(aes(y=ypos)) +
    theme_bw() +
    theme(
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 21),
        title = element_text(size = 25),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill=alpha('white', 0)),
        strip.text = element_text(size = 20)
    ) +
    labs(x = "", y = "", fill = "") +
    facet_grid(.~Type, scales = "free", space = "free") +
    scale_fill_manual(values = c(
        "Bulk Only"="#46A9B5", 
        "Both"="#DDB52C", 
        "scRNA-seq only"="#E74B1C"
    ))
dev.off()


temp = res %>% 
    dplyr::select(cell, Bulk, n, Type) %>% 
    unique() %>%
    group_by(cell) %>% 
    mutate(p = n / sum(n)) %>%
    arrange(cell, desc(Bulk)) %>%
    mutate(labelpos = cumsum(p) - p / 2)


temp = temp[order(temp$p, decreasing = T), ]
cell_order <- temp$cell[temp$Bulk == "Both"]


svg(here("09_bulk/RNA_seq/DEGs_overlap_percentage_bar_DESeq2.svg"), width = 8, height = 6)
ggplot(temp, aes(x=cell, y=p, fill = Bulk)) +
geom_bar(stat = "identity") +
geom_text(aes(label = paste0(100* round(p, 2),"%"),y=labelpos),size = 5) +
theme_bw() +
theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 21),
    title = element_text(size = 25),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
    legend.position = "top",
    strip.text = element_text(size = 20)
) +
labs(x = "", y = "", fill = "") +
facet_grid(.~Type, scales = "free", space = "free") +
scale_fill_manual(values = c(
    "Bulk Only"="#46A9B5", 
    "Both"="#DDB52C", 
    "scRNA-seq only"="#E74B1C"
))
dev.off()


#################################


scrna <- foreach(f = list.files(here("03_each_cells/total/"), pattern="LUAD.xlsx", full.names = T, recursive = T), .combine="rbind") %dopar% {
    temp = read.xlsx(f, rowNames = T)
    temp$cell = basename(dirname(dirname(f)))
    temp$gene = as.character(rownames(temp))
    return(temp)
}

scrna$avg_logFC <- as.numeric(scrna$avg_logFC)
scrna$p_val_adj <- as.numeric(scrna$p_val_adj)

scrna$Type = "Not"
scrna$Type[abs(scrna$avg_logFC) > 0.25 & scrna$p_val_adj < 0.05] <- "LUAD"

table(scrna$Type)



bulk_degs <- read.xlsx("09_bulk/RNA_seq/Disease_DEGs.xlsx", sheet = 2)
bulk_degs <- bulk_degs[bulk_degs$Type != "Not", "symbol"]
overlap_degs <- c(overlap_degs, bulk_degs)

scrna$Bulk = NA
scrna[scrna$Type != "Not", "Bulk"] = "scRNA-seq only"
scrna[is.na(scrna$Bulk), "Bulk"] = "Bulk Only"
scrna[scrna$gene %in% bulk_degs, "Bulk"] = "Both"

degs <- c(bulk_degs, scrna$gene[scrna$Type != "Not"])

scrna = scrna[scrna$gene %in% degs, ]


res = scrna %>% group_by(cell, Bulk) %>% add_tally()
res = res %>% dplyr::select(cell, Bulk, n) %>% unique()

level = c("Bulk Only", "Both", "scRNA-seq only")
res$Bulk = factor(res$Bulk, levels = level)

res = res %>% group_by(cell) %>% arrange(cell, desc(Bulk)) %>% mutate(ypos = cumsum(n))

res = as.data.frame(res)
res$cell = str_replace_all(res$cell, "_", " ")

non_immu = c(
    "Alveolar II",
    "Basal",
    "Ciliated",
    "Club",
    "Epithelial",
    "Fibroblasts",
    "Neuroendocrine"
)

res$Type = "Immu"
res$Type[res$cell %in% non_immu] = "Non-Immu"

res = res[order(res$n, decreasing = T), ]
cell_order <- res$cell[res$Bulk == "Both"]

short_names <- c(
    "APII"="ATII",
    "Alveolar II"="ATII",
    "Basal"="Basal",
    "B cells"="B",
    "CD4+"="CD4",
    "CD4"="CD4",
    "CD8+"="CD8",
    "CD8"="CD8",
    "Ciliated"="Cilia",
    "Club"="Club",
    "DC"="DC",
    "Dendritic"="DC",
    "Endothelial"="EC",
    "Epithelial"="Epi",
    "Exhaust T"="Exh T",
    "Fibroblasts"="Fib",
    "Granulocyte"="Gran",
    "Mast"="Mast",
    "Mascrophages"="Mφ",
    "Macrophages"="Mφ",
    "Monocytes"="Mo",
    "Neuroendocrine"="NE",
    "NK"="NK",
    "Tregs"="Tregs"
)

res$cell = sapply(res$cell, function(x) {
    if (!x %in% names(short_names)) {
        x
    } else {
        short_names[x]
    }
})
res <- res[res$cell != "Mφ", ]

svg(here("09_bulk/RNA_seq/DEGs_overlap_counts_bar_LUAD_DESeq2.svg"), width = 8, height = 6)
ggplot(res, aes(x=cell, y=n, fill = Bulk, label = n)) +
    geom_bar(stat = "identity") +
    geom_text(aes(y=ypos)) +
    theme_bw() +
    theme(
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 21),
        title = element_text(size = 25),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill=alpha('white', 0)),
        strip.text = element_text(size = 20)
    ) +
    labs(x = "", y = "", fill = "") +
    facet_grid(.~Type, scales = "free", space = "free") +
    scale_fill_manual(values = c(
        "Bulk Only"="#46A9B5", 
        "Both"="#DDB52C", 
        "scRNA-seq only"="#E74B1C"
    ))
dev.off()



temp = res %>% 
    dplyr::select(cell, Bulk, n, Type) %>% 
    unique() %>%
    group_by(cell) %>% 
    mutate(p = n / sum(n)) %>%
    arrange(cell, desc(Bulk)) %>%
    mutate(labelpos = cumsum(p) - p / 2)

temp = temp[order(temp$p, decreasing = T), ]
cell_order <- temp$cell[temp$Bulk == "Both"]

temp$cell = factor(temp$cell, levels = cell_order)


svg(here("09_bulk/RNA_seq/DEGs_overlap_percentage_bar_LUAD_DESeq2.svg"), width = 8, height = 6)
ggplot(temp, aes(x=cell, y=p, fill = Bulk)) +
geom_bar(stat = "identity") +
geom_text(aes(label = paste0(100* round(p, 2),"%"),y=labelpos),size = 5) +
theme_bw() +
theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 21),
    title = element_text(size = 25),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
    legend.position = "top",
    strip.text = element_text(size = 20)
) +
labs(x = "", y = "", fill = "") +
facet_grid(.~Type, scales = "free", space = "free") +
scale_fill_manual(values = c(
    "Bulk Only"="#46A9B5", 
    "Both"="#DDB52C", 
    "scRNA-seq only"="#E74B1C"
))
dev.off()



#################################
## 


scrna <- foreach(f = list.files("03_each_cells/total/", pattern="LUSC.xlsx", full.names = T, recursive = T), .combine="rbind") %dopar% {
    temp = read.xlsx(f, rowNames = T)
    temp$cell = basename(dirname(dirname(f)))
    temp$gene = as.character(rownames(temp))
    return(temp)
}

scrna$avg_logFC <- as.numeric(scrna$avg_logFC)
scrna$p_val_adj <- as.numeric(scrna$p_val_adj)

scrna$Type = "Not"
scrna$Type[abs(scrna$avg_logFC) > 0.5 & scrna$p_val_adj < 0.05] <- "LUSC"

table(scrna$Type)


bulk_degs <- read.xlsx("09_bulk/RNA_seq/Disease_DEGs.xlsx", sheet = 3)
bulk_degs <- bulk_degs[bulk_degs$Type != "Not", "symbol"]
overlap_degs <- c(overlap_degs, bulk_degs)

scrna$Bulk = NA
scrna[scrna$Type != "Not", "Bulk"] = "scRNA-seq only"
scrna[is.na(scrna$Bulk), "Bulk"] = "Bulk Only"
scrna[scrna$gene %in% bulk_degs, "Bulk"] = "Both"

degs <- c(bulk_degs, scrna$gene[scrna$Type != "Not"])

scrna = scrna[scrna$gene %in% degs, ]


res = scrna %>% group_by(cell, Bulk) %>% add_tally()
res = res %>% dplyr::select(cell, Bulk, n) %>% unique()

level = c("Bulk Only", "Both", "scRNA-seq only")
res$Bulk = factor(res$Bulk, levels = level)

res = res %>% group_by(cell) %>% arrange(cell, desc(Bulk)) %>% mutate(ypos = cumsum(n))
res = as.data.frame(res)
res$cell = str_replace_all(res$cell, "_", " ")
res$Type = "Immu"
res$Type[res$cell %in% non_immu] = "Non-Immu"

res = res[order(res$n, decreasing = T), ]
cell_order <- res$cell[res$Bulk == "Both"]

short_names <- c(
    "APII"="ATII",
    "Alveolar II"="ATII",
    "Basal"="Basal",
    "B cells"="B",
    "CD4+"="CD4",
    "CD4"="CD4",
    "CD8+"="CD8",
    "CD8"="CD8",
    "Ciliated"="Cilia",
    "Club"="Club",
    "DC"="DC",
    "Dendritic"="DC",
    "Endothelial"="EC",
    "Epithelial"="Epi",
    "Exhaust T"="Exh T",
    "Fibroblasts"="Fib",
    "Granulocyte"="Gran",
    "Mast"="Mast",
    "Mascrophages"="Mφ",
    "Macrophages"="Mφ",
    "Monocytes"="Mo",
    "Neuroendocrine"="NE",
    "NK"="NK",
    "Tregs"="Tregs"
)

res$cell = sapply(res$cell, function(x) {
    if (!x %in% names(short_names)) {
        x
    } else {
        short_names[x]
    }
})
res <- res[res$cell != "Mφ", ]

svg(here("09_bulk/RNA_seq/DEGs_overlap_counts_bar_LUSC_DESeq2.svg"), width = 8, height = 6)
ggplot(res, aes(x=cell, y=n, fill = Bulk, label = n)) +
geom_bar(stat = "identity") +
geom_text(aes(y=ypos)) +
theme_bw() +
theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 21),
    title = element_text(size = 25),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.background = element_rect(fill=alpha('white', 0)),
    strip.text = element_text(size = 20)
) +
labs(x = "", y = "", fill = "") +
facet_grid(.~Type, scales = "free", space = "free") +
scale_fill_manual(values = c(
    "Bulk Only"="#46A9B5", 
    "Both"="#DDB52C", 
    "scRNA-seq only"="#E74B1C"
))
dev.off()



temp = res %>% 
    dplyr::select(cell, Bulk, n, Type) %>% 
    unique() %>%
    group_by(cell) %>% 
    mutate(p = n / sum(n)) %>%
    arrange(cell, desc(Bulk)) %>%
    mutate(labelpos = cumsum(p) - p / 2)


temp = temp[order(temp$p, decreasing = T), ]
cell_order <- temp$cell[temp$Bulk == "Both"]

temp$cell = factor(temp$cell, levels = cell_order)


svg(here("09_bulk/RNA_seq/DEGs_overlap_percentage_bar_LUSC_DESeq2.svg"), width = 8, height = 6)
ggplot(temp, aes(x=cell, y=p, fill = Bulk)) +
geom_bar(stat = "identity") +
geom_text(aes(label = paste0(100* round(p, 2),"%"),y=labelpos),size = 5) +
theme_bw() +
theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 21),
    title = element_text(size = 25),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
    legend.position = "top",
    strip.text = element_text(size = 20)
) +
labs(x = "", y = "", fill = "") +
facet_grid(.~Type, scales = "free", space = "free") +
scale_fill_manual(values = c(
    "Bulk Only"="#46A9B5", 
    "Both"="#DDB52C", 
    "scRNA-seq only"="#E74B1C"
))
dev.off()



##
library(ggpubr)
library(doMC)
library(scales)

bulk <- read.xlsx("09_bulk/RNA_seq/RSEM.xlsx", rowNames = T)
gene_pos <- read.table("09_bulk/RNA_seq/gene_pos.txt", header = F, row.names = 5)
gene_pos$V6 <- make.unique(as.character(gene_pos$V6))
rownames(bulk) <- gene_pos[rownames(bulk), "V6"]


target_genes = c(
    "CYP2B7P", "CKS1B", "IGHG1",
    "MUC1", "PIGR", "SFTPA1", 
    "SFTPA2", "PIGR", "VIM", 
    "SPINK1", "PTMA", "IGFBP7"
)


make_violin_plot <- function(gene) {
    temp = as.data.frame(t(bulk[gene, ]))
    temp$Disease = sapply(rownames(temp), function(x) {
        str_replace_all(x, "\\d+", "")
    })
    colnames(temp)[1] = "value"
    temp = temp[temp$Disease %in% c("LUAD", "LUSC"), ]
    
    
    ggviolin(
        temp, 
        x = "Disease",
        y = "value",
        add = "boxplot",
        fill = "Disease",
        add.params = list(fill = "white"),
        palette = c("LUAD"="#0084D1", "BSPN"="#9CD7D7", "LUSC"="#A0C807", "Not"="grey")
    ) +
        scale_y_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        ) +
        labs(x = "", y = "") +
        theme_bw() +
        theme(
            axis.text = element_text(size = 15),
            legend.position = "none",
            text = element_text(size = 10)
        ) +
        stat_compare_means(comparisons = list(c("LUAD", "LUSC")))
}


registerDoMC(10)
res = foreach(i = target_genes, .errorhandling = "pass") %dopar% {
    ggsave(
        filename = paste0("09_bulk/RNA_violin/", i, ".pdf"),
        plot = make_violin_plot(i),
        width = 3,
        height = 3
    )
}