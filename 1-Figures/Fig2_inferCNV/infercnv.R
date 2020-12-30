options(stringsAsFactors = T)

library(infercnv)
library(Seurat)
library(ggfortify)
library(ggplot2)
library(Rtsne)
library(umap)


setwd("LungCancer10x/")


### ATII
obj <- readRDS("03_each_cells/total/Alveolar_II/seurat.rds")


cells = rownames(obj@meta.data)

temp_file = "anno_ATII"

temp = obj@meta.data[cells, c("Cells", "Disease")]
for (i in 1:ncol(temp)) {
    temp[, i] <- as.character(temp[, i])
}

write.table(temp, temp_file, row.names = F, col.names = F, sep="\t", quote = F)

infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix=obj@raw.data[, cells],
    gene_order_file="12_final_plots/gene_pos.txt",
    annotations_file=temp_file,
    delim="\t",
    ref_group_names=c("LUAD_Normal", "LUSC_Normal")
)



file.remove(temp_file)

infercnv_obj = infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=paste0("11_CNV/ATII"), 
    cluster_by_groups=TRUE, 
    denoise=TRUE,
    HMM=TRUE
)



### Basal
obj <- readRDS("03_each_cells/total/Basal/seurat.rds")


cells = rownames(obj@meta.data)

temp_file = "anno_Basal"

temp = obj@meta.data[cells, c("Cells", "Disease")]
for (i in 1:ncol(temp)) {
    temp[, i] <- as.character(temp[, i])
}

write.table(temp, temp_file, row.names = F, col.names = F, sep="\t", quote = F)

infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix=obj@raw.data[, cells],
    gene_order_file="12_final_plots/gene_pos.txt",
    annotations_file=temp_file,
    delim="\t",
    ref_group_names=c("LUAD_Normal")
)



# file.remove(temp_file)

infercnv_obj = infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=paste0("11_CNV/Basal"), 
    cluster_by_groups=TRUE, 
    denoise=TRUE,
    HMM=TRUE
)


## Compare Non-immune cells of LUAD and LUSC


expr <- readRDS("02_rds/all_cell_expr.rds")
meta <- readRDS("02_rds/meta.rds")


make_infercnv <- function(expr, meta, control="LUAD_Normal", treat="LUAD") {
    meta = meta[meta$Disease %in% c(control, treat), ]
    cells = rownames(meta)
    
    temp_file = paste("anno", treat, sep = "_")
    
    temp = meta[, c("Cells", "Disease")]
    for (i in 1:ncol(temp)) {
        temp[, i] <- as.character(temp[, i])
    }
    
    write.table(temp, temp_file, row.names = F, col.names = F, sep="\t", quote = F)
    
    infercnv_obj = CreateInfercnvObject(
        raw_counts_matrix=expr[, cells],
        gene_order_file="12_final_plots/gene_pos.txt",
        annotations_file=temp_file,
        delim="\t",
        ref_group_names=c(control)
    )
    
    out_dir = paste0("11_CNV/", treat)
    
    dir.create(out_dir, showWarnings = F, recursive = T)
    
    infercnv_obj = infercnv::run(
        infercnv_obj,
        cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
        out_dir=out_dir, 
        cluster_by_groups=TRUE, 
        denoise=TRUE,
        HMM=TRUE
    )
}

make_infercnv(expr, meta, control="LUAD_Normal", treat="LUAD")
make_infercnv(expr, meta, control="LUSC_Normal", treat="LUSC")



read_infercnv_obj <- function(path, filter=FALSE, zscore=TRUE) {
    sample = basename(dirname(path))
    sample = str_split(sample, "_")[[1]]
    
    if (length(sample) == 2) {
        cell_name = sample[1]
    } else {
        cell_name = paste(sample[1], sample[2])
    }
    
    sample = sample[length(sample)]
    
    infercnv_obj <- readRDS(path)
    
    gene_order <- as.data.frame(infercnv_obj@gene_order)
    gene_order$cell = cell_name
    
    if(zscore) {
        expr <- apply(infercnv_obj@expr.data, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))
    } else {
        expr = infercnv_obj@expr.data
    }
    
    if (filter) {
        expr <- expr[apply(expr, 1, function(row) {
            sum(row > 3) > (length(row) / 10)
        }), ]
        
    }
    
    row_orders <- intersect(rownames(expr), rownames(gene_order))
    cbind(expr[row_orders, ], gene_order[row_orders, ])
}



### make heatmap plot of infercnv
make_infercnv_heatmap <- function(path, gene_mark) {
    temp_data <- read_infercnv_obj(path, filter = F, zscore = T)
    temp_data <- temp_data[, !colnames(temp_data) %in% c("chr", "start", "stop", "cell")]
    
    
    temp_meta <- meta[colnames(temp_data), ]
    temp_meta <- temp_meta[order(temp_meta$Disease, temp_meta$Stage), ]
    temp_meta$Disease[temp_meta$Disease == 'LUAD_Normal'] <- "Normal (LUAD)"
    temp_meta$Disease[temp_meta$Disease == 'LUSC_Normal'] <- "Normal (LUSC)"
    
    temp_cells <- rownames(temp_meta)[temp_meta$cell_short == cell]
    
    la <- rowAnnotation(
        Cell = temp_meta[as.character(temp_cells), "cell_short"],
        Disease = temp_meta[as.character(temp_cells), "Disease"],
        Stage = temp_meta[as.character(temp_cells), "Stage"],
        col = list(
            Cell = cell_color,
            Disease = disease_colors,
            Stage = stage_color
        )
    )
    
    temp_pos = data.frame(
        gene = rownames(temp_data),
        chrom = gene_pos[rownames(temp_data), "V1"]
    )
    
    temp_pos <- temp_pos[order(temp_pos$chrom), ]
    
    # gene_mark = which(as.character(temp_pos$gene) %in% as.character(gene_mark))
    # 
    # ba <- HeatmapAnnotation(
    #     # chrom = temp_pos$chrom,
    #     gene = anno_mark(
    #         at=gene_mark,
    #         labels=temp_pos$gene[gene_mark],
    #         side = "bottom"
    #     )
    #     # col = list(
    #     #     chrom = chrom_color
    #     # )
    # )
    # 
    # print(dim(temp_data))
    # print(dim(temp_meta))
    
    pos_level = as.character(c(1:22, "X", "Y"))
    all_temp_chrom = sort(unique(temp_pos$chrom))
    pos_level <- c(pos_level, all_temp_chrom[!all_temp_chrom %in% pos_level])
    
    temp_pos$chrom <- factor(temp_pos$chrom, levels = pos_level)
    temp_pos <- temp_pos[order(temp_pos$chrom), ]
    
    Heatmap(
        t(as.matrix(temp_data[temp_pos$gene, as.character(temp_cells)])),
        name = "Infercnv",
        show_row_names = F,
        # left_annotation = la,
        row_title_rot = 0,
        border = T,
        show_column_names = F,
        column_names_rot = 90,
        # bottom_annotation = ba,
        cluster_columns = F,
        cluster_rows = F,
        column_split = as.character(temp_pos$chrom),
        column_title_rot = 90
    )
}


cell_color = c(
    wes_palette("Darjeeling1"), 
    wes_palette("Moonrise3"),
    wes_palette("Darjeeling2"),
    wes_palette("Moonrise1")
)[1:length(unique(meta$cell_name1))]

names(cell_color) = unique(meta$cell_name1)

disease_colors = c(
    "LUAD" = "#0084D1", 
    "Normal (LUAD)"="#FECC1B", 
    "LUAD_Normal"="#FECC1B", 
    "Normal"="#73BDFF", 
    "Normal (LUSC)"="#778793", 
    "LUSC_Normal"="#778793", 
    "LUSC"="#A0C807", 
    "BSPN"="#6A0019", 
    "BSPN_Normal"="#C5000B"
)


chrom_color <- c(
    wes_palette("Darjeeling1"), 
    wes_palette("Darjeeling2"), 
    wes_palette("Moonrise1"), 
    wes_palette("Moonrise2"), 
    wes_palette("Moonrise3"),
    wes_palette("Royal1"),
    wes_palette("Royal2"),
    wes_palette("Chevalier1"),
    wes_palette("FantasticFox1"),
    wes_palette("GrandBudapest1"),
    wes_palette("GrandBudapest2")
)
chrom_color <- chrom_color[1:length(unique(gene_pos$V1))]
names(chrom_color) <- unique(gene_pos$V1)


stage_color =  c(
    "I"="#65A9A3", 
    "II"="#4A933E", 
    "III"="#EC7A21", 
    "IV"="#D73F47", 
    "LUAD_Normal"="#FECC1B", 
    "LUSC_Normal"="#778793",
    "Normal(LUAD)"="#FECC1B", 
    "Normal(LUSC)"="#778793"
    
)




### Make plots
meta <- readRDS("02_rds/meta.rds")

gene_pos <- read.table("09_bulk/RNA_seq/gene_pos.txt")
rownames(gene_pos) <- make.unique(as.character(gene_pos$V6))


temp_data <- read_infercnv_obj("11_CNV/ATII/run.final.infercnv_obj", filter = F, zscore = T)
temp_data <- temp_data[, !colnames(temp_data) %in% c("chr", "start", "stop", "cell")]


temp_meta <- meta[colnames(temp_data), ]
temp_meta <- temp_meta[order(temp_meta$Disease, temp_meta$Stage), ]
temp_meta$Disease[temp_meta$Disease == 'LUAD_Normal'] <- "Normal (LUAD)"
temp_meta$Disease[temp_meta$Disease == 'LUSC_Normal'] <- "Normal (LUSC)"

# temp_meta <- temp_meta[temp_meta$Disease %in% c("LUAD", "LUSC"), ]

temp_cells <- rownames(temp_meta)[temp_meta$cell_short == cell]

la <- rowAnnotation(
    Disease = temp_meta[as.character(temp_cells), "Disease"],
    Stage = temp_meta[as.character(temp_cells), "Stage"],
    col = list(
        Disease = disease_colors,
        Stage = stage_color
    )
)

temp_pos = data.frame(
    gene = rownames(temp_data),
    chrom = gene_pos[rownames(temp_data), "V1"]
)

temp_pos <- temp_pos[order(temp_pos$chrom), ]

pos_level = as.character(c(1:22, "X", "Y"))
all_temp_chrom = sort(unique(temp_pos$chrom))
pos_level <- c(pos_level, all_temp_chrom[!all_temp_chrom %in% pos_level])

temp_pos$chrom <- factor(temp_pos$chrom, levels = pos_level)
temp_pos <- temp_pos[order(temp_pos$chrom), ]
temp_pos <- temp_pos[!is.na(temp_pos$chrom), ]


pdf("11_CNV/ATII/heatmap.pdf", width = 20, height = 10)
Heatmap(
    t(as.matrix(temp_data[temp_pos$gene, as.character(temp_cells)])),
    col = colorRamp2(c(-1.15, -0.5 , 0, 0.5, 1.15), c("#000074", "white", "white", "white", "#760002")),
    name = "Infercnv",
    show_row_names = F,
    left_annotation = la,
    row_title_rot = 0,
    border = T,
    show_column_names = F,
    column_names_rot = 90,
    # bottom_annotation = ba,
    cluster_columns = F,
    cluster_rows = T,
    row_split = as.character(temp_meta[temp_cells, "Disease"]),
    column_split = as.character(temp_pos$chrom),
    column_title_rot = 90
)
dev.off()


### Correlation between CNV
plot_sd_of_pca <- function(pca) {
    data = data.frame(
        x=1:ncol(pca$x),
        y=apply(pca$x, 2, sd)
    )
    
    ggplot(data, aes(x=x, y=y)) + geom_point()
}


temp_data <- temp_data[order(apply(temp_data, 1, sd), decreasing = T), ]

res.pca = prcomp(t(temp_data[1:500, ]))

autoplot(res.pca, data = meta[colnames(temp_data), ], colour = "Disease", frame = TRUE, frame.type = 'norm')
autoplot(res.pca, data = meta[colnames(temp_data), ], colour = "Stage", frame = TRUE, frame.type = 'norm')

plot_sd_of_pca(res.pca)



res.cor = cor(temp_data)

ggcorrplot::ggcorrplot(res.cor)