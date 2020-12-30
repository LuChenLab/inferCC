library(destiny)
library(Seurat)
library(wesanderson)
library(openxlsx)
library(harmony)
library(dplyr)
library(ggrastr)
library(scatterpie)
library(reticulate)
library(reshape2)
library(tidyr)
library(ggthemes)
library(doMC)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(wesanderson)
library(org.Hs.eg.db)
library(stringr)
library(ComplexHeatmap)
library(circlize)


root.dir = "LungCancer10x/"
full.path <- function(...) { return(paste(root.dir, ..., sep = "/")) }


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

stage_colors = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793")
disease_colors = c("LUAD" = "#0084D1", "LUAD_Normal"="#FECC1B", "Normal"="#73BDFF", "LUSC_Normal"="#778793", "LUSC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B")
tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")

labels_size = 12
title_size = 15

kee = import("kneed")

make_seurat_pipe <- function(expr, meta, n.pcs = 10) {
    all_obj <- CreateSeuratObject(
        expr[, rownames(meta)],
        meta.data = meta
    )
    
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = all_obj@raw.data), value = TRUE)
    percent.mito <- Matrix::colSums(all_obj@raw.data[mito.genes, ])/Matrix::colSums(all_obj@raw.data)
    
    all_obj <- AddMetaData(object =all_obj, metadata = percent.mito, col.name = "percent.mito")
    
    pcs = min(100, min(nrow(all_obj@raw.data), ncol(all_obj@raw.data)) - 1)
    
    all_obj <- NormalizeData(object = all_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    all_obj <- FindVariableGenes(object = all_obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    all_obj <- ScaleData(object = all_obj, vars.to.regress = c("nUMI", "percent.mito"))
    all_obj <- RunPCA(object = all_obj, pc.genes = all_obj@var.genes, do.print = FALSE, pcs.compute = pcs)
    
    all_obj <- RunHarmony(all_obj, c("batch", "PatientID"))
    
    #### 4. make dotplot of SD of PCAs and Harmonys
    tryCatch({
        kee = import("kneed")
        n.pcs = kee$KneeLocator(
            1:ncol(all_obj@dr$pca@cell.embeddings),
            apply(all_obj@dr$pca@cell.embeddings, 2, sd),
            curve='convex',
            direction='decreasing'
        )
        n.pcs = n.pcs$knee
    }, error = function(e) {
        n.pcs = 10
    })
    
    #### 5. tSNE and UMAP
    all_obj <- RunTSNE(
        object = all_obj, 
        dims.use = 1:n.pcs, 
        do.fast = TRUE, 
        reduction.use = "harmony", 
        reduction.name = "tsne",
        perplexity=min(30, floor(pcs / 3) + 1)
    )
    all_obj <- RunUMAP(
        object = all_obj, 
        dims.use = 1:n.pcs, 
        do.fast = TRUE, 
        reduction.use = "harmony",
        n_neighbors = min(30, pcs)
    )
    all_obj <- FindClusters(
        object = all_obj, reduction.type = "harmony", dims.use = 1:10, 
        resolution = c(0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
        print.output = 0, save.SNN = FALSE, force.recalc = TRUE
    )
    
    return (all_obj)
}


calculate_cluster_perc <- function(obj, seed.range = 1:50) {
    meta = obj@meta.data
    
    num = table(meta$Stage)
    
    sf = num / min(num)
    
    # sf = table(meta$ident)
    # sf = sf / min(sf)
    
    temp <- meta %>%
        group_by(ident, Stage) %>%
        add_tally() %>%
        dplyr::select(ident, Stage, n) %>%
        unique() %>%
        as.data.frame()
    
    temp$freq <- apply(temp, 1, function(row) {
        as.numeric(row[3]) / sf[[row[2]]] # / num[[row[2]]]
    })
    
    temp <- temp %>% group_by(ident) %>%
        mutate(freq = freq / sum(freq)) %>%
        as.data.frame()
    
    temp
}


PercentAbove <- function(x, threshold){
    return(length(x = x[x > threshold]) / length(x = x))
}


ModifiedDotPlot <- function(
    object,
    group.by,
    markers,
    cols.use = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA,
    x.axis.rotate = FALSE
) {
    scale.func <- switch(
        EXPR = scale.by,
        'size' = scale_size,
        'radius' = scale_radius,
        stop("'scale.by' must be either 'size' or 'radius'")
    )
    if (!missing(x = group.by)) {
        object <- SetAllIdent(object = object, id = group.by)
    }
    
    genes.plot <- intersect(row.names(object@raw.data), markers$gene)
    
    data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
    colnames(x = data.to.plot) <- genes.plot
    data.to.plot$cell <- rownames(x = data.to.plot)
    data.to.plot$id <- object@ident
    data.to.plot %>% gather(
        key = genes.plot,
        value = expression,
        -c(cell, id)
    ) -> data.to.plot
    data.to.plot %>%
        group_by(id, genes.plot) %>%
        summarize(
            avg.exp = mean(expm1(x = expression)),
            pct.exp = PercentAbove(x = expression, threshold = 0)
        ) -> data.to.plot
    data.to.plot %>%
        ungroup() %>%
        group_by(genes.plot) %>%
        mutate(avg.exp.scale = scale(x = avg.exp)) %>%
        mutate(avg.exp.scale = MinMax(
            data = avg.exp.scale,
            max = col.max,
            min = col.min
        )) ->  data.to.plot
    data.to.plot$genes.plot <- factor(
        x = data.to.plot$genes.plot,
        levels = rev(x = genes.plot)
    )
    # data.to.plot$genes.plot <- factor(
    #   x = data.to.plot$genes.plot,
    #   levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
    # )
    data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
    data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
    data.to.plot <- merge(data.to.plot, markers, by.x = "genes.plot", by.y = "gene")
    
    p <- ggplot(data = data.to.plot, aes(x = id, y = genes.plot, size = pct.exp, color = avg.exp.scale)) + 
        geom_point() + 
        # guides(
        #     color=guide_legend(
        #         ncol = 2,
        #         title.position = "top"
        #     ),
        #     size=guide_legend(
        #         ncol = 1
        #     ),
        #     alpha=guide_legend(
        #         ncol = 1
        #     )
    # ) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
        facet_grid(ident ~ ., scale = "free", space = "free") +
        scale_color_gradient(low = "lightgrey", high = "blue") +
        theme_bw(base_family = "Arial Unicode MS")
    
    if (x.axis.rotate) {
        p = p + theme(
            panel.grid = element_blank(),
            legend.position = "none",
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1),
            strip.text = element_text(size = 15)
        )
    } else {
        p = p + theme(
            panel.grid = element_blank(),
            legend.position = "none",
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.x = element_text(hjust = 0.5, vjust = 1),
            strip.text = element_text(size = 15)
            # axis.title.y = element_text(angle = 90)
        )
    }
    
    
    return(p)
}


### Diffusion map
make_dc_plot <- function(dpt, obj, output, dims = 6, title = "", subtitle = "") {
    
    mtx = as.data.frame(dpt@dm@eigenvectors)
    rownames(mtx) <- rownames(dpt@dm@data_env$data)
    mtx$Stage <- obj@meta.data[rownames(mtx), "Stage"]
    
    for(i in seq(1, dims, 2)) {
        temp <- mtx[, c(paste0("DC", i), paste0("DC", i + 1), "Stage")]
        colnames(temp) <- c("DC1", "DC2", "Stage")
        
        p <- ggplot(temp, aes(x=DC1, y = DC2, color = Stage)) +
            geom_point_rast() +
            scale_color_manual(values = stage_colors) +
            theme_bw(base_family = "Arial Unicode MS") +
            theme(
                aspect.ratio = 1,
                panel.grid = element_blank(),
                # legend.position = "bottom",
                axis.text = element_text(size = 15),
                axis.title = element_text(size = 18),
                legend.text = element_text(size = 15),
                legend.background = element_blank(),
                plot.title = element_text(size = 15),
                plot.subtitle = element_text(size = 10)
            ) +
            guides(color = guide_legend(override.aes = list(size = 5))) +
            labs(
                color = "",
                x = paste0("DC", i), 
                y = paste0("DC", i + 1), 
                title = title,
                subtitle = subtitle
            )
        
        ggsave(
            filename = paste0(output, "_DM_", i, ".pdf"),
            plot = p,
            width = 6,
            height = 7, 
            device = cairo_pdf
        )
    }
}


find_markers <- function(object, group.by = "Stage", n.cores = 4, min.pct = 0.1) {
    
    temp_group <- as.character(sort(unique(object@meta.data[, group.by])))
    
    library(doMC)
    
    registerDoMC(n.cores)
    res =foreach(i = temp_group, .combine = "rbind") %dopar% {
        
        groups = rownames(object@meta.data[object@meta.data[, group.by] == i, ])
        
        if(length(groups) >= 3) {
            temp <- FindMarkers(
                object = object, 
                ident.1 = groups,
                logfc.threshold = 0,
                min.pct = min.pct
            )
            
            temp$ident = i
            temp$gene = rownames(temp)
        } 
        
        temp
    }    
    return(res)
}


make_scatterpie_plot <- function(
    object, 
    color="Disease", 
    alpha=0.8, 
    colors=colors, 
    pt.size = 0.1, 
    guide_ncol=1,
    text_size = 8,
    reduction.use = "tsne",
    r=0.05,
    legend.position="none",
    legend.direction="horizontal",
    legend.breaks=waiver(),
    custom_pos = NULL,
    legend.text = 15
) {
    meta = object@meta.data
    
    coord = object@dr[[reduction.use]]@cell.embeddings
    colnames(coord) = c("x", "y")
    
    meta = merge(meta, coord, by = 0)
    
    scatter = eval(
        parse(
            text=paste0(
                "meta %>% group_by(ident, ", color, ") %>% add_tally() %>% group_by(ident) %>% mutate(x1 = median(x), y1 = median(y), perc = n / sum(n)) %>% dplyr::select(ident, x1, y1, perc, ", color, ") %>% unique()"
            )
        )
    )
    
    colnames(scatter)[colnames(scatter) == color] = "group"
    
    scatter = dcast(scatter, ident + x1 + y1 ~ group, value.var="perc")
    
    scatter[is.na(scatter)] = 0
    scatter$r = min(max(meta$x) - min(meta$x), max(meta$y) - min(meta$y)) * r
    
    clt_color = gg_color_hue(length(unique(meta$ident)))
    names(clt_color) = 1:length(unique(meta$ident))
    
    if (!is.null(custom_pos)) {
        scatter$x1 = custom_pos[scatter$ident, "x1"]
        scatter$y1 = custom_pos[scatter$ident, "y1"]
    }
    
    p <- ggplot() + 
        geom_point_rast(
            aes(x=x, y=y, colour=factor(ident)), 
            data=meta, 
            alpha=alpha, 
            size=pt.size
        ) +
        geom_scatterpie(
            aes(
                x=x1,
                y=y1, 
                group=ident,
                r=r
            ),
            data=scatter,
            cols=as.character(unique(meta[, color])),
            color=NA
        ) +
        coord_equal() +
        geom_text(
            aes(x=x1, y=y1 - r - 0.5, label=ident),
            data=scatter,
            size = text_size
        ) +
        theme_bw(base_family = "Arial Unicode MS") +
        theme(
            panel.grid = element_blank(),
            legend.position = legend.position,
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = legend.text),
            # legend.background = element_blank()
            legend.direction=legend.direction,
            legend.spacing.x = unit(0.1, 'cm')
        ) +
        scale_color_manual(values=clt_color, breaks = c()) + 
        scale_fill_manual(values=colors, breaks = legend.breaks) +
        guides(fill = guide_legend(ncol = guide_ncol))
    
    if (reduction.use == "tsne") {
        p <- p + labs(x="tSNE_1", y="tSNE_2", color="", fill = "")
    } else if (reduction.use == "umap") {
        p <- p + labs(x="UMAP1", y="UMAP2", color="", fill="")
    }
    
    return(p)
}


make_cluster_markers_heatmap <- function(obj, markers, logfc = 0.5) {
    set.seed(42)
    markers <- markers %>%
      filter(p_val_adj < 0.05 & avg_logFC > logfc) %>%
      group_by(ident) %>%
      arrange(desc(avg_logFC)) %>%
      as.data.frame()
    markers$x = 1:nrow(markers)
    
    temp <- markers %>%
      filter(p_val_adj < 0.05 & avg_logFC > 0.25) %>%
      group_by(ident) %>%
      top_n(10, wt=avg_logFC) %>%
      as.data.frame()
    
    obj@meta.data$ident = obj@meta.data$ident
    
    obj@meta.data$Cells = rownames(obj@meta.data)
    cells_order <- obj@meta.data %>%
      group_by(ident) %>%
      arrange(ident) %>%
      sample_n(min(table(obj@meta.data$ident))) %>%
      as.data.frame()
    
    # cells_order <- cells_order[sample(1:nrow(cells_order)), ]
    # cells_order = cells_order[order(as.numeric(cells_order$ident)), ]
    
    ccls <- as.character(wes_palette("GrandBudapest2", length(unique(cells_order$ident)), type = "continuous"))
    names(ccls) <- unique(cells_order$ident)
    
    ta = HeatmapAnnotation(
      Cluster = cells_order$ident,
      col = list(
        Cluster = ccls
      ),
      show_annotation_name = F,
      show_legend = F,
      annotation_name_gp = gpar(fontfamily = "Arial Unicode MS"),
      annotation_legend_param = list(
        labels_gp = gpar(fontfamily = "Arial Unicode MS"),
        title_gp = gpar(fontfamily = "Arial Unicode MS")
      )
    )
    
    ha = rowAnnotation(
      foo = anno_mark(
        at = temp$x, 
        labels = temp$gene,
        side = "left"
      ),
      annotation_name_gp = gpar(fontfamily = "Arial Unicode MS"),
      annotation_legend_param = list(
        labels_gp = gpar(fontfamily = "Arial Unicode MS"),
        title_gp = gpar(fontfamily = "Arial Unicode MS")
      )
    )
    
    
    # pdf(paste0(root.dir, "/cluster_markers.pdf"), width = 5, height = 6)
    Heatmap(
      obj@scale.data[markers$gene, cells_order$Cells],
      name = "Expr",
      cluster_rows = F,
      cluster_columns = F,
      show_row_names = F,
      show_column_names = F,
      show_heatmap_legend = F,
      top_annotation = ta,
      left_annotation = ha,
      column_split = cells_order$ident,
      row_split = markers$ident,
      row_title = NULL,
      heatmap_legend_param = list(
        labels_gp = gpar(fontsize = 15, fontfamily = "Arial Unicode MS"),
        title_gp = gpar(fontsize = 15, fontfamily = "Arial Unicode MS")
      ),
      col = colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow")),
      row_title_gp = gpar(fontfamily = "Arial Unicode MS"),
      column_title_gp = gpar(fontfamily = "Arial Unicode MS"),
      use_raster = T
    )
    # dev.off()
}



meta <- readRDS(full.path("11_CNV/meta.rds"))
expr <- readRDS(full.path("02_rds/all_cell_expr.rds"))


## 1st
for(i in unique(meta$cell_short)) {
    
    if (i %in% c("ATII", "Basal")) {
        next
    }
    
    for (j in c("LUAD", "LUSC")) {
        outdir = full.path("11_CNV/each_cells/", i, "/", j, "/")
        print(outdir)
        temp_meta <- meta[as.character(meta$cell_short) == i & str_detect(as.character(meta$Disease), j) & meta$Malignant, ]
        
        if (nrow(temp_meta) > 100) {
            
            tryCatch({
                dir.create(outdir, showWarnings = F, recursive = T)
                
                if (file.exists(paste0(outdir, "seurat_obj.rds"))) {
                    all_obj <- readRDS(paste0(outdir, "seurat_obj.rds"))
                } else {
                    all_obj <- make_seurat_pipe(expr, temp_meta)
                    
                    saveRDS(all_obj, paste0(outdir, "seurat_obj.rds"))
                }
                
                all_obj@meta.data$ident1 <- as.numeric(all_obj@meta.data$res.0.1) + 1
                all_obj@meta.data$ident <- all_obj@meta.data$ident1
                
                for (c in c(0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)) {
                    p <- DimPlot(
                        all_obj, 
                        reduction.use = "umap", 
                        pt.size = 0.1, 
                        do.label = T, 
                        no.legend = T,
                        do.return = T,
                        label.size = 10,
                        group.by = paste0("res.", c)
                    )
                    
                    p <- p + 
                        theme_bw() +
                        theme(
                            text = element_text(size = 10),
                            legend.position = "none",
                            plot.title = element_text(size = 30, face="bold", hjust = 0.5),
                            axis.title = element_text(size = 25),
                            axis.text = element_text(size = 20),
                            aspect.ratio = 1
                        ) +
                        labs(title = paste0(i, "(", j, ")"))
                    
                    ggsave(
                        filename = paste0(outdir, "umap_cluster_", c, ".pdf"),
                        plot = p,
                        width = 6,
                        height = 6
                    )
                    
                }
  
                p <- make_scatterpie_plot(
                    all_obj, 
                    color = "Stage", 
                    colors = stage_colors, 
                    reduction.use = "umap",
                    legend.position = "bottom",
                    guide_ncol = length(unique(all_obj@meta.data$Stage))
                ) + theme(aspect.ratio = 1)
                
                ggsave(
                    filename = paste0(outdir, "umap_stage.pdf"),
                    plot = p,
                    width = 6,
                    height = 7
                )
                
                
                perc <- calculate_cluster_perc(all_obj)
                
                p <- ggplot(perc, aes(x=factor(ident), y = freq, color = Stage)) + 
                    geom_point_rast() +
                    geom_line(aes(x=ident, y = freq, color = Stage)) +
                    scale_color_manual(values = c(
                        "I"="#65A9A3", 
                        "II"="#4A933E",
                        "III"="#EC7A21",
                        "IV"="#D73F47"
                    )) +
                    labs(x = "", y = "", color = "") +
                    theme_bw() +
                    theme(
                        axis.text = element_text(size = 15),
                        legend.text = element_text(size = 15),
                        legend.position = "bottom",
                        legend.background = element_blank()
                    )
                
                ggsave(
                    filename = paste0(outdir, "stage_perc.pdf"),
                    plot = p,
                    width = 4,
                    height = 4
                )
                
                
                if (!file.exists(paste0(outdir, "cluster_markers.csv"))) {
                    markers <- find_markers(all_obj, group.by = "ident", n.cores = 5, min.pct = 0.1)
                    write.csv(markers, paste0(outdir, "cluster_markers.csv"))
                } else {
                    markers <- read.csv(paste0(outdir, "cluster_markers.csv"), row.names = 1, stringsAsFactors = F)
                }

                temp <- markers %>%
                    filter(p_val_adj < 0.05) %>%
                    group_by(ident) %>%
                    top_n(10, wt=avg_logFC)


                p <- ModifiedDotPlot(all_obj, markers = temp, group.by = "ident") +
                    theme(axis.text = element_text(size = 15), panel.grid = element_blank())

                p

                ggsave(
                    filename = paste0(outdir, "cluster_dotplot.pdf"),
                    plot = p,
                    width = 6,
                    height = 10
                )

                # if (file.exists(paste0(outdir, "DPT.rds"))) {
                #     n.pcs = kee$KneeLocator(
                #         1:ncol(all_obj@dr$pca@cell.embeddings),
                #         apply(all_obj@dr$pca@cell.embeddings, 2, sd),
                #         curve='convex',
                #         direction='decreasing'
                #     )
                #     n.pcs = n.pcs$knee
                #     dpt <- readRDS(paste0(outdir, "DPT.rds"))
                # } else {
                #     
                #     n.pcs = kee$KneeLocator(
                #         1:ncol(all_obj@dr$pca@cell.embeddings),
                #         apply(all_obj@dr$pca@cell.embeddings, 2, sd),
                #         curve='convex',
                #         direction='decreasing'
                #     )
                #     n.pcs = n.pcs$knee
                #     
                #     dmp <- DiffusionMap(all_obj@dr$pca@cell.embeddings[, 1:n.pcs])
                #     dpt <- DPT(dmp)
                #     saveRDS(dpt, paste0(outdir, "DPT.rds"))
                # }
                # 
                # make_dc_plot(dpt, all_obj, paste0(outdir, "DPT"), title = paste0(i, "(", j, ")"), subtitle = paste0("n.pcs=", n.pcs))
            }, error = function(e) {
                print(e)
            })
        }
    }
}


## 2nd

modify = read.xlsx("11_CNV/each_cells/cluster_res.xlsx")


# modify$Cell
for (i in c("ATII", "Basal")) {
    res = modify[modify$Cell == i, "Res"]

    if (i == "NE") {
        res = 0.1
    }
    
    for (j in c("LUAD", "LUSC")) {
        print(paste(i, j))
        
        outdir = paste0(root.dir, "11_CNV/each_cells/", i, "/", j, "/")
        
        tryCatch({
            dir.create(outdir, showWarnings = F, recursive = T)
        
            all_obj <- readRDS(paste0(outdir, "seurat_obj.rds"))
            
            all_obj@meta.data$ident1 <- as.numeric(all_obj@meta.data[, paste0("res.", str_replace(as.character(res), "0+$", ""))]) + 1
            all_obj@meta.data$ident <- all_obj@meta.data$ident1
            
            # if (!file.exists(paste0(outdir, "DPT.rds")) || force) {
            #     n.pcs = 10
            #     dmp <- DiffusionMap(all_obj@dr$pca@cell.embeddings[, 1:n.pcs])
            #     dpt <- DPT(dmp)
            #     saveRDS(dpt, paste0(outdir, "DPT.rds"))
            # } else {
            #     dpt = readRDS(paste0(outdir, "DPT.rds"))
            # }
            
            p <- make_scatterpie_plot(
                all_obj,
                color = "Stage",
                colors = stage_colors,
                reduction.use = "umap",
                legend.position = "bottom",
                guide_ncol = length(unique(all_obj@meta.data$Stage))
            ) + theme(aspect.ratio = 1, panel.grid = element_blank())

            ggsave(
                filename = paste0(outdir, "umap_stage.pdf"),
                plot = p, width = 6, height = 7, device = cairo_pdf
            )


            perc <- calculate_cluster_perc(all_obj)

            p <- ggplot(perc, aes(x=factor(ident), y = freq, color = Stage)) +
              geom_point_rast() +
              geom_line(aes(x=ident, y = freq, color = Stage)) +
              scale_color_manual(values = c(
                  "I"="#65A9A3",
                  "II"="#4A933E",
                  "III"="#EC7A21",
                  "IV"="#D73F47"
              )) +
              labs(x = "", y = "", color = "") +
              theme_bw(base_family = "Arial Unicode MS") +
              theme(
                  panel.grid = element_blank(),
                  axis.text = element_text(size = 15),
                  legend.text = element_text(size = 15),
                  legend.position = "bottom",
                  legend.background = element_blank()
              )

            ggsave(
                filename = paste0(outdir, "stage_perc.pdf"),
                plot = p, width = 4, height = 4
            )
            

            if (!file.exists( paste0(outdir, "cluster_markers.csv")) ) {  # || res != 0.1
                markers <- find_markers(all_obj, group.by = "ident", n.cores = 5, min.pct = 0.1)
                write.csv(markers, paste0(outdir, "cluster_markers.csv"))
            } else {
                markers <- read.csv(paste0(outdir, "cluster_markers.csv"), row.names = 1, stringsAsFactors = F)
            }
            
            temp <- markers %>%
                filter(p_val_adj < 0.05) %>%
                group_by(ident) %>%
                top_n(10, wt=avg_logFC)


            p <- ModifiedDotPlot(all_obj, markers = temp, group.by = "ident") +
                theme(axis.text = element_text(size = 15), panel.grid = element_blank())

            ggsave(
                filename = paste0(outdir, "cluster_dotplot.pdf"),
                plot = p, width = 6, height = 10, device = cairo_pdf
            )
            
            logfc = 0.5
            if (i %in% c("Club", "Cilia")) {
                logfc = 0.25
            }

            h <- make_cluster_markers_heatmap(all_obj, markers, logfc)
            cairo_pdf(paste0(root.dir, "11_CNV/each_cells/", i, "/", j, "/", "cluster_markers.pdf"), width = 5, height = 6)
            draw(h)
            dev.off()
          
        }, error = function(e) {
            print(e)
        })
        
    }
} 



## Stage markers
for (i in modify$Cell) {
  res = modify[modify$Cell == i, "Res"]
  
  for (j in c("LUAD", "LUSC")) {
    print(paste(i, j))
    
    outdir = paste0(root.dir, "11_CNV/each_cells/", i, "/", j, "/")
    
    tryCatch({
      dir.create(outdir, showWarnings = F, recursive = T)
      
      all_obj <- readRDS(paste0(outdir, "seurat_obj.rds"))
      
      if (!file.exists( paste0(outdir, "stage_markers.csv")) ) {  # || res != 0.1
          markers <- find_markers(all_obj, group.by = "Stage", n.cores = 5, min.pct = 0.1)
          write.csv(markers, paste0(outdir, "stage_markers.csv"))
      } else {
          markers <- read.csv(paste0(outdir, "stage_markers.csv"), row.names = 1, stringsAsFactors = F)
      }
      
  
      # h <- make_cluster_markers_heatmap(all_obj, markers, logfc)
      # pdf(paste0(root.dir, "11_CNV/each_cells/", i, "/", j, "/", "cluster_markers.pdf"), width = 5, height = 6)
      # draw(h)
      # dev.off()
      
    }, error = function(e) {
      print(e)
    })
    
  }
} 



## EC LUSC
all_obj <- readRDS(full.path("11_CNV/each_cells/EC/LUSC/seurat_obj.rds"))


all_obj <- FindClusters(
    object = all_obj, reduction.type = "harmony", dims.use = 1:10, 
    resolution = c(0.04), 
    print.output = 0, save.SNN = FALSE, force.recalc = TRUE
)


for (c in c(0.03, 0.04, 0.05)) {
  p <- DimPlot(
      all_obj, 
      reduction.use = "umap", 
      pt.size = 0.1, 
      do.label = T, 
      no.legend = T,
      do.return = T,
      label.size = 10,
      group.by = paste0("res.", c)
  )
  
  p <- p + 
      theme_bw() +
      theme(
          text = element_text(size = 10),
          legend.position = "none",
          plot.title = element_text(size = 30, face="bold", hjust = 0.5),
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          aspect.ratio = 1
      ) +
      labs(title = paste0(i, "(", j, c, ")"))
    
    print(p)
}

all_obj@meta.data$ident = as.character(as.numeric(all_obj@meta.data$res.0.04) + 1)
markers <- find_markers(all_obj, group.by = "ident", n.cores = 5, min.pct = 0.1)


h <- make_cluster_markers_heatmap(all_obj, markers, 0.25)
# pdf(paste0(root.dir, "11_CNV/each_cells/", i, "/", j, "/", "cluster_markers.pdf"), width = 5, height = 6)
draw(h)
# dev.off()



all_obj <- readRDS(full.path("11_CNV/each_cells/Fib/LUAD/seurat_obj.rds"))
all_obj@meta.data$ident = as.character(as.numeric(all_obj@meta.data$res.0.1) + 1)
all_obj <- RunUMAP(all_obj, dims.use = 10, reduction.use = "harmony", n_neighbors = 30, seed.use = 42, reduction.name = "umap8")

DimPlot(all_obj, reduction.use = "umap", group.by = "res.0.1") +
  theme(aspect.ratio = 1)


make_scatterpie_plot(all_obj, "Stage", colors = stage_colors, reduction.use = "umap8", pt.size = 0.5)

DimPlot(all_obj, reduction.use = "tsne", group.by = "res.0.1") +
  theme(aspect.ratio = 1)



## 3rd
# 
# modify = c(
#     "Cilia"=c("LUSC"),
#     "EC"=c("LUSC")
# )
# 
# 
# for (i in names(modify)) {
#     
#     for (j in modify[[i]]) {
#         print(paste(i, j))
#         
#         outdir = paste0(root.dir, "03_each_cells/Batch/", i, "/", j, "/")
#         
#         tryCatch({
#             dir.create(outdir, showWarnings = F, recursive = T)
#             
#             all_obj <- readRDS(paste0(outdir, "seurat.rds"))
#             
#             all_obj@meta.data$ident1 <- as.numeric(all_obj@meta.data$res.0.3) + 1
#             all_obj@meta.data$ident <- all_obj@meta.data$ident1
#             
#             p <- DimPlot(
#                 all_obj, 
#                 reduction.use = "umap", 
#                 pt.size = 0.1, 
#                 do.label = T, 
#                 no.legend = T,
#                 do.return = T,
#                 label.size = 10,
#                 group.by = "ident1"
#             )
#             
#             p <- p + 
#                 theme_bw() +
#                 theme(
#                     text = element_text(size = 10),
#                     legend.position = "none",
#                     plot.title = element_text(size = 30, face="bold", hjust = 0.5),
#                     axis.title = element_text(size = 25),
#                     axis.text = element_text(size = 20),
#                     aspect.ratio = 1
#                 ) +
#                 labs(title = paste0(i, "(", j, ")"))
#             
#             ggsave(
#                 filename = paste0(outdir, "umap_cluster.pdf"),
#                 plot = p,
#                 width = 6,
#                 height = 6
#             )
#             
#             p <- make_scatterpie_plot(
#                 all_obj, 
#                 color = "Stage", 
#                 colors = stage_colors, 
#                 reduction.use = "umap",
#                 legend.position = "bottom",
#                 guide_ncol = length(unique(all_obj@meta.data$Stage))
#             ) + theme(aspect.ratio = 1)
#             
#             ggsave(
#                 filename = paste0(outdir, "umap_stage.pdf"),
#                 plot = p,
#                 width = 6,
#                 height = 7
#             )
#             
#             
#             perc <- calculate_cluster_perc(all_obj)
#             
#             p <- ggplot(perc, aes(x=factor(ident), y = freq, color = Stage)) + 
#                 geom_point() +
#                 geom_line(aes(x=ident, y = freq, color = Stage)) +
#                 scale_color_manual(values = c(
#                     "I"="#65A9A3", 
#                     "II"="#4A933E",
#                     "III"="#EC7A21",
#                     "IV"="#D73F47"
#                 )) +
#                 labs(x = "", y = "", color = "") +
#                 theme_bw() +
#                 theme(
#                     axis.text = element_text(size = 15),
#                     legend.text = element_text(size = 15),
#                     legend.position = "bottom",
#                     legend.background = element_blank()
#                 )
#             
#             ggsave(
#                 filename = paste0(outdir, "stage_perc.pdf"),
#                 plot = p,
#                 width = 4,
#                 height = 4
#             )
#             
#             
#             markers <- find_markers(all_obj, group.by = "ident", n.cores = 5, min.pct = 0.1)
#             write.csv(markers, paste0(outdir, "cluster_markers.csv"))
#             
#             temp <- markers %>%
#                 filter(p_val_adj < 0.05) %>%
#                 group_by(ident) %>%
#                 top_n(10, wt=avg_logFC)
#             
#             
#             p <- ModifiedDotPlot(all_obj, markers = temp, group.by = "ident") +
#                 theme(axis.text = element_text(size = 15))
#             
#             p
#             
#             ggsave(
#                 filename = paste0(outdir, "cluster_dotplot.pdf"),
#                 plot = p,
#                 width = 6,
#                 height = 10
#             )
#             
#         }, error = function(e) {
#             print(e)
#         })
#         
#     }
# } 
# 
# 
# 
# ## final round
# 
# resolution = c(
#     "ATII_LUAD"="res.0.1",
#     "Basal_LUSC"="res.0.03",
#     "EC_LUSC"="res.0.3",
#     "Cilia_LUSC"="res.0.3",
#     "Club_LUAD"="res.0.8"
# )
# 
# for(i in unique(meta$cell_short)) {
#     
#     for (j in c("LUAD", "LUSC")) {
#         outdir = paste0(root.dir, "03_each_cells/Batch/", i, "/", j, "/")
#         # print(outdir)
#         temp_meta <- meta[as.character(meta$cell_short) == i & as.character(meta$Disease) == j, ]
#         
#         if (nrow(temp_meta) > 100) {
#             
#             tryCatch({
#                 dir.create(outdir, showWarnings = F, recursive = T)
#                 # print("obj")
#                 if (file.exists(paste0(outdir, "seurat.rds"))) {
#                     all_obj <- readRDS(paste0(outdir, "seurat.rds"))
#                 } else {
#                     all_obj <- make_seurat_pipe(expr, temp_meta)
#                     
#                     saveRDS(all_obj, paste0(outdir, "seurat.rds"))
#                 }
#                 
#                 # print("res")
#                 if (paste(i, j, sep = "_") %in% names(resolution)) {
#                     res_col = resolution[[paste(i, j, sep = "_")]]
#                 } else {
#                     res_col = "res.0.1"
#                 }
#                 
#                 # print("ident")
#                 if (min(as.numeric(all_obj@meta.data[, res_col])) == 0) {
#                     all_obj@meta.data$ident1 <- as.numeric(all_obj@meta.data[, res_col]) + 1
#                     all_obj@meta.data$ident <- all_obj@meta.data$ident1
#                 } else {
#                     all_obj@meta.data$ident1 <- as.numeric(all_obj@meta.data[, res_col]) - min(as.numeric(all_obj@meta.data[, res_col])) + 1
#                     all_obj@meta.data$ident <- all_obj@meta.data$ident1
#                 }
#                 
#                 # print("uamp")
#                 p <- DimPlot(
#                     all_obj, 
#                     reduction.use = "umap", 
#                     pt.size = 0.1, 
#                     do.label = T, 
#                     no.legend = T,
#                     do.return = T,
#                     label.size = 10,
#                     group.by = "ident1"
#                 )
#                 
#                 p <- p + 
#                     theme_bw() +
#                     theme(
#                         text = element_text(size = 10),
#                         legend.position = "none",
#                         plot.title = element_text(size = 30, face="bold", hjust = 0.5),
#                         axis.title = element_text(size = 25),
#                         axis.text = element_text(size = 20),
#                         aspect.ratio = 1
#                     ) +
#                     labs(title = paste0(i, "(", j, ")"))
#                 
#                 ggsave(
#                     filename = paste0(outdir, "umap_cluster.pdf"),
#                     plot = p,
#                     width = 6,
#                     height = 6
#                 )
#                 
#                 # print("scatter")
#                 p <- make_scatterpie_plot(
#                     all_obj, 
#                     color = "Stage", 
#                     colors = stage_colors, 
#                     reduction.use = "umap",
#                     legend.position = "bottom",
#                     guide_ncol = length(unique(all_obj@meta.data$Stage))
#                 ) + theme(aspect.ratio = 1)
#                 
#                 ggsave(
#                     filename = paste0(outdir, "umap_stage.pdf"),
#                     plot = p,
#                     width = 6,
#                     height = 7
#                 )
#                 
#                 # print("perc")
#                 perc <- calculate_cluster_perc(all_obj)
#                 
#                 p <- ggplot(perc, aes(x=factor(ident), y = freq, color = Stage)) + 
#                     geom_point() +
#                     geom_line(aes(x=ident, y = freq, color = Stage)) +
#                     scale_color_manual(values = c(
#                         "I"="#65A9A3", 
#                         "II"="#4A933E",
#                         "III"="#EC7A21",
#                         "IV"="#D73F47"
#                     )) +
#                     labs(x = "", y = "", color = "") +
#                     theme_bw() +
#                     theme(
#                         axis.text = element_text(size = 15),
#                         legend.text = element_text(size = 15),
#                         legend.position = "bottom",
#                         legend.background = element_blank()
#                     )
#                 
#                 ggsave(
#                     filename = paste0(outdir, "stage_perc.pdf"),
#                     plot = p,
#                     width = 4,
#                     height = 4
#                 )
#                 
#                 # print("markers")
#                 if (!file.exists(paste0(outdir, "cluster_markers.csv"))) {
#                     markers <- find_markers(all_obj, group.by = "ident", n.cores = 5, min.pct = 0.1)
#                     write.csv(markers, paste0(outdir, "cluster_markers.csv"))
#                 } else {
#                     markers <- read.csv(paste0(outdir, "cluster_markers.csv"), row.names = 1, stringsAsFactors = F)
#                 }
#                 
#                 temp <- markers %>%
#                     filter(p_val_adj < 0.05) %>%
#                     group_by(ident) %>%
#                     top_n(10, wt=avg_logFC)
#                 
#                 
#                 p <- ModifiedDotPlot(all_obj, markers = temp, group.by = "ident") +
#                     theme(axis.text = element_text(size = 15))
#                 
#                 ggsave(
#                     filename = paste0(outdir, "cluster_dotplot.pdf"),
#                     plot = p,
#                     width = 6,
#                     height = 10
#                 )
#                 
#                 # print("dpt")
#                 if (file.exists(paste0(outdir, "DPT.rds"))) {
#                     n.pcs = kee$KneeLocator(
#                         1:ncol(all_obj@dr$pca@cell.embeddings),
#                         apply(all_obj@dr$pca@cell.embeddings, 2, sd),
#                         curve='convex',
#                         direction='decreasing'
#                     )
#                     n.pcs = n.pcs$knee
#                     dpt <- readRDS(paste0(outdir, "DPT.rds"))
#                 } else {
#                     
#                     n.pcs = kee$KneeLocator(
#                         1:ncol(all_obj@dr$pca@cell.embeddings),
#                         apply(all_obj@dr$pca@cell.embeddings, 2, sd),
#                         curve='convex',
#                         direction='decreasing'
#                     )
#                     n.pcs = n.pcs$knee
#                     
#                     dmp <- DiffusionMap(all_obj@dr$pca@cell.embeddings[, 1:n.pcs])
#                     dpt <- DPT(dmp)
#                     saveRDS(dpt, paste0(outdir, "DPT.rds"))
#                 }
#                 
#                 # print("dm")
#                 make_dc_plot(dpt, all_obj, paste0(outdir, "DPT"), title = paste0(i, "(", j, ")"), subtitle = paste0("n.pcs=", n.pcs))
#             }, error = function(e) {
#                 print(e)
#             })
#         }
#     }
# }


## Mfuzz 
library(Mfuzz)
library(Seurat)

do_mfuzz <- function(obj, markers) {
    expr = ExpressionSet(
        as.matrix(obj@scale.data[unique(markers$gene),])
    )
    
    m = mestimate(expr)

    tryCatch({
        temp = Dmin(expr, m=m, crange=seq(2,5,1), repeats=3, visu=FALSE)
        
        if (sum(!is.na(temp)) == 0) {
            c = length(unique(obj@meta.data$Stage))
        } else {
            new_temp = temp[!is.na(temp)]
            
            c = kee$KneeLocator(
                1:length(new_temp),
                new_temp, 
                curve='convex', 
                direction='decreasing'
            )
            c = c$knee
            
            c = which(temp == new_temp[c])
        }
        print(paste0("c=", c))
        print(paste0("m=", m))
        
        cl <- mfuzz(expr, c=c, m=m)
        
        if (length(cl$cluster) == 0) {
            res = markers[, c("ident", "gene")]
            colnames(res) <- c("Clt", "gene")
            return(res)
        }
        # extract results
        res = as.data.frame(cl$cluster)
        colnames(res) <- "Clt"
        res$gene = rownames(res)
        
        res = res[order(res$Clt), ]
        return(res)
    }, error = function(e) {
        res = markers[, c("ident", "gene")]
        colnames(res) <- c("Clt", "gene")
        return(res)
    })
}


make_stage_module_heatmap <- function(
    data,
    meta, 
    module, 
    col_fun = NULL, 
    labels_size=15, 
    title_size = 20, 
    cluster=F, 
    do.return=F, 
    show_row_names = F, 
    do.random = F,
    order.by.patient = T,
    mark_genes = NULL,
    cluster_row = T,
    row_split = NULL,
    column_title = "",
    do.balance = F
) {
    if (is.null(col_fun)) {
        col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
    } 
    
    meta$Disease <- as.character(meta$Disease)
    meta$Disease[meta$Disease == "LUAD"] <- "AD"
    meta$Disease[meta$Disease == "LUSC"] <- "SC"
    meta$Disease[meta$Disease == "LUAD_Normal"] <- "NL(AD)"
    meta$Disease[meta$Disease == "LUSC_Normal"] <- "NL(SC)"
    meta$Disease[meta$Disease == "Pleiomorphic Normal"] <- "NL(UPS)"
    meta$Disease[meta$Disease == "Pleiomorphic"] <- "UPS"
    meta$Disease[meta$Disease == "LULC"] <- "LC"
    meta$Disease[meta$Disease == "LULC_Normal"] <- "NL(LC)"
    
    patient_colors = cluster_cols[1:length(unique(meta$PatientID))]
    names(patient_colors) <- sort(as.character(unique(meta$PatientID)))
    
    set.seed(42)
    
    if (do.balance) {
        meta$Cells = rownames(meta)
        
        
        num_cells = sort(table(meta$Stage), decreasing = F)
        
        num_cells = ifelse(num_cells[1] < 20, num_cells[2], num_cells[1])
        
        meta <- meta %>% group_by(Stage) %>% sample_n(num_cells, replace = T) %>% unique() %>% as.data.frame()
        
        rownames(meta) <- meta$Cells
    }
    
    
    if(order.by.patient) {
        meta <- meta[order(meta$Disease, meta$Stage, meta$PatientID), c("Stage", "PatientID")] # "Disease", 
        
        anno_legend_param = list(
            Stage=list(
                direction = "horizontal", 
                nrow = 1,
                labels_gp = gpar(fontsize = labels_size, fontfamily = "Arial Unicode MS"),
                title_gp = gpar(fontsize = title_size, fontfamily = "Arial Unicode MS")
            ),
            PatientID=list(
                direction = "horizontal",
                nrow = 4,
                labels_gp = gpar(fontsize = labels_size, fontfamily = "Arial Unicode MS"),
                title_gp = gpar(fontsize = title_size, fontfamily = "Arial Unicode MS")
            )
            # Disease=list(
            #     direction = "horizontal",
            #     nrow = 1,
            #     labels_gp = gpar(fontsize = labels_size),
            #     title_gp = gpar(fontsize = title_size)
            # )
        )
    } else {
        if (do.random) {
            meta = meta[sample(1:nrow(meta)), ]
        }
        meta <- meta[order(meta$Stage), c("Stage"), drop = F] # meta$Disease, , "Disease"
        
        
        anno_legend_param = list(
            Stage=list(
                direction = "horizontal", 
                nrow = 1,
                labels_gp = gpar(fontsize = labels_size, fontfamily = "Arial Unicode MS"),
                title_gp = gpar(fontsize = title_size, fontfamily = "Arial Unicode MS")
            )
            # Disease=list(
            #     direction = "horizontal",
            #     nrow = 1,
            #     labels_gp = gpar(fontsize = labels_size),
            #     title_gp = gpar(fontsize = title_size)
            # )
        )
    }
    
    
    # top annotation
    ta = HeatmapAnnotation(
        df = meta,
        col = list(
            Stage =  colors[["Stage"]],
            PatientID = patient_colors
            # Disease = colors[["Disease"]]
        ),
        show_legend = T,
        show_annotation_name = T,
        annotation_legend_param = anno_legend_param
    )
    
    if(!is.null(mark_genes)) {
        la = rowAnnotation(
            foo = anno_mark(
                at = which(module %in% mark_genes),
                labels = module[module %in% mark_genes],
                which = "row",
                side = "left"
            ),
            show_legend = F,
            show_annotation_name = F
        )
    } else {
        la = NULL
    }
    
    
    h <- Heatmap(
        name = "Expr",
        data[as.character(module), rownames(meta)],
        cluster_rows = cluster_row,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = show_row_names,
        col = col_fun,
        top_annotation = ta,
        left_annotation = la, 
        column_split = meta$Stage, # paste(meta$Disease, meta$Stage),
        row_split = row_split,
        border = T,
        column_title = column_title,
        column_title_gp = gpar(fontsize=ifelse(column_title == "", 0, 20), fontfamily = "Arial Unicode MS"),
        row_title_gp = gpar(fontsize = 0, fontfamily = "Arial Unicode MS"),
        heatmap_legend_param = list(
            direction = "horizontal",
            labels_gp = gpar(fontsize = labels_size, fontfamily = "Arial Unicode MS"),
            title_gp = gpar(fontsize = title_size, fontfamily = "Arial Unicode MS")
        ),
        use_raster = T
    )
    
    if(do.return) {
        return(h)
    } else {
        draw(
            h, 
            column_title_gp = gpar(fontsize = 20, fontfamily = "Arial Unicode MS"),
            merge_legend = TRUE, 
            heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom"
        )
    }
}





library(ComplexHeatmap)
library(circlize)


colors = list(
    Stage = c(
        "I"="#65A9A3", 
        "II"="#4A933E", 
        "III"="#EC7A21", 
        "IV"="#D73F47", 
        "LUAD_Normal"="#FECC1B", 
        "LUSC_Normal"="#778793"
    ),
    Disease = c(
        "AD" = "#0084D1", 
        "NL(AD)"="#FECC1B", 
        "NL(AD"="#778793", 
        "Normal"="#73BDFF", 
        "NL(SC)"="#778793", 
        "SC"="#A0C807",
        "LC"="#91822B",
        "NL(LC)"="#EDCAB0",
        "UPS"="#EBAAA4",
        "NL(UPS)"="#B43018"
    ),
    Batch = c("1"="#E9B0B7", "2"="#90D0E0", "3"="#769982"),
    Cell = c(
        "ATII"="#FF0000",
        "Basal"="#F2AD00",
        "CD4"="#F98400",
        "CD8"="#5BBCD6",
        "Cilia"="#85D4E3",
        "Club"="#F4B5BD",
        "DC"="#9C964A",
        "EC"="#CDC08C",
        "Epi"="#FAD77B",
        "Fib"="#ECCBAE",
        "Mast"="#D69C4E",
        "Mo"="#ABDDDE",
        "NK"="#000000",
        "NE"="#F3DF6C",
        "Tregs"="#CEAB07",
        "Mφ"="#046C9A"
    )
)

cluster_cols = c(
    as.character(wes_palette("Zissou1")),
    as.character(wes_palette("Royal2")),
    as.character(wes_palette("Royal1")),
    as.character(wes_palette("Darjeeling2")),
    as.character(wes_palette("Darjeeling1")),
    wes_palette("BottleRocket2"),
    wes_palette("Rushmore1")
)



registerDoMC(5)
foreach(i = unique(meta$cell_short), .errorhandling = "pass") %dopar% {
    for (j in c("LUAD", "LUSC")) {
        outdir = paste0(root.dir, "11_CNV/each_cells/", i, "/", j, "/")
        
        print(outdir)
        
        if (!file.exists(paste0(outdir, "seurat_obj.rds"))) {
            next
        }
        print(outdir)
        all_obj <- readRDS(paste0(outdir, "seurat_obj.rds"))

        if (!file.exists(paste0(outdir, "stage_markers.csv"))) {
            markers <- find_markers(all_obj, group.by = "Stage", n.cores = 4, min.pct = 0.1)
            write.csv(markers, paste0(outdir, "stage_markers.csv"))
        } else {
            markers <- read.csv(paste0(outdir, "stage_markers.csv"), row.names = 1, stringsAsFactors = F)
        }

        if (!file.exists(paste0(outdir, "mfuzz.csv")) || nrow(read.csv(paste0(outdir, "mfuzz.csv"))) < 5) {
            temp_markers = markers[markers$avg_logFC > 0.5 & markers$p_val_adj < 0.05, ]
            
            if (nrow(temp_markers) < 25) {
                temp_markers = markers[markers$avg_logFC > 0.25 & markers$p_val_adj < 0.05, ]
            }
            
            res = do_mfuzz(all_obj, temp_markers)

            for(k in unique(res$Clt)) {
                temp = table(temp_markers$ident[temp_markers$gene %in% res$gene[res$Clt == k]])
                res$Clt[res$Clt == k] = paste0("M.", names(temp)[temp == max(temp)])
            }

            write.csv(res, paste0(outdir, "mfuzz.csv"))

        } else {
            res = read.csv(paste0(outdir, "mfuzz.csv"), row.names = 1, stringsAsFactors = F)
        }
    
        temp <- markers %>%
          filter(p_val_adj < 0.05 & avg_logFC > 0.25) %>%
          group_by(ident) %>%
          top_n(5, wt=avg_logFC) %>%
          as.data.frame()
        
        cairo_pdf(paste(outdir, "mfuzz.pdf", sep = "/"), width = 8, height = 6)
        make_stage_module_heatmap(
            as.matrix(all_obj@scale.data),
            all_obj@meta.data,
            res$gene,
            order.by.patient = F,
            do.random = T,
            col_fun = colorRamp2(c(-1, 0, 1), c("purple", "black", "yellow")),
            cluster_row = F,
            row_split = res$Clt,
            mark_genes = temp$gene,
            do.balance = T
        )
        dev.off()
    }
}



### Perform GO on mfuzz module
# files = list.files(paste(root.dir, "03_each_cells/Batch", sep = "/"), pattern = "mfuzz.csv", full.names = T, recursive = T)
# 
# registerDoMC(3)
# res = foreach(f = files, .combine = "rbind") %dopar% {
#     data <- read.csv(f, stringsAsFactors = F, row.names = 1)
#     cell <- basename(dirname(dirname(f)))
#     disease <- basename(dirname(f))
#     res = NULL
#     for (i in unique(data$Clt)) {
#         gene = data[data$Clt == i, "gene"]
#         temp = enrichGO(
#             gene, 
#             OrgDb = org.Hs.eg.db, 
#             keyType       = 'SYMBOL',
#             ont           = "BP",
#             pAdjustMethod = "BH",
#             pvalueCutoff  = 0.01,
#             qvalueCutoff  = 0.05,
#         )
#         
#         temp = as.data.frame(temp)
#         
#         if (nrow(temp) > 0) {
#             temp$ident = i
#             temp$disease = disease
#             temp$cell = cell
#             res = rbind(res, temp)
#         }
#     }
#     res
# }
# 
# res$GeneRatio_perc <- sapply(res$GeneRatio, function(x) {
#     x = str_split(x, "/")[[1]]
#     as.numeric(x[1]) / as.numeric(x[2])
# })
# 
# write.csv(res, paste(root.dir, "03_each_cells/mfuzz_go.csv", sep = "/"))

### GSEA

# c5 <- read.gmt("c5.all.v7.0.symbols.gmt")
# 
# 
# for (i in list.dirs("11_CNV/each_cells/", recursive = T)) {
#     
#     f = paste(i, "cluster_markers.csv", sep = "/")
#     
#     if (file.exists(f)) {
#         tryCatch({
#             markers = read.csv(f, row.names = 1, stringsAsFactors = F)
#             gsea = list()
#             for (i in unique(markers$ident)) {
#                 temp_markers <- markers[markers$ident == i, ]
#                 
#                 geneList = temp_markers$avg_logFC
#                 names(geneList) <- temp_markers$gene
#                 
#                 tryCatch({
#                     gsea[[paste0("C.", i)]] = GSEA(sort(geneList, decreasing = T), TERM2GENE=c5, verbose=FALSE)
#                 }, .error = function(e) {})
#             }
#             
#             saveRDS(gsea, paste(i, "gsea.rds"))
#             markers <- markers[order(markers$avg_logFC, decreasing = T), ]
#         }, .error = function(e) {})
#     }
# }
# 
# 
# list.files("11_CNV/each_cells/", pattern = "gsea.rds", recursive = T, full.names = T)




# ## atii 
# markers <- read.csv("03_each_cells/Batch/ATII/LUAD/cluster_markers.csv", row.names = 1, stringsAsFactors = F)
# markers <- markers[order(markers$avg_logFC, decreasing = T), ]
# 
# gsea = list()
# for (i in unique(markers$ident)) {
#     temp_markers <- markers[markers$ident == i, ]
#     
#     geneList = temp_markers$avg_logFC
#     names(geneList) <- temp_markers$gene
#     
#     gsea[[paste0("C.", i)]] = GSEA(sort(geneList, decreasing = T), TERM2GENE=c5, verbose=FALSE)
# }
# 
# 
# geneSet = c()
# for(i in names(gsea)) {
#     geneSet = c(geneSet, gsea[[i]]@result$Description)
# }
# 
# geneSet = as.data.frame(table(geneSet))
# # geneSet = names(geneSet)[geneSet == length(gsea)]
# 
# target = "GO_ACTIVATION_OF_IMMUNE_RESPONSE"
# 
# plist = list()
# for(i in names(gsea)) {
#     print(i)
# 
#     if (nrow(gsea[[i]]) == 0 ) {
#         next
#     }
#     
#     if (! target %in% gsea[[i]]@result$Description) {
#         next
#     }
#     
#     tryCatch({
#         plist[[i]] <- gseaplot2(
#             gsea[[i]], 
#             geneSetID = which(gsea[[i]]@result$Description == target),
#             pvalue_table = T,
#             title = str_replace_all(i, "^C.", "Cluster ")
#         )
#     }, error = function(e) {
#         print(e)
#     })
# }
# 
# 
# p <- cowplot::plot_grid(plotlist = plist, ncol = 1)
# ggsave(
#     filename = paste0("03_each_cells/Batch/ATII/LUAD/GSEA_clt.pdf"),
#     plot = p,
#     width = 6,
#     height = 8
# )
# 
# 
# 
# ## Basal
# markers <- read.csv("03_each_cells/Batch/Basal/LUSC/cluster_markers.csv", row.names = 1, stringsAsFactors = F)
# markers <- markers[order(markers$avg_logFC, decreasing = T), ]
# 
# gsea = list()
# for (i in unique(markers$ident)) {
#     temp_markers <- markers[markers$ident == i, ]
#     
#     geneList = temp_markers$avg_logFC
#     names(geneList) <- temp_markers$gene
#     
#     gsea[[paste0("C.", i)]] = GSEA(sort(geneList, decreasing = T), TERM2GENE=c5, verbose=FALSE)
# }
# 
# 
# 
# geneSet = c()
# for(i in names(gsea)) {
#     geneSet = c(geneSet, gsea[[i]]@result$Description)
# }
# 
# geneSet = as.data.frame(table(geneSet))
# # geneSet = names(geneSet)[geneSet == length(gsea)]
# 
# target = "GO_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE"
# 
# plist = list()
# for(i in names(gsea)) {
#     print(i)
#     
#     if (nrow(gsea[[i]]) == 0 ) {
#         next
#     }
#     
#     if (! target %in% gsea[[i]]@result$Description) {
#         next
#     }
#     
#     tryCatch({
#         plist[[i]] <- gseaplot2(
#             gsea[[i]], 
#             geneSetID = which(gsea[[i]]@result$Description == target),
#             pvalue_table = T,
#             title = str_replace_all(i, "^C.", "Cluster ")
#         )
#     }, error = function(e) {
#         print(e)
#     })
# }
# 
# 
# p <- cowplot::plot_grid(plotlist = plist, ncol = 1)
# ggsave(
#     filename = paste0("03_each_cells/Batch/Basal/LUSC/GSEA_clt.pdf"),
#     plot = p,
#     width = 6,
#     height = 12
# )


# 
# c5 <- read.gmt("c5.all.v7.0.symbols.gmt")
# 
# files = list.files("03_each_cells/Batch/", pattern = "cluster_markers.csv", full.names = T, recursive = T)
# 
# registerDoMC(10)
# res = foreach (i = files, .errorhandling = "pass") %dopar% {
#     print(i)
#     
#     markers = read.csv(i, row.names = 1, stringsAsFactors = F)
#     markers <- markers[order(markers$avg_logFC, decreasing = T), ]
#     
#     if(file.exists(paste0(dirname(i), "/gsea.rds"))) {
#         gsea = readRDS(paste0(dirname(i), "/gsea.rds"))
#     } else {
#         gsea = list()
#         for (j in unique(markers$ident)) {
#             temp_markers <- markers[markers$ident == j, ]
#             
#             geneList = temp_markers$avg_logFC
#             names(geneList) <- temp_markers$gene
#             
#             gsea[[paste0("C.", j)]] = GSEA(sort(geneList, decreasing = T), TERM2GENE=c5, verbose=FALSE)
#         }
#         
#         saveRDS(gsea, paste0(dirname(i), "/gsea.rds"))
#     }
# 
#     data = NULL
#     for (j in names(gsea)) {
#         temp = as.data.frame(gsea[[j]])
#         
#         if (nrow(temp) == 0) {
#             next
#         }
#         temp$ident = j
#         
#         data = rbind(data, temp)
#     }
#     
#     
#     data = data %>%
#         group_by(ident) %>%
#         top_n(10, wt = enrichmentScore)
#     
#     
#     p <- ggplot(
#         data,
#         aes(x = ident, y = Description, size = enrichmentScore, color = p.adjust)
#     ) + 
#         geom_point() +
#         scale_color_gradientn(colors=rev(wes_palette("Zissou1", 100, type = "continuous"))) +
#         labs(x = "", y="") +
#         theme_bw() +
#         theme(
#             axis.text.x = element_text(size = 15),
#             axis.text.y = element_text(size = 10)
#         )
#     
#     ggsave(
#         filename = paste0(dirname(i), "/gsea_plot.pdf"),
#         plot = p,
#         width = 12,
#         height = 6
#     )
# }



### Scatter plot of markers

## make scatter of markers

make_scatter_plot_of_markers <- function(
    data,
    addition_mark = NULL,
    p_val_cutoff = 0.05,
    logFC_cutoff = 0.25,
    group.by = "ident",
    mark.by = "gene",
    group_width = 0.5,
    group_gap = 0.1,
    text_size = 4,
    dot.size = 0.5,
    mark.size = 5,
    num_mark = 10,
    both_side_mark = FALSE,
    colors = NULL,
    seed = 42,
    x_title = "Cluster",
    y_title = "avgerage logFC",
    axis_text_size = 15,
    axis_title_size = 15,
    point_colors = c("red", "black"),
    legend_position = c(0.9, 0.1),
    legend_text_size = 15,
    legend_point_size = 3,
    aspect.ratio = 1
) {
    library(ggplot2)
    library(wesanderson)
    library(dplyr)
    library(ggrepel)
    set.seed(seed)
    
    data$ident = data[, group.by]
    data$mark = data[, mark.by]
    
    if (is.null(colors)) {
        colors = wes_palette("Zissou1", length(unique(data$ident)), type = "continuous")
    }
    
    rect_x_min = c()
    x_min = 0
    data$x = 0
    for (i in sort(unique(data$ident))) {
        rect_x_min = c(rect_x_min, x_min)
        
        data$x[data$ident == i] = x_min + runif(sum(data$ident == i), min = 0, max = group_width)
        
        x_min = x_min + group_width + group_gap
    }
    
    rect_data = data.frame(
        xmin = rect_x_min,
        xmax = rect_x_min + group_width,
        ident = sort(unique(data$ident)),
        color = as.character(colors)
    )
    
    rect_data$ymin = 0
    rect_data$ymax = 0
    
    for (i in unique(rect_data$ident)) {
        rect_data[rect_data$ident == i, "ymin"] = min(data$avg_logFC[data$ident == i]) 
        rect_data[rect_data$ident == i, "ymax"] = max(data$avg_logFC[data$ident == i]) 
    }
    
    if (both_side_mark) {
        mark = data %>%
            filter(p_val_adj < p_val_cutoff) %>%
            group_by(ident) %>%
            top_n(num_mark, wt = abs(avg_logFC)) %>%
            as.data.frame()
    } else {
        mark = data %>%
            filter(p_val_adj < p_val_cutoff) %>%
            group_by(ident) %>%
            top_n(num_mark, wt = avg_logFC) %>%
            as.data.frame()
    }
    
    if (!is.null(addition_mark)) {
        mark <- mark[!mark$mark %in% addition_mark, ]
    }

    data$col <- paste("adjust P-val <", p_val_cutoff)
    data$col[data$p_val_adj >= p_val_cutoff] <- paste("adjust P-val >=", p_val_cutoff)

    p <- ggplot() +
        geom_rect(   # background
            data = rect_data,
            aes_(
                xmin = ~xmin,
                xmax = ~xmax,
                ymin = ~ymin,
                ymax = ~ymax,
            ),
            fill = "grey75",
            alpha=.5,
            inherit.aes=FALSE
        ) +
        geom_rect(
            data = rect_data,
            aes_(
                xmin = ~xmin,
                xmax = ~xmax,
                ymin = -1 * logFC_cutoff,
                ymax = logFC_cutoff,
                fill = ~I(color)
            ),
            alpha=.9,
            inherit.aes=FALSE
        ) +
        geom_text(
            data=rect_data, 
            aes(x=xmin + (xmax - xmin)/2, y=0, label=ident), 
            size=text_size
        ) +
        geom_point(
            data = data[abs(data$avg_logFC) > logFC_cutoff, ],
            aes(x = x, y = avg_logFC, color = col),
            size = dot.size
        ) +
        geom_text_repel(
            data = mark,
            aes(x = x, y = avg_logFC, label = mark),
            size = mark.size,
            color = "black"
        ) +
        scale_color_manual(values = point_colors) +
        labs(x = x_title, y = y_title, color = "") +
        theme_classic() +
        theme(
            aspect.ratio = aspect.ratio,
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text = element_text(size = axis_text_size),
            axis.title = element_text(size = axis_title_size),
            legend.position = legend_position,
            legend.text = element_text(size = legend_text_size),
            legend.background = element_blank()
        ) +
        guides(color = guide_legend(override.aes = list(size = legend_point_size)))
    
    if (!is.null(addition_mark)) {
        p <- p + 
            geom_label_repel(
                data = data[data$gene %in% addition_mark, ],
                aes(x = x, y = avg_logFC, label = gene),
                size = mark.size,
                color = "red"
            )
    }
    p
}


## collect markers info
registerDoMC(10)
res = foreach(f = list.files("11_CNV/each_cells/", pattern = "cluster_markers.csv", recursive = T, full.names = T), .combine = "rbind") %dopar% {
    indir = dirname(f)
    cell = basename(dirname(indir))
    disease = basename(indir)
    
    res = NULL
    cluster_markers <- paste(indir, "cluster_markers.csv", sep = "/")
    
    if (file.exists(cluster_markers)) {
        temp <- read.csv(cluster_markers, stringsAsFactors = F, row.names = 1)
        
        if (nrow(temp) > 0 && ncol(temp) > 2) {
            temp$Cell = cell
            temp$Disease <- disease
            temp$Type <- "Cluster"
            res = rbind(res, temp)
        }
    }
    
    stage_markers <- paste(indir, "stage_markers.csv", sep = "/")
    if (file.exists(stage_markers)) {
        temp <- read.csv(stage_markers, stringsAsFactors = F, row.names = 1)
        
        if (nrow(temp) > 0 && ncol(temp) > 2) {
            temp$Cell = cell
            temp$Disease <- disease
            temp$Type <- "Stage"
            res = rbind(res, temp)
        }
    }
    
    res
}


res <- res[res$Cell != "each_cells", ]


# res[res$Cell == "ATII" & res$Type == "Cluster" & res$Disease == "LUAD", ]
for(i in unique(res$Cell)) {
    for (j in unique(res$Type)) {
        for (k in unique(res$Disease)) {
            data = res[res$Cell == i & res$Type == j & res$Disease == k, ]
            
            if (nrow(data) == 0) {
                next
            }
            
            p <- make_scatter_plot_of_markers(data, num_mark = 10, text_size = 8, mark.size = 3)
            
            filename = paste(j, "scatter.pdf", sep = "_")
            # print(filename)
            ggsave(
                filename = paste("LungCancer10x/11_CNV/each_cells", i, k, filename, sep = "/"),
                plot = p,
                width = 2 * length(unique(data$ident)),
                height = 6
            )
        }
    }
}

saveRDS(res, "11_CNV/each_cells/markers.rds")


for(i in unique(res$Cell)) {
  for (j in unique(res$Type)) {
    for (k in unique(res$Disease)) {
      data = res[res$Cell == i & res$Type == j & res$Disease == k, ]
      
      if (nrow(data) == 0) {
        next
      }
      
      ego <- enrichGO(gene          = data$gene,
                      keyType       = "SYMBOL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
      
      break
    }
    break
  }
  break
}



data = res[res$Cell == "ATII" & res$Type == "Cluster" & res$Disease == "LUAD", ]
p <- make_scatter_plot_of_markers(
    data, 
    addition_mark = c(
        "ALKBH3", "ALOX5AP", "BRAF", "CASP8", "CCDC130", 
        "EGFR", "EIF2S1", "ERP27", "G6PD", "INSR", 
        "KRAS", "MET", "PSMB4", "RBPJ", "ROS1", 
        "RPS6KA1", "SYT12", "TNFSF9"
    ),
    num_mark = 10, 
    text_size = 8, 
    mark.size = 3
)
p


filename = paste(j, "scatter.pdf", sep = "_")
print(filename)
ggsave(
    filename = paste("LungCancer10x/03_each_cells/Batch", i, k, filename, sep = "/"),
    plot = p,
    width = 8,
    height = 6
)


## make scatter plot of DEGs between LUAD with NLAD and LUSC with NLSC
library(openxlsx)
data = NULL
for (i in list.files("03_each_cells/Batch/", pattern = "normal_cancer", recursive = T, include.dirs = T, full.names = T)) {
    print(i)
    cell = basename(dirname(i))
    path = paste(i, "LUAD.xlsx", sep = "/")
    if (file.exists(path)) {
        temp = read.xlsx(path, rowNames = T)
        temp$Cell = cell
        temp$Disease = "LUAD"
        temp$gene = rownames(temp)
        temp <- temp[temp$avg_logFC < 0, ]
        temp$avg_logFC = abs(temp$avg_logFC)
        data = rbind(data, temp)
    }
    
    path = paste(i, "LUSC.xlsx", sep = "/")
    if (file.exists(path)) {
        temp = read.xlsx(path, rowNames = T)
        temp$Cell = cell
        temp$Disease = "LUSC"
        temp <- temp[temp$avg_logFC < 0, ]
        temp$gene = rownames(temp)
        temp$avg_logFC = -1 * abs(temp$avg_logFC)
        data = rbind(data, temp)
    }
}


p <- make_scatter_plot_of_markers(
    data, 
    num_mark = 10, 
    text_size = 4,
    mark.size = 3, 
    group.by = "Cell", 
    aspect.ratio = 0.5, 
    dot.size = 0,
    x_title = "Cells",
    both_side_mark = T
)

p <- p + 
    scale_y_continuous(breaks = seq(-5, 5, 2.5), labels = abs(seq(-5, 5, 2.5))) +
    annotate("text", x = 0.3, y = -5, label = "LUSC", size = 6) +
    annotate("text", x = 0.3, y = 5, label = "LUAD", size = 6)

ggsave(
    filename = "03_each_cells/Batch/DEGs_all_cells.pdf",
    plot = p,
    width = 10,
    height = 5
)


## GSEA on LUAD & adjacent   LUSC & adjacent

c5 <- read.gmt("h.all.v7.0.symbols.gmt")

files = list.files("03_each_cells/Batch/", pattern = "normal_cancer", full.names = T, recursive = T, include.dirs = T)

registerDoMC(10)
res = foreach (i = files, .errorhandling = "pass") %dopar% {
    print(i)
    
    res = list()
    for (j in c("LUAD", "LUSC")) {
        path = paste(i, "/", j, ".xlsx", sep = "")
        
        if (file.exists(path)) {
            markers = read.xlsx(path, rowNames = T)
            
            geneList = markers$avg_logFC
            names(geneList) <- row.names(markers)
            
            res[[j]] = GSEA(sort(geneList, decreasing = T), TERM2GENE=c5, verbose=FALSE)
        }
    }
    
    saveRDS(res, paste0(i, "/gsea.rds"))
}

## make plots

files = list.files("03_each_cells/Batch/", pattern = "gsea.rds", full.names = T, recursive = T)
files = files[str_detect(files, "normal_cancer")]

gsea <- readRDS(files[1])

as.data.frame(gsea[["LUAD"]])