library(ggplot2)
library(Seurat)
library(openxlsx)
library(fmsb)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(reshape2)
library("wesanderson")
library(ggpubr)


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


### extract temp information
stage_colors = c(
    "Normal"="#73BDFF", 
    "Benign"="#3E6596", 
    "CPI"="#C5000B", 
    "I"="#EC7A21", 
    "II"="#D73F47", 
    "III"="#65A9A3", 
    "IV"="#4A933E",
    "ADC" = "#063373", 
    "LCC"="#FECC1B", 
    "SCC"="#A0C807", 
    "LELC"="#6A0019", 
    "LBT"="#0084D1",
    "Other"="#DA6906"
    )
disease_colors = c(
    "ADC" = "#063373", 
    "LCC"="#FECC1B", 
    "Normal"="#73BDFF", 
    "Other"="#DA6906", 
    "SCC"="#A0C807", 
    "LELC"="#6A0019", 
    "CPI"="#C5000B", 
    "LBT"="#0084D1"
    )
tissue_colors = c(
    "Bronichi"="#FB290F", 
    "Lymph node"="#488F17", 
    "Tumor"="#FC8210"
)
patient_colors = gg_color_hue(50)


make_complex_heatmap <- function(obj, genes.use, order.by = "Stage") {

    # if (!is.null(module_id)) {
    #     module = module[module$ident == module_id, ]
    # }

    # rownames(module) <- module$gene
    # temp_module <- module[, 1, drop=FALSE]
    # temp_module[,1] <- as.character(temp_module[, 1])

    # print(colnames(obj@meta.data))
    meta <- obj@meta.data[, c("Stage", "PatientID")]
    colnames(meta) <- c("Stage", "Patient")
    # print(head(meta))
    # print(order.by)
    meta <- meta[order(meta[, order.by]), ]
    
    mat = as.matrix(obj@scale.data)[genes.use, rownames(meta), drop=F]
    # mat = MinMax(mat, -2.5, 2.5)

    cols_stage = stage_colors[unique(meta$Stage)]
    # print(cols_stage)

    cols_patients = patient_colors[as.numeric(gsub("P", "",unique(meta$Patient)))]
    names(cols_patients) <- unique(meta$Patient)
    # print(cols_patients)

    # cols_cluster = wes_palette("Darjeeling2", length(unique(meta$Cluster)), type = "continuous")
    # names(cols_cluster) <- unique(meta$Cluster)
    # print(cols_cluster)

    ann_colors = list(
        Stage = cols_stage,
        Patient = cols_patients
        # Cluster = cols_cluster
    )

    # print("???")
    ha = HeatmapAnnotation(
        df = meta,
        col = ann_colors
    )

    # print("???")
    ht = Heatmap(
        mat,
        name = "mat",
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        top_annotation = ha,
        # column_title = module_id
    )
    draw(ht)
}



# function to make heatmaps
# :param obj Seurat obj
# :param cluster results from perform_mfuzz
# :param output_prefix: the output file prefix
make_heatmap <- function(obj, genes.use, output_prefix=NULL, group.by = "Stage") {

    height = length(genes.use) / 6
    
    if (height > 45) {
        height = 45
    } else if (height < 5) {
        height = 5
    }
    
    if (is.null(output_prefix)) {
        make_complex_heatmap(obj, genes.use, i) # group.by = group.by)
    } else {
        # ggsave(
        #     paste(output_prefix, i, ".png", sep = ""),
        #     p,
        #     width = 12,
        #     height = height,
        #     dpi = 300,
        #     units = "in"
        # )

        png(output_prefix, width = 12, height = height, res = 600, units = "in")
        make_complex_heatmap(obj, genes.use, order.by = "Stage") # , group.by = group.by)
        dev.off()
    }
    
}



# function to make dotplot
# :param obj Seurat obj
# :param cluster results from perform_mfuzz
# :param output_prefix: the output file prefix
make_dotplot <- function(obj, genes.use, output_prefix=NULL, group.by="Stage", title=NULL) {

    # print(temp_genes)
    p <- DotPlot(
        obj, 
        genes.use,
        group.by = group.by, 
        x.lab.rot = TRUE, 
        do.return = T
    ) + 
    coord_flip()

    if(!is.null(title)) {
        p = p + labs(title = title)
    }
    
    height = length(genes.use) / 6
    
    if (height > 40) {
        height = 40
    } else if (height < 5) {
        height = 5
    }
    
    if (is.null(output_prefix)) {
        print(p)
    } else {
        ggsave(
            output_prefix,
            p,
            width = 12,
            height = height,
            dpi = 300,
            units = "in"
        )
    }
    
}


    
# function to format data for line plots
# :param obj Seurat obj
# :param genes.use is a vector of gene symbol
# :param group.by: which column in object@meta.data
format_data <- function(object, genes.use, group.by) {
    # expr = MinMax(as.matrix(obj@scale.data)[genes.use, , drop=FALSE], -2.5, 2.5)
    expr = as.matrix(obj@scale.data)[genes.use, , drop=FALSE]
    expr = melt(as.matrix(expr))
    colnames(expr) <- c("gene", "cell", "value")

    # temp_mat = melt(as.matrix(object@raw.data[genes.use, , drop = FALSE]))

    # print(head(temp_mat))

    # expr = expr[temp_mat[, 3] != 0, ]
    
    meta = object@meta.data[, group.by, drop = FALSE]
    meta$cell = rownames(meta)
    
    # print(expr)
    # print(meta)
    expr = merge(expr, meta, by = "cell")
    
    expr = as.data.frame(
        expr %>% group_by(cell, Stage) %>% 
            dplyr::summarize(mean = mean(value, na.rm=TRUE))) %>% 
        dplyr::select(cell, Stage, mean) %>% 
        unique()
    
    return(expr)
}


# function to make module line plots
# :param object: Seurat obj
# :param genes.use which genes to use
# :param group.by: column of object@meta.data
make_line_plot_single_module <- function(object, genes.use, stage, group.by = "Stage", title=NULL, output=NULL) {
    expr = format_data(object, genes.use)
    
    p <- ggplot(data = expr, aes(x=cell, y=mean, group = 1)) +
        geom_line() +
        facet_grid(.~Stage, scales = "free_x", space = "free") + 
        theme(axis.text.x = element_blank(), legend.position = "none")

    if (!is.null(title)) {
        p = p + labs(title = title)
    }

    if(is.null(output)) {
        print(p)
    } else {
        ggsave(
            output, plot = p,
            width = 12, height = 6,
            dpi = 600, units = "in"
        )
    }

    # make violin plots to compare
    colnames(expr)[colnames(expr) == "mean"] = "mean_val"
    temp = expr %>% group_by(Stage) %>% mutate(mean_by_stage=mean(mean_val)) %>% dplyr::select(Stage, mean_by_stage) %>% unique()

    # max_stage = temp[temp$mean_by_stage == max(temp$mean_by_stage), "Stage"]

    # print(temp)
    # print(max_stage)
    # my_compare = list()

    # for (i in temp$Stage) {
    #     if (i != stage) {
    #         my_compare[[length(my_compare) + 1]] = c(as.character(i), stage)
    #     }
    # }
    # print(my_compare)
    expr$Stage = factor(expr$Stage, levels = names(stage_colors))

    p <- ggviolin(expr, x = "Stage", y = "mean_val", fill = "Stage",
        palette = stage_colors,
        add = "boxplot", add.params = list(fill = "white"),
        outline=FALSE
    )
        # stat_compare_means(comparisons = my_compare) + # Add significance levels
        # stat_compare_means(label.y = 6 * max(expr$mean_val))

    return(p)
}


# function to make module line plots
# :param object: Seurat obj
# :param genes.use which genes to use
# :param group.by: column of object@meta.data
make_line_plot_multi_module <- function(object, module, stage, group.by = "Stage", title=NULL, output=NULL) {

    expr = NULL
    for(i in unique(module$module_id)) {
        temp = format_data(object, module[module$module_id == i, "gene"])
        temp$module_id = i
        expr = rbind(expr, temp)
    }

    # expr = format_data(object, genes.use)
    
    p <- ggplot(data = expr, aes(x=cell, y=mean, group = 1)) +
        geom_line() +
        facet_grid(module_id~Stage, scales = "free_x", space = "free") + 
        theme(axis.text.x = element_blank(), legend.position = "none")

    if (!is.null(title)) {
        p = p + labs(title = title)
    }

    if(is.null(output)) {
        print(p)
    } else {
        ggsave(
            output, plot = p,
            width = 12, height = 6 * length(unique(module$module_id)),
            dpi = 600, units = "in"
        )
    }

    # make violin plots to compare
    colnames(expr)[colnames(expr) == "mean"] = "mean_val"
    temp = expr %>% group_by(Stage, module_id) %>% mutate(mean_by_stage=mean(mean_val)) %>% dplyr::select(Stage, mean_by_stage, module_id) %>% unique()

    # max_stage = temp[temp$mean_by_stage == max(temp$mean_by_stage), "Stage"]

    # print(temp)
    # print(max_stage)
    # my_compare = list()

    # for (i in temp$Stage) {
    #     if (i != stage) {
    #         my_compare[[length(my_compare) + 1]] = c(as.character(i), stage)
    #     }
    # }

    expr$Stage = factor(expr$Stage, levels = names(stage_colors))

    # print(my_compare)
    p <- ggviolin(expr, x = "Stage", y = "mean_val", fill = "Stage",
        palette = stage_colors,
        add = "boxplot", add.params = list(fill = "white"))+
        # stat_compare_means(comparisons = my_compare) + # Add significance levels
        # stat_compare_means(label.y = 6 * max(expr$mean_val)) +
        facet_grid(module_id~., scales="free_y", space="free_y")

    return(p)
}



# uniform api for makt line plot
# :param object: Seurat obj
# :param genes.use which genes to use
# :param group.by: column of object@meta.data
make_line_plot <- function(object, stage, module=NULL, group.by = "Stage", title=NULL, output=NULL) {
    if (length(unique(module$module_id)) == 1) {
        return (make_line_plot_single_module(
            object=object,
            stage=stage,
            genes.use = module$gene,
            group.by = group.by,
            title = title,
            output = output
        ))
    } else if (!is.null(module)) {
        return (make_line_plot_multi_module(
            object=object,
            stage=stage,
            module = module,
            group.by = group.by,
            title = title,
            output = output
        ))
    }
}



make_radar_plot <- function(object, genes.use, output, module) {
    # expr = MinMax(as.matrix(object@scale.data[genes.use, , drop=FALSE]), -2.5, 2.5)
    # print(expr[1:5, 1:5])
    expr = as.matrix(object@scale.data[genes.use, , drop=FALSE])
    expr = melt(as.matrix(expr))
    colnames(expr) <- c("gene", "cell", "value")

    # print(head(expr))

    expr <- merge(expr, module, by = "gene")
    expr$Stage <- obj@meta.data[expr$cell, "Stage"]
    expr$Patients <- obj@meta.data[expr$cell, "PatientID"]

    expr$module_id <- paste(paste0(expr$module_id), expr$Stage, sep = " -> ")

    expr = as.data.frame(expr %>% dplyr::select(module_id, value, Patients) %>% group_by(module_id, Patients) %>% mutate(mean = mean(value)) %>% dplyr::select(module_id, mean, Patients) %>% unique()) 

    expr = dcast(expr, Patients~module_id, value.var = "mean")
    expr[is.na(expr)] <- 0
    # print(head(expr))

    rownames(expr) <- expr$Patients
    expr = expr[, !colnames(expr) %in% c("Patients"), drop = FALSE]

    # print(head(expr))

    expr <- rbind(
        rep(ceiling(max(expr)), ncol(expr)),
        rep(floor(min(expr)), ncol(expr)),
        expr
    )

    # print(head(expr))

    png(output, width = 12, height = 6, res = 600, units = "in")
    par(
        mfrow=c(1, 2), 
        bty='n', 
        oma = c(0.5,0.5,0,0) + 0.1,
        mar = c(1,0,0,1) + 0.1
    )

    legend_labels = rownames(expr)
    legend_labels = legend_labels[3: length(legend_labels)]
    colors_border = wes_palette("Zissou1", length(legend_labels), type = "continuous")

    radarchart(expr, axistype=1, 
        #custom polygon
        pcol=colors_border, # pfcol=colors_in, plwd=4 , plty=1, 
        #custom the grid
        cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
        caxislabels=seq(expr[2, 1], expr[1, 1], (expr[1, 1] - expr[2, 1]) / 5), 
        #custom labels
        vlcex=0.8 
    )
    plot(0, col="white", cex=0, axes=F, ann=FALSE)

    legend(
        "left",
        legend = legend_labels, 
        bty = "n", pch=20 , 
        col=colors_border , 
        text.col = "grey", 
        cex=1.2, 
        pt.cex=3,
        ncol = 2
    )
    dev.off()
}

# get command line parameters
args = commandArgs(trailingOnly = T)

results = read.xlsx(args[1], sheet=2)
cell_name = args[2]
rds = args[3]
output_dir = paste(args[4], cell_name, sep = "/")
# output_dir = args[4]
cell_name_short = args[5]



dir.create(output_dir, showWarnings = F, recursive = T)
obj <- readRDS(rds)

results = results[results$Cell_name == cell_name, ]


format_genes <- function(results) {
    res = NULL
    for(i in 1:nrow(results)) {
        genes = str_split(results[i, "Genes"], "\\|")[[1]]
        module_id = paste(cell_name_short, results[i, "Stage"], results[i, "Mfuzz_ID"], sep = "_")

        res = rbind(res, data.frame(gene = genes, module_id = module_id))
    }

    return(res)
}


run <- function() {

    disease = str_split(cell_name, "_")[[1]]
    disease = disease[length(disease)]

    patients = obj@meta.data[obj@meta.data$Disease == disease, ]
    patients = unique(patients$PatientID)

    cells = rownames(obj@meta.data[obj@meta.data$PatientID %in% patients, ])

    temp_obj <- CreateSeuratObject(
        obj@raw.data[, cells, drop = F],
        meta = obj@meta.data[cells, , drop = F]
    )

    temp_obj@scale.data = obj@scale.data[, cells, drop = F]

    print(head(temp_obj@meta.data))

    for(i in unique(results$Stage)) {

        print(i)

        temp_results = results[results$Stage == i, , drop=F]

        # print(dim(temp_results))

        # if(length(genes) < 3) {
        #     next
        # }

        module_id = paste(cell_name_short, i, sep = "_") 
        print(module_id)

        module = format_genes(temp_results)

        # print(head(module))

        # output_dir_temp = paste(output_dir, paste(i, results[i, "Msig.ID"], results[i, "Stage"], results[i, "Mfuzz_ID"], sep = "_"), sep = "/")
        # print(output_dir_temp)

        dir.create(output_dir, showWarnings=F, recursive=T)

        print("heatmap")
        # make_heatmap(
        #     temp_obj, 
        #     genes.use = genes, 
        #     output_prefix = paste(
        #         output_dir, 
        #         paste(module_id, "heatmap.png", sep = "_"), 
        #         sep = "/"
        #     ),
        #     group.by="Stage"
        # )

        # p <- DoHeatmap(obj, genes.use = genes, group.by = "Stage", do.plot = F)

        # print("dotplot")
        # make_dotplot(
        #     temp_obj, 
        #     genes.use = genes, 
        #     output_prefix = paste(
        #         output_dir, 
        #         paste(module_id, "dotplot.png", sep = "_"), 
        #         sep = "/"
        #     ),
        #     group.by="Stage"
        # ) 

        print("line")
        p <- make_line_plot(
            temp_obj, 
            module=module, 
            stage=i,
            group.by = "Stage", 
            title=NULL, 
            output=paste(
                output_dir, 
                paste(module_id, "line.png", sep = "_"), 
                sep = "/"
            )
        )

        # save the violin plot
        ggsave(
            paste(
                output_dir, 
                paste(module_id, "violin.png", sep = "_"), 
                sep = "/"
            ),
            width = 6, height = 6 * length(unique(module$module_id)),
            dpi = 600, units = "in"
        ) 

        # print("radar")
        # stages = unique(obj@meta.data$Stage)

        # for (stage in stages) {
        #     temp_patients = unique(temp_obj@meta.data[temp_obj@meta.data$Stage == stage, "PatientID"]) 

        #     temp_meta = unique(temp_obj@meta.data[temp_obj@meta.data$PatientID %in% temp_patients, ])

        #     temp_obj1 = CreateSeuratObject(
        #         temp_obj@data[, rownames(temp_meta), drop=F],
        #         meta = temp_meta
        #     )
        #     temp_obj1@scale.data = temp_obj@scale.data[, rownames(temp_meta), drop=F]

        #     make_radar_plot(
        #         temp_obj1, 
        #         genes.use=genes, 
        #         output=paste(
        #             output_dir, 
        #             paste(module_id, "_radar_", stage, ".png", sep = ""), 
        #             sep = "/"
        #         ), 
        #         module=module
        #     )
        # }
        # make_radar_plot(
        #     temp_obj, 
        #     genes.use=genes, 
        #     output=paste(
        #         output_dir, 
        #         paste(module_id, "_radar.png", sep = ""), 
        #         sep = "/"
        #     ), 
        #     module=module
        # ) 
    }

}

run()
