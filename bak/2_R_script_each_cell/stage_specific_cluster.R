
load_pacakges <- function() {
    library(here)
    library(dplyr)
    library(TCGAbiolinks)
    library(SingleCellExperiment)
    library(openxlsx)
    library(ggpubr)
    library(ggplot2)
    library(doMC)
    library(stringr)
    library(RColorBrewer)
    library(reshape2)
    library(cowplot)
    library(pheatmap)
    library(gridExtra)
}

suppressPackageStartupMessages(load_pacakges())



# function to extract specific markers table
# 只保留高表达的部分
# :param markers: DataFrame of seurat FindMarkers output
# :param stage_specific_cluster: named list, stage -> c(cluster id)
# :param p_val_adj: threshold of adjust p value
# :param logfc: threshold of logFC
extract_stage_specific_markers <- function(
    markers, 
    stage_specific_clusters,
    p_val_adj = 0.05,
    logfc = 0
) {
    
    markers <- markers[markers$avg_logFC > logfc & markers$p_val_adj < p_val_adj, ]
    markers$Stage <- "Not specific"
    
    
    for (i in names(stage_specific_clusters)) {
        markers[markers$ident %in% stage_specific_clusters[[i]], "Stage"] <- i
    }
    
    markers$stage_spec <- paste(markers$ident, markers$Stage, sep = " - ")
    return(markers)
}


# function to calculate the KM-survival info based on TCGA data
# :param markers: selected markers
# :param data: TCGA biolinks Query data
# :param clinical: the clinical data corresponding to data
# :param group.by: only used to kept specific columns
# :rest param: used by TCGAanalyze_SurvivalKM
calculate_km_survival <- function(
    markers, 
    data, 
    clinical,
    group.by = "stage_spec",
    ThreshTop=0.67, 
    ThreshDown=0.33, 
    Survresult=F,
    p.cut = 1
) {
    # get expression value
    expression_raw <- assay(data, "raw_count")
    expression_log <- log2(expression_raw)
    
    # prepare markers gene by TCGA rownames
    genes = as.character(sapply(rownames(expression_log), function(x){strsplit(x, "\\|", perl=T)[[1]][1]}))
    genes <- genes %in% markers$gene
    
    temp_genes = rownames(expression_log)[genes]
    
    tabSurvKM <- TCGAanalyze_SurvivalKM(
        clinical,
        expression_log,
        Genelist = temp_genes,
        Survresult = Survresult,
        ThreshTop=ThreshTop,
        ThreshDown=ThreshDown,
        p.cut = p.cut
    )
    
    tabSurvKM$gene = as.character(sapply(rownames(tabSurvKM), function(x){strsplit(x, "\\|", perl=T)[[1]][1]}))
    tabSurvKM$entrezid = as.character(sapply(rownames(tabSurvKM), function(x){strsplit(x, "\\|", perl=T)[[1]][2]}))
    
    tabSurvKM <- merge(markers[, unique(c("Stage", "ident", "gene", group.by))], tabSurvKM, by = "gene")
    
    return(tabSurvKM)
}






# set colors
stage_colors = c("Benign"="#3E6596", "I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E")
disease_colors = c("ADC" = "#063373", "LCC"="#FECC1B", "Normal"="#73BDFF", "Other"="#DA6906", "SCC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B", "LBT"="#0084D1")
tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")

colors = c(stage_colors, disease_colors, tissue_colors)






# function to make pvalue denstiy plots, based on the KM SA results
# :param markers the resutls of km_survival
make_density_plots <- function(markers) {
    markers$value <- -log10(markers$pvalue)
    markers$ident <- as.character(markers$ident)

    p <- ggdensity(markers, x = "value",
              add = "mean", rug = TRUE,
              color = "ident",
              alpha = 0.1
              )
    
    p <- p + theme(legend.title = element_blank()) +
        labs(x = "-log10(pvalue)")
    
    p <- facet(p, facet.by = "Stage", ncol=1)
    
    return(p)
}


# function to sort list by multiple index
# :param elements: list of c(1, 2)
# :return vector: keys of list
sort_list_by_multi_elements_int <- function(elements) {
    temp = matrix(NA, nrow = 1, ncol=2)
    
   
    for (i in elements) {
        temp <- rbind(temp, sort(i))
    }
    
    temp = na.omit(temp)
    rownames(temp) <- as.character(1:nrow(temp))
    temp = as.data.frame(temp)
    temp$dist = temp[,2] - temp[,1]
    
    res = c()
    for (i in sort(unique(temp[,1]))) {
        temp_data = temp[temp[,1] == i, ]
        res = c(res, rownames(temp_data)[order(temp_data[, 3])])
    }

    res = as.numeric(res)
    return(res)
}




# function to make pvalue violin plots
# :param markers: results of km_survival
# :param colors: color palette
# :param color.by: which column used as color and x, and comparisions
# :param method: used by ggviolin
# :param title: main title
# :param dist: the distance between extra p value text and main vilolin. the times of max(-log10(pvalue))
make_violin_plots <- function(
    markers,
    colors, 
    color.by = "stage_spec",
    method = "wilcox.test", 
    title = "", 
    dist = 2
    ) {
    
    different_levels <- calculate_levels(markers)
    markers$stage_spec <- factor(markers$stage_spec, levels = different_levels)
    
    # construct comparision
    my_comparisons = list()
    stages = as.character(unique(markers[, color.by]))

    for (i in c("Not specific", "Benign", "Normal", "I", "II", "III")) {
        # find the stages that pattern mathced
        matched_stages = stages[str_detect(stages, pattern = paste0("(\\s+)?", i, "(\\s+)?"))]
 
        if (length(matched_stages) > 0) {
            
            # iter over matched stages and make comparison between it and all stages that not matched
            if (length(my_comparisons) == 0) {
                for (j in matched_stages) {
                    for (k in stages) {
                        if (!k %in% matched_stages) {
                            my_comparisons[[length(my_comparisons) + 1]] = c(j, k)
                        } 
                    }
                } 
            }
            
            if (length(matched_stages) >= 2) {
                for (j in 2:length(matched_stages)) {
                    my_comparisons[[length(my_comparisons) + 1]] <- c(matched_stages[j - 1], matched_stages[j])
                }
            }
        }
    }
    
    my_comparisons <- unique(my_comparisons)
    
    my_comparisons_idx = list()
    for(i in my_comparisons) {
        my_comparisons_idx[[length(my_comparisons_idx) + 1]] = c(which(different_levels == i[1]), which(different_levels == i[2]))
    }

    my_comparisons_idx = sort_list_by_multi_elements_int(my_comparisons_idx)
    my_comparisons_new = list()
    for (i in my_comparisons_idx) {
        my_comparisons_new[[length(my_comparisons_new) + 1]] = my_comparisons[[i]]
    }
    
    my_comparisons = my_comparisons_new
    
    markers$pvalue_log10 <- -log10(markers$pvalue)
    
    p1 <- ggviolin(
        markers, 
        x = color.by, 
        y = "pvalue_log10", 
        fill = "Stage",
        add = "boxplot", 
        add.params = list(fill = "white"),
        palette = colors,
        title = "p value"
    ) +
        # add p value into plot
        # label -> p.format (show formatted p value), instead of p.signif (show **)
        # paired -> wilcox
        stat_compare_means(
            comparisons = my_comparisons, 
            label = "p.format", 
            method = method
        ) + # Add significance levels
        stat_compare_means(
            label.y = max(markers$pvalue_log10) * dist,
            label.x = 1.5
        ) +
        theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(
            title = title,
            y = "-log10(pvalue)",
            x = ""
        )
    

    p2 <- ggviolin(
        markers, 
        x = color.by, 
        y = "high", 
        fill = "Stage",
        add = "boxplot", 
        add.params = list(fill = "white"),
        palette = colors,
        title = "p value"
    ) +
        # add p value into plot
        # label -> p.format (show formatted p value), instead of p.signif (show **)
        # paired -> wilcox
        stat_compare_means(
            comparisons = my_comparisons, 
            label = "p.format", 
            method = method
        ) + # Add significance levels
        stat_compare_means(
            label.y = max(markers$pvalue_log10) * dist,
            label.x = 1.5
        ) +
        theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(
            title = title,
            y = "Median survival days of high expression group",
            x = ""
        ) + 
        yscale("log10", .format = TRUE)
    
    # p3 <- ggviolin(
    #     markers, 
    #     x = color.by, 
    #     y = "low", 
    #     fill = "Stage",
    #     add = "boxplot", 
    #     add.params = list(fill = "white"),
    #     palette = colors,
    #     title = "p value"
    # ) +
    #     # add p value into plot
    #     # label -> p.format (show formatted p value), instead of p.signif (show **)
    #     # paired -> wilcox
    #     stat_compare_means(
    #         comparisons = my_comparisons, 
    #         label = "p.format", 
    #         method = method
    #     ) + # Add significance levels
    #     stat_compare_means(
    #         label.y = max(markers$pvalue_log10) * dist,
    #         label.x = 1.5
    #     ) +
    #     theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    #     labs(
    #         title = title,
    #         y = "Median survival days of low expression group",
    #         x = ""
    #     ) + 
    #     yscale("log10", .format = TRUE)
    
    p <- cowplot::plot_grid(p1, p2, nrow = 1)
    return(p)
}






# function to make percentage heatmap of stage componets each cluster
# https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps
# :param obj seurat object
# :param target_dir: output dirctory
make_heatmap <- function(obj, target_dir) {
    
    breaks = seq(1, 100, 0.1)
    colorMap <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks))
    
    stages = sort(unique(obj@meta.data$Stage))
    cluster = sort(unique(obj@meta.data$res.0.6))
    
    # set the height for heatmap
    height = length(cluster)
    if (height < 5) {
        height = 5
    }
    
    plots = list()
    
    # count the stage percentage of each cluster 
    temp <- obj@meta.data %>% 
        dplyr::select(res.0.6, Stage) %>% 
        dplyr::group_by(res.0.6, Stage) %>% 
        dplyr::add_tally() %>%
        unique() %>%
        dplyr::group_by(res.0.6) %>%
        dplyr::mutate(freq = n / sum(n) * 100)
    
    temp <- dcast(temp, res.0.6~Stage, value.var = "freq")
    rownames(temp) <- temp[,1]
    temp <- temp[, 2:ncol(temp)]

    tryCatch({
        temp = pheatmap(
            temp, 
            cluster_rows = F, 
            display_numbers = T,
            color = colorMap,
            breaks = breaks,
            number_format = "%.1f",
            main = "Cluster",
            silent = T
        )
        
        plots[[1]] = temp[[4]]
    }, error = function(e) {

        temp = pheatmap(
            temp, 
            cluster_rows = F, 
            cluster_cols = F,
            display_numbers = T,
            color = colorMap,
            breaks = breaks,
            number_format = "%.1f",
            main = "Cluster",
            silent = T
        )
        
        plots[[1]] = temp[[4]]
    }
    )
    
    
    # count the clsuter percentage of each stage
    temp <- obj@meta.data %>% 
        dplyr::select(res.0.6, Stage) %>% 
        dplyr::group_by(Stage, res.0.6) %>% 
        dplyr::add_tally() %>%
        unique() %>%
        dplyr::group_by(Stage) %>%
        dplyr::mutate(freq = n / sum(n) * 100)
    
    temp <- dcast(temp, res.0.6~Stage, value.var = "freq")
    rownames(temp) <- temp[,1]
    temp <- temp[, 2:ncol(temp)]
    
    tryCatch({
        temp = pheatmap(
            temp, 
            cluster_rows = F, 
            display_numbers = T,
            color = colorMap,
            breaks = breaks,
            number_format = "%.1f",
            main = "Stage",
            silent = T
        )
        
        plots[[2]] = temp[[4]]
    }, error = function(e) {
        temp = pheatmap(
            temp, 
            cluster_rows = F, 
            cluster_cols = F,
            display_numbers = T,
            color = colorMap,
            breaks = breaks,
            number_format = "%.1f",
            main = "Stage",
            silent = T
        )
        
        plots[[2]] = temp[[4]]
    }
    )
    
    
    p <- grid.arrange(arrangeGrob(grobs = plots, nrow = 1))
    
    ggsave(
        filename = paste(target_dir, "heatmap.png", sep = "/"),
        plot = p,
        width = length(stages) * 2 + 1,
        height = height,
        dpi = 600,
        units = "in"
        )
}






# function to calculate the levels for furthor plotting
# :param markers: the results of km_survival
# :param group.by: which column to use
# :param cluster_first: the specific id must be `cluster - stage` or `stage - cluster`, I try to sort clsuter by int
calculate_levels <- function(markers, group.by = "stage_spec", cluster_first = TRUE) {
    idents <- as.character(unique(markers[, group.by]))
    
    stages = c("Not specific", "Normal", "CPI", "Benign", "I", "II", "III", "IV")
    
    order_by_stages = c()
    for(i in stages) {
        # use regexp to detect the qulified str
        pattern = paste0("\\s+", i, "$")
        
        if (!cluster_first) {
            pattern = paste0("^", i, "\\s+")
        }
        
        temp = idents[str_detect(idents, pattern = pattern)]
        
        # this part is kind of complex
        # first -> split the `cluster - stage` or `stage - cluster` by ` - `
        # second -> convert the `cluster` into numberic
        # third -> order these by the cluster
        order_by_stages = c(order_by_stages, temp[order(sapply(temp, function(x) {
            temp = strsplit(x, split = " - ", perl = TRUE)[[1]]
            if(cluster_first) {
                return(as.numeric(temp[1]))
            } else {
                return(as.numeric(temp[2]))
            }
        }))])
    }
    
    return(order_by_stages)
}





# 3. 做KM SA (survival analysis)
# :param rds: path to rds, seurat obj
# :param tcga_dir: the tcga_data from TCGAbiolinks, dump into rds
# :param label: LUAD or LUSC
# :param cluster_info: results from seurat findMarkers
# :param stage_specific_perc: the threshold
# :param marker_p: p_val threshold
# :param marker_logfc: avg_logFC threshold
make_stage_specific_km_sa <- function(
    rds, 
    tcga_dir, 
    label, 
    cluster_info = "annotation_results_by_cluster.xlsx",
    stage_specific_perc = 60,
    marker_p = 0.05,
    marker_logfc = 0,
    skip.heatmap = F
    ) {
    
    setwd(dirname(rds))

    if (!file.exists(cluster_info)) {
        return()
    }
    
    target_dir = "stage_specific"
    
    if (!dir.exists(target_dir)) {
        dir.create(target_dir)
    }
    
    # 计算每个cluster中的细胞，不同时期的构成比例
    obj <- readRDS(rds)
    
    if (!skip.heatmap) {
        make_heatmap(obj, target_dir = target_dir)
    }
    
    print("temp")
    temp <- obj@meta.data %>%
        dplyr::select(res.0.6, Stage) %>%
        dplyr::group_by(res.0.6, Stage) %>%
        dplyr::add_tally() %>%
        unique() %>%
        dplyr::group_by(res.0.6) %>%
        dplyr::mutate(freq = n / sum(n) * 100) %>% 
        dplyr::group_by(res.0.6) %>% 
        dplyr::top_n(1) %>%
        unique()

    print(temp)
    temp <- as.data.frame(temp[temp$freq > stage_specific_perc, ])

    print(temp)

    if (nrow(temp) == 0) {
        return()
    }

    stage_specific_clusters = list()

    for (i in unique(as.character(temp[, "Stage"]))) {
        stage_specific_clusters[[i]] <- c()
    }

    for (i in 1:nrow(temp)) {
        stage_specific_clusters[[temp[i, "Stage"]]] <- c(temp[i, "res.0.6"], stage_specific_clusters[[temp[i, "Stage"]]])
    }

    # print(stage_specific_perc)
    # 提取所需cluster的marker基因
    markers <- read.xlsx(cluster_info, rowNames = T)
    stage_specific_markers <- extract_stage_specific_markers(
        markers,
        stage_specific_clusters,
        p_val_adj = marker_p,
        logfc = marker_logfc
    )

    # 读取所需的TCGA数据
    data <- readRDS(paste(tcga_dir, paste0(label, ".rds"), sep = "/"))
    clinical <- readRDS(paste(tcga_dir, paste0(label, "_clinical.rds"), sep = "/"))
    
    surv <- read.csv(paste(tcga_dir, paste0(label, "_survival.csv"), sep = "/"))
    rownames(surv) <- surv$gene

    print("calculate KM")
    # KM SA
    # res <- calculate_km_survival(stage_specific_markers, data, clinical)
    res <- merge(stage_specific_markers, surv, by = "gene")
    
    if (nrow(res) == 0) {
        return()
    }
    
    # save results

    width = ceiling(length(unique(res$stage_spec)) / 1.5)
    
    if (width < 5) {
        width = 5
    }
    
    p <- make_violin_plots(res, colors = colors, dist = 5)
    ggsave(
        paste(target_dir, paste0("violin_tcga_km_pvalue_", label, "_stage_specific.png"), sep = "/"), 
        width = width * 2,
        height = 12, 
        dpi = 600, 
        units = "in"
        )

    write.csv(res, file = "stage_specific_km_res.csv")
}



