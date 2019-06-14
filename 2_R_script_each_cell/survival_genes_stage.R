# Created at 20190506
# scripts to calculate the KM-survival genes based on Stage specific genes
# 

args = commandArgs(trailingOnly=TRUE)


library(TCGAbiolinks)
library(openxlsx)
library(SingleCellExperiment)
library(ggplot2)
library(doMC)
library(stringr)
library(ggpubr)



# function to calculate the KM-survival info based on TCGA data
# :param markers: selected markers
# :param data: TCGA biolinks Query data
calculate_km_survival <- function(
    markers, 
    data, 
    clinical,
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

    tabSurvKM <- TCGAanalyze_SurvivalKM(clinical,
                                        expression_log,
                                        Genelist = temp_genes,
                                        Survresult = Survresult,
                                        ThreshTop=ThreshTop,
                                        ThreshDown=ThreshDown,
                                        p.cut = p.cut
                                        )
    
    
    tabSurvKM$gene = as.character(sapply(rownames(tabSurvKM), function(x){strsplit(x, "\\|", perl=T)[[1]][1]}))
    tabSurvKM$entrezid = as.character(sapply(rownames(tabSurvKM), function(x){strsplit(x, "\\|", perl=T)[[1]][2]}))
    
    return(tabSurvKM)
}


# function to perform the KM and make plots
# :param input_dir:path to  directory of each cells
# :param label: LUAD or LUSC
# :param tcga_dir: path to tcga_dir
# :param method: t.test, wilcox.test, anova, kruskal.test
# :param dist: distance between p value label and violin, using dist times the max of pvalue_log10
pipeline_km <- function(
    input_dir, 
    label, 
    tcga_dir,
    method="wilcox.test",
    dist=2
    ) {
    if (!label %in% c("LUAD", "LUSC")) {
        stop("label should be LUAD or LUSC")
    }
    
    rds = paste0(label, ".rds")
    
    setwd(input_dir)
    
    input_xlsx ="annotation_results_by_stage.xlsx"
    
    if (!file.exists(input_xlsx)) {
        
    } else {

        data = readRDS(paste(tcga_dir, paste0(label, ".rds"), sep = "/"))
        
        clinical <- paste(tcga_dir, paste0(label, "_clinical.rds"), sep = "/")
        clinical <- readRDS(clinical)
        
        markers <- read.xlsx(input_xlsx, rowNames = T)
        markers <- markers[markers$p_val_adj < 0.05,]
        
        tabSurv <- calculate_km_survival(markers, data, clinical = clinical)
        
        # set factors manully
        markers <- merge(markers, tabSurv, by = "gene")
        markers$ident = factor(markers$ident, levels = c("Normal", "CPI", "Benign", "I", "II", "III", "IV"))
        
        # set colors
        stage_colors = c("Benign"="#3E6596", "I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E")
        disease_colors = c("ADC" = "#063373", "LCC"="#FECC1B", "Normal"="#73BDFF", "Other"="#DA6906", "SCC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B", "LBT"="#0084D1")
        tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")

        colors = c(stage_colors, disease_colors, tissue_colors)
        
        # construct comparision
        my_comparisons = list()
        stages = as.character(unique(markers$ident))
        if ("Benign" %in% stages) {
            for (i in stages) {
                if (i != "Benign") {
                    my_comparisons[[length(my_comparisons) + 1]] = c("Benign", i)
                }
            }
        } else if ("Normal" %in% stages) {
            for (i in stages) {
                if (i != "Normal") {
                    my_comparisons[[length(my_comparisons) + 1]] = c("Normal", i)
                }
            }
        } else {
            for (i in stages) {
                if (i != "I") {
                    my_comparisons[[length(my_comparisons) + 1]] = c("I", i)
                }
            }
        }
        
        p <- ggplot(markers, aes(x = -log10(pvalue), color = ident)) + 
            geom_density() +
            theme(legend.title = element_blank()) +
            scale_color_manual(values = colors)
        ggsave(paste0("density_tcga_km_pvalue_", label, ".png"), width = 6, height = 3, dpi = 600, units = "in")
        
        markers$pvalue_log10 = -log10(markers$pvalue)
        print(my_comparisons)
        p <- ggviolin(
                markers, 
                x = "ident", 
                y = "pvalue_log10", 
                fill = "ident",
                add = "boxplot", 
                add.params = list(fill = "white"),
                palette = colors
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
                label.y = max(markers$pvalue_log10) * 2
            ) +
            theme(legend.title = element_blank()) +
            labs(
                title = gsub(pattern = "_", replacement = " ", x = basename(input_dir)),
                y = "-log10(pvalue)",
                x = ""
            )
        
        # p <- ggplot(markers, aes(x = ident, y = -log10(pvalue), fill = ident)) +
        #     geom_violin() +
        #     theme(legend.title = element_blank())
        ggsave(paste0("boxplot_tcga_km_pvalue_", label, ".png"), width = 6, height = 4, dpi = 600, units = "in")
        
        write.xlsx(markers, paste0("TCGA_km_stage_", label, ".xlsx"))
    }
}


# input_dir = args[1]
# tcga_dir = args[2]
input_dir = getwd()
tcga_dir = "/mnt/raid62/Lung_cancer_10x/TCGA"


## do all SCC
registerDoMC(1)
foreach(i = list.dirs(path = input_dir, full.names = T, recursive = F)) %dopar% {
    print(i)
    
    if (str_detect(basename(i), pattern = "SCC$")) {
        pipeline_km(input_dir = i, label = "LUSC", tcga_dir = tcga_dir)
    } else if (str_detect(basename(i), pattern = "ADC$")) {
        pipeline_km(input_dir = i, label = "LUAD", tcga_dir = tcga_dir)
    } else {
        pipeline_km(input_dir = i, label = "LUSC", tcga_dir = tcga_dir)
        pipeline_km(input_dir = i, label = "LUAD", tcga_dir = tcga_dir)
    }

    # for (j in list.files(i, pattern = "*.rds", full.names = T)) {
    #     if (!basename(j) %in% c("monocle3.rds", "slingshot.rds")) {
    #         print(j)
    #     }
    # }
}
setwd(input_dir)

### calculate the p value between different stages
# function to calculate the p values between different stage of pvalue_log10
calculate_p_value <- function(xlsx, group.by = "ident", target = "pvalue_log10", method = "wilcox.test") {
    meta = read.xlsx(xlsx)
    
    idents = unique(meta[ ,group.by])
    
    # find comparision root
    root = "I"
    if ("Benign" %in% idents) {
        root = "Benign"
    } else if ("Normal" %in% idents) {
        root = "Normal"
    }
    
    pvalues = c()
    pvalue_1 = meta[meta[, group.by] == root, target]
    for (i in idents) {
        if (i != root) {
            pvalue_2 = meta[meta[, group.by] == i, target]
            
            pvalue = eval(parse(text = paste0(method, "(pvalue_1, pvalue_2)$p.value")))
            
            if (is.na(pvalue) || is.null(pvalue)) {pvalue = 0}
            
            pvalues = c(pvalues, pvalue)
        }
    }
    return(pvalues)
}



# function to make violin plots between multiple xlsx files
make_violin_plot <- function(
    files, 
    filename, 
    dist=3, 
    method_pvalue="wilcox.test", 
    method_plot="wilcox.test"
) {
    res = NULL
    for (i in files) {
        print(i)
        
        if(file.exists(i)) {
            label = gsub(pattern = "_(AD|SC)C", x=basename(dirname(i)), replacement = "", perl=T)
            
            temp = data.frame(cell=label, pvalue=calculate_p_value(i, method = method_pvalue))
            
            if (is.null(res)) {
                res = temp
            } else {
                res = rbind(res, temp)
            }
        }
    }
    
    my_comparision = list()
    
    res = res[order(res[, 2]), ]

    root = as.character(res[1, 1])

    for(i in unique(as.character(res[, 1]))) {
        if (i != root) {
            my_comparision[[length(my_comparision) + 1]] <- c(root, i)
        }
    }

    res[, 1] = factor(
        as.character(res[, 1]), 
        levels = c(
            "Alveolar_II", "Basal", "Ciliated", 
            "Club", "Endothelial", "Fibroblasts", 
            "Goblet", "Neuroendocrine",
            
            "B_cells", "CD4", "Dendritic",
            "Erythroid_precursor",
            "Exhaust_T", "Granulocyte", "Monocytes",
            "NK", "T_cells", "Treg"
            )
        )
    
    
    p <- ggviolin(
        res, 
        x = "cell", 
        y = "pvalue", 
        fill = "cell",
        add = "boxplot", 
        add.params = list(fill = "white")
    ) +
        theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(
            title = "ADC",
            y = "",
            x = ""
        ) +
        # add p value into plot
        # label -> p.format (show formatted p value), instead of p.signif (show **)
        # paired -> wilcox
        stat_compare_means(
            comparisons = my_comparision,
            label = "p.format",
            method = method_plot
        ) + # Add significance levels
        stat_compare_means(
            label.y = max(res$pvalue) * dist
        )
    
    ggsave(filename = filename, plot = p, width = 12, height = 6, dpi = 600, units = "in")
    
    return(res)
}

input_dir = "/Volumes/WD/CloudStation/lung_cancer_10X/讨论/2019.05.10/"
single_disease = list.files(path = input_dir, pattern = "*(_ADC|_SCC|_Tumor|_tumor)", full.names = TRUE)

targets = c()
for (i in list.dirs(path = input_dir, full.names = T)) {
    if (!i %in% single_disease) {
        targets <- c(targets, i)
    }
}

res = make_violin_plot(
    files = sapply(targets, function(x){return(paste(x, "TCGA_km_stage_LUSC.xlsx", sep = "/"))}),
    filename = paste(input_dir, "LUSC_pvalue.png", sep = "/"),
    dist = 3
)


res = make_violin_plot(
    files = sapply(targets, function(x){return(paste(x, "TCGA_km_stage_LUAD.xlsx", sep = "/"))}),
    filename = paste(input_dir, "LUAD_pvalue.png", sep = "/"),
    dist = 4
)


make_violin_plot(
    files = sapply(
        list.files(path = input_dir, pattern = "*_ADC", full.names = T), 
        function(x){return(paste(x, "TCGA_km_stage_LUAD.xlsx", sep = "/"))}
        ),
    filename = paste(input_dir, "LUAD_pvalue_without_normal.png", sep = "/"),
    dist = 3
)


make_violin_plot(
    files = sapply(
        list.files(path = input_dir, pattern = "*_SCC", full.names = T), 
        function(x){return(paste(x, "TCGA_km_stage_LUSC.xlsx", sep = "/"))}
        ),
    filename = paste(input_dir, "LUSC_pvalue_without_normal.png", sep = "/"),
    dist = 3
)
 