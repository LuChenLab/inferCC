## TCGA suvival
## make survival analysis to all genes


load_pacakges <- function() {
    library(TCGAbiolinks)
    library(SingleCellExperiment)
    library(doMC)
    library(survival)
    library(survminer)
    library(progress)
    library(stringr)
    library(ddpcr)
    library(doMC)
}

suppressPackageStartupMessages(load_pacakges())


luad = readRDS("LUAD.rds")

lusc = readRDS("LUSC.rds")


make_surviveplot <- function(data, n.cores = 20) {
    expression = t(assay(data, "raw_count"))
    gene_names = as.character(sapply(colnames(expression), function(x){strsplit(x, "\\|", perl=T)[[1]][1]}))
    colnames(expression) <- gene_names
    
    expression = merge(as.data.frame(colData(data)), expression, by = 0)
    
    print("calc high low")
    for (i in colnames(expression)) {
        if (i %in% gene_names) {
            expression[, i] <- ifelse(expression[, i] > median(expression[, i]),'high','low')
        }
    }
    
    pb <- progress_bar$new(total = length(gene_names))
    
    registerDoMC(n.cores)
    foreach(i=gene_names) %dopar% {
         pb$tick()
        if (length(unique(expression[, i])) < 2) {
            return()
        }
        
        surv <- TCGAanalyze_survival(
            expression, 
            clusterCol = i,
            risk.table = F,
            conf.int = F,
            surv.median.line = "hv",
            filename = NULL
        )
        
        temp = surv$plot$layers[[6]]$data
        
        max_x = max(surv$plot$data$time) * 0.8
        p <- surv$plot + geom_text(aes(x = max_x, y = 0.9), label = i, size = 7)
        
        ggplot2::ggsave(filename = paste0(i, ".pdf"), width = 10, height = 5, plot = p, dpi = 600)
        
        col_names <- c("pvalue",  as.character(temp[, "type"]), "gene")
        
        res <- c(
            as.numeric(gsub("p = ", "", surv$plot$layers[[4]]$aes_params$label)),
            as.numeric(temp[, "x1"]),
            i
        )
        
        names(res) <- col_names
        write.table(as.data.frame(res), file = paste0(i, ".txt"), row.names = T, col.names = F, quote = F, sep = "\t")
    }
        
}

root = getwd()
dir.create("Survival_plots/LUAD", recursive = T)
setwd("Survival_plots/LUAD")

make_surviveplot(luad)


dir.create("Survival_plots/LUSC", recursive = T)
setwd("Survival_plots/LUSC")

make_surviveplot(lusc)



