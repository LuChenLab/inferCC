library(pheatmap)
library(stringr)
library(dplyr)
library(reshape2)
library(ggpubr)
library(cowplot)
library(openxlsx)
library(Seurat)
library("wesanderson")


EXHAUST = c(
    "CD274",
    "PDCD1",
    "CTLA4",
    "TIGIT",
    "LAG3",
    "NBL1",
    "CXCL13"
)


STAGE2NUM = c(
    "I"=1,
    "II"=2,
    "III"=3,
    "IV"=4
)


TUMORSIZE2NUM = c(
    "T1"=1,
    "T2"=2,
    "T3"=3,
    "T4"=4
)

N2NUM = c(
    "N0"=0,
    "N1"=1,
    "N2"=2
)


M2NUM = c(
    "M0"=0,
    "M1"=1,
    "M2"=2
)



# Function to make correlation between different genes and meta info
# calculate mean expression by group
# :param object: Seurat2 object
# :param group.by: which meta info
# :param genes.use: genes to test
# :param scale: whether to use scale.data
# :param type: cor.method, see ggscatter. pearson, spearman, kendall
# :param converter: convert character meta info to numeric. If meta info already numeric, leave if null
make_heatmap <- function(
    object, 
    group.by, 
    genes.use, 
    scale=FALSE, 
    type="spearman",
    converter=NULL,
    mean=T,
    set_x_breaks=T,
    label.x = NULL,
    title = ""
) {
    meta = object@meta.data
    
    group.order = sort(unique(meta[, group.by]))
    
    y = sapply(meta[, group.by], function(x) {
        which(group.order == x)
    })
    
    expr = if(scale) object@scale.data else object@raw.data
    
    expr = expr[intersect(genes.use, rownames(expr)), rownames(meta)]
    
    expr = melt(as.matrix(expr))
    expr$group = meta[expr$Var2, group.by]

    expr = expr[!is.na(expr$group),]
        
    if (mean) {
        expr = as.data.frame(
            expr %>% 
                group_by(Var1, group) %>% 
                mutate(m = mean(value)) %>% 
                dplyr::select(Var1, group, m) %>% 
                unique()
            )
        
        colnames(expr)[3] = "value"
    }

    if (!is.null(converter)) {
        expr$group = sapply(as.character(expr$group), function(x) {
            converter[[x]]
        })
    }
    
    expr$group = as.numeric(expr$group)
    
    plist = list()
    
    
    # if(is.null(label.x)) {
    #     label.x = median(temp_expr$group)
    # }
    
    # for(i in sort(unique(expr$Var1))) {
    #     temp_expr = expr[expr$Var1 == i, ]
    #     temp_group = sort(unique(temp_expr$group))
    #     
    #     if(!is.null(converter)) {
    #         temp_x_labels = sapply(sort(unique(temp_expr$group)), function(x) {
    #             names(converter)[converter == x]
    #         })
    #     } else {
    #         temp_x_labels = temp_group
    #     }
    # 
    #     p = ggscatter(
    #         temp_expr,
    #         x = "group",
    #         y = "value",
    #         add = "reg.line",
    #         add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    #         conf.int = TRUE,
    #         cor.coef = TRUE, 
    #         cor.coeff.args = list(method = type, label.x = label.x, label.sep = "\n"),
    #         xlab = "", 
    #         ylab = "",
    #         title = paste(i, title_append)
    #     )
    #     
    #     if(set_x_breaks) {
    #         p = p + scale_x_continuous(
    #             breaks=temp_group,
    #             labels=temp_x_labels
    #         )
    #     }
    #     
    #     plist[[length(plist) + 1]] = p
    # }
    # 
    # plot_grid(plotlist = plist)
    
    colnames(expr)[colnames(expr) == "Var1"] = "Gene"
    
    ggscatter(
        expr,
        x = "group",
        y = "value",
        add = "reg.line",
        palette = "jco",
        xlab = "",
        ylab = "",
        color = "Gene",
        # facet.by = "Gene",
        title = title
    ) +
        stat_cor(aes(color = Gene), method = "spearman") 
}



#####################################################################################
# Try percentage of MetaCell on Correlation

read_lfp <- function(path, meta) {

    if(!file.exists(paste(path, paste0("mc.metacell_mc.Rda"), sep = "/"))) {
        return (NULL)
    }

    load(paste(path, paste0("mc.metacell_mc.Rda"), sep = "/"))
    mc = object
    
    clt = as.data.frame(object@mc)

    load(paste(path, paste0("gstat.metacell.Rda"), sep = "/"))
    gstat = object
    
    
    fp_max = apply(mc@mc_fp, 1, max)
    fp_tot = gstat[intersect(rownames(mc@mc_fp), rownames(gstat)), "tot"]
    
    
    f = fp_max > 0 & fp_tot > 0
    
    lfp = log2(mc@mc_fp[f,])
    
    
    cols = c("N", "Stage", "M", "T")

    colnames(clt) = "Clt"
    clt$cell = rownames(clt)

    res = NULL
    for(i in cols) {
        clt[i] = meta[clt$cell, i]
        
        temp = eval(parse(text = paste0("data.table::dcast(clt, Clt~", i,")")))
        
        rownames(temp) = temp$Clt
        temp = temp[, colnames(temp) != "Clt"]
        
        temp = temp / rowSums(temp)
        
        if(is.null(res)) {
            res = temp
        } else {
            res = cbind(res, temp)
        }
    }
    
    obj <- CreateSeuratObject(
        lfp,
        meta.data = res
    )
    
    obj@scale.data = obj@raw.data
    
    return(obj)
}


args = commandArgs(trailingOnly = T)


# obj <- readRDS("03_each_cells/LUSC/CD8/seurat.rds")
obj <- readRDS(args[1])

obj@meta.data$T <- str_replace_all(obj@meta.data$T, "[abc]", "")

# meta <- read.xlsx("20190707_meta.xlsx")
meta = read.xlsx(args[4])
rownames(meta) <- meta$PatientID

obj@meta.data$TumorSize = meta[obj@meta.data$PatientID, "Diameter.(cm)"]
obj@meta.data$SmokeYear = meta[obj@meta.data$PatientID, "Smoking_year"]
obj@meta.data$SmokePerDay = meta[obj@meta.data$PatientID, "Smoking"]
obj@meta.data$M = meta[obj@meta.data$PatientID, "M"]


# temp = read_lfp("04_metacell/LUSC/CD8/", meta = obj@meta.data)
temp = read_lfp(args[2], meta = obj@meta.data)


dir.create(args[3], showWarnings = F, recursive = T, mode = "0755")


# p1 <- make_heatmap(obj, group.by = "Stage", genes.use = EXHAUST, scale = T, converter = STAGE2NUM, title = "Stage")
# 
# p2 <- make_heatmap(obj, group.by = "T", genes.use = EXHAUST, scale = T, converter = TUMORSIZE2NUM, title = "T")
# 
# p3 <- make_heatmap(obj, group.by = "res.0.8", genes.use = EXHAUST, scale = T, title = "Cluster")
# 
# p4 <- make_heatmap(obj, group.by = "TumorSize", genes.use = EXHAUST, scale = T, title = "Tumor Size (cm)")
# 
# p5 <- make_heatmap(obj, group.by = "N", genes.use = EXHAUST, scale = T, converter = N2NUM, title = "N")
# 
# p6 <- make_heatmap(obj, group.by = "M", genes.use = EXHAUST, scale = T, converter = M2NUM, title = "M")
# 
# 
# ggsave(
#     filename = paste(args[3], "mean.pdf", sep = "/"),
#     plot = plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2),
#     width = 12,
#     height = 12
# )
# 
# 
# p1 <- make_heatmap(obj, group.by = "Stage", genes.use = EXHAUST, scale = T, converter = STAGE2NUM, mean = F, title = "Stage")
# 
# p2 <- make_heatmap(obj, group.by = "T", genes.use = EXHAUST, scale = T, converter = TUMORSIZE2NUM, mean = F, title = "T")
# 
# p3 <- make_heatmap(obj, group.by = "res.0.8", genes.use = EXHAUST, scale = T, mean = F, title = "Cluster")
# 
# p4 <- make_heatmap(obj, group.by = "TumorSize", genes.use = EXHAUST, scale = T, mean = F, title = "Tumor Size (cm)")
# 
# p5 <- make_heatmap(obj, group.by = "N", genes.use = EXHAUST, scale = T, converter = N2NUM, mean = F, title = "N")
# 
# p6 <- make_heatmap(obj, group.by = "M", genes.use = EXHAUST, scale = T, converter = M2NUM, mean = F, title = "M")
# 
# 
# ggsave(
#     filename = paste(args[3], "signle_cell.pdf", sep = "/"),
#     plot = plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2),
#     width = 12,
#     height = 12
# )


if(!is.null(temp)) {
    plist = list()
    # for (i in colnames(temp@meta.data)) {
    #     if(i %in% c("nGene", "nUMI", "orig.ident")) {
    #         next
    #     }
    #     
    #     plist[[length(plist) + 1]] = make_heatmap(temp, group.by = i, genes.use = EXHAUST, scale = T, set_x_breaks = F, label.x = 0.5, title = i)
    # }
    # 
    # ggsave(
    #     filename = paste(args[3], "metacell.pdf", sep = "/"),
    #     plot = plot_grid(plotlist = plist, ncol = 3),
    #     width = 18,
    #     height = 6 * ceiling(length(plist) / 3)
    # )
    
    
    for(i in c("N", "I", "M", "T")) {
        meta = temp@meta.data[, str_starts(colnames(temp@meta.data), i)]
        meta = melt(as.matrix(meta))
        colnames(meta) <- c("Cell", "Group", "Perc")
        
        expr = temp@raw.data[intersect(rownames(temp@raw.data), EXHAUST),]
        expr = melt(as.matrix(expr))
        colnames(expr) <- c("Gene", "Cell", "Expression")
        
        expr = merge(expr, meta, by = "Cell")
        
        if(nrow(expr) < 1) {
            next
        }
        
        p <- ggplot(
            expr,
            aes(x=Expression, y=Perc, color=Group)
        ) +
            geom_point() +
            facet_wrap(.~Gene, scales = "free") +
            stat_cor(aes(color = Group), method = "spearman") +
            geom_smooth(method = "lm", se = T) +
            scale_color_manual(values = wes_palette("Darjeeling1"))
        
        num_genes = length(unique(expr$Gene))
        ggsave(
            filename = paste(args[3], paste0(i, ".pdf"), sep = "/"),
            plot = p,
            width = 6 * ifelse(num_genes > 3, 3, num_genes),
            height = 6 * ceiling(num_genes / 3)
        )
    }
}



