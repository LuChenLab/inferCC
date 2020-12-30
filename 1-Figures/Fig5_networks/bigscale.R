options(stringsAsFactors = F)

library(openxlsx)
library(bigSCale)
library(stringr)
library(R.utils)
library(ggnetwork)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(ggpubr)
library(igraph)
library(wesanderson)
library(doMC)


calculate_network <- function(obj, clustering="recursive") {
    network <- list()
    for(i in unique(obj@meta.data$Disease)) {
        tryCatch({
            sel_cells <- rownames(obj@meta.data[obj@meta.data$Disease == i, ])
            
            mtx = obj@raw.data[, colnames(obj@raw.data) %in% sel_cells]
            
            withTimeout({
                network[[i]] <- compute.network(mtx, gene.names = rownames(mtx), clustering=clustering)
            }, timeout = 18000)
            
            
        }, error = function(x) {
            
        }, TimeoutException = function(ex) {
            message(paste("Timeout. Skipping.", i))
        })
        
    }
        
    return(network)
}

atii <- readRDS("03_each_cells/ATII_Basal/all_atii_cells.rds")
atii_network <- calculate_network(atii, clustering = "direct")
saveRDS(atii_network, "03_each_cells/ATII_Basal/network/atii.rds")

basal <- readRDS("03_each_cells/ATII_Basal/all_basal_cells.rds")
basal_network <- calculate_network(basal, clustering = "direct")
saveRDS(basal_network, "03_each_cells/ATII_Basal/network/basal.rds")



calculate_network_by_stage <- function(obj, disease, clustering="recursive", verbose = F) {
    network <- list()
    for(i in unique(obj@meta.data$Stage)) {
        if (verbose) {
            print(i)
            print(sum(obj@meta.data$Disease == disease & obj@meta.data$Stage == i))
        }
        
        tryCatch({
            sel_cells <- rownames(obj@meta.data[obj@meta.data$Disease == disease & obj@meta.data$Stage == i, ])
            
            mtx = obj@raw.data[, colnames(obj@raw.data) %in% sel_cells]
            
            withTimeout({
                network[[i]] <- compute.network(mtx, gene.names = rownames(mtx), clustering=clustering)
            }, timeout = 18000)
            
            
        }, error = function(x) {
            
        }, TimeoutException = function(ex) {
            message(paste("Timeout. Skipping.", i))
        })
        
    }
    
    return(network)
}


atii_network <- calculate_network_by_stage(atii, "LUAD", clustering = "direct", verbose = T)
sel_cells <- rownames(atii@meta.data[atii@meta.data$Disease == "LUAD" & atii@meta.data$Stage == "II", ])
mtx = atii@raw.data[, colnames(atii@raw.data) %in% sel_cells]
atii_network[["II"]] <- compute.network(mtx, gene.names = rownames(mtx), clustering="normal")
saveRDS(atii_network, "03_each_cells/ATII_Basal/network/atii_stage.rds")

basal_network <- calculate_network_by_stage(basal, "LUSC", clustering = "direct", verbose = T)
saveRDS(basal_network, "03_each_cells/ATII_Basal/network/basal_stage.rds")



## Plot

library(ggnetwork)

atii_network <- readRDS("03_each_cells/ATII_Basal/network/atii.rds")
basal_network <- readRDS("03_each_cells/ATII_Basal/network/basal.rds")

names(atii_network)


make_network_graph <- function(network, top_n = 0.001, title = NULL, genes = NULL, red_label = NULL) {
    
    if (!is.null(genes)) {
        el <- as_edgelist(atii_network_stage[["I"]]$graph, names = TRUE)
        targets <- as.vector(el[el[,1] %in% genes | el[,2] %in% genes, ])
        graph <- induced_subgraph(atii_network_stage[["I"]]$graph, targets)
    } else {
        graph = network$graph
    }
    
    data = fortify(graph)
    for(i in colnames(network$centrality)) {
        data[, i] <- network$centrality[data$vertex.names, i]
    }
    
    data <- data[order(data$Betweenness, decreasing = F), ]
    genes = unique(data$vertex.names[order(data$Betweenness, decreasing = T)])
    
    if (top_n <= 1) {
        genes = genes[1:(round(nrow(data) * top_n))]
    } else {
        genes = genes[1:top_n]
    }
    
    normal_label = data[data$vertex.names %in% setdiff(genes, red_label), ]

    if (!is.null(red_label)) {
        red_label = data[data$vertex.names %in% red_label, ]
    }
   
    p <- ggplot(
        data, 
        aes(x = x, y = y, xend = xend, yend = yend)
    ) +
        geom_edges(color = "grey75") +
        geom_nodes(aes(size = Degree, color = Betweenness)) +
        geom_nodelabel_repel(
            data = normal_label, # function(x) { x[(nrow(x) - round(0.001 * nrow(x))):nrow(x), ]},
            aes(label = vertex.names),
            fontface = "bold", 
            box.padding = unit(1, "lines")
        ) +
        theme_blank() +
        theme(
            plot.title = element_text(size = 18, hjust = 0.5)
        ) +
        scale_colour_gradientn(colors = wes_palette("Zissou1", 100, type = "continuous")) +
        labs(title = title)
    
    if (!is.null(red_label)) {
        p <- p + 
            geom_nodelabel_repel(
                data = red_label, # function(x) { x[(nrow(x) - round(0.001 * nrow(x))):nrow(x), ]},
                aes(label = vertex.names),
                fontface = "bold", 
                box.padding = unit(1, "lines"),
                color = "#DA0118"
            )
    }
    
    p
}


make_network_graph(atii_network[["LUAD"]], title = str_replace_all("LUAD", "_", " "))


for(i in names(atii_network)) {
    print(i)
    tryCatch({
        p <- make_network_graph(atii_network[[i]], title = str_replace_all(i, "_", " "))
        
        ggsave(
            filename = paste0("03_each_cells/ATII_Basal/network/all_atii/ATII_", i, ".pdf"),
            plot = p,
            width = 8,
            height = 6
        )
    }, error = function(x) {})
}


ggsave(
    filename = paste0("03_each_cells/ATII_Basal/network/all_atii/Basal_LUSC.pdf"),
    plot = make_network_graph(basal_network[["LUSC"]], title = str_replace_all("LUSC", "_", " "), top_n = 30),
    width = 8,
    height = 6
)



get_top_markers <- function(network, top=0.001) {
    data = ggnetwork::fortify(network$graph)
    
    for(i in colnames(network$centrality)) {
        data[, i] <- network$centrality[data$vertex.names, i]
    }
    
    data <- data[order(data$Betweenness, decreasing = F), ]
    
    genes = unique(data$vertex.names[order(data$Betweenness, decreasing = T)])
    genes = genes[1:(round(nrow(data) * top))]
    genes
}

gene = as.data.frame(get_top_markers(basal_network[["LUSC"]]))



## Network by Stage

format_gene_value_mtx <- function(base, networks, ident="Betweenness") {
    
    if (!is.null(base)) {
        base = base$centrality[, ident, drop = F]
        colnames(base) <- "Normal"    
    }

    for(i in names(networks)) {
        if (is.null(base)) {
            base = networks[[i]]$centrality[, ident, drop = F]
            colnames(base) <- i
        } else {
            base[, i] = networks[[i]]$centrality[rownames(base), ident]
        }
    }
    
    base[is.na(base)] = 0
    base1 <- do.call(data.frame,lapply(base, function(x) replace(x, is.infinite(x),1)))
    rownames(base1) <- rownames(base)
    base1
}


make_line_plot_of_value <- function(mat, groups) {
    mat <- as.data.frame(mat)
    mat$group = NA
    
    for(i in 1:length(groups)) {
        mat[groups[[i]], "group"] = i
    }
    
    mat$gene = rownames(mat)
    data <- reshape2::melt(mat, id = c("group", "gene"))

    data$variable = factor(data$variable, levels = c("Normal", "I", "II", "III", "IV"))
    # data$variable = as.numeric(data$variable)
    data$group = factor(data$group)

    ggline(
        data, 
        x = "variable", 
        y = "value",
        color = "group",
        add = c("mean_se", "violin"), 
        facet.by = "group"
    ) 
}

atii_network <- readRDS("03_each_cells/ATII_Basal/network/atii.rds")
basal_network <- readRDS("03_each_cells/ATII_Basal/network/basal.rds")
atii_network_stage <- readRDS("03_each_cells/ATII_Basal/network/atii_stage.rds")
basal_network_stage <- readRDS("03_each_cells/ATII_Basal/network/basal_stage.rds")


temp = atii_network_stage[["I"]]$centrality
temp = basal_network_stage[["I"]]$centrality
temp = temp[order(temp$Betweenness, decreasing = T), ]
head(rownames(temp))

registerDoMC(4)
foreach(i = names(atii_network_stage)) %dopar% {
    print(i)
    tryCatch({
        
        if (i == "I") {
            p <- make_network_graph(atii_network_stage[[i]], title = i, top_n=20, red_label = c("PIGR", "SPINK1", "SFTPA2", "MUC1", "CYP2B7P"))
        } else {
            p <- make_network_graph(atii_network_stage[[i]], title = i, top_n=20)
        }
        
        ggsave(
            filename = paste0("03_each_cells/ATII_Basal/network/all_atii/ATII_LUAD_", i, ".pdf"),
            plot = p,
            width = 8,
            height = 6
        )
    }, error = function(x) {})
}


foreach(i = names(basal_network_stage)) {
    print(i)
    tryCatch({
        if (i == "I") {
            p <- make_network_graph(basal_network_stage[[i]], title = i, top_n = 20, red_label = c("VIM", "IGFBP7"))
        } else {
            p <- make_network_graph(atii_network_stage[[i]], title = i, top_n=20)
        }
    
        ggsave(
            filename = paste0("03_each_cells/ATII_Basal/network/all_atii/Basal_LUSC_", i, ".pdf"),
            plot = p,
            width = 8,
            height = 6
        )
    }, error = function(x) {})
}



p <- make_network_graph(atii_network_stage[["I"]], genes = "CYP2B7P", top_n = 1)
ggsave(
    filename = paste0("03_each_cells/ATII_Basal/network/all_atii/ATII_LUAD_CYP2B7P.pdf"),
    plot = p,
    width = 4,
    height = 3
)


p <- make_network_graph(atii_network_stage[["I"]], genes = "PIGR", top_n = 1)
ggsave(
    filename = paste0("03_each_cells/ATII_Basal/network/all_atii/ATII_LUAD_PIGR.pdf"),
    plot = p,
    width = 4,
    height = 3
)


#### cluster ATII
atii_mtx <- format_gene_value_mtx(atii_network[["LUAD_Normal"]], atii_network_stage)
atii_mtx <- apply(atii_mtx, 2, function(col) {
    col / sum(col)
})


atii_mtx_scale <- na.omit(t(scale(t(atii_mtx))))

h <- Heatmap(
    atii_mtx_scale,
    show_row_names = F,
    border = T,
    row_split = 20
)

# draw(h)


row_o <- row_order(h)

atii_markers_group = data.frame(gene = rownames(atii_mtx_scale), Clt = NA)
for(i in 1:length(row_o)) {
    atii_markers_group[row_o[[i]], "Clt"] = i
}

write.csv(atii_markers_group, "03_each_cells/ATII_Basal/network/ATII_network_BC_clt.csv")


p <- make_line_plot_of_value(as.data.frame(atii_mtx_scale), row_o) +
    labs(x = "", y = "") +
    theme(
        legend.position = "none",
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )


ggsave(
    filename = "03_each_cells/ATII_Basal/network/ATII_network_gene_clt.pdf",
    plot = p,
    width = 10,
    height = 8
)


#### clsuter Basal
basal_mtx <- format_gene_value_mtx(NULL, basal_network_stage)
basal_mtx <- apply(basal_mtx, 2, function(col) {
    col / sum(col)
})


basal_mtx_scale <- na.omit(t(scale(t(basal_mtx))))
basal_mtx_scale <- basal_mtx_scale[, colnames(basal_mtx_scale) != "Normal"]

h <- Heatmap(
    basal_mtx_scale,
    show_row_names = F,
    border = T,
    row_split = 10
)

draw(h)


row_o <- row_order(h)


basal_markers_group = data.frame(gene = rownames(basal_mtx_scale), Clt = NA)
for(i in 1:length(row_o)) {
    basal_markers_group[row_o[[i]], "Clt"] = i
}

write.csv(basal_markers_group, "03_each_cells/ATII_Basal/network/Basal_network_BC_clt.csv")


p <- make_line_plot_of_value(as.data.frame(basal_mtx_scale), row_o) +
    labs(x = "", y = "") +
    theme(
        legend.position = "none",
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 15)
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

p

ggsave(
    filename = "03_each_cells/ATII_Basal/network/Basal_network_gene_clt.pdf",
    plot = p,
    width = 10,
    height = 8
)



############################################################################
# 
# Make plot of different BC cluster
# 
############################################################################
atii_markers_group <- read.csv("03_each_cells/ATII_Basal/network/ATII_network_BC_clt.csv", row.names = 1, stringsAsFactors = F)
basal_markers_group <- read.csv("03_each_cells/ATII_Basal/network/Basal_network_BC_clt.csv", row.names = 1, stringsAsFactors = F)


atii <- readRDS("03_each_cells/ATII_Basal/all_atii_cells.rds")
basal <- readRDS("03_each_cells/ATII_Basal/all_basal_cells.rds")

make_heatmap_by_group <- function(markers, obj, disease = "LUAD") {
    
    meta <- obj@meta.data[str_detect(obj@meta.data$Disease, disease), ]
    meta <- meta[order(meta$Disease, meta$Stage), ]
    
    
    colors = list(
        Stage = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793", "BSPN"="#73BDFF"),
        Disease = c("LUAD" = "#0084D1", "NL(AD)"="#FECC1B", "NL(AD"="#778793", "BSPN"="#73BDFF", "NL(SC)"="#778793", "LUSC"="#A0C807"),
        Batch = c("1"="#E9B0B7", "2"="#90D0E0")
    )
    
    ba <- HeatmapAnnotation(
        Disease = sapply(meta$Disease, function(x) {
            disease = c(
                "LUAD"="LUAD",
                "LUAD_Normal"="NL(AD)",
                "LUSC"="LUSC",
                "LUSC_Normal"="NL(SC)"
            )
            
            disease[[x]]
        }),
        Stage = meta$Stage,
        col = colors,
        show_annotation_name = T,
        show_legend = F
    )
    
    Heatmap(
        obj@scale.data[markers$gene, rownames(meta)],
        name = "Expr",
        border = T,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        top_annotation = ba,
        row_split = markers$Clt
    )
}

make_heatmap_by_group(atii_markers_group, atii)

make_violin_by_group <- function(markers, obj, disease = "LUAD") {
    library(ggpubr)
    meta <- obj@meta.data[str_detect(obj@meta.data$Disease, disease), ]

    colors = list(
        Stage = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793", "BSPN"="#73BDFF"),
        Disease = c("LUAD" = "#0084D1", "NL(AD)"="#FECC1B", "NL(AD"="#778793", "BSPN"="#73BDFF", "NL(SC)"="#778793", "LUSC"="#A0C807"),
        Batch = c("1"="#E9B0B7", "2"="#90D0E0")
    )
    
    expr <- melt(as.matrix(obj@scale.data[markers$gene, rownames(meta)]))
    expr$ident <- meta[expr$Var2, "Disease"]
    expr$ident <- sapply(expr$ident, function(x) {
        disease <- c(
            "LUAD"="LUAD",
            "LUSC"="LUSC",
            "LUAD_Normal"="NL(AD)",
            "LUSC_Normal"="NL(SC)"
        )
        
        disease[[x]]
    })
    
    expr[expr$ident %in% c("LUAD", "LUSC"), "ident"] <- meta[expr[expr$ident %in% c("LUAD", "LUSC"), "Var2"], "Stage"]
    expr$group = markers[expr$Var1, "Clt"]
    
    ggviolin(
        data = expr,
        x = "ident",
        y = "value",
        fill = "ident"
        # add = "boxplot"
        # add.params = list(fill = "white")
    ) +
        facet_wrap(.~group, scales = "free")
    
}

p <- make_violin_by_group(atii_markers_group, atii)
ggsave(
    filename = "03_each_cells/ATII_Basal/network/ATII_expr_gene_clt.png",
    plot = p,
    width = 16,
    height = 16,
    dpi = 600
)


p <- make_violin_by_group(basal_markers_group, basal, "LUSC")
ggsave(
    filename = "03_each_cells/ATII_Basal/network/Basal_expr_gene_clt.png",
    plot = p,
    width = 8,
    height = 8,
    dpi = 600
)