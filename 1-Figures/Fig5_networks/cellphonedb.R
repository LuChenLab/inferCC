

library(stringr)

meta <- readRDS(paste0(root.dir, "11_CNV/meta.rds"))
expr <- readRDS(paste0(root.dir, "02_rds/all_cell_expr.rds"))


write_data_for_cellphonedb <- function(meta, expr, outdir) {
    
    dir.create(outdir, showWarnings = F, recursive = T)
    
    temp <- meta[, c("Cells", "cell_short")]
    colnames(temp) <- c("Cell", "cell_type")
    
    temp_expr <- expr[, as.character(temp$Cell)]

    temp_expr$Gene = rownames(temp_expr)
    temp_expr <- temp_expr[, c("Gene", rownames(temp))]
    
    
    write.table(temp, paste0(outdir, "/meta.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    write.table(temp_expr, paste0(outdir, "/counts.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    
    
}


write_data_for_cellphonedb(
    meta[meta$Malignant & str_detect(meta$Disease, "LUAD"), ], 
    expr, "cellphonedb/LUAD/Malignant/"
)

write_data_for_cellphonedb(
    meta[!meta$Malignant & str_detect(meta$Disease, "LUAD"), ], 
    expr, "cellphonedb/LUAD/Non_Malignant/"
)

write_data_for_cellphonedb(
    meta[meta$Malignant & str_detect(meta$Disease, "LUSC"), ], 
    expr, "cellphonedb/LUSC/Malignant/"
)

write_data_for_cellphonedb(
    meta[!meta$Malignant & str_detect(meta$Disease, "LUSC"), ], 
    expr, "cellphonedb/LUSC/Non_Malignant/"
)



for (i in c("LUAD", "LUSC")) {
    temp_meta = meta[meta$Malignant & str_detect(meta$Disease, i), ]
    
    for (j in unique(temp_meta$Stage)) {
        
        print(paste(i, j))
        
        write_data_for_cellphonedb(
            temp_meta[temp_meta$Stage == j, ], 
            expr, paste("cellphonedb", i, "Stage", j, sep = "/")
        )
    }
}


for (i in c("LUAD", "LUSC")) {
    temp_meta = meta[meta$Malignant & str_detect(meta$Disease, i), ]
    
    print(paste(i, "Early"))
    
    write_data_for_cellphonedb(
        temp_meta[temp_meta$Stage %in% c("I", "II"), ], 
        expr, paste("cellphonedb", i, "Phage", "Early", sep = "/")
    )
    
    print(paste(i, "Advanced"))
    
    write_data_for_cellphonedb(
        temp_meta[temp_meta$Stage %in% c("III", "IV"), ], 
        expr, paste("cellphonedb", i, "Phage", "Advanced", sep = "/")
    )
}




## make cell-cell interaction network plot
library(igraph)
library(ggnetwork)
library(RCy3)


make_network_graph <- function(in_dir, save = T, layout = "fruchtermanreingold") {
  
    immu = c(
      "B", "CD4", "CD8", "DC","Gran",
      "Mast", "Mo", "NK", "Tregs", "Mφ"
    )
  
    nodes <- read.table(paste0(in_dir, "interaction_count.txt"), stringsAsFactors = F)
    nodes$Cell = rownames(nodes)
    nodes <- nodes[, c(2, 1)]
    
    nodes$Type = "Stromal"
    nodes$Type[nodes$Cell %in% immu] <- "Immune"
    nodes$Cell[nodes$Cell == "ATII"] = "AT2"
    nodes$Cell[nodes$Cell == "Mo"] = "Mφ"
    
  
    edges <- read.table(paste0(in_dir, "count_network.txt"), header = T, stringsAsFactors = F)
    colnames(edges)[1:2] <- c("from", "to")
    
    for (i in c("from", "to")) {
      edges[edges[, i] == "ATII", i] = "AT2"
      edges[edges[, i] == "Mo", i] = "Mφ"
    }

    g = graph_from_data_frame(edges, directed=TRUE, vertices=nodes)
    
    data <- fortify(g, layout = layout)
  

    # data$Type = "Stromal"
    # data$Type[data$name %in% immu] <- "Immune"
    # data$all_sum <- log10(data$all_sum + 1)
    # 
    # data$name = as.character(data$name)
    # data$name[data$name == "ATII"] = "AT2"
    # data$name[data$name == "Mo"] = "Mφ"
    
    # print(head(data))
    
    p <- ggplot(
        data, 
        aes(x = x, y = y, xend = xend, yend = yend, label = name)
    ) +
        geom_edges(color = "grey25", aes(alpha = count)) +
        geom_nodes(aes(size = all_sum * 5, color = Type)) +
        geom_nodetext_repel(
            size = 6,
            fontface = "bold", 
            box.padding = unit(1, "lines")
        ) +
        theme_blank(base_family = "Arial Unicode MS") +
        theme(
            legend.position = "none",
            legend.text = element_text(size = 10)
        ) +
        scale_color_manual(values = c(
            "Immune"="#E9B0B7",
            "Stromal"="#90D0E0"
        )) +
        labs(
            size = "Total num. of interaction", 
            alpha = "Num. of interaction",
            color = ""
        ) +
        guides(size = F, alpha = F) +
        scale_size(range = c(1, 30)) +
        scale_alpha(range = c(0.01, 1))
    
    if (save) {
        ggsave(
            filename = paste0(in_dir, "network.pdf"),
            plot = p, device = cairo_pdf,
            width = 6,
            height = 6
        )
    } else {
        p
    }

}

layouts = c(
    "fruchtermanreingold", "adj", "circle", 
    "circrand", "eigen", "geodist", "hall", 
    "kamadakawai", "mds", "princoord", "random", 
    "rmds", "segeo", "seham", "spring", "springrepulse", "target"
)
out_dir = "LungCancer/10x/cellphonedb/layout/"
dir.create(out_dir, showWarnings = F, recursive = T)
for (l in layouts) {
    tryCatch({
        print(paste0(out_dir, l, ".pdf"))
        p <- make_network_graph(in_dir, save = F, layout = l)
        ggsave(
            filename = paste0(out_dir, l, ".pdf"),
            plot = p, device = cairo_pdf,
            width = 6,
            height = 6
        )
    }, error=function(e) {
        print(e)
    })

}

temp = fortify(p, layout = "circle")

# disease = c("LUAD", "LUSC", "LUAD_Normal", "LUSC_Normal")
# stages = c("I", "II", "III", "IV")

disease = c("LUAD", "LUSC")  # , "LUAD_Normal", "LUSC_Normal"
stages = c("Malignant", "Non_Malignant") # c("I", "II", "III", "IV")




layouts = c(
  igraph::nicely(), igraph::in_circle(), igraph::on_sphere(),
  igraph::as_bipartite(), igraph::as_star(), igraph::with_dh(),
  igraph::with_fr(), igraph::with_drl(), igraph::with_gem()
)


set.seed(42)
for (i in disease) {
    in_dir = paste0("cellphonedb/", i, "/out/")
    
    if (file.exists(in_dir)) {
        print(in_dir)
        
        p <- make_network_graph(in_dir, save = F, layout = igraph::in_circle())
        print(p)
    }
    
    for (j in stages) {
        in_dir = paste0("cellphonedb/", i, "/", j, "/out/")

        if (file.exists(in_dir)) {
            print(in_dir)

            make_network_graph(in_dir, save = T, layout = igraph::in_circle())
        }
    }
    
    for (j in c("Early", "Advanced")) {
        in_dir = paste0("cellphonedb/", i, "/Phage/", j, "/out/")
        
        if (file.exists(in_dir)) {
          print(in_dir)
          
          make_network_graph(in_dir, save = T, layout = igraph::in_circle())
        }
    }
}

saveRDS(p, "cellphonedb/LUAD/Malignant/out/myIgraph")
createNetworkFromIgraph(p,"cellphonedb/LUAD/Malignant/out/myIgraph")


### make heatmap

draw_heatmap <- function(res, output_file) {
    library(stringr)
    library(reshape2)
    library(ComplexHeatmap)
    library(circlize)
    
    
    data <- melt(res)
    data$Stage <- sapply(data$variable, function(x) {
        str_split(x, "_")[[1]][2]
    })
    data$Disease <- sapply(data$variable, function(x) {
        str_split(x, "_")[[1]][1]
    })
    
    temp_data1 <- dcast(data[data$Disease == "LUAD", ], Cell~Stage, fun.aggregate = mean)
    rownames(temp_data1) <- temp_data1$Cell
    temp_data1 <- temp_data1[, colnames(temp_data1) != "Cell"]
    
    
    temp_data2 <- dcast(data[data$Disease == "LUSC", ], Cell~Stage, fun.aggregate = mean)
    rownames(temp_data2) <- temp_data2$Cell
    temp_data2 <- temp_data2[, colnames(temp_data2) != "Cell"]
    
    or <- intersect(rownames(temp_data1), rownames(temp_data2))
    temp_data <- cbind(
        temp_data1[or, intersect(c("Normal", "I", "II", "III", "IV", "Malignant", "Non"), colnames(temp_data1))], 
        temp_data2[or, intersect(c("Normal", "I", "II", "III", "IV", "Malignant", "Non"), colnames(temp_data2))]
    )
    
    immu = c(
        "B",
        "CD4",
        "CD8",
        "DC",
        #"EC",
        "Gran",
        "Mast",
        "Mo",
        "NK",
        "Tregs"
    )
    
    
    non_immu <- sort(rownames(temp_data)[!rownames(temp_data) %in% immu])
    
    pdf(output_file, width = 6, height = 6)
    h <- Heatmap(
        na.omit(temp_data[c(non_immu, immu), ]), 
        name = "",
        cluster_columns = F,
        cluster_rows = T,
        column_split = c(
            rep("LUAD", ncol(temp_data1)), 
            rep("LUSC", ncol(temp_data2))
        ),
        row_split = c(
            rep("Non-immune", length(non_immu)), 
            rep("Immune", length(immu))
        )
    )
    draw(h)
    dev.off()
}


res = NULL
for (i in disease) {
    in_dir = paste0("cellphonedb/", i, "/out/")
    
    for (j in stages) {
        in_dir = paste0("cellphonedb/", i, "/", j, "/out/")
        
        if (file.exists(in_dir)) {
            print(in_dir)
            
            temp <- read.table(paste0(in_dir, "interaction_count.txt"))
            colnames(temp) = paste(i, j, sep = "_")
            temp$Cell = rownames(temp)
            
            if (is.null(res)) {
                res = temp
            } else {
                res = merge(res, temp, by = "Cell")
            }
        }
    }
}


draw_heatmap(res, "cellphonedb/cell_interaction_heatmap_raw.pdf")

res[, 2:ncol(res)] = apply(res[, 2:ncol(res)], 2, function(col) {
    col / sum(col)
})

draw_heatmap(res, "cellphonedb/cell_interaction_heatmap.pdf")



colnames(temp_data) <- paste(colnames(temp_data), c(rep("LUAD", ncol(temp_data1)), rep("LUSC", ncol(temp_data2))))
temp_data = melt(as.matrix(temp_data))
temp_data$Immu = temp_data$Var1 %in% immu

temp_data$Disease = sapply(temp_data$Var2, function(x) { str_split(x, " ")[[1]][2] })
temp_data$Type = sapply(temp_data$Var2, function(x) { str_split(x, " ")[[1]][1] })


ggplot(temp_data, aes(x=Type, y=Var1, color=value, size=value)) +
    geom_point() +
    facet_grid(Immu~Disease, scale = "free")


genes = c("ATF3", "CTTN", "CYP2B7P", "ANXA1", "ELF3", "KLF6", "S100A13")

cluster_markers <- read.csv("LungCancer10x/11_CNV/each_cells/ATII/LUAD/cluster_markers.csv", row.names = 1, stringsAsFactors = F)
stage_markers <- read.csv("LungCancer10x/11_CNV/each_cells/ATII/LUAD/stage_markers.csv", row.names = 1, stringsAsFactors = F)

temp = cluster_markers[cluster_markers$gene %in% genes, ]
temp1 = stage_markers[stage_markers$gene %in% genes, ]



### interaction pair GO annotation

read_pairs <- function(path) {
  data <- read.table(path, sep = "\t", header = T)
  
  
  temp_pairs = c()
  for (i in 1:nrow(data)) {
    temp = data[i, 13:ncol(data)]
    temp_name = colnames(data)[which(temp > 0) + 12]
    temp_pairs = c(temp_pairs, paste(temp_name, collapse = ","))
  }
  
  data$TargetPair = temp_pairs
  data
}

get_label <- function(path) {
  library(stringr)
  path <- str_split(path, "/")[[1]]
  
  i = 1
  while (i <= length(path)) {
    if (path[i] == "out") {
      
      if (path[i - 2] %in% c("Phage", "Stage")) {
        return (c(path[i - 3], path[i - 1]))
      } else {
        return (c(path[i - 2], path[i - 1]))
      }
    }
    i= i + 1
  }
  return(NULL)
}


files = list.files("LUAD", pattern="significant_means.txt", recursive = T, full.names = T)
files = c(files, list.files("LUSC", pattern="significant_means.txt", recursive = T, full.names = T))


interactions = NULL
for (f in files) {
  print(f)
  data = read_pairs(f) #  by_cell = F
  
  temp = data %>%
    group_by(interacting_pair, TargetPair) %>%
    add_tally() %>%
    dplyr::select(interacting_pair, TargetPair, n) %>%
    unique() %>%
    as.data.frame()
  
  labels = get_label(f)
  temp$ident = labels[2]
  temp$Disease = labels[1]
  
  interactions = rbind(interactions, temp)
}

write.csv(interactions, "interaction.csv")


temp = interactions[str_detect(interactions$TargetPair, "\\.Mo"), ]

temp1 = temp %>%
  dplyr::select(interacting_pair, ident, Disease) %>%
  # group_by(interacting_pair, ident) %>%  # , Disease
  # add_tally() %>%
  unique() %>%
  as.data.frame()

temp1$n = 1

temp2 = dcast(temp1, interacting_pair+ident~Disease, value.var = "n", fun.aggregate = mean, fill = 0)

write.csv(temp2, "MP_pairs.csv")


go = NULL
for (f in files) {
  print(f)
  data = read_pairs(f, by_cell = T)
  
  pairs = colnames(data)[13:ncol(data)]
  pairs = sapply(pairs[pairs != "TargetPair"], function(x) {
    x = str_split(x, "\\.")[[1]]
    x[1]
  })
  pairs = unique(pairs)
  
  for (p in pairs) {
    print(p)
    temp = data[str_detect(data$TargetPair, paste0(",?", p, "[,$]?")), ]
    genes = unique(temp$gene_a, temp$gene_b)
    
    ego <- enrichGO(gene          = genes,
                    keyType       = "SYMBOL",
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = FALSE)
    
    if (is.null(ego)) {
      next
    }
    labels = get_label(f)
    temp = as.data.frame(simplify(ego))
    if (nrow(temp) > 1) {
      temp$cell = p
      temp$Disease = labels[1]
      temp$ident = labels[2]
      go <- rbind(go, temp)
    }
  }
}


write.csv(go, "all_go.csv")


setwd("LungCancer/10x/cellphonedb")
go <- readRDS("all_go.rds")
sel <- read.xlsx("select.xlsx")


make_go_plot <- function(go, count = 0, topn = 10, title = "", subtitle = "") {
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(wesanderson)
  
  go$GeneRatio <- sapply(go$GeneRatio, function(x) {
    x = str_split(x, "/")[[1]]
    as.numeric(x[1]) / as.numeric(x[2])
  })
  
    temp = go %>%
      filter(Count > count & nchar(Description) < 50) %>%
      group_by(ident) %>%
      top_n(topn, wt=-qvalue) %>%
      as.data.frame()
  
  
  if (nrow(temp) == 0) {
    return(NULL)
  }
  
  
  ggplot(temp, aes(x=ident, y=Description, size=GeneRatio, color = qvalue)) +
    geom_point() +
    scale_color_gradientn(colours = rev(wes_palette("Zissou1", 50, type = "continuous"))) +
    theme_classic() + # base_family = "Arial Unicode MS"
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text = element_text(size = 12)
    ) +
    labs(title = title, subtitle = subtitle, x="")
}

plot(density(go$Count))

go$ident <- factor(go$ident, levels = c("Non_Malignant", "Malignant", "Early", "Advanced", "I", "II", "III" ))

for (i in c("LUSC")) { # "LUAD", 
  # ifelse(i == "LUAD", 15, 20)
  for (j in unique(go$cell)) {
    temp = go[go$Disease == i & go$cell == j, ]
      
    # temp = temp[as.character(temp$Description) %in% unique(as.character(sel$Pathway[str_detect(sel$cell, paste(j, i, sep = "_"))])), ]
    temp <- temp[as.character(temp$ident) %in% c("Non_Malignant", "Early", "Advanced"), ]
    
    p <- make_go_plot(
      temp, 
      count = 20, title = j, subtitle = i, 
      topn = 20
    )

    ggsave(
      filename = paste0(i, "/", j, "_go.pdf"),
      plot = p,
      units = "in",
      # dpi = 600,
      width = 8,
      height = 6
    )
  }
}

go$row = 1:nrow(go)

temp = go %>%
  filter(Count > 20 & nchar(Description) < 50) %>%
  group_by(Disease, cell, ident) %>%
  top_n(20, wt=-qvalue) %>%
  as.data.frame()


go$selected <- go$row %in% temp$row
go$GeneRatio <- sapply(go$GeneRatio, function(x) {
  x = str_split(x, "/")[[1]]
  as.numeric(x[1]) / as.numeric(x[2])
})


write.csv(go, "all_go_selected.csv")


for (f in files) {
  print(f)
  
  data = read_pairs(f)
  pairs = colnames(data)[13:ncol(data)]
  pairs = pairs[pairs != "TargetPair"]
  
  go = NULL
  for (p in pairs) {
    print(p)
    temp = data[str_detect(data$TargetPair, paste0(",?", p, "[,$]?")), ]
    genes = unique(temp$gene_a, temp$gene_b)
    
    ego <- enrichGO(gene          = genes,
                    keyType       = "SYMBOL",
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = FALSE)
    
    if (is.null(ego)) {
      next
    }
    
    temp = as.data.frame(simplify(ego))
    if (nrow(temp) > 1) {
      temp$pair = p
      go <- rbind(go, temp)
    }
  }
  
  write.csv(ego, paste(dirname(f), "go_pair.csv", sep = "/"))
}



## Interaction count bar plot
val = NULL
for (i in list.files("cellphonedb/LUAD", pattern = "interaction_count.txt", recursive = T, full.names = T)) {
    temp <- read.table(i)
    lab <- get_label(i)
    temp$ident = lab[2]
    temp$Disease = lab[1]
    temp$cell = rownames(temp)
    
    val <- rbind(val, temp)
}


for (i in list.files("cellphonedb/LUSC", pattern = "interaction_count.txt", recursive = T, full.names = T)) {
  temp <- read.table(i)
  lab <- get_label(i)
  temp$ident = lab[2]
  temp$Disease = lab[1]
  temp$cell = rownames(temp)
  val <- rbind(val, temp)
}


val$Immu <- val$cell %in% immu

val$ident <- str_replace_all(val$ident, "_", " ")
val <- val[val$ident %in% c("Non Malignant", "Early", "Advanced"), ]
val$ident <- factor(val$ident, levels = c("Non Malignant", "Early", "Advanced"))


val$cell[val$cell == "ATII"] = "AT2"
val$cell[val$cell == "Mo"] = "Mφ"


extrafont::font_import()
extrafont::loadfonts()


p <- ggplot(val[val$Immu, ], aes(x=ident, y=all_sum, fill=ident)) +
  geom_bar(stat = "identity") +
  facet_grid(Disease~cell) +
  theme_classic(base_family = "Arial Unicode MS") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = c(0.2, 0.95),
    legend.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) + 
  labs(x="", y="", fill="") +
  guides(fill=guide_legend(nrow=2)) +
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))


ggsave(
  filename = "cellphonedb/immune_barplot.pdf",
  plot = p, width = 8, height = 5, device=cairo_pdf
)


p <- ggplot(val[!val$Immu, ], aes(x=ident, y=all_sum, fill=ident)) +
  geom_bar(stat = "identity") +
  facet_grid(Disease~cell) +
  theme_classic(base_family = "Arial Unicode MS") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = c(0.5, 0.95),
    legend.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) + 
  labs(x="", y="", fill="") +
  guides(fill=guide_legend(nrow=1)) +
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))


ggsave(
  filename = "cellphonedb/non_immune_barplot.pdf",
  plot = p, width = 8, height = 5, device = cairo_pdf
)



## 

temp = reshape2::dcast(val, cell+Disease~ident, value.var = "all_sum", fill = 0, fun.aggregate = mean)
temp$Advanced = temp$Advanced / temp$`Non Malignant`
temp$Early = temp$Early / temp$`Non Malignant`

temp = reshape2::melt(temp[, c("cell", "Disease", "Early", "Advanced")])

temp$variable = factor(temp$variable, levels = c("Early", "Advanced"))



p <- ggplot(temp[temp$cell %in% immu, ], aes(x=variable, y=value, fill=variable)) +
  geom_bar(stat = "identity") +
  facet_grid(Disease~cell) +
  theme_classic(base_family = "Arial Unicode MS") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = c(0.5, 0.95),
    legend.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) + 
  labs(x="", y="", fill="") +
  guides(fill=guide_legend(nrow=1)) +
  scale_fill_manual(values=c("#00AFBB", "#FC4E07"))  # "#E7B800"


ggsave(
  filename = "cellphonedb/immune_barplot_ratio.pdf",
  plot = p, width = 5, height = 3, device = cairo_pdf
)


p <- ggplot(temp[!temp$cell %in% immu, ], aes(x=variable, y=value, fill=variable)) +
  geom_bar(stat = "identity") +
  facet_grid(Disease~cell) +
  theme_classic(base_family = "Arial Unicode MS") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = c(0.5, 0.95),
    legend.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) + 
  labs(x="", y="", fill="") +
  guides(fill=guide_legend(nrow=1)) +
  scale_fill_manual(values=c("#00AFBB", "#FC4E07"))  # "#E7B800", 


ggsave(
  filename = "cellphonedb/non_immune_barplot_ratio.pdf",
  plot = p, width = 5, height = 3, device = cairo_pdf
)



## count the pairs show up in AT2 and Basal and Macrophages