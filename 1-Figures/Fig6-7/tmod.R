load_packages <- function() {
    library(tmod)
    library(stringr)
    library(ggplot2)
    library(wesanderson)
    library(doMC)
    library(dplyr)
}
suppressPackageStartupMessages(load_packages())


make_tmod_plots <- function(data, utest, ztest, pval=0.05, topn=10, N1=0.2, auc = 0.6, scales = "fixed") {

    data = data[data$adj.P.Val < pval, ]
    data = data[data$N1 / max(data$N1) > N1, ]
    data = data[data$AUC > auc, ]
    data = data[order(data$AUC, decreasing = T), ]
    
    data = data %>%
        group_by(ident) %>%
        top_n(topn, wt = AUC)

    ggplot(data, aes(x=ident, y=reorder(Title, adj.P.Val), color=adj.P.Val, size=AUC)) +
        geom_point() +
        scale_color_gradientn(colors=rev(wes_palette("Zissou1", 100, type = "continuous"))) +
        labs(x = "", y="") +
        theme_bw() +
        theme(
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 10)
        )
}


do_tmod <- function(data, output, msig, force=FALSE) {
    if (file.exists(output) && !force) {
        cerno = read.csv(output, row.names = 1, stringsAsFactors = F)
    } else {
        data = data[order(abs(data$avg_logFC), decreasing = T), ]
        cerno = NULL

        if ("ident" %in% colnames(data)) {
            for(i in unique(data$ident)) {
                print(i)

                temp_data = data[data$ident == i, ]
                temp = tmodCERNOtest(temp_data$gene, mset = msig)

                if(nrow(temp) > 0) {
                    temp$ident = i
                    temp$Methods = "CERNO"
                    cerno = rbind(as.data.frame(cerno), as.data.frame(temp))
                }
            }
        } else {
            cerno = tmodCERNOtest(data$gene, mset = msig)
            if (!is.null(cerno) && nrow(cerno) > 0) {
                cerno$Methods = "CERNO"
            }
        }
        
        
        if(!is.null(cerno)) {
            write.csv(cerno, output)
        }
    }
    
    return(cerno)
}



msig <- tmodImportMSigDB("project_scripts/12_final_plots/msigdb_v6.2.xml")
sel <- msig$MODULES$Category == "C5"


files = list.files("03_each_cells/Batch/", pattern = "cluster_markers.csv", full.names = T, recursive = T)

registerDoMC(10)
res = foreach (i = files, .errorhandling = "pass") %dopar% {
    print(i)
    
    data = read.csv(i, row.names = 1, stringsAsFactors = F)
    
    cerno = do_tmod(data, output = paste0(dirname(i), "/cerno.csv"), msig[sel], force = T)
    
    p <- make_tmod_plots(cerno)
    ggsave(
        filename = paste0(dirname(i), "/cerno.pdf"),
        plot = p,
        width = 8,
        height = 6
    )
}


## tmod network
msig <- tmodImportMSigDB("project_scripts/12_final_plots/msigdb_v6.2.xml")
sel <- msig$MODULES$Category == "C5"


atii_network <- readRDS("03_each_cells/ATII_Basal/network/atii.rds")
basal_network <- readRDS("03_each_cells/ATII_Basal/network/basal.rds")
atii_network_stage <- readRDS("03_each_cells/ATII_Basal/network/atii_stage.rds")
basal_network_stage <- readRDS("03_each_cells/ATII_Basal/network/basal_stage.rds")


do_tmod_graph <- function(graph, msig) {
    data = graph$centrality
    
    print(dim(data))
    
    data = data[order(data$Betweenness, decreasing = T), ]
    
    return(tmodCERNOtest(rownames(data), mset = msig))
}


do_tmod_to_network <- function(normal, stages, msig) {
    
    res <- do_tmod_graph(normal, msig)
    res$Methods = "CERNO"
    res$ident = "Normal"
    
    for (i in names(stages)) {
        print(i)
        temp <- do_tmod_graph(stages[[i]], msig)
        temp$Methods = "CERNO"
        temp$ident = i
        
        res = rbind(res, temp)
    }
    
    res
}

atii <- do_tmod_to_network(atii_network[["LUAD_Normal"]], atii_network_stage, msig[sel])
basal <- do_tmod_to_network(basal_network[["LUSC_Normal"]], basal_network_stage, msig[sel])