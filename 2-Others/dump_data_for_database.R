
files = list.files("11_CNV/each_cells", pattern = "seurat_obj.rds", recursive = T, full.names = T)
files

library(doMC)
library(openxlsx)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)



## cells
meta = readRDS("11_CNV/meta.rds")

temp = unique(meta[, c("Cells", "SampleID", "PatientID", "Disease", "cell_short", "Malignant")])


write.csv(temp, "LungCancer10x/03_each_cells/cells.csv", row.names = F)


modify = read.xlsx("11_CNV/each_cells/cluster_res.xlsx")


get_label = function(path) {
    path = str_split(path, "/")[[1]]
    
    for (i in 1:length(path)) {
        if (path[i] %in% c("LUAD", "LUSC")) {
            return(c("cell"=path[i - 1], "disease"=path[i]))
        }
    }
}

isEmpty <- function(x) {
    return(length(x)==0)
}


## Collect coord info
registerDoMC(10)
data = NULL
for (f in files)  {  # data = each, .combine = "rbind" %dopar% 
    if (str_detect(f, "by_stage_M_vs_NM")) {
        next # return(NULL) 
    }
    
    labels = get_label(f)
    
    if (is.null(labels)) {
        next # return(NULL)
    }
    
    res = modify[modify$Cell == labels["cell"], "Res"]
    
    if (labels["cell"] == "NE" || isEmpty(res)) {
        res = 0.1
    }
    print(labels)
    obj <- readRDS(f)
    res = paste("res", res, sep = ".")
    
    umap <- as.data.frame(obj@dr$umap@cell.embeddings)
    umap$Cell_id <- rownames(umap)
    umap <- umap[, c("Cell_id", "UMAP1", "UMAP2")]
    
    tsne = obj@dr$tsne
    if (is.null(tsne)) {
        obj <- RunTSNE(
            object = obj, 
            dims.use = 1:10, 
            do.fast = TRUE, 
            reduction.use = "pca", 
            reduction.name = "tsne",
            perplexity=min(30, floor(ncol(obj@dr$pca@cell.embeddings) / 3) + 1)
        )
        tsne = obj@dr$tsne
        
        saveRDS(obj, f)
    }
    
    tsne <- as.data.frame(tsne@cell.embeddings)
    
    umap <- cbind(umap, tsne[umap$Cell_id, ])
    
    umap$Cell <- labels["cell"]
    umap$Disease <- labels["disease"]
    
    umap$Cluster <- as.numeric(as.character(obj@meta.data[umap$Cell_id, res])) + 1
    data = rbind(data, umap)
    # umap
}



write.csv(data, "LungCancer10x/03_each_cells/coord_for_database.csv")


## collect markers info


res = NULL
for (f in files) {
    
    if (str_detect(f, "by_stage_M_vs_NM")) {
        next
    }
    
    labels = get_label(f)
    
    if (is.null(labels)) {
        next
    }
    
    print(f)
    indir = dirname(f)
    cell = labels[["cell"]]
    disease = labels[["disease"]]
    
    cluster_markers <- paste(indir, "cluster_markers.csv", sep = "/")
    
    if (file.exists(cluster_markers)) {
        temp <- read.csv(cluster_markers, stringsAsFactors = F, row.names = 1)
        
        if (nrow(temp) > 1) {
            temp$Cell = cell
            temp$Disease <- disease
            temp$Type <- "Cluster"
            res = rbind(res, temp)
        }
    }
    
    stage_markers <- paste(indir, "stage_markers.csv", sep = "/")
    if (file.exists(stage_markers)) {
        temp <- read.csv(stage_markers, stringsAsFactors = F, row.names = 1)
        
        if (nrow(temp) > 1) {
            temp$Cell = cell
            temp$Disease <- disease
            temp$Type <- "Stage"
            res = rbind(res, temp)
        }
    }
    
    # res
}
colnames(res) <- str_to_lower(colnames(res))


res1 = read.csv("11_CNV/scRNA_AD_M_vs_NM.csv", row.names = 1, stringsAsFactors = F)
res1$disease = "LUAD"
res1$ident = ifelse(res1$avg_logFC > 0, "Malignant", "Non-Malignant")
res1$type = "Malignant"
colnames(res1) <- str_to_lower(colnames(res1))


res2 = read.csv("11_CNV/scRNA_SC_M_vs_NM.csv", row.names = 1, stringsAsFactors = F)
res2$disease = "LUSC"
res2$ident = ifelse(res2$avg_logFC > 0, "Malignant", "Non-Malignant")
res2$type = "Malignant"
colnames(res2) <- str_to_lower(colnames(res2))

res = rbind(res, res1[, colnames(res)])
res = rbind(res, res2[, colnames(res)])


write.csv(res, "LungCancer10x/03_each_cells/markers_for_database.csv")


## mfuzz

files = list.files("11_CNV/each_cells", pattern = "mfuzz.csv", full.names = T, recursive = T)

registerDoMC(3)
res = foreach(f = files, .combine = "rbind") %dopar% {
    data <- read.csv(f, stringsAsFactors = F, row.names = 1)
    cell <- basename(dirname(dirname(f)))
    disease <- basename(dirname(f))
    
    data$disease = disease
    data$cell = cell

    data
}

write.csv(res, "LungCancer10x/03_each_cells/mfuzz_database.csv")


go = NULL
for (i in unique(res$cell)) {
    print(i)
    for (j in unique(res$disease)) {
        for (k in unique(res$Clt)) {
            temp = res[res$Clt == k & res$disease == j & res$cell == i, ]
            
            if (nrow(temp) > 0) {
                temp = enrichGO(
                    gene = temp$gene, OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL", ont = "BP"
                )
                
                temp = simplify(temp)
                
                temp = as.data.frame(temp)
                
                if (nrow(temp) > 0) {
                    temp$Disease = j
                    temp$cell = i
                    temp$Stage = str_replace_all(k, "^M\\.", "")
                    
                    go = rbind(go, temp)
                }
            }
        }
    }
}

write.csv(go, "LungCancer10x/03_each_cells/mfuzz_go.csv")