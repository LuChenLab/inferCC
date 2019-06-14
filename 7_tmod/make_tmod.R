library(tmod)
library(openxlsx)
library(limma)
# library(GO.db)
# library(org.Hs.eg.db)
# library(KEGGREST)
# library(clusterProfiler)


# make GO tmod
# mtab <- toTable(org.Hs.egGO)
# 
# # mtab <- mtab[ mtab$Ontology == "BP", ]
# temp_m2g <- split(mtab$gene_id, mtab$go_id)
# 
# m2g = list()
# 
# for (i in names(temp_m2g)) {
#     eg = bitr(temp_m2g[[i]], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
#     m2g[[i]] = eg$SYMBOL
# }
# 
# ## remove the rather large object
# rm(mtab)
# ll <- lengths(m2g)
# m2g <- m2g[ ll >= 10 & ll <= 100 ]
# 
# gt <- toTable(GOTERM)
# m <- data.frame(ID=names(m2g))
# m$Title <- gt$Term[ match(m$ID, gt$go_id) ]
# goset <- makeTmod(modules=m, modules2genes=m2g)
goset = readRDS("/mnt/raid62/Lung_cancer_10x/00_data_ingest/00_raw_data/tmod_goSet.rds")

# make KEGG tmod
# pathways <- keggLink("pathway", "hsa")
# ## get pathway Names in addition to IDs
# paths <- sapply(unique(pathways), function(p) keggGet(p)[[1]]$NAME)
# m <- data.frame(ID=unique(pathways), Title=paths)
#
# ## m2g is the mapping from modules (pathways) to genes
# temp_m2g <- split(names(pathways), pathways)

# m2g = list()
# for(i in names(temp_m2g)) {
#     eg = bitr(sapply(temp_m2g[[i]], function(x) {strsplit(x, ":")[[1]][2]}), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
#     m2g[[i]] = eg$SYMBOL
#     
# }

# ## kegg object can now be used with tmod
# kegg <- makeTmod(modules=m, modules2genes=m2g)
kegg = readRDS("/mnt/raid62/Lung_cancer_10x/00_data_ingest/00_raw_data/tmod_kegg.rds")

args = commandArgs(trailingOnly = T)

genes = args[1]
msig = args[2]
expr = args[3]
output_file = args[4]


genes = read.table(genes, header = F)
msig <- tmodImportMSigDB(msig)
sel <- msig$MODULES$Category %in% c("H","C5","C2")

xlsx = read.xlsx(paste(dirname(expr), "annotation_results_by_stage.xlsx", sep = "/"), rowNames = T)
xlsx = xlsx[xlsx$gene %in% genes[, 1],]

# file.remove(args[1])



# function to read and format the normalized_counts from sctransform
read_sctransform <- function(path="paga/normalized_counts.csv.gz") {
    print(path)
    r = gzfile(path)
    data = read.csv(r, row.names = 1)
    
    colnames(data) = gsub("\\.", "-", colnames(data), perl=F)
    colnames(data) = gsub("^X", "", colnames(data), perl=T)
    return(data)
}

# expr = read_sctransform(expr)
obj = readRDS(expr)

E = as.matrix(obj@scale.data)

# get design
stagesToNum <- c("I"=1, "II"=2, "III"=3, "IV"=4)

design <- obj@meta.data[, "Stage", drop = F]
design$Stage = sapply(design$Stage, function(x){stagesToNum[x]})

fit <- eBayes(lmFit(E, design))


# get msd
x <- tmodLimmaTopTable(fit, coef="Stage")

msd = x[rev(order(x$msd.Stage)), ]
temp_genes = rownames(x)[rownames(x) %in% genes[, 1]]

print("msd")
res = tmodUtest(temp_genes, mset=msig[sel], qval = Inf)
res$ident = "Msig"

temp = tmodUtest(temp_genes, mset=goset, qval = Inf)
temp$ident = "GO"
res = rbind(res, temp)

temp = tmodUtest(temp_genes, mset=kegg, qval = Inf)

if (!is.null(temp) && nrow(temp) > 0) {
    temp$ident = "KEGG"
    res = rbind(res, temp)
}


print("logfc")
res_logfc = tmodUtest(unique(xlsx[rev(order(xlsx$avg_logFC)), "gene"]), mset=msig[sel], qval=Inf)

if (!is.null(res_logfc) && nrow(res_logfc) > 0) {
    res_logfc$ident = "Msig"
}


temp = tmodUtest(unique(xlsx[rev(order(xlsx$avg_logFC)), "gene"]), mset=goset, qval = Inf)
if (!is.null(temp) && nrow(temp) > 0) {
    temp$ident = "GO"
    res_logfc = rbind(res_logfc, temp)
}

temp = tmodUtest(unique(xlsx[rev(order(xlsx$avg_logFC)), "gene"]), mset=kegg, qval = Inf)
if (!is.null(temp) && nrow(temp) > 0) {
    temp$ident = "KEGG"
    res_logfc = rbind(res_logfc, temp)
}

print("adj_p")
res_adj_p = tmodUtest(unique(xlsx[order(xlsx$p_val_adj), "gene"]), mset=msig[sel], qval=Inf)
if (!is.null(res_adj_p) && nrow(res_adj_p) > 0) {
    res_adj_p$ident = "Msig"
}

temp = tmodUtest(unique(xlsx[order(xlsx$p_val_adj), "gene"]), mset=goset, qval = Inf)
if (!is.null(temp) && nrow(temp) > 0) {
    temp$ident = "GO"
    res_adj_p = rbind(res_adj_p, temp)
}

temp = tmodUtest(unique(xlsx[order(xlsx$p_val_adj), "gene"]), mset=kegg, qval = Inf)
if (!is.null(temp) && nrow(temp) > 0) {
    temp$ident = "KEGG"
    res_adj_p = rbind(res_adj_p, temp)
}


print("pca")
x = E[xlsx$gene, ]
pca <- prcomp(t(x), scale.=T)
eigen <- pca$x[,1]
cors <- t(cor(eigen, t(E)))
ord <- order(abs(cors), decreasing=TRUE)

res_pc = tmodUtest(xlsx$gene[ ord ], mset = msig)
if (!is.null(res_pc) && nrow(res_pc) > 0) {
    res_pc$ident = "Msig"
}

temp = tmodUtest(xlsx$gene[ ord ], mset=goset, qval = Inf)
if (!is.null(temp) && nrow(temp) > 0) {
    temp$ident = "GO"
    res_pc = rbind(res_pc, temp)
}

temp = tmodUtest(xlsx$gene[ ord ], mset=kegg, qval = Inf)
if (!is.null(temp) && nrow(temp) > 0) {
    temp$ident = "KEGG"
    res_pc = rbind(res_pc, temp)
}


wb = createWorkbook()
addWorksheet(wb, "gene_module")
writeData(wb, 1, genes)

addWorksheet(wb, "msd")
writeData(wb, 2, res[res$N1 > 0, ])

addWorksheet(wb, "avg_logFC")
writeData(wb, 3, res_logfc[res_logfc$N1 > 0, ])

addWorksheet(wb, "p_val_adj")
writeData(wb, 4, res_adj_p[res_adj_p$N1 > 0, ])


addWorksheet(wb, "PC")
writeData(wb, 5, res_pc[res_pc$N1 > 0, ])

saveWorkbook(wb, output_file, overwrite = T)

