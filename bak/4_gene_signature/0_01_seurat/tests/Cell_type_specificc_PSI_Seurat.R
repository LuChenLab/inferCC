
load(file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/All_Cell_specific_SJ_intron_centric_PSI2.RData")

library(Seurat)
library(dplyr)
library(Matrix)

psi_sub <- te_sub[!duplicated(te_sub), ]
psi_sub[is.na(psi_sub)] <- 0
# 1. Create Seurat object
ucb <- CreateSeuratObject(raw.data = psi_sub, min.cells = 3, min.genes=10, is.expr = 0, project="UCB_cluster")

# 3. Normalizing the data
ucb <- NormalizeData(object = ucb)

# 4. Detection of variable genes across the single cells
ucb <- FindVariableGenes(object = ucb, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0, x.high.cutoff = 4.5, y.cutoff = -5)
length(x = ucb@var.genes)

# 5. Scaling the data
ucb <- ScaleData(ucb)  #add 'do.cpp=F' if an error occurs

# 6. Perform linear dimensional reduction (PCA)
ucb <- RunPCA(ucb, pc.genes = ucb@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = ucb, pcs.use = 1:9)
PCAPlot(ucb, dim.1 = 1, dim.2 = 2)
ucb = ProjectPCA(object = ucb, do.print = FALSE)
PCHeatmap(object = ucb, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
# good: PC1-11

# 7. Determine statistically significant principal components
ucb <- JackStraw(object = ucb, num.replicate = 100)
JackStrawPlot(object = ucb, PCs = 1:20)
# good: PC1-20
PCElbowPlot(object = ucb)
## good: PC1-10

# 8. Cluster the cells
ucb <- FindClusters(ucb, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE, temp.file.location="./")
PrintFindClustersParams(ucb)

# 9. Run Non-linear dimensional reduction (tSNE)
ucb <- RunTSNE(object = ucb, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
TSNEPlot(object = ucb, do.label = FALSE, pt.size = 1.5, label.size=6)

load("/mnt/data5/BGI/UCB/ExpMat_NewID/Seurat.RData")
tsne <- as.data.frame(ucb@dr$tsne@cell.embeddings)
dim(tsne)
dim(Cell_type)
tsne$CellType <- as.factor(Cell_type$CellType)
library(ggplot2)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()

#pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Tsne_Cell_specific_SJ.pdf")
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne.pdf")
TSNEPlot(object = ucb, do.label = TRUE, pt.size = 1.5, label.size=6)
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, colour = CellType))+
  geom_point()
dev.off()


# . Finding differentially expressed genes (cluster biomarkers)
## If an error occurs, please exit and re-load the saved data
ucb.markers <- FindAllMarkers(object = ucb, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
ucb.markers <- ucb.markers[ucb.markers$p_val_adj<0.05,]
ucb.markers %>% group_by(cluster) %>% top_n(10) -> top10
ucb.markers %>% group_by(cluster) %>% top_n(2) -> top2
write.table(ucb.markers,file ="/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_New_markers.txt",sep = "\t",quote = F)
#pdf("New_markers_heatmap.pdf",width = 15,height = 20)
DoHeatmap(ucb, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
DoHeatmap(ucb, genes.use = top2$gene, slim.col.label = TRUE, remove.key = FALSE)
#dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_New_markers.pdf",width = 20,height = 20)
FeaturePlot(object = ucb, features.plot = top2$gene, cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()



pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_New_markers_cluster8.pdf",width = 20,height = 10)
FeaturePlot(object = ucb, features.plot = with(ucb.markers,gene[cluster==8])[1:12], cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()

save(ucb, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_Seurat_final.Rda")














#### 
head(ucb.markers)
load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
dim(sj[sj$sj %in% ucb.markers[ucb.markers$cluster==8,]$gene,c("SJ_start_exon", "SJ_end_exon")])
[1] 30  2
gene_cluster8 <- union(sj[sj$sj %in% ucb.markers[ucb.markers$cluster==8,]$gene,"SJ_start_exon"], sj[sj$sj %in% ucb.markers[ucb.markers$cluster==8,]$gene, "SJ_end_exon"])

library(org.Hs.eg.db)##### data base target to human-geneIDs
library(AnnotationDbi)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)

columns(org.Hs.eg.db)
mapIds(org.Hs.eg.db, keys = gene_cluster8, column = "ENTREZID", keytype = "ENSEMBL")
mapIds(org.Hs.eg.db, keys = gene_cluster8, column = "SYMBOL", keytype = "ENSEMBL")
#ENSG00000116824 ENSG00000136051 ENSG00000110876 ENSG00000184007 ENSG00000189067
#          "CD2"        "WASHC4"        "SELPLG"        "PTP4A2"         "LITAF"
#ENSG00000162434 ENSG00000141293 ENSG00000142875 ENSG00000100100 ENSG00000115687
#         "JAK1"         "SKAP1"        "PRKACB"       "PIK3IP1"          "PASK"
#ENSG00000163519 ENSG00000113810 ENSG00000138795 ENSG00000074966 ENSG00000131507
#        "TRAT1"          "SMC4"          "LEF1"           "TXK"        "NDFIP1"
#ENSG00000168685 ENSG00000068796 ENSG00000133872 ENSG00000078668 ENSG00000135046
#         "IL7R"         "KIF2A"         "SARAF"         "VDAC3"         "ANXA1"

BP_cluster8 <- enrichGO(gene = gene_cluster8,
                        OrgDb  = org.Hs.eg.db,
                        keyType= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

print(dotplot(BP_cluster8,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(BP_cluster8,firstSigNodes = 15, useInfo = "all"))
enrichMap(BP_cluster8,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

BP_cluster8 <- setReadable(BP_cluster8, OrgDb = org.Hs.eg.db)
gofilter(BP_cluster8, level = 3)@result

gene8 <- as.character(mapIds(org.Hs.eg.db, keys = gene_cluster8, column = "ENTREZID", keytype = "ENSEMBL"))
kk <- enrichKEGG(gene         = gene8,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
nrow(kk@result)
[1] 0

mkk <- enrichMKEGG(gene 	= gene8,
                   organism = 'hsa')
nrow(mkk@result)
[1] 0




gene_cluster7 <- union(sj[sj$sj %in% ucb.markers[ucb.markers$cluster==7,]$gene,"SJ_start_exon"], sj[sj$sj %in% ucb.markers[ucb.markers$cluster==7,]$gene, "SJ_end_exon"])
gene_cluster7 <- as.character(na.omit(unlist(strsplit(gene_cluster7,split="\\|"))))

BP_cluster7 <- enrichGO(gene = gene_cluster7,
                        OrgDb  = org.Hs.eg.db,
                        keyType= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

print(dotplot(BP_cluster7,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(BP_cluster7,firstSigNodes = 15, useInfo = "all"))
enrichMap(BP_cluster7,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

BP_cluster7 <- setReadable(BP_cluster7, OrgDb = org.Hs.eg.db)
gofilter(BP_cluster7, level = 3)@result





gene_cluster3 <- union(sj[sj$sj %in% ucb.markers[ucb.markers$cluster==3,]$gene,"SJ_start_exon"], sj[sj$sj %in% ucb.markers[ucb.markers$cluster==3,]$gene, "SJ_end_exon"])
gene_cluster3 <- as.character(na.omit(unlist(strsplit(gene_cluster3,split="\\|"))))

BP_cluster3 <- enrichGO(gene = gene_cluster3,
                        OrgDb  = org.Hs.eg.db,
                        keyType= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

print(dotplot(BP_cluster3,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(BP_cluster3,firstSigNodes = 15, useInfo = "all"))
enrichMap(BP_cluster3,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

BP_cluster3 <- setReadable(BP_cluster3, OrgDb = org.Hs.eg.db)
gofilter(BP_cluster3, level = 3)@result




gene_cluster2 <- union(sj[sj$sj %in% ucb.markers[ucb.markers$cluster==2,]$gene,"SJ_start_exon"], sj[sj$sj %in% ucb.markers[ucb.markers$cluster==2,]$gene, "SJ_end_exon"])
gene_cluster2 <- as.character(na.omit(unlist(strsplit(gene_cluster2,split="\\|"))))

BP_cluster2 <- enrichGO(gene = gene_cluster2,
                        OrgDb  = org.Hs.eg.db,
                        keyType= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

print(dotplot(BP_cluster3,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(BP_cluster3,firstSigNodes = 15, useInfo = "all"))
enrichMap(BP_cluster3,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

BP_cluster2 <- setReadable(BP_cluster2, OrgDb = org.Hs.eg.db)
gofilter(BP_cluster2, level = 3)@result





gene_cluster10 <- union(sj[sj$sj %in% ucb.markers[ucb.markers$cluster==10,]$gene,"SJ_start_exon"], sj[sj$sj %in% ucb.markers[ucb.markers$cluster==10,]$gene, "SJ_end_exon"])
gene_cluster10 <- as.character(na.omit(unlist(strsplit(gene_cluster10,split="\\|"))))

BP_cluster10 <- enrichGO(gene = gene_cluster10,
                        OrgDb  = org.Hs.eg.db,
                        keyType= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

print(dotplot(BP_cluster10,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(BP_cluster10,firstSigNodes = 15, useInfo = "all"))
enrichMap(BP_cluster10,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

BP_cluster10 <- setReadable(BP_cluster10, OrgDb = org.Hs.eg.db)
gofilter(BP_cluster10, level = 3)@result


#### compareCluster enrichGO
genelist <- list(cluster2 = gene_cluster2, cluster3 = gene_cluster3, cluster7 = gene_cluster7, cluster8 = gene_cluster8, cluster10 = gene_cluster10)
GOcompare <- compareCluster(genelist, 
                      		fun="enrichGO", 
                      		OrgDb  = org.Hs.eg.db, 
                      		keyType= 'ENSEMBL', 
                      		ont  = "BP", 
                      		pAdjustMethod = "BH")


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_cluster378_compareCluster_enrichGO.pdf", height = 20, width = 20)
print(dotplot(GOcompare,showCategory = 10,font.size = 16, title = "top10"))
print(dotplot(GOcompare,showCategory = 20,font.size = 16, title = "top20"))
dev.off()



#### compareCluster enrichDO
genelist_entr <- list(cluster2 = as.character(mapIds(org.Hs.eg.db, keys = gene_cluster2, column = "ENTREZID", keytype = "ENSEMBL")),
					  cluster3 = as.character(mapIds(org.Hs.eg.db, keys = gene_cluster3, column = "ENTREZID", keytype = "ENSEMBL")), 
					  cluster7 = as.character(mapIds(org.Hs.eg.db, keys = gene_cluster7, column = "ENTREZID", keytype = "ENSEMBL")), 
					  cluster8 = as.character(mapIds(org.Hs.eg.db, keys = gene_cluster8, column = "ENTREZID", keytype = "ENSEMBL")),
					  cluster10 = as.character(mapIds(org.Hs.eg.db, keys = gene_cluster10, column = "ENTREZID", keytype = "ENSEMBL")))

DOcompare <- compareCluster(genelist_entr, 
							fun="enrichDO", 
							pAdjustMethod = "BH", 
							ont = "DO")

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_cluster378_compareCluster_enrichDO.pdf", height = 20, width = 20)
print(dotplot(DOcompare,showCategory = 10,font.size = 16, title = "top10"))
print(dotplot(DOcompare,showCategory = 20,font.size = 16, title = "top20"))
dev.off()


mapIds(org.Hs.eg.db, keys = gene_cluster3, column = "SYMBOL", keytype = "ENSEMBL")




#### heatmap of cluster 2,3,7,8,10 marker SJ PSI

with(ucb.markers, gene[cluster %in% c(2,3,7,8,10)])
load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
psi_sj <- psi_sj_same_start_table[with(ucb.markers, gene[cluster %in% c(2,3,7,8,10)]), ]
dim(psi_sj)
[1] 1393 3574

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

identical(row.names(as.data.frame(ucb@ident)),row.names(Cell_type))
[1] TRUE
Cell_type$Cluster <- ucb@ident
Cell_type <- Cell_type[order(Cell_type$Cluster),]

psi_sj[is.na(psi_sj)] <- 0

library(pheatmap)
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_cluster378_heatmap.pdf", height = 20, width = 20)
pheatmap(psi_sj[,row.names(Cell_type)], annotation_col = Cell_type, cluster_rows = F, cluster_cols = F, show_colnames = F, show_rownames = F)
dev.off()

library(RColorBrewer)
colors = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100)

pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_cluster378_heatmap2.pdf", height = 20, width = 20)
pheatmap(psi_sj[,row.names(Cell_type[Cell_type$Cluster %in% c(2,3,7,8,10),])], color = colors, annotation_col = Cell_type, cluster_rows = F, cluster_cols = F, show_colnames = F, show_rownames = F)
dev.off()



#### CD2

mapIds(org.Hs.eg.db, keys = gene_cluster8, column = "SYMBOL", keytype = "ENSEMBL")
# ENSG00000116824
which(grepl(sj$SJ_start_exon, pattern="ENSG00000116824") | grepl(sj$SJ_end_exon, pattern="ENSG00000116824"))
sj[which(grepl(sj$SJ_start_exon, pattern="ENSG00000116824") | grepl(sj$SJ_end_exon, pattern="ENSG00000116824")),]$sj

#load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
dim(psi_sj_same_start_table[row.names(psi_sj_same_start_table) %in% sj[which(grepl(sj$SJ_start_exon, pattern="ENSG00000116824") | grepl(sj$SJ_end_exon, pattern="ENSG00000116824")),]$sj,])
[1]   24 3574
dim(psi_sj[row.names(psi_sj) %in% sj[which(grepl(sj$SJ_start_exon, pattern="ENSG00000116824") | grepl(sj$SJ_end_exon, pattern="ENSG00000116824")),]$sj,])

CD2_sj <- sj[which(grepl(sj$SJ_start_exon, pattern="ENSG00000116824") | grepl(sj$SJ_end_exon, pattern="ENSG00000116824")),]$sj
CD2_psi <- psi_sj_same_start_table[row.names(psi_sj_same_start_table) %in% CD2_sj,]


pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_CD2_FeaturePlot.pdf",width = 10,height = 10)
FeaturePlot(object = ucb, features.plot = CD2_sj[CD2_sj %in% row.names(ucb@data)], cols.use = c("grey", "blue"), reduction.use = "tsne", no.legend = F)
dev.off()



pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_CD2_My_FeaturePlot.pdf", height = 6,width = 6)
for (i in 1:nrow(CD2_psi)){
	tmp_data <- data.frame(tsne, PSI = as.numeric(CD2_psi[i,row.names(tsne)]))
	tmp_data <- tmp_data[rev(order(tmp_data$PSI, decreasing=T)),]
	tmp_data[is.na(tmp_data$PSI),]$PSI <- 0
	tmp_data[tmp_data$PSI < 0.15,]$PSI <- 0
	print(ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
  		geom_point(aes(colour = PSI), size = 2)+
  		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "PSI"))+
  		ggtitle(row.names(CD2_psi)[i])+
  		theme(panel.background = element_blank())
)
}#statements
dev.off()


TPM <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/TPM_mat.txt", header = T, row.names = 1, stringsAsFactors = F)
TPM["ENSG00000116824",]

tmp_data <- data.frame(tsne, GE = log10(as.numeric(TPM["ENSG00000116824",row.names(tsne)])+1))
tmp_data <- tmp_data[rev(order(tmp_data$GE, decreasing=T)),]
#tmp_data[tmp_data$GE < 0.15,]$GE <- 0
pdf("/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_specific_PSI_tsne_CD2_Gene_Expression_My_FeaturePlot.pdf", height = 6,width = 6)
ggplot(data = tmp_data, aes(x = tSNE_1, y = tSNE_2))+
 		geom_point(aes(colour = GE), size = 2)+
 		scale_colour_gradient(low = "grey80", high = "#EF3B2C", na.value = "grey90", guide = guide_legend(title = "GE"))+
 		ggtitle("ENSG00000116824")+
 		theme(panel.background = element_blank())
dev.off()












