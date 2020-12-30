library(Seurat)
library(Matrix)
library(dplyr)

#args <- commandArgs(T)
#if(length(args) != 4 ){stop("Usage: Rscript UCB_Seurat_analysis.R  <ExpMat.txt>  <Gene_info.txt>  <batch_list>  <Outdir>")}
#setwd("RC")

# *** 0. Load the data ***
mat_all = read.table( "RC_mat.rmSuspB.txt",header = T,row.names = 1)

## Convert GeneID to GeneName in matrix
gname = read.table( "Gene_info.txt",header=T,row.names=1,sep = "\t")
gname = gname[!grepl("ERCC-",rownames(gname)), ]
#gname = gname[grepl("protein_coding|^TR_|^IG_",gname$Biotype), ]
gname = gname[grepl("protein_coding",gname$Biotype), ]
gname = gname[!duplicated(gname$GeneName), ]
flag = rownames(mat_all) %in% rownames(gname)
mat_all = mat_all[flag,]
rownames(mat_all) <- gname[rownames(mat_all),1]

## Split matrix by each batch
batch = read.table( "batch.rmSuspB.txt", header=F, stringsAsFactors=F, sep = "\t")
b1 = batch[batch$V2==1, "V1"] 
b2 = batch[batch$V2==2, "V1"] 
b3 = batch[batch$V2==3, "V1"] 
b4 = batch[batch$V2==4, "V1"] 
mat_b1 = mat_all[, b1]
mat_b2 = mat_all[, b2]
mat_b3 = mat_all[, b3]
mat_b4 = mat_all[, b4]

# *** 1. Setup the Seurat objects ***
s1 = CreateSeuratObject( raw.data=mat_b1, min.cells = 3, min.genes=500, is.expr=1, project="Batch1")
s2 = CreateSeuratObject( raw.data=mat_b2, min.cells = 3, min.genes=500, is.expr=1, project="Batch2")
s3 = CreateSeuratObject( raw.data=mat_b3, min.cells = 3, min.genes=500, is.expr=1, project="Batch3")
s4 = CreateSeuratObject( raw.data=mat_b4, min.cells = 3, min.genes=500, is.expr=1, project="Batch4")
s1@meta.data$batch = "UCB1"
s2@meta.data$batch = "UCB3"
s3@meta.data$batch = "UCB3S"
s4@meta.data$batch = "UCB4"

## Normalization and Scaling 
s1 = NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = 1000000)
s2 = NormalizeData(s2, normalization.method = "LogNormalize", scale.factor = 1000000)
s3 = NormalizeData(s3, normalization.method = "LogNormalize", scale.factor = 1000000)
s4 = NormalizeData(s4, normalization.method = "LogNormalize", scale.factor = 1000000)
s1 = ScaleData(s1)
s2 = ScaleData(s2)
s3 = ScaleData(s3)
s4 = ScaleData(s4)

## Determine genes to use for CCA, must be highly variable in at least 2 datasets
s1 = FindVariableGenes(s1, do.plot = F)
s2 = FindVariableGenes(s2, do.plot = F)
s3 = FindVariableGenes(s3, do.plot = F)
s4 = FindVariableGenes(s4, do.plot = F)
ob.list = list( s1,s2,s3,s4 )
genes.use = c()
for (i in 1:length(ob.list)) {
  genes.use = c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 2000))
}
genes.use = names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use = genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}
length(genes.use) #1958 genes


# *** 2. Perform a canonical correlation analysis (CCA) ***
Combine = RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 30)
pdf("Raw_CC.pdf",18,8)
p1 = DimPlot(object = Combine, reduction.use = "cca", group.by = "batch", pt.size = 0.5, do.return = TRUE)
p2 = VlnPlot(object = Combine, features.plot = "CC1", group.by = "batch", do.return = TRUE)
plot_grid(p1, p2)
dev.off()

## CC Selection
pdf("CC_selection1.pdf",10,8)
MetageneBicorPlot(Combine, grouping.var = "batch", dims.eval = 1:30)  #1-20
dev.off()
pdf("CC_selection2.pdf",10,20)
DimHeatmap(object = Combine, reduction.type = "cca", cells.use = 500, dim.use = 1:30, do.balanced = TRUE)  #1-20
dev.off()

## Run rare non-overlapping filtering
Combine = CalcVarExpRatio(object = Combine, reduction.type = "pca", grouping.var = "batch", dims.use = 1:20)
Combine = SubsetData(object = Combine, subset.name = "var.ratio.pca", accept.low = 0.5)


# *** 3. Align the CCA subspaces ***
Combine = AlignSubspace(object = Combine, reduction.type = "cca", grouping.var = "batch", dims.align = 1:20)
pdf("Align_CC.pdf",10,8)
p1 = VlnPlot(object = Combine, features.plot = "ACC1", group.by = "batch", do.return = TRUE)
p2 = VlnPlot(object = Combine, features.plot = "ACC2", group.by = "batch", do.return = TRUE)
plot_grid(p1, p2)
dev.off()


# *** 4. t-SNE and Clustering ***
Combine = FindClusters(Combine, reduction.type = "cca.aligned", dims.use = 1:20, save.SNN = T, resolution = 0.8)
Combine = RunTSNE(Combine, reduction.use = "cca.aligned", dims.use = 1:20)
pdf("tSNE.pdf",10,8)
TSNEPlot(Combine, do.return = T, pt.size = 1.5, group.by = "batch")
TSNEPlot(Combine, do.label = T, do.return = T, pt.size = 1.5)
dev.off()

save(Combine, file = "save.Rda")
#load("save.Rda")
write.table(Combine@dr$tsne@cell.embeddings, file = "tSNE_pos.txt",sep = "\t", quote = F, col.names = NA, row.names=T)

# *** 5. Known Markers Plots ***
known_markers <- c("GYPA", "ALAS2", "HBB", "MS4A1","CD79A","CD3D","CD44","CD4","CD8A", "FOXP3", "NKG7", "LYZ", "PTPRC", "PF4", "CD34", "ENG", "NT5E", "ITGB1","MME","GP9")
known_markers <- known_markers[ known_markers %in% rownames(Combine@raw.data) ]
# Erythrocyte: GYPA ALAS2 HBB
# B cell: MS4A1 CD79A
# T cell: CD3D CD44 CD4(Th) CD8A(Tc) FOXP3(Treg)
# NK: NKG7
# Monocyte: LYZ
# All leukocyte: PTPRC
# Platelet: PF4 GP9
# Hematopietic stem cell: CD34
# ? Embryonic stem cell: NT5E ENG ITGB1
pdf("Known_markers_heatmap.pdf",width = 10,height = 10)
DoHeatmap(Combine, genes.use = known_markers, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("Known_markers.pdf",width = 15,height = 16)
FeaturePlot(object = Combine, features.plot = known_markers, cols.use = c("grey", "blue"))
dev.off()


### Hemocytes
# Hematopoietic stem cells: CD34 PROM1(CD133)
# T cells: CD3E CD3D CD3G CD5 CD6 CD7 TRAV2 TRBV2 TRGV1 TRDV1
# B cells: CD19 MS4A1(CD20) CD79A
# NK cells: FCGR3A(CD16A) FCGR3B(CD16B) NCAM1(CD56) KLRD1(CD94) KIR2DL1(CD158A) KLRB1(CD161,NK1.1) NKG7
# Dendritic cells: ITGAX(CD11C) CD83 CD209 IL3RA(CD123) HLA-DRA 
# Monocyte/Macrophagocyte: CD14 CD33 FCGR1A(CD64) CD68 LYZ
# Granulocyte: ANPEP(CD13) FUT4(CD15) CD33
# Erythrocyte: SLC4A1(CD233) GYPA(CD235A) HBB ALAS2
# Megakaryocyte/Platelet: ITGA2B(CD41) GP1BA(CD42A) GP1BB(CD42B) ITGAV(CD51) ITGB3(CD61) PF4 GP9
HTS_markers = c("CD34", "PROM1")
T_markers = c("CD3E", "CD3D", "CD3G", "CD5", "CD6", "CD7")
B_markers = c("CD19", "MS4A1", "CD79A")
NK_markers = c("FCGR3A", "FCGR3B", "NCAM1", "KLRD1", "KIR2DL1", "KLRB1", "NKG7")
DC_markers = c("ITGAX", "CD83", "CD209", "IL3RA", "HLA-DRA")
Mono_markers = c("CD14", "CD33", "FCGR1A", "CD68", "LYZ")
Gran_markers = c("ANPEP", "FUT4", "CD33")
Eryt_markers = c("SLC4A1", "GYPA", "HBB", "ALAS2")
Mega_markers = c("ITGA2B", "GP1BA", "GP1BB", "ITGAV", "ITGB3", "PF4", "GP9")
All_Hemocytes_Markers = c(HTS_markers,T_markers,B_markers,NK_markers,DC_markers,Mono_markers,Gran_markers,Eryt_markers,Mega_markers)
All_Hemocytes_Markers = All_Hemocytes_Markers[ All_Hemocytes_Markers %in% rownames(Combine@raw.data) ]
pdf("Hemocyte_markers_heatmap.pdf",width = 15,height = 20)
DoHeatmap(Combine, genes.use = All_Hemocytes_Markers, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("Hemocyte_markers.pdf",width = 16,height = 32)
FeaturePlot(object = Combine, features.plot = All_Hemocytes_Markers, cols.use = c("grey", "blue") )
dev.off()

### T cells
# Th/Treg: CD4
# Tc: CD8A CD8B
# Th1: IFNG CCR5
# Th2: GATA3 CXCR4 CCR3
# Th17: KLRB1 CCR4 CCR6
# Treg: CTLA4 FOXP3 IL2RA
T_cd4cd8 = c("CD4", "CD8A", "CD8B")
Th1_markers = c("IFNG", "CCR5")
Th2_markers = c("GATA3", "CXCR4", "CCR3")
Th17_markers = c("KLRB1","CCR4","CCR6")
Treg_markers = c("CTLA4", "FOXP3", "IL2RA")
All_Tcells_Markers = c(T_cd4cd8, Th1_markers,Th2_markers,Th17_markers,Treg_markers)
All_Tcells_Markers = All_Tcells_Markers[ All_Tcells_Markers %in% rownames(Combine@raw.data) ]
pdf("Tcells_markers_heatmap.pdf",width = 10,height = 10)
DoHeatmap(Combine, genes.use = All_Tcells_Markers, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("Tcells_markers.pdf",width = 16,height = 14)
FeaturePlot(object = Combine, features.plot = All_Tcells_Markers, cols.use = c("grey", "blue") )
dev.off()


# *** 6. Annotation for cell type ***
current.cluster.ids <- c(0, 1, 2, 3, 4, 5)
new.cluster.ids <- c("CD4+ T Cells", "B Cells", "CD8+ T Cells", "NK cells", "Monocytes", "Megakaryocytes")
Combine@ident <- plyr::mapvalues(x = Combine@ident, from = current.cluster.ids, to = new.cluster.ids)
pdf("tSNE_anno.pdf",width = 10,height = 7)
TSNEPlot(object = Combine, do.label = TRUE, pt.size = 1.5, label.size=5)
dev.off()
write.table(Combine@ident,file = "Cell_type.txt",sep = "\t", quote = F, col.names = FALSE)
write.table(table(Combine@ident),file = "Cell_type_stat.txt",sep = "\t", quote = F, col.names = FALSE, row.names = FALSE)


# *** 7. Find New Markers ***
#C.marker0 = FindConservedMarkers( object = Combine, ident.1=0, grouping.var = "batch" )  #bad result
All.markers = FindAllMarkers(object = Combine, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)
All.top10.markers = All.markers  %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf("New_markers_heatmap.pdf",width = 15,height = 20)
DoHeatmap(Combine, genes.use = All.top10.markers$gene, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("New_markers.pdf",width = 16,height = 40)
FeaturePlot(Combine, features.plot = All.top10.markers$gene, cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
write.table(All.markers, file = "New_markers.txt",sep = "\t", quote = F, col.names = TRUE, row.names = FALSE)

T.markers = FindMarkers( object = Combine, ident.1=c(0,2), only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1 )
pdf("NewT_markers_heatmap.pdf",width = 15,height = 10)
DoHeatmap(Combine, genes.use = rownames(T.markers), slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("NewT_markers.pdf",width = 16,height = 22)
FeaturePlot(Combine, features.plot = rownames(T.markers), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
write.table(T.markers, file = "NewT_markers.txt",sep = "\t", quote = F, col.names = TRUE, row.names = TRUE)



########## CD45 isoforms
load("save.Rda")
cd45 = read.table("UCB_CD45.xls", header=T, row.names=1, sep="\t")
cd45 = cd45[names(Combine@ident), ]
Combine@meta.data = cbind( Combine@meta.data, cd45)
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="RO", colors=c("gray","red") )
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="RA", colors=c("gray","red") )
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="RB", colors=c("gray","red") )
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="RC", colors=c("gray","red") )
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="RAB", colors=c("gray","red") )
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="RBC", colors=c("gray","red") )
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="RAC", colors=c("gray","red") )
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="RABC", colors=c("gray","red") )
library(RColorBrewer)
color8 = brewer.pal(8,"Set1")
TSNEPlot(object = Combine, do.label = FALSE, pt.size = 1.5, group.by="Main", colors=c("gray",color8) )

######### AS number (this will destroy raw data)
#as = read.table("Raw_Cell_stat.addRead.newid.xls",header=T,row.names = 1)
#as = as[rownames(Combine@dr$pca@cell.embeddings),]
#Combine@dr$pca@cell.embeddings = cbind( Combine@dr$pca@cell.embeddings,  100 )
#Combine@dr$pca@cell.embeddings[,"PC1"] = log10(as$Total)
#FeaturePlot(object = Combine, features.plot = 'PC1', cols.use = c("green","red"), reduction.use = "tsne")

###### New signiture 
load("save.Rda")
current.cluster.ids <- c(0, 1, 2, 3, 4, 5)
new.cluster.ids <- c("CD4+ T Cells", "B Cells", "CD8+ T Cells", "NK cells", "Monocytes", "Megakaryocytes")
Combine@ident <- plyr::mapvalues(x = Combine@ident, from = current.cluster.ids, to = new.cluster.ids)

sign_Treg = c("NT5E","CD3D","CD3G","CD3E","CD4","CD5","ENTPD1","CTLA4","IZUMO1R","TNFRSF18","IL2RA","ITGAE","LAG3","TGFB1","LRRC32","TNFRSF4","SELL","FOXP3","STAT5A","STAT5B","LGALS1","IL10","IL12A","EBI3","TGFB1")
sign_CD8T_Act = c("CD69","CCR7","CD27","BTLA","CD40LG","IL2RA","CD3E","CD47","EOMES","GNLY","GZMA","GZMB","PRF1","IFNG","CD8A","CD8B","FASLG","LAMP1","LAG3","CTLA4","HLA-DRA","TNFRSF4","ICOS","TNFRSF9","TNFRSF18")
sign_T_TermDiff = c("TIGIT","PDCD1","CD274","CTLA4","LAG3","HAVCR2","CD244","CD160")
sign_G1_S = c("BRCA1","BRCA2","CCNE1","CCNE2","CCNG2","CDC25A","CDC45","CDC6","CDKN1A","CDKN2C","CDKN3","DHFR","E2F1","E2F5","H1F0","H1FNT","H1FOO","H1FX","H2AFB1","H2AFB2","H2AFB3","H2AFJ","H2AFV","H2AFVP1","H2AFX","H2AFY","H2AFY2","H2AFZ","H2AFZP1","H2AFZP2","H2AFZP3","H2AFZP4","H2AFZP5","H2AFZP6","H2BFM","H2BFS","H2BFWT","H2BFXP","H3F3A","H3F3AP1","H3F3AP2","H3F3B","H3F3C","HIST1H1A","HIST1H1B","HIST1H1C","HIST1H1D","HIST1H1E","HIST1H1PS1","HIST1H1PS2","HIST1H1T","HIST1H2AA","HIST1H2AB","HIST1H2AC","HIST1H2AD","HIST1H2AE","HIST1H2AG","HIST1H2AH","HIST1H2AI","HIST1H2AJ","HIST1H2AK","HIST1H2AL","HIST1H2AM","HIST1H2APS1","HIST1H2APS2","HIST1H2APS3","HIST1H2APS4","HIST1H2APS5","HIST1H2BA","HIST1H2BB","HIST1H2BC","HIST1H2BD","HIST1H2BE","HIST1H2BF","HIST1H2BG","HIST1H2BH","HIST1H2BI","HIST1H2BJ","HIST1H2BK","HIST1H2BL","HIST1H2BM","HIST1H2BN","HIST1H2BO","HIST1H2BPS1","HIST1H2BPS2","HIST1H3A","HIST1H3B","HIST1H3C","HIST1H3D","HIST1H3E","HIST1H3F","HIST1H3G","HIST1H3H","HIST1H3I","HIST1H3J","HIST1H3PS1","HIST1H4A","HIST1H4B","HIST1H4C","HIST1H4D","HIST1H4E","HIST1H4F","HIST1H4G","HIST1H4H","HIST1H4I","HIST1H4J","HIST1H4K","HIST1H4L","HIST1H4PS1","HIST2H2AA3","HIST2H2AA4","HIST2H2AB","HIST2H2AC","HIST2H2BA","HIST2H2BB","HIST2H2BC","HIST2H2BD","HIST2H2BE","HIST2H2BF","HIST2H3A","HIST2H3C","HIST2H3D","HIST2H3PS2","HIST2H4A","HIST2H4B","HIST3H2A","HIST3H2BA","HIST3H2BB","HIST3H3","HIST4H4","MCM2","MCM6","MSH2","NASP","NPAT","PCNA","RRM1","RRM2","SLBP","TYMS")
sign_G2_M = c("AURKA","BIRC5","BUB1","BUB1B","CCNA2","CCNB1","CCNB2","CCNF","CDC20","CDC25B","CDC25C","CDK1","CDKN2D","CENPA","CENPF","CKS2","KIF20A","PLK1","RACGAP1","TOP2A")

sign_Treg = sign_Treg[ sign_Treg %in% rownames(Combine@raw.data) ]
sign_CD8T_Act = sign_CD8T_Act[ sign_CD8T_Act %in% rownames(Combine@raw.data) ]
sign_T_TermDiff = sign_T_TermDiff[ sign_T_TermDiff %in% rownames(Combine@raw.data) ]
sign_G1_S = sign_G1_S[ sign_G1_S %in% rownames(Combine@raw.data) ]
sign_G2_M = sign_G2_M[ sign_G2_M %in% rownames(Combine@raw.data) ]

tcell=names(Combine@ident[ Combine@ident=="CD4+ T Cells"|Combine@ident=="CD8+ T Cells" ])

pdf("sign_Treg_heatmap.pdf",width = 10,height = 5)
DoHeatmap(Combine, cells.use=tcell, genes.use = sign_Treg, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("sign_CD8T_Act_heatmap.pdf",width = 10,height = 5)
DoHeatmap(Combine, cells.use=tcell, genes.use = sign_CD8T_Act, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("sign_T_TermDiff_heatmap.pdf",width = 10,height = 2)
DoHeatmap(Combine, cells.use=tcell, genes.use = sign_T_TermDiff, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("sign_G1_S_heatmap.pdf",width = 10,height = 28)
DoHeatmap(Combine, cells.use=tcell, genes.use = sign_G1_S, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
pdf("sign_G2_M_heatmap.pdf",width = 10,height = 5)
DoHeatmap(Combine, cells.use=tcell, genes.use = sign_G2_M, slim.col.label = TRUE, remove.key = FALSE)
dev.off()
