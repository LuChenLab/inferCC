### Created by Zhang at 2019.01.23
### Plot velocyto

library(velocyto.R)
library(pagoda2)


ldat <- read.loom.matrices(input_loom)

emat <- ldat$spliced
mito.genes <- grep(pattern = "^mt-", x = rownames(emat), value = TRUE)
hba.genes <- grep(pattern = "^Hba", x = rownames(emat), value = TRUE)
hbb.genes <- grep(pattern = "^Hbb", x = rownames(emat), value = TRUE)
discard.genes <- c(mito.genes,hba.genes,hbb.genes)
# this dataset has already been pre-filtered, but this is where one would do some filtering
emat <- emat[,colSums(emat)>=1e3]
emat <- emat[unique(rownames(emat)),]
emat <- emat[!(rownames(emat) %in% discard.genes),]
r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
pdf('2018jz14.adjuct.pdf',12,8)
r$adjustVariance(plot=T,do.par=T,gam.k=10)
dev.off()

r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')
r$getKnnClusters(method=igraph::multilevel.community,type='PCA',name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)



pdf('2018jz14.tSNE.filter.pdf',12,8)
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"Xist"],main='Xist')  
dev.off()


emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,rownames(r$counts)]
nmat <- nmat[,rownames(r$counts)] # restrict to cells that passed p2 filter
# take cluster labels
cluster.label <- r$clusters$PCA[[1]]
cell.colors <- pagoda2:::fac2col(cluster.label)
# take embedding
emb <- r$embeddings$PCA$tSNE
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

pdf('2018jz14.trace.filter.pdf')
show.velocity.on.embedding.cor(
    emb,
    rvel.cd,
    n=300,
    scale='sqrt',
    cell.colors=ac(cell.colors,alpha=0.5),
    cex=0.8,
    arrow.scale=5,
    show.grid.flow=TRUE,
    min.grid.cell.mass=0.5,
    grid.n=40,
    arrow.lwd=1,
    do.par=F,
    cell.border.alpha = 0.1
)
dev.off()


gene <- "Olig2"
pdf('svzm.Olig2.pdf',15,8)
gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 20,kGenes=1,fit.quantile=fit.quantile,cell.emb=emb,cell.colors=cell.colors,cell.dist=cell.dist,show.gene=gene,old.fit=rvel.cd,do.par=T)
dev.off()

grep(pattern = "^Olig", x = rownames(nmat), value = TRUE)