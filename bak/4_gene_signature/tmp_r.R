
    random_state = 1
    data = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Alveolar_II/stage_gene_sig_cluster_1_sctransform/RENN.csv.gz"
    tmp_genes = c("/mnt/raid62/Lung_cancer_10x/scripts/tmp.txt")
    output = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Alveolar_II/stage_gene_sig_cluster_1_sctransform/RENN_WGCNA"
    rds = "/mnt/raid62/Lung_cancer_10x/02_figures_each_cell/Alveolar_II/Alveolar_II.rds"
    group_by = "res.0.6"
    
    library(WGCNA)
    library(ggplot2)
    library(Seurat)

    enableWGCNAThreads()
     
    set.seed(random_state)
            
    r = gzfile(data)
    data = read.csv(r, row.names = 1)
    
    genes = read.table(tmp_genes, header = F)
    genes = genes[, 1]
    # remove(tmp_genes)
    
    mtx = data[intersect(as.character(genes), rownames(data)),]
   
    sampleTree = hclust(dist(mtx), method = "average");
    png(file = paste0(output, "_sampleClustering.png"), width = 12, height = 9, res = 600, units = "in")
    par(cex = 0.6);
    par(mar = c(0,4,2,0))
    plot(sampleTree, 
         main = "Sample clustering to detect outliers", 
         sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
    dev.off()
    
    # 作图选择power
    # Call the network topology analysis function
    sft = pickSoftThreshold(t(mtx), verbose = 5)
    # Plot the results:
    
    power = sft$powerEstimate
    
    if (is.null(power) | is.na(power)) {
        power = 3
    } else {
        if (power < 1) {
            power = 1
        }
        if (power > 30) {
            power = 30
        }
        
        tryCatch({
            png(file = paste0(output, "_softthreshold.png"), width = 12, height = 9, res = 600, units = "in")
            par(mfrow = c(1,2));
            cex1 = 0.9;
            # Scale-free topology fit index as a function of the soft-thresholding power
            plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                 xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
                 main = paste("Scale independence"));
            text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                 labels=power,cex=cex1,col="red");
            # this line corresponds to using an R^2 cut-off of h
            abline(h=0.90,col="red")
            # Mean connectivity as a function of the soft-thresholding power
            plot(sft$fitIndices[,1], sft$fitIndices[,5],
                 xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
                 main = paste("Mean connectivity"))
            text(sft$fitIndices[,1], sft$fitIndices[,5], labels=sft$powerEstimate, cex=cex1,col="red")
            dev.off()  
        }, finally = {})
    }
    
    ## make network
    net = blockwiseModules(t(mtx), power = power,
                           TOMType = "unsigned", minModuleSize = 30,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = T,
                           saveTOMFileBase = paste0(output, "_wcgna"),
                           verbose = 3)
    
    
    res = as.data.frame(net$colors)
    write.table(res, file = paste0(output, "_data.txt"), quote=F, col.name=F, sep = "	")
    
    
    # open a graphics window
    png(file = paste0(output, "_clustering.png"), width = 12, height = 9, res = 600, units = "in")
    # Convert labels to colors for plotting
    mergedColors = labels2colors(net$colors)
    # Plot the dendrogram and the module colors underneath
    plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
                        
    dev.off()
    
    
    obj <- readRDS(rds)
    
    groups = unique(res[, 1])
    
    image_dir = paste(output, "heatmap", sep = "_")
    dir.create(image_dir, showWarnings = F)
    
    for(i in groups) {
        temp = res[res[, 1] == i, , drop=F]
        p <- DoHeatmap(
            obj, 
            group.by = group_by, 
            genes.use = rownames(temp), 
            cex.col=0,
            slim.col.label = TRUE, 
            remove.key = TRUE,
            do.plot = F
        )
        
        height = nrow(temp) / 8
        
        if(height < 5) {
            height = 5
        } else if(height > 40){
            height = 40
        }
        
        ggsave(filename = paste(image_dir, paste0(i, ".png"), sep = "/"),
               plot = p,
               width = 6,
               height = height,
               limitsize = F,
               units = "in",
               dpi = 600)
    }
    