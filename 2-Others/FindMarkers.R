library(Seurat)


## AD vs SC
markers = NULL
for (i in unique(meta$cell_short)) {
    print(i)
    
    ident1 = meta$Cells[meta$cell_short == i & meta$Disease == "LUAD"]
    ident2 = meta$Cells[meta$cell_short == i & meta$Disease == "LUSC"]
    
    if (length(ident1) > 3 && length(ident2) > 3) {
        temp = FindMarkers(
            obj, 
            ident.1 = ident1,
            ident.2 = ident2,
            logfc.threshold = 0
        )
        
        temp$Cell = i
        temp$gene = rownames(temp)
        
        markers = rbind(markers, temp)
    }
}



markers = NULL
for (i in unique(meta$cell_short)) {
    print(i)
    
    ident1 = meta$Cells[meta$cell_short == i & meta$Disease == "LUAD"]
    ident2 = meta$Cells[meta$cell_short == i & meta$Disease == "LUAD_Normal"]
    
    if (length(ident1) > 3 && length(ident2) >= 3) {
        temp = FindMarkers(
            obj, 
            ident.1 = ident1,
            ident.2 = ident2,
            logfc.threshold = 0,
            min.cells.group = 0
        )
        
        temp$Cell = i
        temp$gene = rownames(temp)
        
        markers = rbind(markers, temp)
    }
}



markers = NULL
for (i in unique(meta$cell_short)) {
    print(i)
    
    ident1 = meta$Cells[meta$cell_short == i & meta$Disease == "LUSC"]
    ident2 = meta$Cells[meta$cell_short == i & meta$Disease == "LUSC_Normal"]
    
    if (length(ident1) > 3 && length(ident2) > 0) {
        temp = FindMarkers(
            obj, 
            ident.1 = ident1,
            ident.2 = ident2,
            min.cells.group = 0,
            logfc.threshold = 0
        )
        
        temp$Cell = i
        temp$gene = rownames(temp)
        
        markers = rbind(markers, temp)
    }
}