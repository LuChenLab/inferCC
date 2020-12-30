## Extract all marker genes


load_pacages <- function() {
    if(!require("argparse")) {
        install.packages("argparse")
        
        library("argparse")
    }
    
    if(!require("Seurat")) {
        install.packages("Seurat")
        
        library("Seurat")
    }
    
    if(!require("openxlsx")) {
        install.packages("openxlsx")
        library("openxlsx")
    }
    
    if(!require("logging")) {
        install.packages("logging")
        
        library("logging")
    }
}

load_pacages()
basicConfig()
## all  "SCGB1A1", "KRT5",  "DNAI1", "FCER1A", "FCGR3A", "NGFR", "KRT5", "MUC5AC",  "P63", "IGTA6",  "KRT14"
## 2018jz10 "SCGB1A1", "KRT5", "DNAI1", "MUC5AC", "FCER1A", "FCGR3A", "NGFR", "KRT5"，
## 2018jz33  "SCGB1A1", "KRT5", "DNAI1", "FCER1A", "FCGR3A", "NGFR", "KRT5"
## 2018jz34 "SCGB1A1", "KRT5", "DNAI1", "FCER1A", "FCGR3A","NGFR", "KRT5", "MUC5AC",  "KRT14"
## 2018jz35  "SCGB1A1", "KRT5", "DNAI1", "FCER1A", "FCGR3A","NGFR", "KRT5", "MUC5AC",  "KRT14"
## 2018jz7  "SCGB1A1", "KRT5",  "FCER1A", "FCGR3A", "NGFR", "KRT5"
## 2018jz8 "SCGB1A1", "KRT5",  "DNAI1", "FCER1A", "FCGR3A", "NGFR", "KRT5"
## 2018jz4  "SCGB1A1", "KRT5",  "FCER1A", "FCGR3A", "NGFR", "KRT5", "KRT14"
## 2018jz5  "SCGB1A1", "KRT5",  "DNAI1", "FCER1A", "FCGR3A", "NGFR", "KRT5",   "KRT14"

# png("feature.png", width = 1200, height = 1200)
# FeaturePlot(
#     object = data, 
#     features.plot = c(
#        
#         ),
#     cols.use = c("grey", "blue"), 
#     reduction.use = "tsne"
# )
# dev.off()

# "MUC5AC", "P63", "IGTA6", "KRT14"

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument(
    "-i", 
    "--input", 
    help="Path to input directoy"
)

parser$add_argument(
    "-o", 
    "--output", 
    help="Path to output xlsx"
)

args = parser$parse_args()

input_dir = args$input
output_xlsx = args$output


# input_dir = "/Volumes/WD/Seurat/results/Seurat_results/cluster"

wb = createWorkbook()

sheet_index = 1
for (i in list.dirs(input_dir)) {
    if (i == input_dir) {
        next
    }
    
    loginfo(i)
    # print(i)
    
    rds <- paste(i, "2_clustered.rds", sep = "/")

    data <- readRDS(rds)
    
    markers <- FindAllMarkers(data, min.pct = 0.75)
    
    ws = addWorksheet(wb, sheetName = basename(i))
    writeData(wb, sheet_index, markers)
    
    sheet_index = sheet_index + 1
}

saveWorkbook(wb, file = output, overwrite = TRUE )


xlsx = "/Volumes/WD/Seurat/results/Seurat_results/gene.xlsx"


get_top_genes <- function(sheet, num=50) {
    genes = c()
    
    for (i in unique(sheet$cluster)) {
        tmp = sheet[sheet$cluster == i, ]
        
        tmp = tmp[order(tmp$p_val_adj), ]
        
        genes = c(genes, tmp$gene[1:num])
    }
    
    return(sheet[sheet$gene %in% genes, c("gene", "cluster")])
}

test = get_top_genes(data)

wb = loadWorkbook(xlsx)

sheet_names = sheets(wb)

res = NA
sheet_index = 1
for (i in sheet_names) {
    data = read.xlsx(xlsx, sheet = sheet_index)
    
    sheet_index = sheet_index + 1
    
    tmp = get_top_genes(data)
    
    colnames(tmp) <- c("gene", i)
    
    if (is.na(res)) {
        res = tmp
    } else {
        res = merge(res, tmp, by = "gene", all=TRUE)
    }
}


gene_matrix <- read.xlsx("/Volumes/WD/Seurat/results/Seurat_results/gene2.xlsx", rowNames = T)

guoguoji <- read.xlsx("/Volumes/WD/Seurat/results/Guoguiji/mouse_cell_altas.xlsx")

guoguoji <- guoguoji[regexpr("lung", guoguoji$Tissue, ignore.case = T, perl = T) != -1, ]
guoguoji <- guoguoji[, c("alias", "gene")]
guoguoji$gene <- toupper(guoguoji$gene)


data <- unique(merge(gene_matrix, guoguoji, by = "gene", all.x=T)[, c("gene", "alias")])

wb = loadWorkbook("/Volumes/WD/Seurat/results/Seurat_results/gene2.xlsx")
addWorksheet(wb, sheetName = "MCA")
writeData(wb, 2, data)
saveWorkbook(wb, "/Volumes/WD/Seurat/results/Seurat_results/gene2.xlsx", overwrite = T)
