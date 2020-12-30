### Created by Zhang at 2019.01.23
### This script is used to format marker genes from Guoguoji


library(openxlsx)

data <- read.xlsx("/Users/zhangyiming/Library/Mobile Documents/com~apple~CloudDocs/gene list/mouse_cell_altas.xlsx")


format_cell_name <- function(data) {
    data <- data[4]
    return(strsplit(data, " ?\\(", perl = T)[[1]][1])
}


data$alias <- apply(data, 1, format_cell_name)


data <- data[regexpr("(lung|lymph|Bronichi)", data$Tissue, perl = TRUE, ignore.case = TRUE) != -1,]

### Only select the top 10 genes of each cells
top <- 10

res = data.frame(cell_name=NA, gene=NA)

for (i in unique(data$alias)) {
    print(i)
    tmp = data[data$alias == i, ]
    tmp = tmp[order(tmp$p_val), ]
    tmp = tmp[1:top, c("alias", "gene")]
    colnames(tmp) <- c("cell_name", "gene")
    
    res <- rbind(res, tmp)
}

res$gene <- toupper(res$gene)

res <- na.omit(res)


wb = loadWorkbook("/Users/zhangyiming/Library/Mobile Documents/com~apple~CloudDocs/gene list/mouse_cell_altas.xlsx")
# addWorksheet(wb, "top")
writeData(wb, 2, res)
saveWorkbook(wb, file = "/Users/zhangyiming/Library/Mobile Documents/com~apple~CloudDocs/gene list/mouse_cell_altas.xlsx", overwrite = T)
