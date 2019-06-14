library(ggplot2)
library(openxlsx)
library(dplyr)


args = commandArgs(trailingOnly = T)

setwd(args[1])

read_data <- function(input_file, sheet=2) {
    data = read.xlsx(input_file, sheet = sheet)
    
    data["GeneRatio"] = sapply(data["GeneRatio"], function(x) {
        temp = strsplit(x, "/")[[1]]
        return (as.numeric(temp[1]) / as.numeric(temp[2]))
    })
    
    data["ident"] = sapply(data["ident"], as.character)
    
    return(data)
}


res = "results.xlsx"


kk = read_data(res, sheet = 2)

# if (nrow(kk) > 60) {
#     kk <- kk %>% group_by(ident) %>% top_n(-20, wt=p.adjust)
# }

p <- ggplot(
    data = kk, 
    aes(x=ident, y=Description, size=GeneRatio, color = p.adjust)
) + 
    geom_point() + 
    scale_color_gradientn(colours = rainbow(5)) +
    labs(x="", title = "KEGG")


height = nrow(kk) / 5
if (height < 5) {
    height = 5
} else if (height > 45) {
    height = 45
}

ggsave(
    filename = "kegg.png", 
    plot = p, 
    width = 12, 
    height = height, 
    dpi = 600, 
    units = "in",
    limitsize = F
)


go = read_data(res, sheet = 3)

for (i in unique(go$ONTOLOGY)) {
    temp <- go[go$ONTOLOGY == i, ]

    p <- ggplot(
        data = temp,
        aes(y=ident, x=Description, size=GeneRatio, color = p.adjust)
    ) +
        geom_point() +
        scale_color_gradientn(colours = rainbow(5)) +
        labs(y="", title = i) +
        # facet_wrap(~ONTOLOGY) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


    height = nrow(temp) / 5
    if (height < 5) {
        height = 5
    } else if (height > 45) {
        height = 45
    }

    ggsave(
        filename = paste0(i, ".png"),
        plot = p,
        width = height,
        height = 12,
        dpi = 600,
        units = "in",
        limitsize = F
    )
}


do = read_data(res, sheet = 4)

# if (nrow(do) > 60) {
#     do <- do %>% group_by(ident) %>% top_n(-20, wt=p.adjust)
# }

p <- ggplot(
    data = do, 
    aes(x=ident, y=Description, size=GeneRatio, color = p.adjust)
) + 
    geom_point() + 
    scale_color_gradientn(colours = rainbow(5)) +
    labs(x="", title = "DOSE")


height = nrow(do) / 5
if (height < 5) {
    height = 5
} else if (height > 45) {
    height = 45
}

ggsave(
    filename = "do.png", 
    plot = p, 
    width = 12, 
    height = height, 
    dpi = 600, 
    units = "in",
    limitsize = F
)


