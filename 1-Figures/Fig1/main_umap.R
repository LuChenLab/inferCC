options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(wesanderson)
library(scales)
library(ggrepel)
library(ggrastr)
library(ggthemes)

extrafont::loadfonts()

root.dir = "LungCancer10x/"

full.path <- function(...) { return(paste(root.dir, ..., sep = "/")) }


umap <- read.csv(full.path("02_rds/umap.csv"), row.names = 1, stringsAsFactors = F)
meta = readRDS("11_CNV/meta.rds")
# meta <- read.csv(full.path("02_rds/meta_after_singleR.csv"), row.names = 1, stringsAsFactors = F)
# 
# nm_meta <- read.xlsx(full.path("02_rds/NM_meta.xlsx"))
# rownames(nm_meta) <- nm_meta$SampleID
# 
# meta[meta$Batch == 3, "Stage"] <- nm_meta[meta$SampleID[meta$Batch == 3], "Stage"]
# meta[meta$Batch == 3, "Disease"] <- nm_meta[meta$SampleID[meta$Batch == 3], "Disease"]
# 
meta$Disease = sapply(meta$Disease, function(x) {
    disease = c(
        "LUAD"="LUAD",
        "LUAD_Normal"="NL(AD)",
        "LUSC"="LUSC",
        "LUSC_Normal"="NL(SC)",
        "Normal"="Normal",
        "Tumor"="Tumor",
        "LULC"="LULC",
        "LULC_Normal"="NL(LC)",
        "Pleiomorphic"="UPS",
        "Pleiomorphic Normal"="NL(UPS)"
    )

    disease[x]
})
meta$cell_short[meta$cell_name == "Mo"] = "Mφ"
# unique(meta$Disease)


### set colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


stage_colors = c("I"="#65A9A3", "II"="#4A933E", "III"="#EC7A21", "IV"="#D73F47", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#778793")
disease_colors = c(
    "LUAD" = "#0084D1", 
    "Normal (LUAD)"="#FECC1B", 
    "Normal"="#73BDFF", 
    "Normal (LUSC)"="#778793", 
    "LUSC"="#A0C807", 
    "Tumor"="#6A0019", 
    "CPI"="#C5000B",
    "NL(AD)"="#FECC1B", 
    "NL(AD"="#778793", 
    "Normal"="#73BDFF", 
    "NL(SC)"="#778793", 
    "LULC"="#91822B",
    "NL(LC)"="#EDCAB0",
    "UPS"="#EBAAA4",
    "NL(UPS)"="#B43018"
)
tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")
patient_colors = gg_color_hue(45)
names(patient_colors) <- c(
    sapply(1:23, function(x){sprintf("PA%02d", x)}),
    sapply(1:12, function(x){sprintf("PS%02d", x)}),
    sapply(1:10, function(x){sprintf("DL%01d", x)})
)


### Bar plots of immu cells issues
# immu_cells = c(
#     "B cells",
#     "CD4",
#     "CD8",
#     "DC",
#     "Dendritic",
#     "Exhaust T",
#     "Endothelial",
#     "Erythroid precursor",
#     "Granulocyte",
#     "Mast",
#     "Macrophages",
#     "Monocytes",
#     "NK",
#     "T_cells",
#     "Treg"
# )
# 
# short_names <- c(
#     "APII"="ATII",
#     "Alveolar II"="ATII",
#     "Basal"="Basal",
#     "B cells"="B",
#     "CD4"="CD4",
#     "CD8"="CD8",
#     "Ciliated"="Cilia",
#     "Club"="Club",
#     "DC"="DC",
#     "Dendritic"="DC",
#     "Endothelial"="EC",
#     "Epithelial"="Epi",
#     "Fibroblasts"="Fib",
#     "Granulocyte"="Gran",
#     "Mast"="Mast",
#     "Mascrophages"="Mφ",
#     "Macrophages"="Mφ",
#     "Monocytes"="Mo",
#     "Neuroendocrine"="NE",
#     "NK"="NK",
#     "Tregs"="Tregs"
# )
# 
# short_names <- as.data.frame(short_names)
# 
# temp_meta = meta[meta$cell_name %in% immu_cells, ]
# temp_meta = temp_meta[temp_meta$Batch != 3,]
# 
# unique(temp_meta$cell_name)


### Make umap plot

make_cell_types_plot3 <- function(
    coord,
    pt.size = 0.1, 
    alpha = 0.5,
    text_size = 10,
    axis_size = 18,
    axis_title_size = 20,
    title_size = 25,
    group_by = "",
    label = "cell_name",
    title = "Cell type",
    legend.position = "none",
    colors = NULL
) {
    data = coord[, c("UMAP1", "UMAP2", group_by, label)]
    
    if (group_by == "") {
        group_by = "cell_name"
    }
    
    p <- eval(
        parse(text=paste0(
                "ggplot(data = data, aes(x=UMAP1, y=UMAP2, color=", group_by, "))"
            )
        )
    )
    
    if (!is.null(label)) {

        if (label == "PatientID") {
            
            immu <- c( 
                "B", "CD4", "CD8",
                "DC", "Mast", "Mφ",
                "NK", "Tregs", "Gran"
            )
            
            temp = coord %>%
                group_by(cell_name, PatientID) %>%
                add_tally() %>%
                dplyr::select(cell_name, PatientID, n) %>%
                unique() %>%
                as.data.frame()
            
            temp$Immu = ifelse(temp$cell_name %in% immu, "Immu", "Not")
            
            temp1 = dcast(temp, PatientID~Immu, value.var = "n", fun.aggregate = mean, fill = 0)
            temp1$p = temp1$Immu / temp1$Not
            
            temp = rbind(
                temp[temp$cell_name %in% immu & temp$PatientID %in% temp1$PatientID[temp1$p >= 5], ],
                temp[!temp$cell_name %in% immu & temp$PatientID %in% temp1$PatientID[temp1$p < 5], ]
            ) %>%
                group_by(PatientID, Immu) %>%
                
            
                group_by(PatientID) %>%
                top_n(1, wt = n) %>%
                as.data.frame()
            
            temp$clt = paste(temp$cell_name, temp$PatientID)
            
            text_loc = coord[paste(coord$cell_name, coord$PatientID) %in% temp$clt, ] %>%
                group_by(PatientID) %>%
                mutate(x = median(UMAP1), y = median(UMAP2)) %>%
                dplyr::select(x, y, PatientID) %>%
                unique() %>%
                as.data.frame()
            
            
        } else {
            text_loc = eval(
                parse(text=paste0("data %>% group_by(", label, ") %>% mutate(x = median(UMAP1), y = median(UMAP2)) %>% dplyr::select(x, y, ", label, ")"))
            )
            
            text_loc = text_loc %>% unique() %>% as.data.frame()
        }
    }
    
    p <- p + geom_point_rast(size = pt.size, alpha = alpha)
    
    if (!is.null(label)) {
        p <- p + eval(parse(text=paste0("geom_text_repel(data = text_loc,  aes(x=x, y=y, label = ", label, "), color = \"black\", size = text_size, alpha = 1)")))
    }
    
    p <- p + 
        theme_bw(base_family = "Arial Unicode MS") +
        labs(
            title = title
        ) +
        theme(
            plot.title = element_text(hjust = 0.5, face="bold", size = title_size),
            legend.position = legend.position,
            axis.text = element_text(size = axis_size),
            axis.title = element_text(size = axis_title_size),
            panel.grid = element_blank()
        )
    
    if (!is.null(colors)) {
        p <- p + scale_color_manual(values = colors)
    }
    p
}


#### prepare data
coord = read.csv("02_rds/umap.csv", row.names = 1)
colnames(coord) <- c("UMAP1", "UMAP2")

coord$Disease = meta[rownames(coord), "Disease"]
coord$cell_name = meta[rownames(coord), "cell_short"]
coord$PatientID = meta[rownames(coord), "PatientID"]


##### cell types
cell_colors = c(
    '#e6194B', '#3cb44b', '#ffe119', 
    '#4363d8', '#f58231', '#911eb4', 
    '#42d4f4', '#f032e6', '#bfef45', 
    '#fabebe', '#469990', '#e6beff', 
    '#9A6324', '#fffac8', '#800000',
    '#aaffc3', '#808000', '#ffd8b1',
    '#000075', '#a9a9a9', '#ffffff', 
    '#000000'
)
set.seed(1)
# cell_colors = wes_palette("Zissou1", length(unique(meta$cell_name1)), type = "continuous")

cell_colors = tableau_color_pal(palette = "Tableau 20")(length(unique(meta$cell_short)))

# cell_colors = cell_colors[1:length(unique(meta$cell_name1))]
names(cell_colors) <- sample(as.character(unique(meta$cell_short)))

p <- make_cell_types_plot3(
    coord, 
    group_by="cell_name", 
    title_size = 25, 
    text_size = 8,
    colors = cell_colors
)


p

cairo_pdf(filename = full.path("01_first_plots/cell_types_umap.pdf"), width = 6, height = 6)
print(p)
dev.off()





##### patients
# coord1 = rbind(
#     coord[str_detect(coord$PatientID, "^P(A|S)"),],
#     coord[str_detect(coord$PatientID, "^DL"),]
# )


immu <- c( 
    "B", "CD4", "CD8",
    "DC", "Mast", "Mφ",
    "NK", "Tregs", "Gran"
)


temp = coord %>%
    group_by(cell_name, PatientID) %>%
    add_tally() %>%
    dplyr::select(cell_name, PatientID, n) %>%
    unique() %>%
    # group_by(PatientID) %>%
    # top_n(1, wt = n) %>%
    as.data.frame()

temp = dcast(temp, PatientID~cell_name, value.var = "n", fun.aggregate = mean, fill = 0)


coord$PatientID <- meta[rownames(coord), "PatientID"]

p <- make_cell_types_plot3(
    coord, 
    title_size = 25, 
    text_size = 6, 
    group_by="PatientID", 
    label = "PatientID", 
    title = "Patients",
    colors = patient_colors
)

p

ggsave(
    filename = full.path("01_first_plots/patient_umap.pdf"),
    plot = p, width = 6, height = 6, device = cairo_pdf
)



##### Disease
coord$Disease = meta[rownames(coord), "Disease"]

coord <- coord[rownames(meta[order(meta$Batch, decreasing = T), ]), ]


p <- make_cell_types_plot3(
    coord, 
    title_size = 25, 
    text_size = 8, 
    group_by = "Disease",
    label = NULL, 
    title = "Disease",
    legend.position = c(0.5, 0.05),
    colors = disease_colors
)

p <- p + theme(
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    legend.spacing = unit(0.05, 'cm'),
    # legend.background = element_blank()
) + labs(color = "") +
    guides(
        color = guide_legend(
            override.aes = list(size = 5, fill=NA),
            nrow = 2
        )
    ) + ylim(-15, 15)

#                             

ggsave(
    filename = full.path("01_first_plots/disease_umap.pdf"),
    plot = p, width = 6, height = 6, device = cairo_pdf
)


##### Transcripts counts
coord$nUMI = log10(meta[rownames(coord), "nUMI"] + 1)


p <- make_cell_types_plot3(
    coord, 
    title_size = 25, 
    text_size = 8, 
    group_by = "nUMI",
    label = NULL, 
    title = "Transcripts counts",
    legend.position = c(0.8, 0.05)
)

p <- p + theme(
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
) + 
    labs(color = "log10(counts + 1)") +
    scale_colour_gradient(
        low = "lightgrey",
        high = "blue"
    ) + 
    guides(color = guide_colorbar(
        title.position = "bottom",
        barwidth = 10
    ))

# p

ggsave(
    filename = full.path("01_first_plots/transcript_umap.pdf"),
    plot = p, width = 6, height = 6, device = cairo_pdf
)


### 
aucell <- readRDS("02_rds/02_aucell/cells_AUC.rds")
# data <- readRDS("02_rds/02_aucell/cells_assignment.rds")
# selectedThresholds <- getThresholdSelected(data)
temp = as.data.frame(t(as.data.frame(aucell@assays$data$AUC)))

short_names <- c(
    "APII"="AT2",
    "Alveolar II"="AT2",
    "Basal"="Basal",
    "B cells"="B",
    "CD4+"="CD4",
    "CD4"="CD4",
    "CD8+"="CD8",
    "CD8"="CD8",
    "Ciliated"="Cilia",
    "Club"="Club",
    "DC"="DC",
    "Dendritic"="DC",
    "Endothelial"="EC",
    "Epithelial"="Epi",
    # "Exhaust T"="Exh T",
    "Fibroblasts"="Fib",
    "Granulocyte"="Gran",
    "Mast"="Mast",
    "Mascrophages"="Mφ",
    "Macrophages"="Mφ",
    "Monocytes"="Mφ",  # rename to macrophages
    "Neuroendocrine"="NE",
    "NK"="NK",
    "Tregs"="Tregs",
    "Treg"="Tregs"
)

auc_thr <- c(
    "AT2"=0.18, "B"=0.1, "Basal"=0.23, "CD4"=0.2, "CD8"=0.18, "Cilia"=0.07,
    "Mast"=0.09, "NE"=0.06, "Tregs"=0.1, "NK"=0.1, "Mφ"=0.09, "Gran"=0.14, 
    "Fib"=0.1, "Epi"=0.38, "EC"=0.1, "DC"=0.1, "Club"=0.2
)

for (i in colnames(temp)) {
    
    if (!i %in% names(short_names) || !short_names[[i]] %in% names(auc_thr)) {
        next
    }
    
    print(i)
    coord$AUC = temp[rownames(coord), i]
    
    thr = auc_thr[[short_names[[i]]]]
    
    coord$AUC[coord$AUC < thr] = 0
    
    p <- make_cell_types_plot3(
        coord, 
        title_size = 25, 
        text_size = 8, 
        group_by = "AUC",
        label = NULL, 
        title = i,
        legend.position = "none"
    ) 
    
    p <- p + theme(
            legend.direction = "horizontal",
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)
        ) + 
        labs(title = paste0(short_names[[i]], " (AUC >= ", thr, ")")) +
        scale_colour_gradient(
            low = "lightgrey",
            high = "blue"
        ) + 
        guides(color = guide_colorbar(
            title.position = "bottom",
            barwidth = 10
        ))
    
    ggsave(
        filename = paste0("01_first_plots/aucell/", paste0(short_names[[i]], ".pdf")),
        plot = p, width = 6, height = 6, device = cairo_pdf
    )
}