
load_pacakges <- function() {
    library(Seurat)
    library(gridExtra)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(harmony)
    library(cowplot)
    library(ggthemes)
    library(dplyr)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(openxlsx)
    library(UBL)
    
    library(stackoverflow)
    library(scales)
    library(DOSE)
    library(data.table)
    library(RColorBrewer)

    library(SingleCellExperiment)
    library(tidyr)
    library(reticulate)
    library(stringr)
    library(ggrastr)
}


suppressPackageStartupMessages(load_pacakges())


custom_colors <- c(
    "#1B62A3", "#9CB8DD", "#EE6719", "#F8A966", "#278F36", 
    "#8CC873", "#C70F1F", "#EF7E82", "#7D509A", "#B59BC7",
    "#76423A", "#B48880", "#CB60A1", "#F1A3C5", "#7E7F83", "#BDC3C4", 
    "#ABB126", "#D1D57A", "#1DAFC3", "#8BCFDA"
)


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


### extract temp information
stage_colors = c("I"="#EC7A21", "II"="#D73F47", "III"="#65A9A3", "IV"="#4A933E", "LUAD_Normal"="#FECC1B", "LUSC_Normal"="#DA6906")
disease_colors = c("LUAD" = "#0084D1", "LUAD_Normal"="#FECC1B", "Normal"="#73BDFF", "LUSC_Normal"="#DA6906", "LUSC"="#A0C807", "LELC"="#6A0019", "CPI"="#C5000B")
tissue_colors = c("Bronichi"="#FB290F", "Lymph node"="#488F17", "Tumor"="#FC8210")
patient_colors = gg_color_hue(35)
names(patient_colors) <- c(
    sapply(1:23, function(x){sprintf("PA%02d", x)}),
    sapply(1:12, function(x){sprintf("PS%02d", x)})
)


colors = c(stage_colors, disease_colors, tissue_colors, patient_colors)



get_smote_match <- function(meta, num.cells=NULL) {
    # Why Am I using R
    # I hate R for real
    # 现在问题比较有意思，通过这个理论上确实是能划分
    # 统计每个Stage中其他数据转化为数字，然后通过他来算一下看看
    convert_int <- function(meta, column="") {
        # meta = apply(meta, 2, as.character)
        for (i in colnames(meta)) {
            if (i != column) {

                temp = sort(unique(as.character(meta[, i])))
                idx = 1

                for (j in temp) {
                    temp = as.character(meta[, i]) == j
                    meta[, i][temp] <- idx
                    idx = idx + 1
                }
            }
        }
        
        return(as.data.frame(meta))
    }
        
    meta1 = convert_int(meta[, c("PatientID", "Gender", "Age", "Disease", "Stage")])  
    # rownames(meta1) = rownames(meta)
    # meta1 = na.omit(meta1)

    temp = table(meta1$Stage)
    
    temp = summary(as.numeric(temp))
    
    if (is.null(num.cells)) {
        num.cells = max(100, temp[3])
    }
    
    test =  1 / (table(meta1$Stage) / num.cells)
    test[test > 1] = 1
    
    test = split(unname(test),names(test))

    for(i in 1:ncol(meta1)) {
        meta1[,i] <- as.factor(meta1[,i])
        print(unique(meta1[, i]))
    }

    newData <- SmoteClassif(Stage ~ PatientID + Disease + Gender + Age, meta1, dist = "HEOM", C.perc = test)
    
    return(newData)
}


dump_seurat_for_scanpy <- function(object, path) {
    
    dir.create(path, showWarnings = F)
    
    gz = gzfile(paste(path, "raw.csv.gz", sep = "/"), "w+")
    write.csv(as.matrix(object@raw.data), file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "scale.csv.gz", sep = "/"), "w+")
    write.csv(object@scale.data, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "meta.csv.gz", sep = "/"), "w+")
    write.csv(object@meta.data, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "pca.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$pca@cell.embeddings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "pca_gene.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$pca@gene.loadings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "harmony.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$harmony@cell.embeddings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "harnomy_gene.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$harmony@gene.loadings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "tsne.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$tsne@cell.embeddings, file = gz)
    close(gz)
    
    gz = gzfile(paste(path, "umap.csv.gz", sep = "/"), "w+")
    write.csv(object@dr$umap@cell.embeddings, file = gz)
    close(gz)
}



### function to make modified DimPlot
myDimPlot <- function(
    object,
    reduction.use = "pca",
    dim.1 = 1,
    dim.2 = 2,
    cells.use = NULL,
    pt.size = 1,
    do.return = FALSE,
    do.bare = FALSE,
    cols.use = NULL,
    group.by = "ident",
    pt.shape = NULL,
    do.hover = FALSE,
    data.hover = 'ident',
    do.identify = FALSE,
    do.label = FALSE,
    label.size = 4,
    no.legend = FALSE,
    coord.fixed = FALSE,
    no.axes = FALSE,
    dark.theme = FALSE,
    plot.order = NULL,
    cells.highlight = NULL,
    cols.highlight = 'red',
    sizes.highlight = 1,
    plot.title = NULL,
    vector.friendly = FALSE,
    png.file = NULL,
    png.arguments = c(10,10, 100),
    na.value = 'grey50',
    ...
) {
    #first, consider vector friendly case
    if (vector.friendly) {
        previous_call <- blank_call <- png_call <-  match.call()
        blank_call$pt.size <- -1
        blank_call$do.return <- TRUE
        blank_call$vector.friendly <- FALSE
        png_call$no.axes <- TRUE
        png_call$no.legend <- TRUE
        png_call$do.return <- TRUE
        png_call$vector.friendly <- FALSE
        png_call$plot.title <- NULL
        blank_plot <- eval(blank_call, sys.frame(sys.parent()))
        png_plot <- eval(png_call, sys.frame(sys.parent()))
        png.file <- SetIfNull(x = png.file, default = paste0(tempfile(), ".png"))
        ggsave(
            filename = png.file,
            plot = png_plot,
            width = png.arguments[1],
            height = png.arguments[2],
            dpi = png.arguments[3]
        )
        to_return <- AugmentPlot(plot1 = blank_plot, imgFile = png.file)
        file.remove(png.file)
        if (do.return) {
            return(to_return)
        } else {
            print(to_return)
        }
    }
    embeddings.use <- GetDimReduction(
        object = object,
        reduction.type = reduction.use,
        slot = "cell.embeddings"
    )
    if (length(x = embeddings.use) == 0) {
        stop(paste(reduction.use, "has not been run for this object yet."))
    }
    cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
    dim.code <- GetDimReduction(
        object = object,
        reduction.type = reduction.use,
        slot = "key"
    )
    dim.codes <- paste0(dim.code, c(dim.1, dim.2))
    data.plot <- as.data.frame(x = embeddings.use)
    # data.plot <- as.data.frame(GetDimReduction(object, reduction.type = reduction.use, slot = ""))
    cells.use <- intersect(x = cells.use, y = rownames(x = data.plot))
    data.plot <- data.plot[cells.use, dim.codes]
    ident.use <- as.factor(x = object@ident[cells.use])
    if (group.by != "ident") {
        ident.use <- as.factor(x = FetchData(
            object = object,
            vars.all = group.by
        )[cells.use, 1])
    }
    data.plot$ident <- ident.use
    data.plot$x <- data.plot[, dim.codes[1]]
    data.plot$y <- data.plot[, dim.codes[2]]
    data.plot$pt.size <- pt.size
    if (!is.null(x = cells.highlight)) {
        # Ensure that cells.highlight are in our data.frame
        if (is.character(x = cells.highlight)) {
            cells.highlight <- list(cells.highlight)
        } else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
            cells.highlight <- as.list(x = cells.highlight)
        }
        cells.highlight <- lapply(
            X = cells.highlight,
            FUN = function(cells) {
                cells.return <- if (is.character(x = cells)) {
                    cells[cells %in% rownames(x = data.plot)]
                } else {
                    cells <- as.numeric(x = cells)
                    cells <- cells[cells <= nrow(x = data.plot)]
                    rownames(x = data.plot)[cells]
                }
                return(cells.return)
            }
        )
        # Remove groups that had no cells in our dataframe
        cells.highlight <- Filter(f = length, x = cells.highlight)
        if (length(x = cells.highlight) > 0) {
            if (!no.legend) {
                no.legend <- is.null(x = names(x = cells.highlight))
            }
            names.highlight <- if (is.null(x = names(x = cells.highlight))) {
                paste0('Group_', 1L:length(x = cells.highlight))
            } else {
                names(x = cells.highlight)
            }
            sizes.highlight <- rep_len(
                x = sizes.highlight,
                length.out = length(x = cells.highlight)
            )
            cols.highlight <- rep_len(
                x = cols.highlight,
                length.out = length(x = cells.highlight)
            )
            highlight <- rep_len(x = NA_character_, length.out = nrow(x = data.plot))
            if (is.null(x = cols.use)) {
                cols.use <- 'grey'
            }
            cols.use <- c(cols.use[1], cols.highlight)
            size <- rep_len(x = pt.size, length.out = nrow(x = data.plot))
            for (i in 1:length(x = cells.highlight)) {
                cells.check <- cells.highlight[[i]]
                index.check <- match(x = cells.check, rownames(x = data.plot))
                highlight[index.check] <- names.highlight[i]
                size[index.check] <- sizes.highlight[i]
            }
            plot.order <- sort(x = unique(x = highlight), na.last = TRUE)
            plot.order[is.na(x = plot.order)] <- 'Unselected'
            highlight[is.na(x = highlight)] <- 'Unselected'
            highlight <- as.factor(x = highlight)
            data.plot$ident <- highlight
            data.plot$pt.size <- size
            if (dark.theme) {
                cols.use[1] <- 'white'
            }
        }
    }
    if (!is.null(x = plot.order)) {
        if (any(!plot.order %in% data.plot$ident)) {
            stop("invalid ident in plot.order")
        }
        plot.order <- rev(x = c(
            plot.order,
            setdiff(x = unique(x = data.plot$ident), y = plot.order)
        ))
        data.plot$ident <- factor(x = data.plot$ident, levels = plot.order)
        data.plot <- data.plot[order(data.plot$ident), ]
    }
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) +
        geom_point_rast(mapping = aes(colour = factor(x = ident), size = pt.size))
    if (!is.null(x = pt.shape)) {
        shape.val <- FetchData(object = object, vars.all = pt.shape)[cells.use, 1]
        if (is.numeric(shape.val)) {
            shape.val <- cut(x = shape.val, breaks = 5)
        }
        data.plot[, "pt.shape"] <- shape.val
        p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) +
            geom_point_rast(mapping = aes(
                colour = factor(x = ident),
                shape = factor(x = pt.shape),
                size = pt.size
            ))
    }
    if (!is.null(x = cols.use)) {
        p <- p + scale_colour_manual(values = cols.use, na.value=na.value)
    }
    if(coord.fixed){
        p <- p + coord_fixed()
    }
    p <- p + guides(size = FALSE)
    p2 <- p +
        xlab(label = dim.codes[[1]]) +
        ylab(label = dim.codes[[2]]) +
        scale_size(range = c(min(data.plot$pt.size), max(data.plot$pt.size)))
    p3 <- p2 +
        SetXAxisGG() +
        SetYAxisGG() +
        SetLegendPointsGG(x = 6) +
        SetLegendTextGG(x = 12) +
        no.legend.title +
        theme_bw() +
        NoGrid()
    if (dark.theme) {
        p <- p + DarkTheme()
        p3 <- p3 + DarkTheme()
    }
    p3 <- p3 + theme(legend.title = element_blank())
    if (!is.null(plot.title)) {
        p3 <- p3 + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
    }
    if (do.label) {
        data.plot %>%
            dplyr::group_by(ident) %>%
            summarize(x = median(x = x), y = median(x = y)) -> centers
        p3 <- p3 +
            geom_point_rast(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
            geom_text(data = centers, mapping = aes(label = ident), size = label.size)
    }
    if (no.legend) {
        p3 <- p3 + theme(legend.position = "none")
    }
    if (no.axes) {
        p3 <- p3 + theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank()
        )
    }
    if (do.identify || do.hover) {
        if (do.bare) {
            plot.use <- p
        } else {
            plot.use <- p3
        }
        if (do.hover) {
            if (is.null(x = data.hover)) {
                features.info <- NULL
            } else {
                features.info <- FetchData(object = object, vars.all = data.hover)
            }
            return(HoverLocator(
                plot = plot.use,
                data.plot = data.plot,
                features.info = features.info,
                dark.theme = dark.theme
            ))
        } else if (do.identify) {
            return(FeatureLocator(
                plot = plot.use,
                data.plot = data.plot,
                dark.theme = dark.theme,
                ...
            ))
        }
    }
    if (do.return) {
        if (do.bare) {
            return(p)
        } else {
            return(p3)
        }
    } else {
        if (do.bare) {
            print(p)
        } else {
            print(p3)
        }
    }
}

environment(myDimPlot) <- asNamespace('Seurat')




make_tsne_plot <- function(object, color="Disease", alpha=0.8, colors=colors, pt.size = 0.1, guide_ncol=1) {
    meta = object@meta.data
    
    meta = merge(meta, object@dr$tsne@cell.embeddings, by = 0)
    
    command <- paste('ggplot(meta, aes(x=tSNE_1, y=tSNE_2, color=', color, '))', sep = '')
    
    p = eval(parse(text = command))
    
    p <- p + 
        geom_point_rast(alpha=alpha, size=pt.size) +
        # SetXAxisGG() +
        # SetYAxisGG() +
        # SetLegendPointsGG(x = 6) +
        # SetLegendTextGG(x = 12) +
        no.legend.title +
        theme_bw() +
        # NoGrid() +
        theme(
            legend.position = "none",
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 20)
        ) +
        guides(color = guide_legend(ncol=guide_ncol, override.aes = list(size = 10, alpha = 1))) +
        scale_color_manual(values=colors)
    
    return(p)
}
environment(make_tsne_plot) <- asNamespace('Seurat')





SmoteClassif1 <- function(form, dat, C.perc = "balance",
                          k = 5, repl = FALSE, dist = "Euclidean", p = 2)
    # Args:
    # form    a model formula
    # dat    the original training set (with the unbalanced distribution)
    # C.perc  named list containing each class percentage of under- or 
    #         over-sampling to apply between 0 and 1. The user may provide
    #         only a subset of the existing classes where sampling is to
    #         be applied. Alternatively it may be "balance" (the default) or
    #         "extreme", cases where the sampling percentages are automatically
    #         estimated.
    # k       is the number of neighbors to consider as the pool from where
    #         the new examples are generated
# repl    is it allowed to perform sampling with replacement, when 
#         performing under-sampling
# dist    is the distance measure to be used (defaults to "Euclidean")
# p       is a parameter used when a p-norm is computed
#
# Returns: a new data frame modified through the smote algorithm

{
    if (any(is.na(dat))) {
        stop("The data set provided contains NA values!")
    }
    # the column where the target variable is
    tgt <- which(names(dat) == as.character(form[[2]]))
    names <- sort(unique(dat[, tgt]))
    li <- class.freq(dat, tgt)
    if (tgt < ncol(dat)) {
        orig.order <- colnames(dat)
        cols <- 1:ncol(dat)
        cols[c(tgt, ncol(dat))] <- cols[c(ncol(dat), tgt)]
        dat <- dat[, cols]
    }
    
    if (is.list(C.perc)) {
        names.und <- names(which(C.perc < 1))
        names.ove <- names(which(C.perc > 1))
        names.same <- setdiff(names, union(names.und, names.ove))
        
        # include examples from classes unchanged
        newdata <- dat[which(dat[, ncol(dat)] %in% names.same), ]
        
        if (length(names.und)) {  # perform under-sampling
            for (i in 1:length(names.und)) {
                Exs <- which(dat[, ncol(dat)] == names.und[i])
                sel <- sample(Exs,
                              as.integer(C.perc[[names.und[i]]] * length(Exs)),
                              replace = repl)
                newdata <- rbind(newdata, dat[sel, ])
            }
        }
        if (length(names.ove)) { # perform over-sampling
            for (i in 1:length(names.ove)) {
                if(length(which(dat[, ncol(dat)] == names.ove[i])) == 1){
                    warning(paste("SmoteClassif :: Unable to use SmoteClassif in a bump with 1 example.
                                  Introducing replicas of the example."), call.=FALSE)
                    newdata <- rbind(newdata, dat[rep(which(dat[, ncol(dat)] == names.ove[i]),C.perc[names.ove[i]]),])
                } else if (length(which(dat[, ncol(dat)] == names.ove[i])) <= k){
                    warning(paste("SmoteClassif :: Nr of examples is less or equal to k.\n Using k =",
                                  length(which(dat[, ncol(dat)] == names.ove[i]))-1, 
                                  "in the nearest neighbours computation in this bump."), call.=FALSE)
                    Origk <- k
                    k <- length(which(dat[, ncol(dat)] == names.ove[i]))-1
                    newExs <- Smote.exsClassif1(dat[which(dat[, ncol(dat)] == names.ove[i]), ],
                                                ncol(dat),
                                                li[[3]][ove[i]]/li[[2]][ove[i]] + 1,
                                                k,
                                                dist,
                                                p)
                    # add original rare examples and synthetic generated examples
                    newdata <- rbind(newdata, newExs, 
                                     dat[which(dat[,ncol(dat)] == names.ove[i]),])
                    k <- Origk
                } else {
                    newExs <- Smote.exsClassif1(dat[which(dat[, ncol(dat)] == names.ove[i]), ],
                                                ncol(dat),
                                                C.perc[[names.ove[i]]],
                                                k,
                                                dist,
                                                p)
                    # add original rare examples and synthetic generated examples
                    newdata <- rbind(newdata, newExs,
                                     dat[which(dat[, ncol(dat)] == names.ove[i]), ])
                }
            }
        }
    } else {
        if (C.perc == "balance") {  
            li[[3]] <- round(sum(li[[2]])/length(li[[2]]), 0) - li[[2]]
        } else if (C.perc == "extreme") {
            med <- sum(li[[2]])/length(li[[2]])
            li[[3]] <- round(med^2/li[[2]] * sum(li[[2]])/sum(med^2/li[[2]]), 0) - li[[2]]
        } else {
            stop("Please provide a list with classes to under-/over-sample
                 or alternatively indicate 'balance' or 'extreme'.")
        }
        und <- which(li[[3]] < 0) # classes to under-sample
        ove <- which(li[[3]] > 0) #classes to over-sample
        same <- which(li[[3]] == 0) # unchanged classes
        
        # include examples from classes unchanged
        newdata <- dat[which(dat[, ncol(dat)] %in% li[[1]][same]), ]
        
        if (length(und)) { #perform under-sampling
            for (i in 1:length(und)) { 
                Exs <- which(dat[, ncol(dat)] == li[[1]][und[i]])
                sel <- sample(Exs,
                              as.integer(li[[2]][und[i]] + li[[3]][und[i]]),
                              replace = repl)
                newdata <- rbind(newdata, dat[sel, ])
            }
        }
        
        if (length(ove)) { #perform over-sampling
            for (i in 1:length(ove)) {
                if(length(which(dat[, ncol(dat)] == li[[1]][ove[i]])) == 1){
                    warning(paste("SmoteClassif :: Unable to use SmoteClassif in a bump with 1 example.
                                  Introducing replicas of the example."), call.=FALSE)
                    newdata <- rbind(newdata, dat[rep(which(dat[, ncol(dat)] == li[[1]][ove[i]]), li[[3]][ove[i]]),])
                } else if(length(which(dat[, ncol(dat)] == li[[1]][ove[i]]))<= k){
                    warning(paste("SmoteClassif :: Nr of examples is less or equal to k.\n Using k =",
                                  length(which(dat[, ncol(dat)] == li[[1]][ove[i]]))-1, 
                                  "in the nearest neighbours computation in this bump."), call.=FALSE)
                    Origk <- k
                    k <- length(which(dat[, ncol(dat)] == li[[1]][ove[i]]))-1
                    newExs <- Smote.exsClassif1(dat[which(dat[, ncol(dat)] == li[[1]][ove[i]]), ],
                                                ncol(dat),
                                                li[[3]][ove[i]]/li[[2]][ove[i]] + 1,
                                                k,
                                                dist,
                                                p)
                    # add original rare examples and synthetic generated examples
                    newdata <- rbind(newdata, newExs, 
                                     dat[which(dat[,ncol(dat)] == li[[1]][ove[i]]),])
                    k <- Origk
                } else {
                    newExs <- Smote.exsClassif1(dat[which(dat[, ncol(dat)] == li[[1]][ove[i]]), ],
                                                ncol(dat),
                                                li[[3]][ove[i]]/li[[2]][ove[i]] + 1,
                                                k,
                                                dist,
                                                p)
                    # add original rare examples and synthetic generated examples
                    newdata <- rbind(newdata, newExs, 
                                     dat[which(dat[,ncol(dat)] == li[[1]][ove[i]]),])
                }
            } 
        }
        
        }
    
    if (tgt < ncol(dat)) {
        newdata <- newdata[,cols]
        dat <- dat[,cols]
    }
    
    newdata
}


# ===================================================
# Obtain a set of smoted examples for a set of rare cases.
# L. Torgo, Feb 2010
# P.Branco, Mar,Apr 2015
# ---------------------------------------------------
Smote.exsClassif1 <- function(dat, tgt, N, k, dist, p)
    # INPUTS:
    # dat   are the rare cases (the minority class cases)
    # tgt    is the name of the target variable
    # N      is the percentage of over-sampling to carry out;
    # k      is the number of nearest neighbors to use for the generation
    # dist   is the distance function to use for the neighbors computation
    # p      is an integer used when a "p-norm" distance is selected
    # OUTPUTS:
    # The result of the function is a (N-1)*nrow(dat) set of generated
    # examples with rare class on the target
{
    nomatr <- c()
    T <- matrix(nrow = dim(dat)[1], ncol = dim(dat)[2] - 1)
    for (col in seq.int(dim(T)[2])) { 
        if (class(dat[, col]) %in% c('factor', 'character')) {
            T[, col] <- as.integer(dat[, col])
            nomatr <- c(nomatr, col)
        } else {
            T[, col] <- dat[, col]
        }
    }
    
    print(nomatr)
    nC <- dim(T)[2]
    nT <- dim(T)[1]
    
    # check if there is enough data to determine the k neighbors
    if (nT <= k) {
        stop("Trying to determine ",k, " neighbors for a subset with only ",
             nT, " examples")
    }
    
    kNNs <- neighbours(tgt, dat, dist, p, k)
    
    nexs <-  as.integer(N - 1) # nr of examples to generate for each rare case
    extra <- as.integer(nT * (N - 1 - nexs)) # the extra examples to generate
    idx <- sample(1:nT, extra)
    newM <- matrix(nrow = nexs * nT + extra, ncol = nC)    # the new cases
    if (nexs) {
        for (i in 1:nT) {
            for (n in 1:nexs) {
                # select randomly one of the k NNs
                neig <- sample(1:k, 1)
                
                # the attribute values of the generated case
                difs <- T[kNNs[i, neig], ] - T[i, ]
                newM[(i - 1) * nexs + n, ] <- T[i, ] + runif(1) * difs
                for (a in nomatr) {
                    # nominal attributes are randomly selected among the existing values
                    # of seed and the selected neighbor 
                    newM[(i - 1) * nexs + n, a] <- c(T[kNNs[i, neig], a], 
                                                     T[i, a])[1 + round(runif(1), 0)]
                }
            }
        }
    }
    if (extra) {
        count <- 1
        for (i in idx) {    
            # select randomly one of the k NNs
            neig <- sample(1:k, 1)
            # the attribute values of the generated case
            difs <- T[kNNs[i, neig], ] - T[i, ]
            newM[nexs * nT + count, ] <- T[i, ] + runif(1) * difs
            for (a in nomatr) {
                newM[nexs * nT + count, a] <- c(T[kNNs[i, neig], a], 
                                                T[i, a])[1 + round(runif(1), 0)]
            }
            count <- count + 1
        }
    }
    newCases <- data.frame(newM)
    
    # for (a in nomatr){
    # 
    #     newCases[, a] <- factor(newCases[, a],
    #                             levels = 1:nlevels(dat[, a]),
    #                             labels = unique(dat[, a]))
    # }
    newCases[, tgt] <- factor(rep(dat[1, tgt], nrow(newCases)),
                              levels = levels(dat[, tgt]))
    colnames(newCases) <- colnames(dat)
    newCases
}

environment(SmoteClassif1) <- asNamespace('UBL')
environment(Smote.exsClassif1) <- asNamespace('UBL')



# Function to make combined bar plots
# 1. percentage plot by Stage
# 2. percentage plot by Disease
# 3. percentage plot by PatientID
# 4. barplot of number of cells
# 5. boxplots of number of transcripts
# Note: custom_colors needs to exists in namespace, cause this function using customed color palette
# :param object: Seurat object 2.3+
# :param group.by: how to make groups
# :param group.type: if int is set, sort the plot order by numeric, else by str
# :param rel_height: this function separate legend and main plot into two grid (up and down), so using this to adjust the relative position of legend and main plot
# :param with_disease: whethter plot the disease percentage
# :param ncol: the legend ncol of patient plot
# :return ggplot2 object
make_combined_barplot <- function(
    object, 
    group.by = "res.0.8", 
    group.type = "int", 
    rel_height = c(0.2, 1),
    with_disease=TRUE,
    ncol=3
) {
    ## make percentage bar plot with stage
    meta <- object@meta.data
    meta$group <- as.character(meta[, group.by])
    # meta$newTissue <- as.character(meta$Tissue)
    # meta$newTissue[meta$newTissue == "Tumor" & meta$Disease != ""] <- as.character(meta$Disease[meta$newTissue == "Tumor" & meta$Disease != ""])

    plot_order = c(unique(meta$group), "Total")
    
    if(group.type == "int") {
        plot_order = as.character(c(sort(as.numeric(unique(meta$group))), "Total"))
    }
    
    ### make percentage plot of stage
    temp <- as.data.frame(
        meta %>% 
            dplyr::select(group, Stage) %>% 
            dplyr::filter(!Stage %in% c("Normal", "CPI")) %>% 
            dplyr::group_by(group, Stage) %>% 
            dplyr::add_tally() %>% 
            unique() %>% 
            dplyr::group_by(group) %>% 
            dplyr::mutate(freq = n / sum(n)) %>% unique()
    )
    
    temp1 <- as.data.frame(table(meta$Stage))
    temp1 <- temp1[!temp1[,1] %in% c("Normal", "CPI"), ]
    colnames(temp1) <- c("Stage", "n")
    
    temp1$group = "Total"
    temp1$freq = temp1$n / sum(as.numeric(temp1$n))
    
    temp1 <- temp1[, colnames(temp)]
    
    temp <- rbind(temp, temp1)
    
    temp$group <- factor(temp$group, levels = plot_order)
    
    p1 <- ggplot(data = temp, aes(x=group, y=freq, fill=Stage)) + 
        geom_bar(stat = "identity", width = .75) + 
        coord_flip() + 
        theme(legend.title = element_blank(), legend.justification = "center") +
        labs(y="Fraction of cells", x = "") +
        # scale_x_discrete(breaks=plot_order) + 
        guides(fill = guide_legend(ncol=4)) + 
        scale_fill_manual(values=colors)
    
    p1_legend <- get_legend(p1)
    
    p1 <- p1 + theme(legend.position = "none")
    
    
    ### make percentage of tissue
    temp <- as.data.frame(
        meta %>% 
            dplyr::select(group, Disease) %>% 
            dplyr::group_by(group, Disease) %>% 
            dplyr::add_tally() %>% 
            unique() %>% 
            dplyr::group_by(group) %>% 
            dplyr::mutate(freq = n / sum(n)) %>% 
            unique()
    )
    
    temp1 <- as.data.frame(table(meta$Disease))
    
    colnames(temp1) <- c("Disease", "n")
    temp1$group = "Total"
    temp1$freq = temp1$n / sum(as.numeric(temp1$n))
    
    temp1 <- temp1[, colnames(temp)]
    temp <- rbind(temp, temp1)
    
    temp$group <- factor(temp$group, levels = plot_order)
    
    ## combine colors of tissue and disease
    p2 <- ggplot(data = temp, aes(x=group, y=freq, fill=Disease)) + 
        geom_bar(stat = "identity", width = .75) + 
        coord_flip() + 
        theme(legend.position = "right", legend.title = element_blank(), legend.justification = "center") +
        labs(y="Fraction of cells", x = "") +
        # scale_x_discrete(breaks=plot_order) + 
        guides(fill = guide_legend(ncol=1)) + 
        scale_fill_manual(values=colors)
    
    p2_legend <- get_legend(p2)
    
    p2 <- p2 + theme(legend.position = "none")
    
    
    ### make percentage barplot of patients
    temp <- as.data.frame(
        meta %>% 
            dplyr::select(group, PatientID) %>% 
            dplyr::group_by(group, PatientID) %>% 
            dplyr::add_tally() %>% 
            unique() %>% 
            dplyr::group_by(group) %>% 
            dplyr::mutate(freq = n / sum(n)) %>% 
            unique()
    )
    
    temp1 <- as.data.frame(table(object@meta.data$PatientID))
    colnames(temp1) <- c("PatientID", "n")
    temp1$group = "Total"
    temp1$freq = temp1$n / sum(temp1$n)
    
    temp <- rbind(temp, temp1)
    temp$group <- factor(temp$group, levels = plot_order)
    p3 <- ggplot(data = temp, aes(x=group, y=freq, fill=PatientID)) + 
        geom_bar(stat = "identity", width = .75) + 
        coord_flip() + 
        theme(legend.position = "right", legend.title = element_blank(), legend.justification = "center") +
        labs(y="Fraction of cells", x = "") +
        # scale_x_discrete(breaks=plot_order) + 
        guides(fill = guide_legend(ncol=ncol)) + 
        scale_fill_manual(values=colors)
    p3_legend <- get_legend(p3)
    
    p3 <- p3 + theme(legend.position = "none")
    
    
    ### make bar plot of number of cells
    temp <- as.data.frame(
        meta %>% 
            dplyr::select(group) %>% 
            dplyr::group_by(group) %>% 
            dplyr::add_tally() %>% 
            unique()
    )
    
    temp1 <- data.frame(group="Total", n=nrow(meta))
    temp <- rbind(temp, temp1)
    temp$group <- factor(temp$group, levels = plot_order)
    p4 <- ggplot(temp, aes(x=group, y=n, fill=stata_pal("s1color")(1))) + 
        geom_bar(stat = "identity") +
        coord_flip() + 
        theme(legend.position = "none") +
        labs(y="Number of cells", x = "") +
        # scale_x_discrete(breaks=plot_order) +
        scale_y_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        )
    
    ### make boxplot of transcript number
    meta1 = meta
    meta1$group <- "Total"
    meta <- rbind(meta, meta1)
    meta$group <- factor(meta$group, levels = plot_order)

    p5 <- ggplot(meta, aes(x=group, group=group, y=nUMI)) + 
        geom_boxplot(fill = hc_pal()(1)) + 
        coord_flip() + 
        # scale_x_discrete(breaks=plot_order[1:(length(plot_order) - 1)]) +
        labs(y="Number of transcripts", x = "") + 
        scale_y_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        )
    
    plots = list(p1_legend, p2_legend, p3_legend, ggplot() + theme_void(), ggplot() + theme_void() , p1, p2, p3, p4, p5)
    
    if (!with_disease) {
        plots = list(p1_legend, p3_legend, ggplot() + theme_void(), ggplot() + theme_void(), p1, p3, p4, p5)
    }
    
    p <- plot_grid(plotlist=plots, nrow = 2, rel_heights = rel_height)
    return(p)
}



find_markers <- function(object, n.cores=10, group.by = "Stage") {
    
    temp_group <- sort(unique(object@meta.data[, group.by]))
    
    res = NULL
    for(i in temp_group){

        groups = rownames(object@meta.data[object@meta.data[, group.by] == i, ])
        
        if(length(groups) >= 3) {
            temp <- FindMarkers(
                object = object, 
                ident.1 = groups,
                logfc.threshold = 0
            )
            
            temp$ident = i
            temp$gene = rownames(temp)
            res = rbind(res, temp)
        } 
    }    
    return(res)
}



do_kegg <- function(eg, cluster=NA, pvalueCutoff = 0.05) {
    kk <- enrichKEGG(gene     = eg$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = pvalueCutoff)
    kk <- as.data.frame(kk)
    
    if (nrow(kk) == 0) {
        return(NULL)
    }
    kk$cluster <- cluster
    return(kk)
}


do_go <- function(eg, cluster = NA, pvalueCutoff = 0.01, qvalueCutoff = 0.05, cutoff=0.7) {
    
    res = NULL
    for(i in c("BP", "CC", "MF")) {
        ego <- enrichGO(gene      = eg$ENTREZID,
                        keyType       = 'ENTREZID',
                        OrgDb         = org.Hs.eg.db,
                        ont           = i,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = pvalueCutoff,
                        qvalueCutoff  = qvalueCutoff,
                        readable      = FALSE)
        
        if(is.null(ego)) {
            return(ego)
        }
        ego <- clusterProfiler::simplify(ego, cutoff=cutoff, by="p.adjust", select_fun=min)
        
        ego = as.data.frame(ego)
        
        if (nrow(ego) > 0) {
            ego$ONTOLOGY = i
        }
        
        if (is.null(res)) {
            res = ego
        } else {
            res = rbind(res, ego)
        }
    }
    
    if (nrow(res) == 0) {
        return(NULL)
    }
    
    res$cluster <- cluster
    
    return(res)
}


do_do <- function(eg, cluster = NA, pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 500) {
    do <- enrichDO(gene     = eg$ENTREZID,
                   ont           = "DO",
                   pvalueCutoff  = pvalueCutoff,
                   pAdjustMethod = "BH",
                   minGSSize     = minGSSize,
                   maxGSSize     = maxGSSize,
                   qvalueCutoff  = qvalueCutoff,
                   readable      = TRUE)
    
    do = as.data.frame(do)
    
    if (nrow(do) == 0) {
        return(NULL)
    }
    do$cluster = cluster
    
    return(do)
} 


get_entrzid <- function(markers) {
    res = NULL
    for(i in unique(markers$ident)) {
        print(i)
        
        temp <- markers[markers$ident == i & markers$p_val_adj < 0.05 & markers$avg_logFC > 0.5, ]
        eg <- bitr(unique(temp$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
        
        if(!is.null(eg) && nrow(eg) > 0) {
            eg$ident = i
            
            if(is.null(res)) {
                res = eg
            } else {
                res = rbind(res, eg)
            }
        }
    }
    
    return(res)
}

