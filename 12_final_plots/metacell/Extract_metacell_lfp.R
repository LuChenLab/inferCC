args = commandArgs(trailingOnly = T)

files = list.files("04_metacell/LUAD/", pattern = "mc.metacell_mc.Rda", recursive=T, full.name = T)
files = c(files, list.files("04_metacell/LUSC/", pattern = "mc.metacell_mc.Rda", recursive=T, full.name = T))

read_lfp <- function(path) {

    load(paste(path, paste0("mc.metacell_mc.Rda"), sep = "/"))
    mc = object
    
    load(paste(path, paste0("gstat.metacell.Rda"), sep = "/"))
    gstat = object
    
    
    fp_max = apply(mc@mc_fp, 1, max)
    fp_tot = gstat[intersect(rownames(mc@mc_fp), rownames(gstat)), "tot"]
    
    
    f = fp_max > 0 & fp_tot > 0
    
    lfp = round(log2(mc@mc_fp[f,]), 2)
    
    return(lfp)
}


for(i in files) {
    print(i)

    data = read_lfp(dirname(i))

    write.csv(t(data), paste(dirname(i), "lfp.csv", sep = "/"))
}


