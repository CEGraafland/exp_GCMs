#############################################################################
# Propgation of evidence
#############################################################################
#setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list=ls())
library(bnlearn)
library(RColorBrewer)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(visualizeR)
########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"

cubicsets5 <- tas_historical_10d_akima_cubic_corrected
cubicsets5extra <- tas_historical_cmip5_extra_10d_akima_cubic
cubicsets5left <- tas_historical_cmip5_left_10d_akima_cubic
namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)

ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")
ex.subset <- c("CMIP5_EC.EARTH_r1i1p1","CMIP5_MIROC5_r1i1p1","CMIP5_MRI.CGCM3_r1i1p1","CMIP5_CSIRO.Mk3.6.0_r1i1p1","CMIP5_CNRM.CM5_r1i1p1","CMIP5_HadGEM2.CC_r1i1p1","CMIP5_HadGEM2.ES_r1i1p1")
cubicsets <- c(listinterim,cubicsets5,cubicsets5extra,cubicsets5left)[ex.subset]



# cubicsets <- c(tas_historical_10d_akima_cubic_corrected,listinterim)

namescubics <- names(cubicsets)
shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
                   pattern = "_r12i1p1", replacement = "")
shortnames5 <- gsub(gsub(names(cubicsets5), pattern = "_r1i1p1", replacement = ""),
                    pattern = "_r12i1p1", replacement = "")
############################################################################
# Functions to load hciterations of models in a list 
############################################################################
loadIterations <- function(pattern,permused, hctype = NULL, ncep = FALSE, it = NULL){
  if(is.null(hctype)){hctype <- ""}
  if(ncep == FALSE){fullpattern <- paste0("data/hciterations/",pattern)} else {fullpattern <- pattern}
  
  hc_interim_list <- list.files(paste0(fullpattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0(fullpattern,"/perm",permused,hctype))
  hc_interim_names <- gsub(".rda", "", hc_interim_names)
  if(!is.null(it)){
    hc_interim_list <- hc_interim_list[grep(it,hc_interim_list)]
    hc_interim_names <- hc_interim_names[grep(it,hc_interim_names)]
  }
  
  hc_interim_networks <- list()
  
  hc_interim_networks <- lapply(hc_interim_list, function(x){get(load(x))})
  names(hc_interim_networks) <- hc_interim_names
  interimsizes <- sapply(hc_interim_networks,narcs)
  hc_interims <- hc_interim_networks[order(interimsizes)]
  return(hc_interims)
}
################################################################################
# How do the models of various gcms explain the data of various gcms?
################################################################################
it <- "1700_1800"
permk <- 3
hc_gcms5 <- lapply(paste0(namescubics5), loadIterations, permused = permk, it = it)
hc_gcms5extra <- lapply(paste0("CMIP5_extra/",namescubics5extra), loadIterations, permused = permk, it = it)
hc_gcms5left<- lapply(paste0("CMIP5_left/",namescubics5left), loadIterations, permused = permk, it = it)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk, it = it) 

ex.subset.sel <- c("CMIP5_EC.EARTH_r1i1p1","CMIP5_MIROC5_r1i1p1","CMIP5_MRI.CGCM3_r1i1p1","CMIP5_CSIRO.Mk3.6.0_r1i1p1","CMIP5_CNRM.CM5","CMIP5_HadGEM2.CC_r1i1p1","CMIP5_HadGEM2.ES")
selection_hc_gcms <- (unlist(c(hc_gcms5,hc_gcms5extra,hc_gcms5left),recursive = FALSE))[paste0(ex.subset.sel,"_",permk,"_",it,"i")]

#selection_hc_gcms <- (unlist(c(hc_interim,hc_gcms5,hc_gcms5extra,hc_gcms5left),recursive = FALSE))[paste0(ex.subset.sel,"_",permk,"_",it,"i")]



# hc_gcms<- lapply(namescubics, loadIterations, permused = 3)
# names(hc_gcms) <- shortnames
data_gcms <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE)))

modelsize <- 18

# selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
# hc_gcms <- NULL
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)
data_gcms <- NULL

hc_gcms <- NULL

loglik_selection_datasets <- matrix(data = NA, nrow = length(selection_fits_gcms), ncol = length(data_gcms), dimnames = list(names(selection_fits_gcms),names(data_gcms)))

for(i in 1:length(selection_fits_gcms)){
  lo <- sapply(X = data_gcms, logLik, object = selection_fits_gcms[[i]])
  loglik_selection_datasets[i,] <- lo
}


abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")

#######################################################################################
# Automatic V477
#######################################################################################

nodenumber <- 477
whichnode <- paste0("V",nodenumber)
whichdata <- "int_"

if (length(nodenumber) >=1){whichnode <- paste0(whichnode, collapse = "")}

numberBN1 <- 17
numberBN1 <- 18
# numberBN2 <- 83
# numberBN2 <- 87
# numberBN <- t(as.vector(c(numberBN1,numberBN2)))
numberBN <- t(as.vector(c(numberBN1)))
iterationslimits <- apply(numberBN, MARGIN = 2, function(x) c((x-1)*100,x*100))
whichiterations <- apply(iterationslimits, MARGIN = 2, FUN = function(x) paste0(x[[1]],"_",x[[2]]))

whichsize <- nedges_int_hc[c(numberBN1,numberBN2)]
whichloglik <- logliks_int_hc[c(numberBN1,numberBN2)]

files <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}


posposmix <- list()
posposmix <- lapply(namescubics, function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posposmix[[1]]
# for(i in 1:length(shortnames)){
  # posposmix[[i]] <- eval(parse(text = paste0("prop_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
# }
names(posposmix) <- propfiles


prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

files <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
# propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}

posnegmix <- lapply(namescubics, function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))

# for(i in 1:length(shortnames)){
  # posnegmix[[i]] <- eval(parse(text = paste0("propneg_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
# }
names(posnegmix) <- propfiles

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


prop_dif_pos_and_neg <- mapply(function(griddata,proppos,propneg) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),griddata = cubicsets,proppos = posposmix,propneg = posnegmix, SIMPLIFY = FALSE)
#prop_dif_pos_and_neg <- prop_dif_pos_and_neg[sort(names(prop_dif_pos_and_neg))]

enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposandnegclims$Members <- namescubics

col <- c(rev(col.b),col.r)
plotname <- paste0("figs/evidence_propagation/interim_explore/interim_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")
plotname <- paste0("figs/evidence_propagation/interim_explore/subset_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")
plotname <- paste0("figs/evidence_propagation/interim_explore/subset2_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")
pdf(plotname, width = 10, height = 3)
spatialPlot(enspropdifsposandnegclims,names.attr = namescubics,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,layout = c(4,2),
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
#######################################################################################
# Automatic V478
#######################################################################################

nodenumber <- 478
whichnode <- paste0("V",nodenumber)
whichdata <- "int_"

if (length(nodenumber) >=1){whichnode <- paste0(whichnode, collapse = "")}

numberBN1 <- 17
numberBN1 <- 18
# numberBN2 <- 83
# numberBN2 <- 87
# numberBN <- t(as.vector(c(numberBN1,numberBN2)))
numberBN <- t(as.vector(c(numberBN1)))
iterationslimits <- apply(numberBN, MARGIN = 2, function(x) c((x-1)*100,x*100))
whichiterations <- apply(iterationslimits, MARGIN = 2, FUN = function(x) paste0(x[[1]],"_",x[[2]]))

whichsize <- nedges_int_hc[c(numberBN1,numberBN2)]
whichloglik <- logliks_int_hc[c(numberBN1,numberBN2)]

files <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}


posposmix <- list()
posposmix <- lapply(namescubics, function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
# posposmix[[1]]
# for(i in 1:length(shortnames)){
  # posposmix[[i]] <- eval(parse(text = paste0("prop_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
# }
names(posposmix) <- propfiles


prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

files <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}

posnegmix <- lapply(namescubics, function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))

# for(i in 1:length(shortnames)){
  # posnegmix[[i]] <- eval(parse(text = paste0("propneg_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
# }
names(posnegmix) <- propfiles

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


prop_dif_pos_and_neg <- mapply(function(griddata,proppos,propneg) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),griddata = cubicsets,proppos = posposmix,propneg = posnegmix, SIMPLIFY = FALSE)
#prop_dif_pos_and_neg <- prop_dif_pos_and_neg[sort(names(prop_dif_pos_and_neg))]

enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposandnegclims$Members <- namescubics

col <- c(rev(col.b),col.r)

plotname <- paste0("figs/evidence_propagation/interim_explore/interim_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")
plotname <- paste0("figs/evidence_propagation/interim_explore/subset_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")
plotname <- paste0("figs/evidence_propagation/interim_explore/subset2_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")

pdf(plotname, width = 10, height = 3)
spatialPlot(enspropdifsposandnegclims,names.attr = namescubics,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col, layout = c(4,2),
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

#######################################################################################
# Automatic V95
#######################################################################################

nodenumber <- 95
whichnode <- paste0("V",nodenumber)
whichdata <- "int_"

if (length(nodenumber) >=1){whichnode <- paste0(whichnode, collapse = "")}

numberBN1 <- 17
numberBN1 <- 18
# numberBN2 <- 83
# numberBN2 <- 87
# numberBN <- t(as.vector(c(numberBN1,numberBN2)))
numberBN <- t(as.vector(c(numberBN1)))
iterationslimits <- apply(numberBN, MARGIN = 2, function(x) c((x-1)*100,x*100))
whichiterations <- apply(iterationslimits, MARGIN = 2, FUN = function(x) paste0(x[[1]],"_",x[[2]]))

whichsize <- nedges_int_hc[c(numberBN1,numberBN2)]
whichloglik <- logliks_int_hc[c(numberBN1,numberBN2)]

files <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}


posposmix <- list()
posposmix <- lapply(namescubics, function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posposmix[[1]]
for(i in 1:length(namescubics)){
  posposmix[[i]] <- eval(parse(text = paste0("prop_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
}
names(posposmix) <- propfiles


prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

files <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}

posnegmix <- lapply(namescubics, function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))

# for(i in 1:length(shortnames)){
  # posnegmix[[i]] <- eval(parse(text = paste0("propneg_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
# }
names(posnegmix) <- propfiles

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline",
            lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop_dif_pos_and_neg <- mapply(function(griddata,proppos,propneg) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),griddata = cubicsets,proppos = posposmix,propneg = posnegmix, SIMPLIFY = FALSE)
#prop_dif_pos_and_neg <- prop_dif_pos_and_neg[sort(names(prop_dif_pos_and_neg))]

enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposandnegclims$Members <- namescubics

col <- c(rev(col.b),col.r)

plotname <- paste0("figs/evidence_propagation/interim_explore/interim_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")
plotname <- paste0("figs/evidence_propagation/interim_explore/subset_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")
plotname <- paste0("figs/evidence_propagation/interim_explore/subset2_",whichiterations,"_plus",whichnode,"pos_and_neg_perm3.pdf")

pdf(plotname, width = 10, height = 3)
spatialPlot(enspropdifsposandnegclims,names.attr = namescubics,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,layout = c(4,2),
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()

#######################################################################################
# Automatic V337V645
#######################################################################################

nodenumber <- c(337,645)
whichnode <- paste0("V",nodenumber)
whichdata <- "int_"
whichnodecol <- paste0(whichnode[1],"cubic",whichnode[2])
if (length(nodenumber) >=1){whichnodecol <- paste0(whichnode, collapse = "cubic")}

numberBN1 <- 18
numberBN <- t(as.vector(c(numberBN1)))
iterationslimits <- apply(numberBN, MARGIN = 2, function(x) c((x-1)*100,x*100))
whichiterations <- apply(iterationslimits, MARGIN = 2, FUN = function(x) paste0(x[[1]],"_",x[[2]]))


###########################################################################
# V337cubicV645 + Positive
###########################################################################

files <- paste0("results/propagation/perm3/plus",whichnodecol,"pos/prop_",shortnames,"_hc_",whichiterations,"i_",whichnodecol,"_equal2.rda")
posposmix <- lapply(files,function(x)get(load(file = x)))
names(posposmix)<- shortnames

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[1]]), SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)

spatialPlot(enspropdifsclims, as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

##############################################################################
# V337cubicV645 + Negative
##############################################################################
files <- paste0("results/propagation/perm3/plus",whichnodecol,"neg/propneg_",shortnames,"_hc_",whichiterations,"i_",whichnodecol,"_equal2.rda")
posnegmix <- lapply(files,function(x)get(load(file = x)))
names(posnegmix)<- shortnames

posnegmix <- eval(parse(text = paste0("propneg_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2")))

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[1]]), SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)


spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


prop_dif_pos_and_neg <- mapply(function(griddata,proppos,propneg) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),griddata = cubicsets,proppos = posposmix,propneg = posnegmix, SIMPLIFY = FALSE)
#prop_dif_pos_and_neg <- prop_dif_pos_and_neg[sort(names(prop_dif_pos_and_neg))]

enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposandnegclims$Members <- namescubics

col <- c(rev(col.b),col.r)

plotname <- paste0("figs/evidence_propagation/interim_explore/interim_",whichiterations,"_plus",whichnodecol,"pos_and_neg_perm3.pdf")
pdf(plotname)
spatialPlot(enspropdifsposandnegclims,names.attr = namescubics,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
