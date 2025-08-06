

#############################################################################
# Propgation of evidence
#############################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list=ls())
library(visualizeR)
library(RColorBrewer)
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"
cubicsets <- c(tas_historical_10d_akima_cubic_corrected,listinterim)

namescubics <- names(cubicsets)
shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
                   pattern = "_r12i1p1", replacement = "")
shortnames
############################################################################
# Functions to load hciterations of models in a list 
############################################################################
loadIterations <- function(pattern,permused, hctype = NULL){
  if(is.null(hctype)){hctype <- ""}
  hc_interim_list <- list.files(paste0("data/hciterations/",pattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0("data/hciterations/",pattern,"/perm",permused,hctype))
  hc_interim_names <- gsub(".rda", "", hc_interim_names)
  
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
hc_gcms <- lapply(namescubics, loadIterations, permused = 3)
names(hc_gcms) <- shortnames
data_gcms <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE)))

modelsize <- 18

selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)

hc_gcms <- NULL

loglik_selection_datasets <- matrix(data = NA, nrow = length(selection_fits_gcms), ncol = length(data_gcms), dimnames = list(names(selection_fits_gcms),names(data_gcms)))

for(i in 1:length(selection_fits_gcms)){
  lo <- sapply(X = data_gcms, logLik, object = selection_fits_gcms[[i]])
  loglik_selection_datasets[i,] <- lo
}


abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")

#######################################################################################
# Automatic V81
#######################################################################################

nodenumber <- 81
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

files <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/pos",whichnode,"pos/prop_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}


posposmix <- list()
posposmix <- lapply(shortnames, function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posposmix[[1]]
for(i in 1:length(shortnames)){
  posposmix[[i]] <- eval(parse(text = paste0("prop_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
}
names(posposmix) <- propfiles


prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

files <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/pos",whichnode,"neg/propneg_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}
all.equal(propneg_CMIP5_CanESM2_hc_1700_1800i_V81_equal2,propneg_CMIP5_CNRM.CM5_hc_1700_1800i_V81_equal2)

posnegmix <- lapply(shortnames, function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))

for(i in 1:length(shortnames)){
  posnegmix[[i]] <- eval(parse(text = paste0("propneg_",shortnames,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
}
names(posnegmix) <- propfiles

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

#######################################################################################
# Automatic V81V227
#######################################################################################

nodenumber <- c(81,227)
whichnode <- paste0("V",nodenumber)
whichdata <- "int_"
whichnodecol <- whichnode
if (length(nodenumber) >=1){whichnodecol <- paste0(whichnode, collapse = "")}

numberBN1 <- 18
numberBN <- t(as.vector(c(numberBN1)))
iterationslimits <- apply(numberBN, MARGIN = 2, function(x) c((x-1)*100,x*100))
whichiterations <- apply(iterationslimits, MARGIN = 2, FUN = function(x) paste0(x[[1]],"_",x[[2]]))

###########################################################################
# V81 + V227 + Positive
###########################################################################

files <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/plus",whichnode[1],"plus",whichnode[2],"pos/prop_",shortnames,"_hc_",whichiterations,"i_",whichnodecol,"_equal22.rda")
propfiles <- lapply(files,function(x)get(load(file = x)))
names(propfiles)<- shortnames

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[6]]), SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)

plotname <- paste0("figs/",whichiterations,"_BNgcms_plus",whichnode[1],"plus",whichnode[2],"pos_perm3.pdf")
pdf(plotname)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
##############################################################################
# V81 + V227 + Negative
##############################################################################
files <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/plus",whichnode[1],"plus",whichnode[2],"neg/propneg_",shortnames,"_hc_",whichiterations,"i_",whichnodecol,"_equal22.rda")
propfiles <- lapply(files,function(x)get(load(file = x)))
names(propfiles)<- shortnames


prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[6]]), SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)

plotname <- paste0("figs/",whichiterations,"_BNgcms_plus",whichnode[1],"plus",whichnode[2],"neg_perm3.pdf")
pdf(plotname)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()

###############################################################################
# V81 + V227 - Positive
###############################################################################
files <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/plus",whichnode[1],"min",whichnode[2],"pos/prop_",shortnames,"_hc_",whichiterations,"i_",whichnodecol,"_equal2min2.rda")
propfiles <- lapply(files,function(x)get(load(file = x)))
names(propfiles)<- shortnames

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[6]]), SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)

plotname <- paste0("figs/",whichiterations,"_BNgcms_plus",whichnode[1],"min",whichnode[2],"pos_perm3.pdf")
pdf(plotname)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
###############################################################################
# V81 + V227 - Negative 
###############################################################################
files <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/plus",whichnode[1],"min",whichnode[2],"neg/propneg_",shortnames,"_hc_",whichiterations,"i_",whichnodecol,"_equal2min2.rda")
whichf <- c(1,2,3,4,5,6,7,8,9,10)
propfiles <- lapply(files[whichf],function(x)get(load(file = x)))
names(propfiles)<- shortnames[whichf]

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[6]]), SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev[whichf]
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)

plotname <- paste0("figs/",whichiterations,"_BNgcms_plus",whichnode[1],"min",whichnode[2],"neg_perm3.pdf")
pdf(plotname)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
#######################################################################################
# Automatic V81V171
#######################################################################################

nodenumber <- c(81,171)
whichnode <- paste0("V",nodenumber)
whichdata <- "int_"
whichnodecol <- whichnode
if (length(nodenumber) >=1){whichnodecol <- paste0(whichnode, collapse = "")}

numberBN1 <- 18
numberBN <- t(as.vector(c(numberBN1)))
iterationslimits <- apply(numberBN, MARGIN = 2, function(x) c((x-1)*100,x*100))
whichiterations <- apply(iterationslimits, MARGIN = 2, FUN = function(x) paste0(x[[1]],"_",x[[2]]))

###############################################################################
# V81 + V171 0 Positive
###############################################################################
files <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/plus",whichnode[1],"neu",whichnode[2],"pos/prop_",shortnames,"_hc_",whichiterations,"i_",whichnodecol,"_equal20.rda")
propfiles <- lapply(files,function(x)get(load(file = x)))
names(propfiles)<- shortnames

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[6]]), SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)

plotname <- paste0("figs/",whichiterations,"_BNgcms_plus",whichnode[1],"neu",whichnode[2],"pos_perm3.pdf")
pdf(plotname)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()


###############################################################################
# V81 + V171 0 Negative 
###############################################################################
files <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/plus",whichnode[1],"neu",whichnode[2],"neg/propneg_",shortnames,"_hc_",whichiterations,"i_",whichnodecol,"_equal20.rda")
propfiles <- lapply(files,function(x)get(load(file = x)))
names(propfiles)<- shortnames

prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[6]]), SIMPLIFY = FALSE)
enspropdifsclims <- bindGrid(prop_dif, dimension = c("member"))
enspropdifsclims$Members <- abrev
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)

plotname <- paste0("figs/",whichiterations,"_BNgcms_plus",whichnode[1],"neu",whichnode[2],"neg_perm3.pdf")
pdf(plotname)
spatialPlot(enspropdifsclims,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
