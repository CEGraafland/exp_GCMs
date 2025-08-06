#############################################################################
# Propgation of evidence CSIRO runs
#############################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
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
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models CSIRO
###########################################################################################
cubicsets <- get(load("data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda"))
# cubicsets <- cubicsets[c(1,4,5,6,7,9)]
namescubics <- gsub("-",".",gsub(".ncml","",gsub("https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip5/historical/historical_day_","",sapply(cubicsets,function(x)attr(x,"dataset")))))
############################################################################
# Functions to load hciterations of models in a list 
############################################################################
loadIterations <- function(pattern,permused, hctype = NULL,it = NULL){
  if(is.null(hctype)){hctype <- ""}
  hc_interim_list <- list.files(paste0("data/hciterations/CMIP5_CSIRO/",pattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0("data/hciterations/CMIP5_CSIRO/",pattern,"/perm",permused,hctype))
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
###########################################################################
#
###########################################################################
permk <- 3
it <- "1700_1800"
hc_gcmsCSIRO <- lapply(paste0(namescubics), loadIterations, permused = permk, it = it)
data_gcmsCSIRO <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE)))

modelsize <- 18

selection_hc_gcms <-  mapply(function(x,y) x[[grep((y),x)]], x = hc_gcmsCSIRO, y = modelsize, SIMPLIFY = FALSE)
hc_gcmsCSIRO <- NULL
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcmsCSIRO, SIMPLIFY = FALSE)
data_gcms <- NULL
#######################################################################################
# Automatic V81
#######################################################################################
nodenumber <- 81
whichnode <- paste0("V",nodenumber)
whichdata <- "int_"

if (length(nodenumber) >=1){whichnode <- paste0(whichnode, collapse = "")}

numberBN1 <- 17
numberBN1 <- 18
numberBN <- t(as.vector(c(numberBN1)))
iterationslimits <- apply(numberBN, MARGIN = 2, function(x) c((x-1)*100,x*100))
whichiterations <- apply(iterationslimits, MARGIN = 2, FUN = function(x) paste0(x[[1]],"_",x[[2]]))

files <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP5_CSIRO/prop_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}

posposmix <- list()
posposmix <- lapply(namescubics, function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posposmix[[1]]
for(i in 1:length(namescubics)){
  posposmix[[i]] <- eval(parse(text = paste0("prop_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2")))
}
names(posposmix) <- propfiles


prop_dif_pos <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsposclims <- bindGrid(prop_dif_pos, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposclims$Members <- namescubics
enspropdifsposclims$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsposclims, names.attr = namescubics,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

files <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP5_CSIRO/propneg_",namescubics,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files,function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(files)){load(file = files[i])}

posnegmix <- lapply(namescubics, function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posnegmix
names(posnegmix) <- propfiles
prop_dif_neg <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsnegclims <- bindGrid(prop_dif_neg, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsnegclims$Members <- namescubics
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsnegclims,names.attr = namescubics,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_neg[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop_dif_pos_and_neg <- mapply(function(griddata,proppos,propneg) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")), griddata),griddata = cubicsets,proppos = posposmix,propneg = posnegmix, SIMPLIFY = FALSE)


enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg[c(2:10,1)], dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposandnegclims$Members <- namescubics[c(2:10,1)]
col <- c(rev(col.b),col.r)

plotname <-  paste0("figs/evidence_propagation/CMIP5_CSIRO/CMIP5_CSIRO_",whichiterations,"_plus",whichnode,"pos.pdf")
pdf(plotname,width= 12, height= 3)
plotname <-  paste0("figs/evidence_propagation/CMIP5_CSIRO/CMIP5_CSIRO_",whichiterations,"_plus",whichnode,"pos.png")
png(plotname,units = 'in', width= 12, height= 3,res = 180)

spatialPlot(enspropdifsposandnegclims,names.attr = namescubics[c(2:10,1)],nrow = 2,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
