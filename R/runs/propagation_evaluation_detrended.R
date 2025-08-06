#############################################################################
# Propgation of evidence
#############################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list=ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"

cubicsets5 <- tas_historical_10d_akima_cubic_corrected
namescubics5 <- names(cubicsets5)
shortnames5 <- gsub(gsub(namescubics5, pattern = "_r1i1p1", replacement = ""),
                    pattern = "_r12i1p1", replacement = "")

cubicsets5.rcp85 <- tas_rcp85_cmip5_10d_akima_cubic_corrected
cubicsets.future <- c(cubicsets5.rcp85)
institutions.overlap <- c("IPSL.CM5A","CESM1","EC.EARTH_r12i1p1","Can","MRI","MIROC.ES","MIROC5","CNRM","MPI","GFDL.ESM","Nor")
institutions <- rep(institutions.overlap,each =1)
CMIPS.5 <- rep(c("CMIP5"),length(institutions.overlap))
combinations.5 <- cbind(CMIPS.5,institutions)
namesort.5 <- unlist(apply(combinations.5,MARGIN = 1,function(x) grep(paste0(x[2]),names(cubicsets.future))))

# For overlap # without overlap
overlap = TRUE
if(overlap == TRUE){cubicsets.future <- cubicsets.future[namesort.5]
} else if (overlap == FALSE){cubicsets.future <- cubicsets.future[-namesort.5]}

namescubics.future <-names(cubicsets.future)
namescubics5.rcp85 <- names(cubicsets5.rcp85)


cubicsets <- c(cubicsets5,cubicsets.future)

namescubics <- names(cubicsets)
shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
                   pattern = "_r12i1p1", replacement = "")
shortnames
############################################################################
# Functions to load hciterations of models in a list 
############################################################################
loadIterations <- function(pattern,permused, hctype = NULL){
  if(is.null(hctype)){hctype <- ""}
  hc_interim_list <- list.files(paste0("data/hciterations_detrended/",pattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0("data/hciterations_detrended/",pattern,"/perm",permused,hctype))
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
data_gcms <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(detrendGrid(x), rms = TRUE)))
hc_gcms$CMIP5_NorESM1.M_rcp85
modelsize <- 18

hc_gcms$CMIP5_GFDL.ESM2M_rcp85

selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
hc_gcms <- NULL
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)
data_gcms <- NULL
################################################################################
#
################################################################################
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

# nedges_int_hc <- lapply(selection)
# whichsize <- nedges_int_hc[c(numberBN1,numberBN2)]
# whichloglik <- logliks_int_hc[c(numberBN1,numberBN2)]

# filesncep <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",names(listncep),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# filesjra55 <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",names(listjra55),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files5extra <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",namescubics5extra,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files5left <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",namescubics5left,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5.detrended <- paste0("results/propagation/perm3/pos",whichnode,"pos_detrended/prop_",shortnames5,"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended.rda")
files5.rcp85.detrended <- paste0("results/propagation/perm3/pos",whichnode,"pos_detrended/CMIP5_rcp85_detrended/prop_",namescubics5.rcp85,"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended.rda")
files5.future.detrended <- paste0("results/propagation/perm3/pos",whichnode,"pos_detrended/CMIP5_rcp85_detrended/prop_",namescubics.future,"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended.rda")
# files5left.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP5_rcp85/prop_",namescubics5left.rcp85,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files6.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6_ssp585/prop_",namescubics6.ssp585,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files6left.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6_ssp585/prop_",names(cubicsets6left.ssp585),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")



# files6 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6/prop_",namescubics6,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files6left <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6/prop_",names(cubicsets6left),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
filesint.detrended <- paste0("results/propagation/perm3/pos",whichnode,"pos_detrended/prop_",names(listinterim),"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended.rda")
propfiles <- lapply(c(files5.detrended,files5.future.detrended),function(x)load(file = x))
for (i in 1:length(c(files5.detrended,files5.future.detrended))){load(file = c(files5.detrended,files5.future.detrended)[i])}

posposmix <- list()
combinames <- c(shortnames5,namescubics.future)
posposmix <- lapply(combinames, function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended"))))

all.equal(posposmix[[1]],posposmix[[2]])
for(i in 1:length(combinames)){
  posposmix[[i]] <- eval(parse(text = paste0("prop_",combinames[i],"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended")))
}
names(posposmix) <- propfiles
names(cubicsets)
names(posposmix)
length(posposmix)

namescubics5.rcp85
# namescubics5left.rcp85
#ready.prop.ind <- c(namescubics5,namescubics5extra,namescubics5left,namescubics5.rcp85,namescubics6,namescubics6left,namescubics6all.ssp585,names(listinterim),names(listncep),names(listjra55))
# gsub("_hc_1700_1800i_V81_equal2", "", gsub("prop_","",names(posposmix)))%in%names(cubicsets)
length(cubicsets)
length(posposmix)
cubicsets$CMIP5_CanESM2_r1i1p1
cubicsets$CMIP5_IPSL.CM5A.LR_r1i1p1_rcp85

#names(cubicsets[ready.prop.ind])
prop_dif_pos <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, 
                       griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsposclims <- bindGrid(prop_dif_pos, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposclims$Members <- combinames
enspropdifsposclims$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsposclims, names.attr = combinames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


files5.detrended <- paste0("results/propagation/perm3/pos",whichnode,"neg_detrended/propneg_",shortnames5,"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended.rda")
)
# files5extra <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",namescubics5extra,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files5left <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",namescubics5left,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files5.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP5_rcp85/propneg_",namescubics5.rcp85,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5.future.detrended <- paste0("results/propagation/perm3/pos",whichnode,"neg_detrended/CMIP5_rcp85_detrended/propneg_",namescubics.future,"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended.rda")
# files5left.rcp85<- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP5_rcp85/propneg_",namescubics5left.rcp85,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files6.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP6_ssp585/propneg_",namescubics6.ssp585,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files6left.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP6_ssp585/propneg_",names(cubicsets6left.ssp585),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files6 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP6/propneg_",namescubics6,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# files6left <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP6/propneg_",names(cubicsets6left),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# filesint <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listinterim),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# filesncep <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listncep),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
# filesjra55 <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listjra55),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(c(files5.detrended,files5.future.detrended),function(x)load(file = x))
for (i in 1:length(c(files5.detrended,files5.future.detrended))){load(file = c(files5.detrended,files5.future.detrended)[i])}

posnegmix <- lapply(combinames, function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2_detrended"))))
posnegmix
names(posnegmix) <- propfiles


prop_dif_neg <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, 
                       griddata = cubicsets, SIMPLIFY = FALSE)

enspropdifsnegclims <- bindGrid(prop_dif_neg, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsnegclims$Members <- combinames
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsnegclims,names.attr = combinames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_neg[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


names(posposmix)
names(posnegmix)
names(cubicsets)

cubicsets$CMIP5_CanESM2_r1i1p1,
cubicsets$CMIP5_CNRM.CM5_r1i1p1,

cubicsets$CMIP5_IPSL.CM5A.LR_r1i1p1_rcp85,
cubicsets$CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85,
cubicsets$CMIP5_EC.EARTH_r12i1p1_rcp85,
cubicsets$CMIP5_CanESM2_r1i1p1_rcp85
prop_dif_pos_and_neg <- mapply(function(griddata,proppos,propneg) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")), griddata),
                               griddata = cubicsets,
                               proppos = posposmix,propneg = posnegmix, SIMPLIFY = FALSE)
prop_dif_pos_and_neg <- prop_dif_pos_and_neg[sort(names(prop_dif_pos_and_neg))]


# ind6 <- grep("CMIP6",names(prop_dif_pos_and_neg))
# indearth <- grep("EARTH_r1i|EARTH_r2i",names(prop_dif_pos_and_neg))

enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposandnegclims$Members <- sort(combinames)

col <- c(rev(col.b),col.r)
spatialPlot(enspropdifsposandnegclims,names.attr = sort(combinames),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
