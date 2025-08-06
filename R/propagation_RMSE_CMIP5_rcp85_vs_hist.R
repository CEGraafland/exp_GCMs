#############################################################################
# Propgation of evidence CMIP5 rcp85 and CMIP6 ssp585
#############################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
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
load("../Data/tas_ncep_10d.rda")
load("data/tas_JRA55_10d_akima_cubic.rda")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")
load("data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda")
load("data/tas_rcp85_cmip5_left_10d_akima_cubic.rda")
#load("data/tas_historical_cmip6_left_10d_akima_cubic_corrected.rda")
#load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")

###########################################################################################
# Create namelist with models + interim data
###########################################################################################
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"

listjra55 <- list(tas_JRA55_10d_akima_cubic)
names(listjra55) <- "JRA55_10d_akima_cubic"

listncep <- list(tas_ncep_10d)
names(listncep) <- "ncep_10d"

cubicsets5 <- tas_historical_10d_akima_cubic_corrected
cubicsets5extra <- tas_historical_cmip5_extra_10d_akima_cubic
cubicsetsearth <- tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic
cubicsets5left <- tas_historical_cmip5_left_10d_akima_cubic

cubicsets5.rcp85 <- tas_rcp85_cmip5_10d_akima_cubic_corrected
# Leave out bcc "CMIP5_bcc.csm1.1.m_r1i1p1_rcp85" because does not reach 2100 (one year atrasado)
xout <- which(names(tas_rcp85_cmip5_left_10d_akima_cubic) == "CMIP5_bcc.csm1.1.m_r1i1p1_rcp85")
cubicsets5left.rcp85 <- tas_rcp85_cmip5_left_10d_akima_cubic[-xout]


cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsets5.rcp85,cubicsets5left.rcp85,listinterim,listjra55,listncep)
names(cubicsets)

cubicsets.future <- c(cubicsets5.rcp85,cubicsets5left.rcp85)
names(cubicsets.future)


namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics5.rcp85 <- names(cubicsets5.rcp85)
namescubics5left.rcp85 <- names(cubicsets5left.rcp85)


namescubicsearth <- names(cubicsetsearth)
namescubics <- names(cubicsets)
shortnames <- #gsub(
  gsub(names(cubicsets), pattern = "_r1i1p1", replacement = "")
#,
#pattern = "_r12i1p1", replacement = "")
shortnames2 <- names(cubicsets)
abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")

############################################################################
# Function to load hciterations of models in a list 
############################################################################
# pattern <- "CMIP5_CanESM2_r1i1p1"
# permused <- 3
# it= "1700_1800|1800_1900"
# hctype <- NULL
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

# loadIterations("CMIP5_CanESM2_r1i1p1",3,it= "1700_1800")

##################################################################
# Make al datasets standardized
##################################################################
permk <- 3
# Anomalies and standarization over seasons (networks are learned with this set)
data_gcms_out <- lapply(cubicsets,function(x) TimeCoordsAnom_from_Grid_rms(x, rms = TRUE))
data_gcms_out_df <- lapply(data_gcms_out, function(x) as.data.frame(x))
data_gcms_anom_out <- lapply(data_gcms_out, function(x) mat2Dto3Darray(x,attr(x,"Xcoords"), attr(x,"Ycoords")))
grid_gcms_anom <- mapply(function(x,y) {x$Data <- y ;return(x)}, x = cubicsets, y = data_gcms_anom_out,SIMPLIFY = FALSE)
# Global standarazation over the above: FIELD
grid_gcms_anom_scaled <- lapply(grid_gcms_anom, 
                                function(x) scaleGrid(x,
                                                      spatial.frame ="field",
                                                      type = "standardize"))

# Extract the FIELD standarized data
data_gcms_anom_scaled <- lapply(grid_gcms_anom_scaled,
                                function(x) { x <- redim(x,drop = TRUE)
                                xdata <- array3Dto2Dmat(x$Data) 
                                return(xdata)}) 
data_gcms_anom_scaled <- lapply(data_gcms_anom_scaled, function(x) as.data.frame(x))
# adapt to k-th permutation
assign(paste0("data_gcms_anom_scaled_",permk),lapply(data_gcms_anom_scaled, function(x) x[,permutations[[permk]]]))
length(data_gcms_anom_scaled)

it <- "1700_1800"
hc_gcms5 <- lapply(paste0(namescubics5), loadIterations, permused = permk, it = it)
hc_gcms5extra <- lapply(paste0("CMIP5_extra/",namescubics5extra), loadIterations, permused = permk, it = it)
hc_gcms5left<- lapply(paste0("CMIP5_left/",namescubics5left), loadIterations, permused = permk, it = it)
hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk, it = it)
hc_gcms5.rcp85 <- lapply(paste0("FUTURE_",namescubics5.rcp85), loadIterations, permused = permk, it = it)
hc_gcms5left.rcp85<- lapply(paste0("FUTURE_CMIP5_left/FUTURE_",namescubics5left.rcp85), loadIterations, permused = permk, it = it)

# hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk, it = it)
# hc_gcms6left <- lapply(paste0("CMIP6_left/",namescubics6left), loadIterations, permused = permk, it = it)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_jra55 <- lapply(c("JRA55_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_ncep <- lapply(c("../Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 

# Without cmip6 without earth
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcms5.rcp85,hc_gcms5left.rcp85,hc_interim,hc_jra55,hc_ncep)
length(hc_gcms)

abrev
names(hc_gcms) <- shortnames2


# Choose between 'own optimums' or constant magnitude
# modelsize <- own_optimums
modelsize <- 18
selection_hc_gcms <-  mapply(function(x,y) x[[grep((y),x)]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
names(selection_fits_gcms_scaled)


load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_hist_vs_fut.rda")

#############################################################################################
# Emergent constraints: propagation.
#############################################################################################

futures <- grep("rcp85",names(cubicsets))
reans <- length(cubicsets)-3+ 1:3
hists <- (1:length(cubicsets))[-c(futures,reans)]

ids1 <- sapply(names(cubicsets)[hists], function(x)which(x == gsub("_rcp85","",names(cubicsets)[futures])))
ids2 <- sapply(names(cubicsets)[futures], function(x) which(gsub("_rcp85","",x) == names(cubicsets)[hists]))

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

sets.reans <- cubicsets[reans]
sets.hists <- cubicsets[hists][!sapply(ids1,is.integer0)]
names(sets.hists)
sets.futures <- cubicsets[futures][!sapply(ids2,is.integer0)]
names(sets.futures)
all.equal(sort(names(sets.hists)),gsub("_rcp85","",sort(names(sets.futures))))
sets.hists <- sets.hists[sort(names(sets.hists))]
sets.futures <- sets.futures[sort(names(sets.futures))]

whichnode <- "V81"
whichiterations <- it
filesint <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listinterim),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
filesncep <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listncep),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
filesjra55 <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listjra55),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files.reans <- lapply(c(filesint,filesncep,filesjra55),function(x)load(file = x))
for (i in 1:length(c(filesint,filesncep,filesjra55))){load(file = c(filesint,filesncep,filesjra55)[i])}

# posnegmix.reans <- lapply(names(sets.reans), function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
# names(posnegmix.reans)
# names(posnegmix.reans) <- names(files.reans)
# all.equal(posnegmix$propneg_CMIP5_CanESM2_hc_1700_1800i_V81_equal2,
#           posnegmix$propneg_CMIP5_CNRM.CM5_hc_1700_1800i_V81_equal2)
# 
# prop_dif_neg <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata = cubicsets, SIMPLIFY = FALSE)
# enspropdifsnegclims <- bindGrid(prop_dif_neg, dimension = c("member"))
# enspropdifsnegclims$Members <- combinames
# col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
# spatialPlot(enspropdifsnegclims,names.attr = combinames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# 
# 
# 
# prop_dif_pos_and_neg <- mapply(function(proppos,propneg,griddata) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")), griddata),proppos = posposmix,propneg = posnegmix, griddata = cubicsets, SIMPLIFY = FALSE)
# prop_dif_pos_and_neg <- prop_dif_pos_and_neg[sort(names(prop_dif_pos_and_neg))]
# ind6 <- grep("CMIP6",names(prop_dif_pos_and_neg))
# indearth <- grep("EARTH_r1i|EARTH_r2i",names(prop_dif_pos_and_neg))
# 
# enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg[-c(ind6)], dimension = c("member"))
# enspropdifsposandnegclims$Members <- sort(combinames)[-c(ind6)]
# enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg[-c(ind6,indearth)], dimension = c("member"))
# enspropdifsposandnegclims$Members <- sort(combinames)[-c(ind6,indearth)]
# col.b
# sort(combinames[-c(ind6,indearth)])
# col <- c(rev(col.b),col.r)
# spatialPlot(enspropdifsposandnegclims,names.attr = sort(combinames)[-c(ind6,indearth)],as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
#             col.regions= col,
#             set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# 
# prop_dif_pos_and_neg$prop_CMIP5_ACCESS1.0_r1i1p1_hc_1700_1800i_V81_equal2$Data
# prop_dif_pos_and_neg$prop_CMIP5_ACCESS1.3_r1i1p1_hc_1700_1800i_V81_equal2$Data
# 
# 
# 





#######################################################################################
# Automatic V81 All CMIP5 
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


files5.all.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP5_rcp85/prop_",names(sets.futures),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles.all.rcp85 <- lapply(files5.all.rcp85,function(x)load(file = x))
for (i in 1:length(propfiles.all.rcp85)){load(file = files5.all.rcp85[i])}

files5.all.hist <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",names(sets.hists),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles.all.hist <- lapply(files5.all.hist,function(x)load(file = x))
for(i in 1:length(files5.all.hist)){load(file = files5.all.hist[i])}

files.reans <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",names(sets.reans),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles.reans <- lapply(files.reans,function(x)load(file = x))
for(i in 1:length(propfiles.reans)){load(file = files.reans[i])}



posposmix.rcp85 <- list()
posposmix.rcp85 <- lapply(names(sets.futures), function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posposmix.hists <- list()
posposmix.hists <- lapply(propfiles.all.hist, function(x) eval(parse(text = paste0(x))))
posposmix.reans <- list()
posposmix.reans <- lapply(propfiles.reans, function(x) eval(parse(text = paste0(x))))



# for(i in 1:length(combinames)){
#   posposmix[[i]] <- eval(parse(text = paste0("prop_",combinames[i],"_hc_",whichiterations,"i_",whichnode,"_equal2")))
# }
names(posposmix.rcp85) <- propfiles.all.rcp85
names(posposmix.hists) <- propfiles.all.hist
names(posposmix.reans) <- propfiles.reans

prop_dif_pos.rcp85 <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix.rcp85, griddata = sets.futures, SIMPLIFY = FALSE)
enspropdifsposclims.rcp85 <- bindGrid(prop_dif_pos.rcp85, dimension = c("member"))
enspropdifsposclims.rcp85$Members <-names(sets.futures)
enspropdifsposclims.rcp85$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsposclims.rcp85, names.attr = names(sets.futures),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.rcp85[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop_dif_pos.hists <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix.hists, griddata = sets.hists, SIMPLIFY = FALSE)
enspropdifsposclims.hists <- bindGrid(prop_dif_pos.hists, dimension = c("member"))
enspropdifsposclims.hists$Members <-names(sets.hists)
enspropdifsposclims.hists$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsposclims.hists, names.attr = names(sets.hists),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.hists[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop_dif_pos.reans <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix.reans, griddata = sets.reans, SIMPLIFY = FALSE)
enspropdifsposclims.reans <- bindGrid(prop_dif_pos.reans, dimension = c("member"))
enspropdifsposclims.reans$Members <-names(sets.reans)
enspropdifsposclims.reans$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsposclims.reans, names.attr = names(sets.reans),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.reans[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))




files5.all.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP5_rcp85/propneg_",names(sets.futures),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles.all.rcp85 <- lapply(files5.all.rcp85,function(x)load(file = x))
for (i in 1:length(propfiles.all.rcp85)){load(file = files5.all.rcp85[i])}

files5.all.hists <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(sets.hists),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles.all.hist <- lapply(files5.all.hists,function(x)load(file = x))
for (i in 1:length(propfiles.all.hist)){load(file = files5.all.hists[i])}

files.reans <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(sets.reans),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles.reans <- lapply(files.reans,function(x)load(file = x))
for (i in 1:length(propfiles.reans)){load(file = files.reans[i])}


posnegmix.rcp85 <- lapply(names(sets.futures), function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
names(posnegmix.rcp85) <- propfiles.all.rcp85

posnegmix.hists <- lapply(propfiles.all.hist, function(x) eval(parse(text = paste0(x))))
names(posnegmix.hists) <- propfiles.all.hist

posnegmix.reans <- lapply(propfiles.reans, function(x) eval(parse(text = paste0(x))))
names(posnegmix.reans) <- propfiles.reans

prop_dif_neg.rcp85 <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix.rcp85, griddata = sets.futures, SIMPLIFY = FALSE)
enspropdifsnegclims.rcp85 <- bindGrid(prop_dif_neg.rcp85, dimension = c("member"))
enspropdifsnegclims.rcp85$Members <- names(sets.futures)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsnegclims.rcp85,names.attr = names(sets.futures),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_neg.rcp85[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop_dif_neg.hists <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix.hists, griddata = sets.hists, SIMPLIFY = FALSE)
enspropdifsnegclims.hists <- bindGrid(prop_dif_neg.hists, dimension = c("member"))
enspropdifsnegclims.hists$Members <- names(sets.hists)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsnegclims.hists,names.attr = names(sets.hists),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_neg.hists[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop_dif_neg.reans <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix.reans, griddata = sets.reans, SIMPLIFY = FALSE)
enspropdifsnegclims.reans <- bindGrid(prop_dif_neg.reans, dimension = c("member"))
enspropdifsnegclims.reans$Members <- names(sets.reans)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsnegclims.reans,names.attr = names(sets.reans),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_neg.reans[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))



prop_dif_pos_and_neg.rcp85 <- mapply(function(proppos,propneg,griddata) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),proppos = posposmix.rcp85,propneg = posnegmix.rcp85, griddata = sets.futures, SIMPLIFY = FALSE)
enspropdifsposandnegclims.rcp85 <- bindGrid(prop_dif_pos_and_neg.rcp85, dimension = c("member"))

prop_dif_pos_and_neg.hists <- mapply(function(proppos,propneg,griddata) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),proppos = posposmix.hists,propneg = posnegmix.hists, griddata = sets.hists, SIMPLIFY = FALSE)
enspropdifsposandnegclims.hists <- bindGrid(prop_dif_pos_and_neg.hists, dimension = c("member"))

prop_dif_pos_and_neg.reans <- mapply(function(proppos,propneg,griddata) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),proppos = posposmix.reans,propneg = posnegmix.reans, griddata = sets.reans, SIMPLIFY = FALSE)
enspropdifsposandnegclims.reans <- bindGrid(prop_dif_pos_and_neg.reans, dimension = c("member"))

enspropdifsposandnegclims.rcp85$Members <- names(sets.futures)
enspropdifsposandnegclims.hists$Members <- names(sets.hists)
enspropdifsposandnegclims.reans$Members <- names(sets.reans)
col.b

col <- c(rev(col.b),col.r)


#cmip5
spatialPlot(enspropdifsposandnegclims.rcp85,names.attr =names(sets.futures),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.rcp85[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg.rcp85[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
spatialPlot(enspropdifsposandnegclims.hists,names.attr =names(sets.hists),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.hists[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg.hists[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

#############################################################################################
# Check if hellinger coefficients are diagnostic for RMSE between propagation fields
#############################################################################################
load("results/hellinger_coefficient/hellinger_coefficients.rda")
load("results/hellinger_coefficient/hellinger_coefficients.rcp85.rda")
load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_hist_vs_fut.rda")
load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_CMIP6_hist_vs_fut.rda")

bd.cmips.fut.vs.hist <- diag(-log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(sets.hists),names(sets.futures)])
names(bd.cmips.fut.vs.hist) <- names(sets.hists)

bd.cmips.fut.vs.ncep <- -log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(sets.futures),"ncep_10d"]
bd.cmips.hist.vs.int <- -log(hellinger_coefficients)[names(sets.hists),"interim_10d_akima_cubic"]
bd.cmips.hist.fut.vs.int <- -log(hellinger_coefficients)[c(names(sets.hists),names(sets.futures)),"interim_10d_akima_cubic"]
bd.cmips.hist.vs.ncep <- -log(hellinger_coefficients)[names(sets.hists),"ncep_10d"]
bd.cmips.hist.vs.JRA55 <- -log(hellinger_coefficients)[names(sets.hists),"JRA55_10d_akima_cubic"]
bd.cmips.reans <- -log(hellinger_coefficients)[c(names(sets.hists),names(sets.futures),names(sets.reans)),c(names(sets.hists),names(sets.futures),names(sets.reans))]

rownames(hellinger_coefficients)
keltoC <-function(x)(x-273)%%273
gtimes.grid <- function(x){
  # x <- sets.hists$CMIP5_ACCESS1.0_r1i1p1
  m <- x$Data
  v <- as.vector(cos(x$xyCoords$y/(180)*pi))
  for(i in 1:dim(x$Data)[2]){
    m[,i,] <- x$Data[,i,]*v[i]}
  weighted.gtimeseries <- apply(m, MARGIN= c(1),FUN =sum)/sum(v*dim(x$Data)[3])
  x$Data <- m
  return(x)
}

weighted.props.cmip5s.rcp85 <- lapply(prop_dif_pos_and_neg.rcp85,gtimes.grid)
weighted.props.cmip5s.hists <- lapply(prop_dif_pos_and_neg.hists,gtimes.grid)
weighted.props.cmip5s.reans <- lapply(prop_dif_pos_and_neg.reans,gtimes.grid)

RMSE.deltaprop.w <- numeric(length(weighted.props.cmip5s.rcp85))
for (i in 1:length(weighted.props.cmip5s.rcp85)){
  x1 <- redim(weighted.props.cmip5s.rcp85[[i]],drop = TRUE)$Data
  x2 <- redim(weighted.props.cmip5s.hists[[i]],drop = TRUE)$Data
  
  RMSE.deltaprop.w[i] <- sqrt(sum((x1-x2)^2)/length(x1))
}

RMSE.intprop.w <- numeric(length(weighted.props.cmip5s.hists))
for (i in 1:length(weighted.props.cmip5s.hists)){
  x1 <- redim(weighted.props.cmip5s.hists[[i]],drop = TRUE)$Data
  x2 <- redim(weighted.props.cmip5s.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2,drop = TRUE)$Data
  
  RMSE.intprop.w[i] <- sqrt(sum((x1-x2)^2)/length(x1))
}
names(RMSE.intprop.w)<- names(sets.hists)

cor.deltaprop.w <- numeric(length(weighted.props.cmip5s.rcp85))
for (i in 1:length(weighted.props.cmip5s.rcp85)){
  x1 <- redim(weighted.props.cmip5s.rcp85[[i]],drop = TRUE)$Data
  x2 <- redim(weighted.props.cmip5s.hists[[i]],drop = TRUE)$Data
  
  cor.deltaprop.w[i] <- cor(as.vector(x1),as.vector(x2))
}
names(cor.deltaprop.w)<- names(sets.hists)

weighted.props <- c(weighted.props.cmip5s.hists,weighted.props.cmip5s.rcp85,weighted.props.cmip5s.reans)
weighted.props <- c(prop_dif_pos_and_neg.hists,prop_dif_pos_and_neg.rcp85,prop_dif_pos_and_neg.reans)
weighted.props <- c(prop_dif_pos.hists,prop_dif_pos.rcp85,prop_dif_pos.reans)
weighted.props <- c(prop_dif_neg.hists,prop_dif_neg.rcp85,prop_dif_neg.reans)
cor.prop.w <- matrix(nrow = length(weighted.props),ncol = length(weighted.props))
for(j in 1:length(weighted.props)){
for (i in 1:length(weighted.props)){
  x1 <- redim(weighted.props[[i]],drop = TRUE)$Data
  x2 <- redim(weighted.props[[j]],drop = TRUE)$Data
  
  cor.prop.w[i,j] <- cor(as.vector(x1),as.vector(x2),method = "pearson")
}
}
rownames(cor.prop.w)<- colnames(cor.prop.w)<-c(names(sets.hists),names(sets.futures),names(sets.reans))


cor.intprop.w <- numeric(length(weighted.props.cmip5s.hists))
for (i in 1:length(weighted.props.cmip5s.hists)){
  x1 <- redim(weighted.props.cmip5s.hists[[i]],drop = TRUE)$Data
  x2 <- redim(weighted.props.cmip5s.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2,drop = TRUE)$Data
  
  cor.intprop.w[i] <- cor(as.vector(x1),as.vector(x2))
}
names(cor.intprop.w )<- names(sets.hists)

weighted.props.cmip5s <- c(weighted.props.cmip5s.hists,weighted.props.cmip5s.rcp85)
weighted.props.cmip5s <- c(prop_dif_pos_and_neg.hists,prop_dif_pos_and_neg.rcp85)
weighted.props.cmip5s <- c(prop_dif_pos.hists,prop_dif_pos.rcp85)
weighted.props.cmip5s <- c(prop_dif_neg.hists,prop_dif_neg.rcp85)

x2 <- redim(weighted.props.cmip5s.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2,drop = TRUE)$Data
x2 <- redim(prop_dif_pos_and_neg.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2,drop = TRUE)$Data
x2 <- redim(prop_dif_pos.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2,drop = TRUE)$Data
x2 <- redim(prop_dif_neg.reans$propneg_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2,drop = TRUE)$Data

cor.intprop.w <- numeric(length(weighted.props.cmip5s))
for (i in 1:length(weighted.props.cmip5s)){
  x1 <- redim(weighted.props.cmip5s[[i]],drop = TRUE)$Data
  cor.intprop.w[i] <- cor(as.vector(x1),as.vector(x2))
}

names(cor.intprop.w )<- c(names(sets.hists),names(sets.futures))

############################################################################
# diagnostic (1): RMSE Relationship Bhatacharrya distance fut - hist
# and probability propagation fields fut - hist.
############################################################################
df1 <- data.frame("RMSE.change" = RMSE.deltaprop.w,"Bchange" = bd.cmips.fut.vs.hist)
lm1 <-lm(RMSE.change~Bchange,df1)
summary(lm1)
Bvalues <- seq(20, 40, 1)
predictedBvalues <- predict(lm1,list(Bchange=Bvalues))

plot(bd.cmips.fut.vs.hist,RMSE.deltaprop.w,xlab = expression("Bscore: d"[B]*"(BN"[fut]*",BN"[hist]*")"),pch = 16,
     ylab = expression("RMSE(P(X|X"[e]*")"[fut]*",P(X|X"[e]*")"[hist]*")"), main = "Relationship Bhatacharrya distance and probability propagation fields: Future vs history")

lines(Bvalues,predictedBvalues, col = "red", lwd = 2)
lines(x=c(20,40) , y = c(lm1$coefficients[2]*20+lm1$coefficients[1],lm1$coefficients[2]*40+lm1$coefficients[1]),col = 'red')
text(20,0.30, paste("R^2 = ",round(summary(lm1)$r.squared,3)))
text(20,0.28,paste("cor = ",round(cor(bd.cmips.fut.vs.hist,RMSE.deltaprop.w,method = "pearson"),3)))

cor(bd.cmips.fut.vs.hist,RMSE.deltaprop.w,method = "pearson")
############################################################################
# diagnostic (2): RMSE Relationship Bhatacharrya distance hist - era interim
# and probability propagation fields hist - era interim
############################################################################
df2 <- data.frame("RMSE.score" = RMSE.intprop.w,"Bscore" = bd.cmips.hist.vs.int)
lm2 <-lm(RMSE.score~Bscore,df2)
summary(lm2)

plot(bd.cmips.hist.vs.int,RMSE.intprop.w, xlab = expression("Bscore: d"[B]*"(BN"[hist]*",BN"[interim]*")"),pch = 16,
     ylab = expression("RMSE(P(X|X"[e]*")"[hist]*",P(X|X"[e]*")"[interim]*")"), main = "Relationship Bhatacharrya distance and probability propagation fields: history vs ERA-Interim")
lines(x=c(30,55) , y = c(lm2$coefficients[2]*30+lm2$coefficients[1],lm2$coefficients[2]*55+lm2$coefficients[1]),col = 'red')
text(30,0.18, paste("R^2 = ",round(summary(lm2)$r.squared,3)))
text(30,0.16,paste("cor = ",round(cor(bd.cmips.hist.vs.int,RMSE.intprop.w,method = "pearson"),3)))

cor(bd.cmips.hist.vs.int,RMSE.intprop.w)
text(bd.cmips.hist.vs.int,RMSE.intprop.w,names(RMSE.intprop.w))

############################################################################
# diagnostic (4): Correlation Bhatacharrya distance hist,  - era interim
# and probability propagation fields hist - era interim
############################################################################
df4 <- data.frame("cor.score" = cor.intprop.w[names(sets.hists)],"Bscore" = bd.cmips.hist.vs.int)
lm4 <-lm(cor.score~Bscore,df4)
summary(lm4)

plot(bd.cmips.hist.vs.int,cor.intprop.w[names(sets.hists)], xlab = expression("Bscore: d"[B]*"(BN"[hist]*",BN"[interim]*")"),pch = 16,
     ylab = expression("cor(P(X|X"[e]*")"[hist]*",P(X|X"[e]*")"[interim]*")"), main = "Relationship Bhatacharrya distance and probability propagation fields: history vs ERA-Interim")
lines(x=c(30,80) , y = c(lm6$coefficients[2]*30+lm6$coefficients[1],lm6$coefficients[2]*80+lm6$coefficients[1]),col = 'red')
text(30,0.75, paste("R^2 = ",round(summary(lm4)$r.squared,3)))
text(30,0.74,paste("cor = ",round(cor(bd.cmips.hist.vs.int,cor.intprop.w[names(sets.hists)],method = "pearson"),3)))

cor(bd.cmips.hist.vs.int,cor.intprop.w[names(sets.hists)], method = "pearson")
text(bd.cmips.hist.fut.vs.int,cor.intprop.w,names(cor.intprop.w))

############################################################################
# diagnostic (5): CORRELATION Bhatacharrya distance fut - hist
# and probability propagation fields fut - hist.
############################################################################
df5 <- data.frame("cor.change" = cor.deltaprop.w,"Bchange" = bd.cmips.fut.vs.hist)
lm5 <-lm(cor.change~Bchange,df5)
summary(lm5)
Bvalues <- seq(20, 40, 1)
predictedBvalues <- predict(lm5,list(Bchange=Bvalues))

plot(bd.cmips.fut.vs.hist,cor.deltaprop.w,xlab = expression("Bscore: d"[B]*"(BN"[fut]*",BN"[hist]*")"),pch = 16,
     ylab = expression("cor(P(X|X"[e]*")"[fut]*",P(X|X"[e]*")"[hist]*")"), main = "Relationship Bhatacharrya distance and probability propagation fields: Future vs history")

lines(Bvalues,predictedBvalues, col = "red", lwd = 2)
lines(x=c(20,40) , y = c(lm5$coefficients[2]*20+lm5$coefficients[1],lm5$coefficients[2]*40+lm5$coefficients[1]),col = 'red')
text(20,0.7, paste("R^2 = ",round(summary(lm5)$r.squared,3)))
text(20,0.6,paste("cor = ",round(cor(bd.cmips.fut.vs.hist,cor.deltaprop.w,method = "pearson"),3)))

cor(bd.cmips.fut.vs.hist,cor.deltaprop.w,method = "pearson")

############################################################################
# diagnostic (6): Correlation Bhatacharrya distance hist, fut - era interim
# and probability propagation fields hist - era interim
############################################################################
df6 <- data.frame("RMSE.score" = cor.intprop.w,"Bscore" = bd.cmips.hist.fut.vs.int)
lm6 <-lm(RMSE.score~Bscore,df6)
summary(lm6)

plot(bd.cmips.hist.fut.vs.int,cor.intprop.w, xlab = expression("Bscore: d"[B]*"(BN"[hist]*",BN"[interim]*")"),pch = 16,
     ylab = expression("cor(P(X|X"[e]*")"[hist]*",P(X|X"[e]*")"[interim]*")"), main = "Relationship Bhatacharrya distance and probability propagation fields: history vs ERA-Interim")
lines(x=c(30,80) , y = c(lm6$coefficients[2]*30+lm6$coefficients[1],lm6$coefficients[2]*80+lm6$coefficients[1]),col = 'red')
text(30,0.75, paste("R^2 = ",round(summary(lm6)$r.squared,3)))
text(30,0.74,paste("cor = ",round(cor(bd.cmips.hist.fut.vs.int,cor.intprop.w,method = "pearson"),3)))

cor(bd.cmips.hist.fut.vs.int,cor.intprop.w, method = "pearson")
text(bd.cmips.hist.fut.vs.int,cor.intprop.w,names(cor.intprop.w))
############################################################################
# diagnostic (7): Correlation Bhatacharrya distance hist,fut 
# and probability propagation fields hist,fut 
############################################################################
cor.prop.w.d <- cor.prop.w
bhat.gcms.d <-bd.cmips.reans
diag(bhat.gcms.d) <- diag(cor.prop.w.d)<-NA
bhat.gcms.d[which(bhat.gcms.d== Inf, arr.ind = TRUE)]<- NA

upper.tri()
df7 <- data.frame("cor.scores" = as.vector(cor.prop.w.d),"Bscore" = as.vector(bhat.gcms.d))
lm7 <-lm(cor.scores~Bscore,df7)
summary(lm7)
# bd

plot(bhat.gcms.d,cor.prop.w.d,  xlab = expression("Bscore: d"[B]*"(BN"[modelA]*",BN"[modelB]*")"),pch = 16,
     ylab = expression("cor(P(X|X"[e]*")"[modelA]*",P(X|X"[e]*")"[modelB]*")"), main = "Relationship Bhatacharrya distance and probability propagation fields: history vs ERA-Interim")
lines(x=c(10,100) , y = c(lm7$coefficients[2]*10+lm7$coefficients[1],lm7$coefficients[2]*100+lm7$coefficients[1]),col = 'red')
text(10,0.75, paste("R^2 = ",round(summary(lm7)$r.squared,3)))
text(10,0.74,paste("cor = ",round(cor(df7$cor.scores, df7$Bscore, method = "pearson", use = "pairwise.complete.obs"),3)))

points(as.vector(bd.cmips.reans["interim_10d_akima_cubic" ,c(names(sets.hists))]),as.vector(cor.prop.w["interim_10d_akima_cubic" ,c(names(sets.hists))])
       ,pch = 16, col = "red")
cor(as.vector(bd.cmips.reans["interim_10d_akima_cubic" ,c(names(sets.hists))]),as.vector(cor.prop.w["interim_10d_akima_cubic" ,c(names(sets.hists))]),use = "pairwise.complete.obs")

points(as.vector(bd.cmips.reans["interim_10d_akima_cubic" ,c(names(sets.futures))]),as.vector(cor.prop.w["interim_10d_akima_cubic" ,c(names(sets.futures))])
       ,pch = 16, col = "blue")


cor(bhat.gcms.d[upper.tri(bhat.gcms.d)],cor.prop.w[upper.tri(cor.prop.w)],use = "pairwise.complete.obs",method = "pearson")
cor(as.vector(bd.cmips.reans["interim_10d_akima_cubic" ,c(names(sets.hists))]),as.vector(cor.prop.w["interim_10d_akima_cubic" ,c(names(sets.hists))]),use = "pairwise.complete.obs")
cor(as.vector(bd.cmips.reans["interim_10d_akima_cubic",c(names(sets.hists),names(sets.futures))]),as.vector(cor.prop.w["interim_10d_akima_cubic",c(names(sets.hists),names(sets.futures))]),use = "pairwise.complete.obs")
names(sets.reans)

cor(as.vector(bd.cmips.reans),as.vector(cor.prop.w), method = "pearson")
text(bd.cmips.hist.fut.vs.int,cor.intprop.w,names(cor.intprop.w))

##############################################################################
# diagnostics (3) en (4)
##############################################################################
plot(RMSE.deltaprop.w,RMSE.intprop.w)
cor(RMSE.deltaprop.w,RMSE.intprop.w)

plot(bd.cmips.hist.vs.int,RMSE.deltaprop.w)
cor(bd.cmips.hist.vs.int,RMSE.deltaprop.w)
##############################################################################
# Order w.r.t. interim
##############################################################################
fut.int.order <- order(bd.cmips.reans["interim_10d_akima_cubic",names(sets.futures)])
### weighted probability climatologies 
ens.weighted.props.cmip5s.rcp85 <- bindGrid(weighted.props.cmip5s.rcp85[fut.int.order], dimension = c("member"))

names.fut <- gsub("CMIP5_","",gsub("r1i1p1_rcp85","",names(sets.futures)))
spatialPlot(ens.weighted.props.cmip5s.rcp85,names.attr =names.fut[fut.int.order],as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.rcp85[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg.rcp85[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 1)))
### unweighted positive probability climatologies 
ens.propdifsposclims.rcp85 <- bindGrid(prop_dif_pos.rcp85[fut.int.order], dimension = c("member"))
ens.propdifsposclims.rcp85$Members <-names(sets.futures)[fut.int.order]
ens.propdifsposclims.rcp85$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(ens.propdifsposclims.rcp85, names.attr = names(sets.futures)[fut.int.order],as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.rcp85[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

### unweighted positive and negative probability climatologies 
ens.propdifsposandnegclims.rcp85 <- bindGrid(prop_dif_pos_and_neg.rcp85[fut.int.order], dimension = c("member"))
ens.propdifsposandnegclims.rcp85$Members <- names(sets.futures)[fut.int.order]

plotname <-  "figs/weights/orderpropagation/propagation_fut_vs_interim_wrt_bhat.pdf"
pdf(plotname,width= 10, height= 20)

spatialPlot(ens.propdifsposandnegclims.rcp85,names.attr =names(sets.futures)[fut.int.order],as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.rcp85[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg.rcp85[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
######################################################################################################################################
hist.int.order <- order(bd.cmips.reans["interim_10d_akima_cubic",names(sets.hists)])
### weighted probability climatologies
ens.weighted.props.cmip5s.hists <- bindGrid(c(list(weighted.props.cmip5s.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2),weighted.props.cmip5s.hists[hist.int.order]), dimension = c("member"))

names.hists<- gsub("CMIP5_","",gsub("r1i1p1","",names(sets.hists)))
spatialPlot(ens.weighted.props.cmip5s.hists,names.attr =c("interim",names.hists[hist.int.order]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.rcp85[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg.rcp85[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 1)))

### unweighted positive probability climatologies 
ens.propdifsposclims.hists <- bindGrid(c(list(prop_dif_pos.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2),prop_dif_pos.hists[hist.int.order]), dimension = c("member"))
ens.propdifsposclims.hists$Members <-names(sets.hists)[hist.int.order]
ens.propdifsposclims.hists$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(ens.propdifsposclims.hists, names.attr = c("interim",names(sets.hists)[fut.int.order]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.hists[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

### unweighted positive and negative probability climatologies 
ens.propdifsposandnegclims.hists <- bindGrid(c(list(prop_dif_pos_and_neg.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2),prop_dif_pos_and_neg.hists[hist.int.order]), dimension = c("member"))
ens.propdifsposandnegclims.hists$Members <-names(sets.hists)[hist.int.order]
ens.propdifsposandnegclims.hists$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)

plotname <-  "figs/weights/orderpropagation/propagation_hist_vs_interim_wrt_bhat.pdf"
pdf(plotname,width= 10, height= 20)
spatialPlot(ens.propdifsposandnegclims.hists,names.attr =c("interim",names(sets.hists)[hist.int.order]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos.hists[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg.hists[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
dev.off()
##############################################################################################################
# 5 per cent 
##############################################################################################################

col.blue <- rev(brewer.pal(5,"Blues"))
col.blue
col.red <- brewer.pal(5,"Reds")

col.div <- c(col.blue,"white","white", col.red)
col.div[4:5] <- "white"
col.b <- c(col.blue,col.red)

dataRMS2_hist <- data_gcms_anom_scaled[c("interim_10d_akima_cubic",names(sets.hists[hist.int.order]))]

CMIP5_data_5<- lapply(dataRMS2_hist, function(x) quantile(x$V81, c(0.95)))

Big2ind_hist <- mapply(function(x,y)which(x$V81 >= y), x = dataRMS2_hist, y = CMIP5_data_5, SIMPLIFY = FALSE)

extremes_hist <- mapply(function(x,z) colMeans(x[z,]), x = dataRMS2_hist, z = Big2ind_hist, SIMPLIFY = FALSE)

dataRMS2_hist$CMIP5_CNRM.CM5_r1i1p1[Big2ind_hist$CMIP5_CNRM.CM5_r1i1p1,]

climBig2_hist <- mapply(function(x,w) quantity2clim(w, what = "up 95",x),
                        x = cubicsets[c("interim_10d_akima_cubic",names(sets.hists[hist.int.order]))],w = extremes_hist, SIMPLIFY = FALSE)
climBig2_hist.ens <- bindGrid(climBig2_hist, dimension = "member")
climBig2_hist.ens$Members <-c("interim",names(sets.hists)[hist.int.order])

plotname <-  "figs/weights/orderpropagation/climatology_hist_interim_V81_Q95.pdf"
pdf(plotname,width= 10, height= 20)
spatialPlot(climBig2_hist.ens, lonCenter = 180,
            names.attr =c("interim",names(sets.hists)[hist.int.order]), 
            main = "climatology 5% warmest months V81",
            backdrop.theme = "coastline", rev.colors = TRUE, as.table = TRUE,
            col.regions = col.div,set.min = -3, set.max = 3,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)),
            at = seq(-3,3,0.5))
dev.off()
############################################################################################################################
# 
############################################################################################################################
dataRMS2_fut <- data_gcms_anom_scaled[names(sets.futures[fut.int.order])]

CMIP5_data_fut_5<- lapply(dataRMS2_fut, function(x) quantile(x$V81, c(0.95)))

Big2ind_fut <- mapply(function(x,y)which(x$V81 >= y), x = dataRMS2_fut, y = CMIP5_data_fut_5, SIMPLIFY = FALSE)

extremes_fut <- mapply(function(x,z) colMeans(x[z,]), x = dataRMS2_fut, z = Big2ind_fut, SIMPLIFY = FALSE)

climBig2_fut <- mapply(function(x,w) quantity2clim(w, what = "up 95",x),
                        x = cubicsets[names(sets.futures[fut.int.order])],w = extremes_fut, SIMPLIFY = FALSE)
climBig2_fut.ens <- bindGrid(climBig2_fut, dimension = "member")
climBig2_fut.ens$Members <-names(sets.futures)[fut.int.order]

plotname <-  "figs/weights/orderpropagation/climatology_fut_interim_V81_Q95.pdf"
pdf(plotname,width= 10, height= 20)
spatialPlot(climBig2_fut.ens, lonCenter = 180,
            names.attr =names(sets.futures)[fut.int.order], 
            main = "climatology 5% warmest months V81 RCP85",
            backdrop.theme = "coastline", rev.colors = TRUE, as.table = TRUE,
            col.regions = col.div,set.min = -3, set.max = 3,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)),
            at = seq(-3,3,0.5))
dev.off()
############################################################################################################################
