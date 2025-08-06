#############################################################################
# Prin Comp Evaluation: Evaluation CMIP5/CMIP6
#############################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(sparsebn)
library(RColorBrewer)
library(gaussDiff)
########################################################################################################
# Model Evaluation CMIP5/CMIP6
########################################################################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("data/tas_JRA55_10d_akima_cubic.rda")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")
load("data/tas_historical_cmip6_left_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
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
cubicsets6 <- tas_historical_cmip6_10d_akima_cubic_corrected
cubicsets6left <- tas_historical_cmip6_left_10d_akima_cubic_corrected
# With CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,cubicsets6,listinterim,listjra55,listncep)
names(cubicsets)

# Without CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,listinterim,listjra55,listncep)
names(cubicsets)
# Without earth
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,listinterim,listjra55,listncep)
names(cubicsets)
# Only CMIP6
cubicsets <- c(cubicsets6,cubicsets6left,listinterim,listjra55,listncep)
names(cubicsets)


namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
namescubics6left <- names(cubicsets6left)
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
hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk, it = it)
hc_gcms6left <- lapply(paste0("CMIP6_left/",namescubics6left), loadIterations, permused = permk, it = it)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_jra55 <- lapply(c("JRA55_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_ncep <- lapply(c("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 
# With cmip6
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_gcms6,hc_interim,hc_jra55,hc_ncep)
# Without cmip6
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_interim,hc_jra55,hc_ncep)
length(hc_gcms)
# Without cmip6 without earth
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_interim,hc_jra55,hc_ncep)
length(hc_gcms)
# only cmip6 
hc_gcms <- c(hc_gcms6,hc_gcms6left,hc_interim,hc_jra55,hc_ncep)
length(hc_gcms)
# names(hc_gcms) <- shortnames
abrev
names(hc_gcms) <- shortnames2
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
##########################################################################
# Global dataset evaluation from BN 
##########################################################################
whichBNFits <- selection_fits_gcms_scaled
whichData <- data_gcms_anom_scaled
whichSet <- "interim_10d_akima_cubic"
whichSet <- "ncep_10d"
whichSet <- "JRA55_10d_akima_cubic"
#########################################################################
# pickedEOFsgrids is prepared under EOFs logliks preparation. 
#####################################################################
EOFS_preparation <- function(whichPCs,EOFs,PCs){
  selected_EOFs <- EOFs[,whichPCs]
  selected_PCs <- PCs[,whichPCs]
  if(is.null(dim(selected_PCs))){
    dim(selected_PCs) <- c(length(selected_PCs),1)
    dim(selected_EOFs) <- c(length(selected_EOFs),1)
  }
  
  projection_EOFs <- matrix(data = 0,nrow = nrow(selected_PCs),ncol = nrow(selected_EOFs))
  for(i in 1:nrow(selected_PCs)){
    for(j in 1:ncol(selected_PCs)){
      projection_EOFs[i,] <- projection_EOFs[i,] + selected_EOFs[,j]*selected_PCs[i,j]
    }
  }
  return(projection_EOFs)
}

gridforEOF <- grid_gcms_anom_scaled
prinCompSet <- prinComp(gridforEOF[[whichSet]], n.eofs = 360)
expvars <- attr(prinCompSet[[1]][[1]],"explained_variance")
EOFs <- prinCompSet[[1]][[1]]$EOFs
PCs <- prinCompSet[[1]][[1]]$PCs

gridforEOF <- grid_gcms_anom_scaled

prinCompSets <- lapply(grid_gcms_anom_scaled, function(x) prinComp(x, n.eofs = 360))
prinCompSets$CMIP5_CanESM2_r1i1p1$tas
prinCompSets$interim_10d_akima_cubic
expvarsSets <- lapply(prinCompSets, function(x)attr(x[[1]][[1]],"explained_variance"))
EOFsSets <- lapply(prinCompSets, function(x) x[[1]][[1]]$EOFs)
PCsSets <- lapply(prinCompSets, function(x) x[[1]][[1]]$PCs)



varpercentages <- c(0.25,0.50,0.75,0.90,0.95,0.99)
# varpercentages <- c(0.50,0.75,0.825,0.90,0.95,0.99)
pickedEOFsSets$CMIP5_CanESM2_r1i1p1
pickedEOFsSets <- lapply(expvarsSets, FUN = function(y)lapply(varpercentages, function(x) which(y > x)[[1]]))
pickedEOFsdataSets <- mapply(function(x,y,z) lapply(x,function(x) EOFS_preparation(1:x,y,z)),x = pickedEOFs,y = EOFsSets, z=PCsSets,SIMPLIFY = FALSE)

pickedEOFsdataSets$CMIP5_CanESM2_r1i1p1

# templatedataSets <- grid_gcms_anom_scaled
templatesdataSets <- lapply(grid_gcms_anom_scaled, function(x) rep(list(x),length(varpercentages)))
templatesdataSets$CMIP5_CanESM2_r1i1p1[[1]]
pickedEOFsgridsSets <- mapply(function(a,b) mapply(function(x,y) {z<- x; z$Data[1,,,] <- y; return(z)},x = a,y = b, SIMPLIFY = FALSE),a= templatesdataSets, b= pickedEOFsdataSets, SIMPLIFY  = FALSE)

pickedEOFsgridsSets$CMIP5_CanESM2_r1i1p1
# pickedEOFsnames <- paste0(whichSet,"_v.exp_",as.character(varpercentages))
# names(pickedEOFsgrids) <- names(pickedEOFsdata) <- pickedEOFsnames

data_pickedEOFsSets <- lapply(pickedEOFsdataSets, function(y)lapply(y, function(x) as.data.frame(x)))
data_pickedEOFsSets$
whichData <- c(data_gcms_anom_scaled,data_pickedEOFs)

ordersets <- TRUE

loglik_selection_datasets_sets <- matrix(data = NA, nrow = length(whichBNFits), ncol = length(varpercentages), dimnames = list(names(whichBNFits),as.character(varpercentages)))
for(i in 1:length(whichBNFits)){
  lo <- sapply(X = data_pickedEOFsSets[[i]], logLik, object = whichBNFits[[i]])
  loglik_selection_datasets_sets[i,] <- lo
}


###############################################################################
# Visualize logliks example set versus amount of variance own dataset explained
###############################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/exampleclustering/logliks_ex_i_vs_variance_ex_i.pdf")
pdf(plotname, width = 14,height = 7)

ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")

par(mfrow = c(1,2))
plot(varpercentages*100,loglik_selection_datasets_sets[ex.subset[1],],
     ylim = c(-150000,40000),col = rainbow(5)[1],type = "l",cex.lab = 1.5,
     xlab = "% variance in data",
     ylab =  "log P(% variance in data| BN from data)",
     )
for(i in 2:length(ex.subset)){
  lines(varpercentages*100,loglik_selection_datasets_sets[ex.subset[i],],col = rainbow(5)[[i]])}
legend("bottomleft",legend = gsub("_r1i1p1","",gsub("CMIP5_","",ex.subset)),pch = 1,col = rainbow(5),cex = 1,bty = "n",title= "data",title.adj = 0.1)

###############################################################################
# Visualize Amount of EOFs needed to explain amount of variance versus amount of variance
###############################################################################

plot(varpercentages*100, pickedEOFsSets[[1]],ylim = c(0,250),col = rainbow(5)[1],type = "l", 
     xlab = "% variance in data",
     ylab =  "|EOFs| needed to explain % variance in data",
      cex.lab = 1.5)
for(i in 2:length(ex.subset)){
  lines(varpercentages*100,pickedEOFsSets[[i]],col = rainbow(5)[i])}
legend("topleft",legend = gsub("_r1i1p1","",gsub("CMIP5_","",ex.subset)),pch = 1,col = rainbow(5),cex = 1,bty = "n",title= "data",title.adj = 0.1)
dev.off()




