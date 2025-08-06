##################################################################################################
# Propagation of evidence ssp585 2071-2100
##################################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list=ls())
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
load("data/tas_ssp585_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_ssp585_cmip6_left_10d_akima_cubic_corrected.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models
###########################################################################################
cubicsets6.ssp585 <- tas_ssp585_cmip6_10d_akima_cubic_corrected
cubicsets6left.ssp585 <- tas_ssp585_cmip6_left_10d_akima_cubic_corrected

cubicsets.future <- c(cubicsets6.ssp585,cubicsets6left.ssp585)


names(cubicsets.future)

institutions.overlap <- c("IPSL","CNRM","CESM","EC","CanESM5","MRI","MIROC","MPI","GFDL","Nor")
institutions <- rep(institutions.overlap,each =1)
CMIPS.6 <- rep(c("CMIP6"),length(institutions.overlap))
combinations.6 <- cbind(CMIPS.6,institutions)
namesort.6 <- unlist(apply(combinations.6,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),names(cubicsets.future))))
namesort.6.remaining <- (1:length(cubicsets.future))[-namesort.6]

# TAKE 1
# cubicsets.future <- cubicsets.future[namesort.6]
# TAKE 2 
cubicsets.future <- cubicsets.future[namesort.6.remaining]


namescubics6.ssp585 <- names(cubicsets6.ssp585)
namescubics6left.ssp585 <- names(cubicsets6left.ssp585)

shortnames.future <- gsub("CMIP6_",gsub(names(cubicsets.future), pattern = "_r1i1p1", replacement = ""),replacement = "")

#,
#pattern = "_r12i1p1", replacement = "")
# shortnames2 <- names(cubicsets)
abrev <-  gsub("_akima_cubic",gsub("CMIP6_",shortnames.future,replacement = ""),replacement = "")


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
data.future_gcms_out <- lapply(cubicsets.future,function(x) TimeCoordsAnom_from_Grid_rms(x, rms = TRUE))
data.future_gcms_out_df <- lapply(data.future_gcms_out, function(x) as.data.frame(x))
data.future_gcms_anom_out <- lapply(data.future_gcms_out, function(x) mat2Dto3Darray(x,attr(x,"Xcoords"), attr(x,"Ycoords")))
grid.future_gcms_anom <- mapply(function(x,y) {x$Data <- y ;return(x)}, x = cubicsets.future, y = data.future_gcms_anom_out,SIMPLIFY = FALSE)
# Global standarazation over the above: FIELD
grid.future_gcms_anom_scaled <- lapply(grid.future_gcms_anom, 
                                       function(x) scaleGrid(x,
                                                             spatial.frame ="field",
                                                             type = "standardize"))

# Extract the FIELD standarized data
data.future_gcms_anom_scaled <- lapply(grid.future_gcms_anom_scaled,
                                       function(x) { x <- redim(x,drop = TRUE)
                                       xdata <- array3Dto2Dmat(x$Data) 
                                       return(xdata)}) 
data.future_gcms_anom_scaled <- lapply(data.future_gcms_anom_scaled, function(x) as.data.frame(x))


it <- "1700_1800"
# hc_gcms5.rcp85 <- lapply(paste0("FUTURE_",namescubics5.rcp85), loadIterations, permused = permk, it = it)
# hc_gcms5left.rcp85<- lapply(paste0("FUTURE_CMIP5_left/FUTURE_",namescubics5left.rcp85), loadIterations, permused = permk, it = it)
hc_gcms6.ssp585 <- lapply(paste0("FUTURE_CMIP6/FUTURE_",namescubics6.ssp585), loadIterations, permused = permk, it = it)
hc_gcms6left.ssp585 <- lapply(paste0("FUTURE_CMIP6_left/FUTURE_",namescubics6left.ssp585), loadIterations, permused = permk, it = it)
# hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk, it = it)
# hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk, it = it)




hc_gcms.future <- c(hc_gcms6.ssp585,hc_gcms6left.ssp585)
names(hc_gcms.future) <- c(namescubics6.ssp585,namescubics6left.ssp585)


# TAKE 1
# namesort.hc_gcms.6 <- unlist(apply(combinations.6,MARGIN = 1,function(x) grep(paste0(x[2]),names(hc_gcms.future))))
# namesort.hc_gcms.6 == namesort.6
# hc_gcms.future <- hc_gcms.future[namesort.hc_gcms.6]

# TAKE 2
hc_gcms.future <- hc_gcms.future[namesort.6.remaining]


# Choose between 'own optimums' or constant magnitude
# modelsize <- own_optimums
modelsize <- 18
selection_hc_gcms.future <-  mapply(function(x,y) x[[grep((y),x)]], x = hc_gcms.future, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms.future <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms.future, y = data.future_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms.future_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms.future, y = data.future_gcms_anom_scaled, SIMPLIFY = FALSE)

####################################################################################
# Propagation V81
####################################################################################
#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms.future_scaled)){
  
  assign(paste0("prop_",names(hc_gcms.future) [i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms.future_scaled[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",names(hc_gcms.future) [i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV81pos/CMIP6_ssp585/prop_",names(hc_gcms.future) [i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2.rda"))
  
}
#################################################################################
# Negative evidence. V81 (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms.future)){
  
  assign(paste0("propneg_",names(hc_gcms.future)[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms.future[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",names(hc_gcms.future)[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV81neg/CMIP6_ssp585/propneg_",names(hc_gcms.future)[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2.rda"))
  
}