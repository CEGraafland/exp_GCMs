#############################################################################
# Model Evaluation CMIP5 CMIP6 Treeseg
#############################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
#############################################################################
rm(list = ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
#install.packages("sparsebnUtils")
#install.packages("sparsebn")
library(sparsebnUtils)
library(sparsebn)
library(RColorBrewer)
library(gaussDiff)
library(visualizeR)
########################################################################################################
# Model Evaluation CMIP5/CMIP6
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
load("../Data/tas_ncep_10d.rda")
load("data/tas_JRA55_10d_akima_cubic.rda")
load("data/tas_interim_10d_akima_cubic.rda")

load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")

load("data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda")
load("data/tas_rcp85_cmip5_left_10d_akima_cubic.rda")

load("data/tas_historical_cmip6_left_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")

load("data/tas_ssp585_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_ssp585_cmip6_left_10d_akima_cubic_corrected.rda")

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


cubicsets6 <- tas_historical_cmip6_10d_akima_cubic_corrected
cubicsets6left <- tas_historical_cmip6_left_10d_akima_cubic_corrected
cubicsets6all<- c(cubicsets6,cubicsets6left)

cubicsets6.ssp585 <- tas_ssp585_cmip6_10d_akima_cubic_corrected
cubicsets6left.ssp585 <- tas_ssp585_cmip6_left_10d_akima_cubic_corrected
cubicsets6all.ssp585 <- c(cubicsets6.ssp585,cubicsets6left.ssp585)

ssp585.out <-gsub("_ssp585","",names(cubicsets6all.ssp585)) %in% gsub("_historical","",names(cubicsets6all))
cubicsets6all.ssp585 <-cubicsets6all.ssp585[ssp585.out]
hist6.out<-gsub("_historical","",names(cubicsets6all)) %in% gsub("_ssp585","",names(cubicsets6all.ssp585))
cubicsets6all <- cubicsets6all[hist6.out]
all.equal(sort(gsub("_historical","",names(cubicsets6all))),sort(gsub("_ssp585","",names(cubicsets6all.ssp585))))

# # With CMIP6
# cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,cubicsets6,listinterim,listjra55,listncep)
# names(cubicsets)

names(cubicsets)
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsets5.rcp85,cubicsets5left.rcp85,cubicsets6,cubicsets6left, cubicsets6.ssp585,cubicsets6left.ssp585,listinterim,listjra55,listncep)
names(cubicsets)<- gsub(names(cubicsets), pattern = "_historical",replacement ="") 


namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics5.rcp85 <- names(cubicsets5.rcp85)
namescubics5left.rcp85 <- names(cubicsets5left.rcp85)

namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
namescubics6left <- names(cubicsets6left)
namescubics6all <- gsub(names(cubicsets6all), pattern = "_historical",replacement ="")
namescubics6.ssp585 <- names(cubicsets6.ssp585)
namescubics6left.ssp585 <- names(cubicsets6left.ssp585)
namescubics6all.ssp585 <- names(cubicsets6all.ssp585)


namescubicsearth <- names(cubicsetsearth)
namescubics <- names(cubicsets)
#shortnames <- #gsub(
#  gsub(names(cubicsets), pattern = "_r1i1p1", replacement = "")
#,
#pattern = "_r12i1p1", replacement = "")
#shortnames2 <- names(cubicsets)
#abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")




#############################################################################################
# Emergent constraints: (1) Hellinger distance vs global mean temperature increase
#############################################################################################
names(cubicsets)


futures.rcp85 <- grep("rcp85",names(cubicsets))
futures.ssp585 <- grep("ssp585",names(cubicsets))
futures <- c(futures.rcp85,futures.ssp585)
reans <- c(grep("interim",names(cubicsets)),grep("JRA55",names(cubicsets)),grep("ncep",names(cubicsets)))
hists.5 <- grep("CMIP5",names(cubicsets)[-futures.rcp85])
hists.6 <- grep("CMIP6",names(cubicsets)[-futures.ssp585])
hists <- c(hists.5,hists.6)

ids1.5 <- sapply(names(cubicsets)[hists.5], function(x)which(x == gsub("_rcp85","",names(cubicsets)[futures.rcp85])))
ids2.5 <- sapply(names(cubicsets)[futures.rcp85], function(x) which(gsub("_rcp85","",x) == names(cubicsets)[hists.5]))
ids1.6 <- sapply(names(cubicsets)[hists.6], function(x)which(x == gsub("_ssp585","",names(cubicsets)[futures.ssp585])))
ids2.6 <- sapply(names(cubicsets)[futures.ssp585], function(x) which(gsub("_ssp585","",x) == names(cubicsets)[hists.6]))



is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

sets.reans <- cubicsets[reans]
sets.hists <- cubicsets[hists][!sapply(c(ids1.5,ids1.6),is.integer0)]

names(sets.hists)
sets.futures <- cubicsets[futures][!sapply(c(ids2.5,ids2.6),is.integer0)]
names(sets.futures)
all.equal(sort(names(sets.hists)),gsub("_ssp585","",gsub("_rcp85","",sort(names(sets.futures)))))
sets.hists <- sets.hists[sort(names(sets.hists))]
sets.futures <- sets.futures[sort(names(sets.futures))]
names(sets.hists)
names(sets.futures)

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
hc_ncep <- lapply(c("../Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 

hc_gcms5.rcp85 <- lapply(paste0("FUTURE_",namescubics5.rcp85), loadIterations, permused = permk, it = it)
hc_gcms5left.rcp85<- lapply(paste0("FUTURE_CMIP5_left/FUTURE_",namescubics5left.rcp85), loadIterations, permused = permk, it = it)
hc_gcms6.ssp585 <- lapply(paste0("FUTURE_CMIP6/FUTURE_",namescubics6.ssp585), loadIterations, permused = permk, it = it)
hc_gcms6left.ssp585<- lapply(paste0("FUTURE_CMIP6_left/FUTURE_",namescubics6left.ssp585), loadIterations, permused = permk, it = it)



# With cmip6
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcms5.rcp85,hc_gcms5left.rcp85,hc_gcms6,hc_gcms6left,hc_gcms6.ssp585,hc_gcms6left.ssp585,hc_interim,hc_jra55,hc_ncep)
names(hc_gcms)<- c(namescubics5,namescubics5extra,namescubics5left,namescubics5.rcp85,namescubics5left.rcp85,namescubics6,namescubics6left,namescubics6.ssp585,namescubics6left.ssp585,names(listinterim),names(listjra55),names(listncep))
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
names(selection_hc_gcms)
names(data_gcms_anom_scaled)
########################################################################################
# Calculate hellinger distance bayesian networks
# KL distances give NaN (because there are x such that P(x) = 0 does not imply Q(x) = 0)
########################################################################################
names(selection_fits_gcms_scaled)
fits_NEL <- lapply(selection_fits_gcms_scaled,bnlearn:::as.graphNEL)
edgelistssparse <- lapply(fits_NEL, edgeList)
data_gcms_anom_scaled_sparse <- lapply(data_gcms_anom_scaled_3,sparsebnData, type = "continuous")
fits_COVS <- mapply(get.covariance, x = edgelistssparse, data = data_gcms_anom_scaled_sparse)
fits_COVmats <- lapply(fits_COVS, as.matrix)
fits_Means <- lapply(data_gcms_anom_scaled, function(x) colMeans(x[permutations[[3]]]) )

hellinger_coefficients <- matrix(data = NA, nrow = length(fits_COVmats), ncol = length(fits_COVmats), dimnames = list(names(fits_COVmats),names(fits_COVmats)))
i <- 1
for(i in 1:length(fits_COVmats)){
  kl <- mapply(function(b,d) normdiff(mu1 =fits_Means[[i]], mu2 = b, sigma1 = fits_COVmats[[i]], sigma2 = d, method = "Hellinger"),b = fits_Means, d = fits_COVmats)
  hellinger_coefficients[i,] <- kl
}

save(hellinger_coefficients, file ="results/hellinger_coefficient/hellinger_coefficients_CMIP5_CMIP6_hist_vs_fut.rda")

#load(file ="results/hellinger_coefficient/hellinger_coefficients_CMIP5_CMIP6_hist_vs_fut.rda")
#dimnames(hellinger_coefficients) <- list(names(fits_COVmats),names(fits_COVmats))

hellinger_coefficients
