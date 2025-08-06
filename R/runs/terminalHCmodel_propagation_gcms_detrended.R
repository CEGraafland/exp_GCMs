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
########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")
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

modelsize <- 18

selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
hc_gcms <- NULL
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)
data_gcms <- NULL
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
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended"),
       file = paste0("results/propagation/perm",whichperm,"/posV81pos_detrended/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended.rda"))
  
}

#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 8:10){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV81pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2.rda"))
  
}


#################################################################################
# Negative evidence. V81 (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended"),
       file = paste0("results/propagation/perm",whichperm,"/posV81neg_detrended/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended.rda"))
  
}
#################################################################################
# Negative evidence. V81 (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 8:10){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV81neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2.rda"))
  
}
