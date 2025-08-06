#############################################################################
# Propgation of evidence
#############################################################################
#setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
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
# load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"

cubicsets <- listinterim

# cubicsets <- c(tas_historical_10d_akima_cubic_corrected,listinterim)

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
hc_gcms<- lapply(namescubics, loadIterations, permused = 3)
names(hc_gcms) <- shortnames
data_gcms <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE)))

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

####################################################################################
# Propagation V477
####################################################################################
#################################################################################
# Single evidence. V477 postive (+ +)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V477_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(477),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V477_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV477pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V477_equal2.rda"))
  
}
#################################################################################
# Negative evidence. V477 (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V477_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(477),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V477_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV477neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V477_equal2.rda"))
  
}

####################################################################################
# Propagation V478
####################################################################################
#################################################################################
# Single evidence. V478 postive (+ +)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V478_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(478),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V478_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV478pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V478_equal2.rda"))
  
}
#################################################################################
# Negative evidence. V478 (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V478_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(478),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V478_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV478neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V478_equal2.rda"))
  
}

####################################################################################
# Propagation V95
####################################################################################
#################################################################################
# Single evidence. V95 postive (+ +)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V95_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(95),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V95_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV95pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V95_equal2.rda"))
  
}



#########################################################################
# Multiple evidence
#########################################################################
data <- TimeCoordsAnom_from_Grid_rms(cubicsets$interim_10d_akima_cubic, rms = TRUE)
VCs <- attr(data, "VertexCoords")
EurAsia <- which(VCs$x >= 0 & VCs$x <= 180 & VCs$y >= 25 & VCs$y <= 60)
#################################################################################
# Double evidence. V337 tot en met V645 postive (+ +)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3

  for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V337cubicV645_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = EurAsia,
                                     valueEvidence = rep(2,length(EurAsia)),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V337cubicV645_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/plusV337cubicV645pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V337cubicV645_equal2.rda"))
  
}
#################################################################################
# Double evidence. V337 tot en met V645 negative (+ +)
#################################################################################
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3

for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V337cubicV645_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = EurAsia,
                                     valueEvidence = rep(2,length(EurAsia)),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V337cubicV645_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/plusV337cubicV645neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V337cubicV645_equal2.rda"))
  
}


#################################################################################
# Double evidence. V81 V227 postive (+ +)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 6:10){
#for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81plusV227pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22.rda"))
  
}

#################################################################################
# Double evidence. V81 V227 negative (+ +)
#################################################################################
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 6:10){
# for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81plusV227neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22.rda"))
  
}
#################################################################################
# Double evidence. V81 V227 postive (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,-2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81minV227pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2.rda"))
  
}
#################################################################################
# Double evidence. V81 V227 postive (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,-2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81minV227neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2.rda"))
  
}


#################################################################################
# Double evidence. V81 V171 postive (+ n)
#################################################################################

x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81neuV171pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20.rda"))
  
}
#################################################################################
# Double evidence. V81 V171 negative (+ n)
#################################################################################

x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81neuV171neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20.rda"))
  
}

#################################################################################
# Single evidence. V171 postive (+ +)
#################################################################################

x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(171),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV171pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2.rda"))
  
}
#################################################################################
# Negative evidence. V171 (+ -)
#################################################################################

x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(171),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV171neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2.rda"))
  
}


#############################################################################
# Propgation of evidence CMIP6
#############################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list=ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
########################################################################################################
# Model Evaluation CMIP6
########################################################################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"
cubicsets6 <- tas_historical_cmip6_10d_akima_cubic_corrected
namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
# namescubics <- names(cubicsets)
# shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
#                    pattern = "_r12i1p1", replacement = "")
# shortnames
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
# hc_gcms <- lapply(namescubics6, loadIterations, permused = 3)
hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = 3)
names(hc_gcms6) <- namescubics6
data_gcms <- lapply(cubicsets6,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE)))

modelsize <- 18

selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms6, y = modelsize, SIMPLIFY = FALSE)
hc_gcms6 <- NULL
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)
data_gcms <- NULL
####################################################################################
# Propagation V81
####################################################################################
#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",namescubics6[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",namescubics6[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV81pos/CMIP6/prop_",namescubics6[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2.rda"))
  
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
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",namescubics6[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",namescubics6[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV81neg/CMIP6/propneg_",namescubics6[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2.rda"))
  
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

#################################################################################
# Double evidence. V81 V227 postive (+ +)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 6:10){
  #for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81plusV227pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22.rda"))
  
}

#################################################################################
# Double evidence. V81 V227 negative (+ +)
#################################################################################
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 6:10){
  # for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81plusV227neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal22.rda"))
  
}
#################################################################################
# Double evidence. V81 V227 postive (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,-2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81minV227pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2.rda"))
  
}
#################################################################################
# Double evidence. V81 V227 postive (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,-2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81minV227neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V227_equal2min2.rda"))
  
}


#################################################################################
# Double evidence. V81 V171 postive (+ n)
#################################################################################

x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81neuV171pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20.rda"))
  
}
#################################################################################
# Double evidence. V81 V171 negative (+ n)
#################################################################################

x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20"),
       file = paste0("results/propagation/perm",whichperm,"/plusV81neuV171neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81V171_equal20.rda"))
  
}

#################################################################################
# Single evidence. V171 postive (+ +)
#################################################################################

x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(171),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV171pos/prop_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2.rda"))
  
}
#################################################################################
# Negative evidence. V171 (+ -)
#################################################################################

x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms)){
  
  assign(paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(171),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2"),
       file = paste0("results/propagation/perm",whichperm,"/posV171neg/propneg_",shortnames[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V171_equal2.rda"))
  
}


###############
# V477 V95
###############

