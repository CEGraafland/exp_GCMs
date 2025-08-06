###################################################################################
# CMIP5_subset_5d
###################################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_cmip5_exsubset_5d_akima_cubic.rda")
nlat <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[2]
nlon <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[3]
set.seed(5)
perm1 <- sample(1:(nlat*nlon))

cubicsets <- tas_historical_cmip5_exsubset_5d_akima_cubic
namescubics <- names(cubicsets)

# #dir.create(paste0("data/hciterations/CMIP5_5d/"))
# directorynames <- names(tas_historical_cmip5_exsubset_5d_akima_cubic)
# # 
# for(i in 1:length(directorynames)){
# dir.create(paste0("data/hciterations/CMIP5_5d/",directorynames[i]))
# }
#  k <- 1
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_5d/",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/CMIP5_5d/",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/CMIP5_5d/",directorynames[i],"/perm",k,"test"))
# }
####################################################################################
# Datasets cmip5 5d, permutations k, full data
####################################################################################
k <- 1
s <- 5
for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[perm1]
  
  start <- NULL
  steps <- 100
  last <- 5000
  m <- 0
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 5000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[perm1]
  
  load(paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_4900_5000i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_4900_5000i")))
  steps <- 100
  last <- 10000
  
  for (m in 80:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP5_5d / TRAIN / indTrain to be defined
# 1 permutation. 
# in two parts: up to 5000 
########################################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_cmip5_exsubset_5d_akima_cubic.rda")
ntime <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[1]
nlat <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[2]
nlon <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[3]

set.seed(5)
perm1 <- sample(1:(nlat*nlon))

cubicsets <- tas_historical_cmip5_exsubset_5d_akima_cubic
namescubics <- names(cubicsets)

set.seed(5)
indTRAIN1 <- sample(1:ntime, size =ntime/2)
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]




k <- 1
s <- 5
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[perm1]
  
  start <- NULL
  steps <- 100
  last <- 5000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 5000
for(s in 1:length(cubicsets)){
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[perm1]

  load(paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_5000_5100i.rda"))
  start <- eval(get(paste0("train1_",namescubics[set],"_",k,"_5000_5100i")))
  steps <- 100
  last <- 10000

  for (m in 51:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)

    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)

    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"),
         file = paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))


    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP5_5d / TEST /
# First permutation. 
# in two parts: up to 5000 
########################################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_cmip5_exsubset_5d_akima_cubic.rda")
ntime <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[1]
nlat <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[2]
nlon <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[3]

set.seed(5)
perm1 <- sample(1:(nlat*nlon))

cubicsets <- tas_historical_cmip5_exsubset_5d_akima_cubic
namescubics <- names(cubicsets)

set.seed(5)
indTRAIN1 <- sample(1:ntime, size =ntime/2)
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]




k <- 1
s <- 5
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[perm1]
  
  start <- NULL
  steps <- 100
  last <- 5000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# for(s in 1:length(cubicsets)){
#   
#   set <- s
#   # make data TEST and TRAIN
#   dataset <- cubicsets[[set]]
#   df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
#   testdata <- df[indTEST1,]
#   data <- testdata[permutations[[k]]]
#   
#   load(paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_2000_2100i.rda"))
#   start <- eval(get(paste0("test1_",namescubics[set],"_",k,"_2000_2100i")))
#   steps <- 100
#   last <- 10000
#   
#   for (m in 21:(last/steps)) {
#     i <- m*steps
#     j <- i+steps
#     berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
#     
#     assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
#     
#     save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
#          file = paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
#     
#     
#     if(m==0){
#       start <- berekening
#     } else if(narcs(berekening) == narcs(start)){
#       break
#     } else {start <- berekening}
#   }
# }



###################################################################################
# CMIP5_subset_20d
###################################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_cmip5_exsubset_20d_akima_cubic.rda")
nlat <- dim(tas_historical_cmip5_exsubset_20d_akima_cubic$tas_2T_Interim$Data)[2]
nlon <- dim(tas_historical_cmip5_exsubset_20d_akima_cubic$tas_2T_Interim$Data)[3]
set.seed(5)
perm1 <- sample(1:(nlat*nlon))

cubicsets <- tas_historical_cmip5_exsubset_20d_akima_cubic
namescubics <- names(cubicsets)

# dir.create(paste0("data/hciterations/CMIP5_20d/"))
# directorynames <- names(tas_historical_cmip5_exsubset_20d_akima_cubic)
# #
# for(i in 1:length(directorynames)){
# dir.create(paste0("data/hciterations/CMIP5_20d/",directorynames[i]))
# }
#  k <- 1
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_20d/",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/CMIP5_20d/",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/CMIP5_20d/",directorynames[i],"/perm",k,"test"))
# }
####################################################################################
# Datasets cmip5 20d, permutations k, full data
####################################################################################
k <- 1
s <- 1
for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[,perm1]
  
  start <- NULL
  steps <- 100
  last <- 5000

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_20d/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 5000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[perm1]
  
  load(paste0("data/hciterations/CMIP5_20d/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_4900_5000i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_4900_5000i")))
  steps <- 100
  last <- 10000
  
  for (m in 80:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_20d/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP5_20d / TRAIN / indTrain to be defined
# 1 permutation. 
# in two parts: up to 5000 
########################################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_cmip5_exsubset_20d_akima_cubic.rda")
ntime <- dim(tas_historical_cmip5_exsubset_20d_akima_cubic$tas_2T_Interim$Data)[1]
nlat <- dim(tas_historical_cmip5_exsubset_20d_akima_cubic$tas_2T_Interim$Data)[2]
nlon <- dim(tas_historical_cmip5_exsubset_20d_akima_cubic$tas_2T_Interim$Data)[3]

set.seed(5)
perm1 <- sample(1:(nlat*nlon))

cubicsets <- tas_historical_cmip5_exsubset_20d_akima_cubic
namescubics <- names(cubicsets)

set.seed(5)
indTRAIN1 <- sample(1:ntime, size =ntime/2)
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]




k <- 1
s <- 5
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[perm1]
  
  start <- NULL
  steps <- 100
  last <- 5000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_20d/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 5000
for(s in 1:length(cubicsets)){
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[perm1]
  
  load(paste0("data/hciterations/CMIP5_20d/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_5000_5100i.rda"))
  start <- eval(get(paste0("train1_",namescubics[set],"_",k,"_5000_5100i")))
  steps <- 100
  last <- 10000
  
  for (m in 51:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"),
         file = paste0("data/hciterations/CMIP5_20d/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# ########################################################################################
# # CMIP5_5d / TEST /
# # First permutation. 
# # in two parts: up to 5000 
# ########################################################################################
# setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
# rm(list = ls())
# library(bnlearn)
# library(transformeR)
# library(magrittr)
# library(igraph)
# source("../R/Functions/BasicNetworkFunctions.R")
# load("data/tas_historical_cmip5_exsubset_5d_akima_cubic.rda")
# ntime <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[1]
# nlat <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[2]
# nlon <- dim(tas_historical_cmip5_exsubset_5d_akima_cubic$tas_2T_Interim$Data)[3]
# 
# set.seed(5)
# perm1 <- sample(1:(nlat*nlon))
# 
# cubicsets <- tas_historical_cmip5_exsubset_5d_akima_cubic
# namescubics <- names(cubicsets)
# 
# set.seed(5)
# indTRAIN1 <- sample(1:ntime, size =ntime/2)
# # samplesize is amount of months (t in grid) nvars is 18 * 36
# samplesize <- dim(cubicsets[[1]]$Data)[1]
# nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# # indices for test is what is left
# indTEST1 <- (1:samplesize)[-indTRAIN1]
# 
# 
# 
# 
# k <- 1
# s <- 5
# for(s in 1:length(cubicsets)){
#   
#   set <- s
#   # make data TEST and TRAIN
#   dataset <- cubicsets[[set]]
#   df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
#   testdata <- df[indTEST1,]
#   data <- testdata[perm1]
#   
#   start <- NULL
#   steps <- 100
#   last <- 5000
#   
#   for (m in 0:(last/steps)) {
#     i <- m*steps
#     j <- i+steps
#     berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
#     
#     assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
#     
#     save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
#          file = paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
#     
#     
#     if(m==0){
#       start <- berekening
#     } else if(narcs(berekening) == narcs(start)){
#       break
#     } else {start <- berekening}
#   }
# }
# 
# # for(s in 1:length(cubicsets)){
# #   
# #   set <- s
# #   # make data TEST and TRAIN
# #   dataset <- cubicsets[[set]]
# #   df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
# #   testdata <- df[indTEST1,]
# #   data <- testdata[permutations[[k]]]
# #   
# #   load(paste0("data/hciterations/CMIP5_5d/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_2000_2100i.rda"))
# #   start <- eval(get(paste0("test1_",namescubics[set],"_",k,"_2000_2100i")))
# #   steps <- 100
# #   last <- 10000
# #   
# #   for (m in 21:(last/steps)) {
# #     i <- m*steps
# #     j <- i+steps
# #     berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
# #     
# #     assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
# #     
# #     save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
# #          file = paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
# #     
# #     
# #     if(m==0){
# #       start <- berekening
# #     } else if(narcs(berekening) == narcs(start)){
# #       break
# #     } else {start <- berekening}
# #   }
# # }