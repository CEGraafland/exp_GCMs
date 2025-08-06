###################################################################################
# HC iterations for GCMs and Interim_cubic_akima
###################################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"
cubicsets <- c(tas_historical_10d_akima_cubic_corrected,listinterim)

namescubics <- names(cubicsets)
shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
                   pattern = "_r12i1p1", replacement = "")
shortnames
#####################################################################################
# create directories
#####################################################################################
k <- 2


 # for(i in 1:length(namescubics)){
 #  dir.create(paste0("data/hciterations/",namescubics[i],"/perm",k))
 #   dir.create(paste0("data/hciterations/",namescubics[i],"/perm",k,"train"))
 #   dir.create(paste0("data/hciterations/",namescubics[i],"/perm",k,"test"))
 # }
 # 
 # 
 # dir.create(paste0("data/interim_10d_akima_cubic"))
 # dir.create(paste0("data/interim_10d_akima_cubic/perm3"))
 # dir.create(paste0("data/interim_10d_akima_cubic/perm3train"))
 # dir.create(paste0("data/interim_10d_akima_cubic/perm3test"))



####################################################################################
# Different datasets, permutations k, full data
####################################################################################
# k <- 3
k <- 2
# for(s in 1:5){
# for(s in 6:10){
# for(s in c(5,4)){
# for(s in c(3)){ MISSCHIEN 8 NOG DOEN !!!
# for(s in c(10,9)){
for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]

  start <- NULL
  steps <- 100
  last <- 10000

  for (m in 0:(last/steps)) {
      i <- m*steps
      j <- i+steps
      berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
      assign(paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
      save(list = paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
      if(m==0){
        start <- berekening
      } else if(narcs(berekening) == narcs(start)){
        break
      } else {start <- berekening}
  }
}

# Until 2000
   for(s in 1:length(cubicsets)){
   set <- s
   
   dataset <- cubicsets[[set]]
   df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
   data <- df[permutations[[k]]]
   
   start <- NULL
   steps <- 100
   last <- 2000
   
   for (m in 0:(last/steps)) {
     i <- m*steps
     j <- i+steps
     berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
     assign(paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
     
     save(list = paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), 
          file = paste0("data/hciterations/",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
     
     if(m==0){
       start <- berekening
     } else if(narcs(berekening) == narcs(start)){
       break
     } else {start <- berekening}
   }
 }

# Up from 2100

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0(shortnames[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}


########################################################################################
# Hill Climbing iterations.
# kth permutation. 
# Train  
# Use indTrain1 which is indices of tas_ncep_10d.
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 2

# for(s in 1:5){
# for(s in 6:10){
# for(s in c(5)){
# for(s in c(10,9)){
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]

  start <- NULL
  steps <- 100
  last <- 10000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# until 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# up from 2100

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  load(paste0("data/hciterations/",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("train1_",shortnames[set],"_",k,"_2000_2100i")))
  
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# Hill Climbing iterations.
# kth permutation. 
# Test. 
# Use indTrain1 which is indices of tas_ncep_10d.
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 2

# for(s in 1:5){
# for(s in 6:10){
# for(s in c(5,4)){
for(s in c(10,9)){
#for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 10000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# until 2000
for(s in 1:length(cubicsets)){

set <- s
# make data TEST and TRAIN
dataset <- cubicsets[[set]]
df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
testdata <- df[indTEST1,]
data <- testdata[permutations[[k]]]

start <- NULL
steps <- 100
last <- 2000

for (m in 0:(last/steps)) {
  i <- m*steps
  j <- i+steps
  berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
  
  assign(paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
  
  save(list = paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
       file = paste0("data/hciterations/",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
  
  
  if(m==0){
    start <- berekening
  } else if(narcs(berekening) == narcs(start)){
    break
  } else {start <- berekening}
}
}

# Up from 2100
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  load(paste0("data/hciterations/",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("test1_",shortnames[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}


# Add interim? or maybe let solo
#####################################################################
# Make initial data two times or one time like before.
# initialdataRMS_int <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_interim_10dnew, rms = TRUE))
# initialdataRMS_int2 <- as.data.frame(TimeCoordsAnom_from_Grid_stand1_2times(tas_interim_10dnew))
# initialdataRMS_int3 <- as.data.frame(TimeCoordsAnom_from_Grid_stand1_2times(tas_interim_10dnew,twotimes = TRUE))
# 
# berekening <- hc(initialdataRMS_int, max.iter = 1800, score = "bic-g",start = NULL)
# berekening2 <- hc(initialdataRMS_int2, max.iter = 1800, score = "bic-g",start = NULL)
# berekening3 <- hc(initialdataRMS_int3, max.iter = 1800, score = "bic-g",start = NULL)
# par(mfrow = c(1,2))
# dev.off()
# plot.Meteodag(TimeCoordsAnom_from_Grid_rms(tas_interim_10dnew, rms = TRUE),berekening,lis = TRUE)
# plot.Meteodag(TimeCoordsAnom_from_Grid_stand1_2times(tas_interim_10dnew),berekening2,lis = TRUE)
# plot.Meteodag(TimeCoordsAnom_from_Grid_stand1_2times(tas_interim_10dnew,twotimes = TRUE),berekening3,lis = TRUE)
# 
# c1 <- cor(initialdataRMS_int)
# cv1 <- cov(initialdataRMS_int)
# c2 <- cor(initialdataRMS_int2)
# cv2 <- cov(initialdataRMS_int2)
# c3 <- cor(initialdataRMS_int3)
# cv3 <- cov(initialdataRMS_int3)


# all.equal(cv1,cv3)

###################################################################################
# HC iterations for JRA55
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_JRA55_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


listJRA55 <- list(tas_JRA55_10d_akima_cubic)
names(listJRA55) <- "JRA55_10d_akima_cubic"
cubicsets <- listJRA55
namescubics <- names(cubicsets)


#####################################################################################
# create directories JRA55
#####################################################################################
dir.create(paste0("data/hciterations/JRA55_10d_akima_cubic"))
dir.create(paste0("data/hciterations/JRA55_10d_akima_cubic/perm3"))
dir.create(paste0("data/hciterations/JRA55_10d_akima_cubic/perm3train"))
dir.create(paste0("data/hciterations/JRA55_10d_akima_cubic/perm3test"))

####################################################################################
# permutations k, full JRA55 data
####################################################################################
k <- 3
for(s in namescubics){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 10000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(set,"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(set,"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",set,"/perm",k,"/",set,"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}
########################################################################################
# JRA55 / TRAIN / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 10000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}
########################################################################################
# JRA55 / TEST / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 10000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

###################################################################################
# HC iterations for CMIP5 extra GCMs and CMIP5 earth and CMIP6 
###################################################################################
###################################################################################
# CMIP5 extra 
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


# listinterim <- list(tas_interim_10d_akima_cubic)
# names(listinterim) <- "interim_10d_akima_cubic"
cubicsets <- tas_historical_cmip5_extra_10d_akima_cubic
# shortnames <- names(tas_historical_cmip5_extra_10d_akima_cubic)
namescubics <- names(cubicsets)
# no removal of pattern r1i1p1 because we want to differ e.g. earth_ r12i1p1 (as used in first set)
# from earth_r1i1p1 and earth_r2i1p1 (as used in CMIP5_extra)
# shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
#                   pattern = "_r12i1p1", replacement = "")

#####################################################################################
# create directories CMIP5 extra / plus folder CMIP6 / plus folder CMIP5 earth. 
#####################################################################################
# # datasets_cmip6: created in load_and_prepare_data_cmip6_cmip5extra_cmip5earth.R
# load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
# directorynames <- gsub(names(tas_historical_cmip6_10d_akima_cubic_corrected), pattern = "_historical", replacement = "") 
# # gsub(gsub(datasets_cmip6, pattern = "_historical", replacement = ""),pattern = "-",replacement =".")
# k <- 2
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP6/",directorynames[i]))
# }
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP6/",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/CMIP6/",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/CMIP6/",directorynames[i],"/perm",k,"test"))
# }

# # datasets_cmip5_earth: created in load_and_prepare_data_cmip6_cmip5extra_cmip5earth.R
# directorynames <- gsub(gsub(names(tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic), pattern = "_historical", replacement = "_HIST"),pattern = "-",replacement =".")
# # directorynames <- gsub(gsub(datasets_cmip5_earth, pattern = "_historical", replacement = "_HIST"),pattern = "-",replacement =".")
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",directorynames[i]))
# }
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",directorynames[i],"/perm",k,"test"))
# }

# # datasets_cmip5_extra: created in load_and_prepare_data_cmip6_cmip5extra_cmip5earth.R
# directorynames <- names(tas_historical_cmip5_extra_10d_akima_cubic)
# # directorynames <- gsub(gsub(datasets_cmip5_extra, pattern = "_historical", replacement = ""),pattern = "-",replacement =".")
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_extra/",directorynames[i]))
# }
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_extra/",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/CMIP5_extra/",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/CMIP5_extra/",directorynames[i],"/perm",k,"test"))
# }

# The difference between all datasets_cmip5 and datasets_cmip5_extra are the ones who already were used 
# in the first analisis: We manually removed "CMIP5_MPI-ESM-MR_r1i1p1_historical" because we already had
# "CMIP5_MPI-ESM-LR_r1i1p1_historical":
# setdiff(datasets_cmip5,datasets_cmip5_extra)

####################################################################################
# Different datasets, permutations k, full data
# divide in two parts: until 1800 and after 1800
# Adapt name earth and CMIP6: Add or left out historical?
# replace "cubicsets"
# replace last
# other for loop with rest
####################################################################################
k <- 2

# for(s in 1:5){
# for(s in 6:10){
# for(s in c(10,9)){
  for(s in namescubics){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 10000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(set,"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(set,"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_extra/",set,"/perm",k,"/",set,"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Until 2100
# for(s in 1:5){
# for(s in 6:10){
# for(s in c(10,9)){

for(s in namescubics){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(set,"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(set,"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_extra/",set,"/perm",k,"/",set,"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2100 
for(s in namescubics){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_extra/",set,"/perm",k,"/",set,"_",k,"_2000_2100i.rda"))
  
  start <- eval(get(paste0(set,"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(set,"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(set,"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_extra/",set,"/perm",k,"/",set,"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}


########################################################################################
# CMIP5_extra / TRAIN / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_extra/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_extra/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("train1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_extra/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP5_extra / TEST / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_extra/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_extra/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("test1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_extra/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

###################################################################################
# CMIP5_left
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


cubicsets <- tas_historical_cmip5_left_10d_akima_cubic
namescubics <- names(cubicsets)

#  # directory cmip5_left: created in load_and_prepare_data_cmip6_cmip5extra_cmip5earth.R
#  directorynames <- names(tas_historical_cmip5_left_10d_akima_cubic)
# 
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_left/",directorynames[i]))
# }
# k <- 3
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_left/",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/CMIP5_left/",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/CMIP5_left/",directorynames[i],"/perm",k,"test"))
# }
####################################################################################
# Datasets cmip5 left, permutations k, full data
####################################################################################
k <- 3

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 2000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_1900_2000i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_1900_2000i")))
  steps <- 100
  last <- 10000
  
  for (m in 20:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}


########################################################################################
# CMIP5_left / TRAIN / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("train1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP5_left / TEST / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("test1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_left/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}




###################################################################################
# CMIP5 earth 
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

# listinterim <- list(tas_interim_10d_akima_cubic)
# names(listinterim) <- "interim_10d_akima_cubic"
cubicsets <- tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic
# shortnames <- names(tas_historical_cmip6_10d_akima_cubic)
namescubics <- names(cubicsets)
foldernamescubics <- gsub(names(cubicsets),pattern = "historical",replacement = "HIST")
# no removal of pattern r1i1p1 because we want to differ e.g. earth_ r12i1p1 (as used in first set)
# from earth_r1i1p1 and earth_r2i1p1 (as used in CMIP5_extra)

####################################################################################
# Datasets cmip5 earth, permutations k, full data
####################################################################################
k <- 2

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",foldernamescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 2000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",foldernamescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_2000_2100i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 20:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",foldernamescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP5 earth / TRAIN / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",foldernamescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_ONLY_HIST/",foldernamescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("train1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",foldernamescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP5 earth / TEST / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",foldernamescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",foldernamescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("test1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_EARTH_ONLYHIST/",foldernamescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

###################################################################################
# CMIP5_CSIRO
###################################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
load("../Data/interim/tas_interim_10dnew.rda")
load("data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")


cubicsets <- tas_historical_cmip5_CSIRO_10d_akima_cubic
namescubics <- gsub("-",".",gsub(".ncml","",gsub("https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip5/historical/historical_day_","",sapply(cubicsets,function(x)attr(x,"dataset")))))

# # directory cmip5_CSIRO:
#directorynames <- namescubics
# # dir.create(paste0("data/hciterations/CMIP5_CSIRO"))
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_CSIRO/",directorynames[i]))
# }
# k <- 3
#  for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP5_CSIRO/",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/CMIP5_CSIRO/",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/CMIP5_CSIRO/",directorynames[i],"/perm",k,"test"))
#  }
####################################################################################
# Datasets cmip5 CSIRO, permutations k, full data
####################################################################################
k <- 5
for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 1700
  last <- 1700
  m <- 0
  
   for (m in 0:0) {
  # for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_CSIRO/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 1700

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP5_CSIRO/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_0_1700i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_0_1700i")))
  steps <- 100
  last <- 1800
  m <- 17
  #for (m in 17:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP5_CSIRO/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  #}
}



############################################################################################
# FUTURE CMIP5 rcp 85
############################################################################################

###################################################################################
# HC iterations for GCMs and Interim_cubic_akima
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

cubicsets <- tas_rcp85_cmip5_10d_akima_cubic_corrected
namescubics <- names(cubicsets)
shortnames <- gsub(names(cubicsets), pattern = "_r1i1p1", replacement = "")
#####################################################################################
# create directories
#####################################################################################
# k <- 3
# for(i in 1:length(namescubics)){
#   dir.create(paste0("data/hciterations/FUTURE_",namescubics[i]))
# }
# 
# for(i in 1:length(namescubics)){
#  dir.create(paste0("data/hciterations/FUTURE_",namescubics[i],"/perm",k))
#   dir.create(paste0("data/hciterations/FUTURE_",namescubics[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/FUTURE_",namescubics[i],"/perm",k,"test"))
# }



####################################################################################
# Different datasets, permutations k, full data
####################################################################################
k <- 3

# Until 2000
for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2100

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0(shortnames[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}


########################################################################################
# Hill Climbing iterations.
# kth permutation. 
# Train  
# Use indTrain1 which is indices of tas_ncep_10d.
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3

# until 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# up from 2100

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  load(paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("train1_",shortnames[set],"_",k,"_2000_2100i")))
  
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# Hill Climbing iterations.
# kth permutation. 
# Test. 
# Use indTrain1 which is indices of tas_ncep_10d.
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3

# until 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2100
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  load(paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("test1_",shortnames[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}



###################################################################################
# FUTURE CMIP5_left rcp85
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("data/tas_rcp85_cmip5_left_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


cubicsets <- tas_rcp85_cmip5_left_10d_akima_cubic
namescubics <- names(cubicsets)

#  # directory cmip5_left FUTURE:
#  directorynames <- names(tas_rcp85_cmip5_left_10d_akima_cubic)
# 
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/FUTURE_CMIP5_left/FUTURE_",directorynames[i]))
# }
# 
# k <- 3
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/FUTURE_CMIP5_left/FUTURE_",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/FUTURE_CMIP5_left/FUTURE_",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/FUTURE_CMIP5_left/FUTURE_",directorynames[i],"/perm",k,"test"))
# }
####################################################################################
# Datasets FUTURE cmip5 left, permutations k, full data
####################################################################################
k <- 3

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP5_left/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 2000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/FUTURE_CMIP5_left/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_1900_2000i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_1900_2000i")))
  steps <- 100
  last <- 10000
  
  for (m in 20:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP5_left/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

###################################################################################
# CMIP6 
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


# listinterim <- list(tas_interim_10d_akima_cubic)
# names(listinterim) <- "interim_10d_akima_cubic"
cubicsets <- tas_historical_cmip6_10d_akima_cubic_corrected
# shortnames <- names(tas_historical_cmip6_10d_akima_cubic)
namescubics <- gsub(names(cubicsets),pattern = "_historical", replacement = "")
# no removal of pattern r1i1p1 because we want to differ e.g. earth_ r12i1p1 (as used in first set)
# from earth_r1i1p1 and earth_r2i1p1 (as used in CMIP5_extra)
shortnamesa <- gsub(gsub(gsub(gsub(names(cubicsets), 
                                   pattern = "Amon", replacement = ""),
                              pattern = "_historical", replacement = ""),
                         pattern = "_r1i1p1f2", replacement = ""),
                    pattern = "_r1i1p1f3", replacement = "")
shortnamesb <- gsub(shortnamesa[2:length(shortnamesa)],pattern = "_r1i1p1f1", replacement = "")
shortnames <- c(shortnamesa[1],shortnamesb)

####################################################################################
# Datasets cmip6, permutations k, full data
####################################################################################
k <- 2

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 2000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_1900_2000i.rda"))
  
  start <- eval(get(paste0(shortnames[set],"_",k,"_1900_2000i")))
  steps <- 100
  last <- 10000
  
  for (m in 20:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"/",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP6 / TRAIN / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("train1_",shortnames[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"train/train1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP6 / TEST / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("test1_",shortnames[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6/",namescubics[set],"/perm",k,"test/test1_",shortnames[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

###################################################################################
# CMIP6_left
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_historical_cmip6_left_10d_akima_cubic_corrected.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


cubicsets <- tas_historical_cmip6_left_10d_akima_cubic_corrected
namescubics <- names(cubicsets)

# # directory cmip6_left:
# directorynames <- names(tas_historical_cmip6_left_10d_akima_cubic_corrected)
# 
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP6_left/",directorynames[i]))
# }
# k <- 3
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/CMIP6_left/",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/CMIP6_left/",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/CMIP6_left/",directorynames[i],"/perm",k,"test"))
# }
####################################################################################
# Datasets cmip6 left, permutations k, full data
####################################################################################
k <- 3

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 2000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_1900_2000i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_1900_2000i")))
  steps <- 100
  last <- 10000
  
  for (m in 20:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}


########################################################################################
# CMIP6_left / TRAIN / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("train1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP6_left / TEST / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  load(paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("test1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/CMIP6_left/",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}






####################################################################################
# FUTURE CMIP6 ssp545
####################################################################################

setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/tas_ssp585_cmip6_10d_akima_cubic_corrected.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


# listinterim <- list(tas_interim_10d_akima_cubic)
# names(listinterim) <- "interim_10d_akima_cubic"
cubicsets <- tas_ssp585_cmip6_10d_akima_cubic_corrected
namescubics <- names(cubicsets)
# shortnames <- names(tas_historical_cmip6_10d_akima_cubic)
# namescubics <- gsub(names(cubicsets),pattern = "_historical", replacement = "")
# no removal of pattern r1i1p1 because we want to differ e.g. earth_ r12i1p1 (as used in first set)
# from earth_r1i1p1 and earth_r2i1p1 (as used in CMIP5_extra)
#shortnamesa <- gsub(gsub(gsub(gsub(names(cubicsets), 
 #                                  pattern = "Amon", replacement = ""),
  #                            pattern = "_ssp585", replacement = ""),
   #                      pattern = "_r1i1p1f2", replacement = ""),
    #                pattern = "_r1i1p1f3", replacement = "")
#shortnamesb <- gsub(shortnamesa[2:length(shortnamesa)],pattern = "_r1i1p1f1", replacement = "")
#shortnames <- c(shortnamesa[1],shortnamesb)

# # directory cmip6 FUTURE:
#   directorynames <- names(tas_ssp585_cmip6_10d_akima_cubic_corrected)
#  
#  for(i in 1:length(directorynames)){
#    dir.create(paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",directorynames[i]))
#  }
#  
#  k <- 3
#  for(i in 1:length(directorynames)){
#    dir.create(paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",directorynames[i],"/perm",k))
#    dir.create(paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",directorynames[i],"/perm",k,"train"))
#    dir.create(paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",directorynames[i],"/perm",k,"test"))
# }

####################################################################################
# Datasets cmip6 FUTURE, permutations k, full data
####################################################################################
k <- 3
for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 2000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_1900_2000i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_1900_2000i")))
  steps <- 100
  last <- 10000
  
  for (m in 20:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP6 FUTURE / TRAIN / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (NOT YET DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# Up from 2000
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  learndata <- df[indTRAIN1,]
  # testdata <- df[indTEST1,]
  data <- learndata[permutations[[k]]]
  
  load(paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("train1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("train1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"train/train1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

########################################################################################
# CMIP6 FUTURE / TEST / Use indTrain1 which is indices of tas_ncep_10d.
# Third permutation. 
# in two parts: up to 2000 (NOT YET DONE)
# from 2100 (NOT YET DONE)
########################################################################################
# load indTRAIN1 random 180 of 360
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]

k <- 3
for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

for(s in 1:length(cubicsets)){
  
  set <- s
  # make data TEST and TRAIN
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  testdata <- df[indTEST1,]
  data <- testdata[permutations[[k]]]
  
  load(paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_2000_2100i.rda"))
  start <- eval(get(paste0("test1_",namescubics[set],"_",k,"_2000_2100i")))
  steps <- 100
  last <- 10000
  
  for (m in 21:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    
    assign(paste0("test1_",shortnames[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("test1_",namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP6/FUTURE_",namescubics[set],"/perm",k,"test/test1_",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}



###################################################################################
# FUTURE CMIP6_left ssp585
###################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("data/tas_ssp585_cmip6_left_10d_akima_cubic_corrected.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


cubicsets <- tas_ssp585_cmip6_left_10d_akima_cubic_corrected
namescubics <- names(cubicsets)

#  # directory cmip6_left FUTURE:
#  directorynames <- names(tas_ssp585_cmip6_left_10d_akima_cubic_corrected)
# 
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/FUTURE_CMIP6_left/FUTURE_",directorynames[i]))
# }
# 
# k <- 3
# for(i in 1:length(directorynames)){
#   dir.create(paste0("data/hciterations/FUTURE_CMIP6_left/FUTURE_",directorynames[i],"/perm",k))
#   dir.create(paste0("data/hciterations/FUTURE_CMIP6_left/FUTURE_",directorynames[i],"/perm",k,"train"))
#   dir.create(paste0("data/hciterations/FUTURE_CMIP6_left/FUTURE_",directorynames[i],"/perm",k,"test"))
# }
####################################################################################
# Datasets FUTURE ssp585 cmip6 left, permutations k, full data
####################################################################################
k <- 3

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 100
  last <- 2000
  
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP6_left/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}

# HIER VERDER GAAN VANAF 2000

for(s in 1:length(cubicsets)){
  set <- s
  
  dataset <- cubicsets[[set]]
  df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(dataset,rms = TRUE))
  data <- df[permutations[[k]]]
  
  load(paste0("data/hciterations/FUTURE_CMIP6_left/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_1900_2000i.rda"))
  
  start <- eval(get(paste0(namescubics[set],"_",k,"_1900_2000i")))
  steps <- 100
  last <- 10000
  
  for (m in 20:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0(namescubics[set],"_",k,"_",i,"_",j,"i"), 
         file = paste0("data/hciterations/FUTURE_CMIP6_left/FUTURE_",namescubics[set],"/perm",k,"/",namescubics[set],"_",k,"_",i,"_",j,"i.rda"))
    
    if(m==0){
      start <- berekening
    } else if(narcs(berekening) == narcs(start)){
      break
    } else {start <- berekening}
  }
}





