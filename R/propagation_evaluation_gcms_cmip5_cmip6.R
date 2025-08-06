#############################################################################
# Model Evaluation CMIP5 and CMIP6 hist vs fut
#############################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
#############################################################################
rm(list = ls())
library(bnlearn)
 library(magrittr)
 library(reshape2)
 # install.packages("ggplot2")
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
 source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
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
 
 
 cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsets5.rcp85,cubicsets5left.rcp85,cubicsets6,cubicsets6left, cubicsets6.ssp585,cubicsets6left.ssp585,listinterim,listjra55,listncep)
 names(cubicsets)<- gsub(names(cubicsets), pattern = "_historical",replacement ="") 
 
 
 namescubics5 <- names(cubicsets5)
 shortnames5 <- gsub(gsub(names(cubicsets5), pattern = "_r1i1p1", replacement = ""),pattern = "_r12i1p1", replacement = "")
 namescubics5extra <- names(cubicsets5extra)
 namescubics5left <- names(cubicsets5left)
 namescubics5.rcp85 <- names(cubicsets5.rcp85)
 namescubics5left.rcp85 <- gsub(names(cubicsets5left.rcp85), pattern = "_historical",replacement ="") 
 
 #namescubics6 <- names(cubicsets6)
 namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
namescubics6left <- names(cubicsets6left)
namescubics6left <- gsub(names(cubicsets6left), pattern = "_historical",replacement ="") 
 namescubics6all <- gsub(names(cubicsets6all), pattern = "_historical",replacement ="")
 #namescubics6all <- names(cubicsets6all)
 namescubics6.ssp585 <- names(cubicsets6.ssp585)
 namescubics6left.ssp585 <- names(cubicsets6left.ssp585)
 namescubics6all.ssp585 <- names(cubicsets6all.ssp585)
 
 
 namescubicsearth <- names(cubicsetsearth)
 namescubics <- names(cubicsets)
 
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
 
 #  hellinger_coefficients[names(sets.hists),]
 #  names(sets.hists) %in% rownames(hellinger_coefficients)
 #   names(sets.hists)[!names(sets.hists) %in% rownames(hellinger_coefficients)]
 #   rownames(hellinger_coefficients)[!rownames(hellinger_coefficients)%in% c(names(sets.hists),names(sets.futures),names(sets.reans))]
 # ############################################################################
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
 
 #########################################################################
 # Load BNs
 #########################################################################
 it <- "1700_1800"
 hc_gcms5 <- lapply(paste0(namescubics5), loadIterations, permused = permk, it = it)
 hc_gcms5extra <- lapply(paste0("CMIP5_extra/",namescubics5extra), loadIterations, permused = permk, it = it)
 hc_gcms5left<- lapply(paste0("CMIP5_left/",namescubics5left), loadIterations, permused = permk, it = it)
 hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk, it = it)
 hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk, it = it)
 hc_gcms6left <- lapply(paste0("CMIP6_left/",names(cubicsets6left)), loadIterations, permused = permk, it = it)
 hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
 hc_jra55 <- lapply(c("JRA55_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
 hc_ncep <- lapply(c("../Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 

 hc_gcms5.rcp85 <- lapply(paste0("FUTURE_",namescubics5.rcp85), loadIterations, permused = permk, it = it)
 hc_gcms5left.rcp85<- lapply(paste0("FUTURE_CMIP5_left/FUTURE_",names(cubicsets5left.rcp85)), loadIterations, permused = permk, it = it)
 hc_gcms6.ssp585 <- lapply(paste0("FUTURE_CMIP6/FUTURE_",namescubics6.ssp585), loadIterations, permused = permk, it = it)
 hc_gcms6left.ssp585<- lapply(paste0("FUTURE_CMIP6_left/FUTURE_",namescubics6left.ssp585), loadIterations, permused = permk, it = it)
 
 
 
 # Join BNs With cmip6
 hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcms5.rcp85,hc_gcms5left.rcp85,hc_gcms6,hc_gcms6left,hc_gcms6.ssp585,hc_gcms6left.ssp585,hc_interim,hc_jra55,hc_ncep)
 names(hc_gcms)<- c(namescubics5,namescubics5extra,namescubics5left,namescubics5.rcp85,namescubics5left.rcp85,namescubics6,namescubics6left,namescubics6.ssp585,namescubics6left.ssp585,names(listinterim),names(listjra55),names(listncep))
 
 
 
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
 ##########################################################################
 # Global dataset evaluation from BN 
 ##########################################################################
 whichBNFits <- selection_fits_gcms_scaled
whichData <- data_gcms_anom_scaled
ordersets <- TRUE

loglik_selection_datasets <- matrix(data = NA, nrow = length(whichBNFits), ncol = length(whichData), dimnames = list(names(whichBNFits),names(whichData)))
for(i in 1:length(whichBNFits)){
  lo <- sapply(X = whichData, logLik, object = whichBNFits[[i]])
  loglik_selection_datasets[i,] <- lo
}

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

filesncep <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",names(listncep),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
filesjra55 <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",names(listjra55),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5extra <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",namescubics5extra,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5left <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",namescubics5left,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5 <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",namescubics5,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP5_rcp85/prop_",namescubics5.rcp85,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5left.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP5_rcp85/prop_",namescubics5left.rcp85,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files6.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6_ssp585/prop_",namescubics6.ssp585,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files6left.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6_ssp585/prop_",names(cubicsets6left.ssp585),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")


files6 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6/prop_",namescubics6,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files6left <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6/prop_",names(cubicsets6left),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
filesint <- paste0("results/propagation/perm3/pos",whichnode,"pos/prop_",names(listinterim),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(c(files5,files5extra,files5left,files5.rcp85,files5left.rcp85,files6,files6left,files6.ssp585,files6left.ssp585,filesint,filesncep,filesjra55),function(x)load(file = x))
for (i in 1:length(c(files5,files5extra,files5left,files5.rcp85,files5left.rcp85,files6,files6left,files6.ssp585,files6left.ssp585,filesint,filesncep,filesjra55))){load(file = c(files5,files5extra,files5left,files5.rcp85,files5left.rcp85,files6,files6left,files6.ssp585,files6left.ssp585,filesint,filesncep,filesjra55)[i])}


all.equal(c(namescubics6,namescubics6left),namescubics6all)

posposmix <- list()
combinames <- c(shortnames5,namescubics5extra,namescubics5left,namescubics5.rcp85,namescubics5left.rcp85,namescubics6,names(cubicsets6left),namescubics6.ssp585,namescubics6left.ssp585,names(listinterim),names(listncep),names(listjra55))
posposmix <- lapply(combinames, function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))

prop
all.equal(posposmix[[1]],posposmix[[2]])
for(i in 1:length(combinames)){
  posposmix[[i]] <- eval(parse(text = paste0("prop_",combinames[i],"_hc_",whichiterations,"i_",whichnode,"_equal2")))
}
names(posposmix) <- propfiles
names(cubicsets)
names(posposmix)
length(posposmix)

namescubics5.rcp85
namescubics5left.rcp85
#ready.prop.ind <- c(namescubics5,namescubics5extra,namescubics5left,namescubics5.rcp85,namescubics6,namescubics6left,namescubics6all.ssp585,names(listinterim),names(listncep),names(listjra55))
# gsub("_hc_1700_1800i_V81_equal2", "", gsub("prop_","",names(posposmix)))%in%names(cubicsets)
length(cubicsets)
length(posposmix)

#names(cubicsets[ready.prop.ind])
prop_dif_pos <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, griddata = cubicsets, SIMPLIFY = FALSE)
enspropdifsposclims <- bindGrid(prop_dif_pos, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposclims$Members <- combinames
enspropdifsposclims$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsposclims, names.attr = combinames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


files5 <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",shortnames5,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
names(files5)
files5extra <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",namescubics5extra,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5left <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",namescubics5left,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP5_rcp85/propneg_",namescubics5.rcp85,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files5left.rcp85<- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP5_rcp85/propneg_",namescubics5left.rcp85,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files6.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP6_ssp585/propneg_",namescubics6.ssp585,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files6left.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP6_ssp585/propneg_",names(cubicsets6left.ssp585),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files6 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP6/propneg_",namescubics6,"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
files6left <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP6/propneg_",names(cubicsets6left),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
filesint <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listinterim),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
filesncep <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listncep),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
filesjra55 <- paste0("results/propagation/perm3/pos",whichnode,"neg/propneg_",names(listjra55),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(c(files5,files5extra,files5left,files5.rcp85,files5left.rcp85,files6,files6left,files6.ssp585,files6left.ssp585,filesint,filesncep,filesjra55),function(x)load(file = x))
propfiles[[1]]
for (i in 1:length(c(files5,files5extra,files5left,files5.rcp85,files5left.rcp85,files6,files6left,files6.ssp585,files6left.ssp585,filesint,filesncep,filesjra55))){load(file = c(files5,files5extra,files5left,files5.rcp85,files5left.rcp85,files6,files6left,files6.ssp585,files6left.ssp585,filesint,filesncep,filesjra55)[i])}

posnegmix <- lapply(combinames, function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posnegmix
names(posnegmix) <- propfiles
all.equal(posnegmix$propneg_CMIP5_CanESM2_hc_1700_1800i_V81_equal2,
posnegmix$propneg_CMIP5_CNRM.CM5_hc_1700_1800i_V81_equal2)

prop_dif_neg <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata = cubicsets, SIMPLIFY = FALSE)

enspropdifsnegclims <- bindGrid(prop_dif_neg, dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsnegclims$Members <- combinames
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsnegclims,names.attr = combinames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_neg[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


names(posposmix)
names(posnegmix)
names(cubicsets)
prop_dif_pos_and_neg <- mapply(function(griddata,proppos,propneg) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")), griddata),griddata = cubicsets,proppos = posposmix,propneg = posnegmix, SIMPLIFY = FALSE)
prop_dif_pos_and_neg <- prop_dif_pos_and_neg[sort(names(prop_dif_pos_and_neg))]


ind6 <- grep("CMIP6",names(prop_dif_pos_and_neg))
indearth <- grep("EARTH_r1i|EARTH_r2i",names(prop_dif_pos_and_neg))

enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg[-c(ind6)], dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposandnegclims$Members <- sort(combinames)[-c(ind6)]
enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg[-c(ind6,indearth)], dimension = c("member"),skip.temporal.check = TRUE)
enspropdifsposandnegclims$Members <- sort(combinames)[-c(ind6,indearth)]
col.b
sort(combinames[-c(ind6,indearth)])
col <- c(rev(col.b),col.r)
spatialPlot(enspropdifsposandnegclims,names.attr = sort(combinames)[-c(ind6,indearth)],as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


############################################################################
#
############################################################################
##############################################################################################################
# 
##############################################################################################################
IPweight.split <- function(bdists,ref.perf,sD,sS,scaleDS = TRUE,meth = c("exponential","hyperbolic","expstep","extreme")){
  refID <- which(rownames(bdists) == ref.perf)
  D <- bdists[ref.perf,-refID]
  S <- bdists[-refID,-refID]
  if(isTRUE(scaleDS)){D <- D/(max(D));S <- S/max(S)}

  W <- numeric(length=length(D))
  W.ind <- numeric(length=length(D))
  W.perf <- numeric(length=length(D))
  names(W.ind)<-names(W.perf)<- names(D)
  
  for(i in 1:length(names(D))){
    name.i <-names(D)[i]
    nameID.S <- which(rownames(S) == name.i)
    nameID.D <- which(names(D) == name.i)
    Di <- D[nameID.D]
    Sij <- S[nameID.S,][-nameID.S]
    
    if(meth == "exponential"){
      W.perf[i] <- exp(-Di^2/sD^2)
      sumS <- sum(exp(-Sij^2/sS^2))
      W.ind[i] <- 1/(1+sumS)
    } else if (meth == "hyperbolic"){
      W.perf[i] <- 1/Di
      sumS <- sum(1/Sij)
      W.ind[i] <- 1/(1+sumS)
    } else if (meth== "expstep"){  
      Step.Sij <-Sij
      Step.Sij[Sij>sS]<- 0
      Step.Sij[Sij<=sS]<- 1
      
      W.perf[i] <- exp(-Di^2/sD^2)
      sumS <- sum(Step.Sij)
      W.ind[i] <- 1/(1+sumS)
  } else if (meth == "extreme"){
    W.perf[i] <- 1-(Di-min(D))/(max(D)-min(D))
    sumS <- sum(1/Sij)
    W.ind[i] <- 1/(1+sumS)
  }
  }
  W <- W.perf*W.ind
  W <- W/sum(W)
  W.perf<- W.perf/sum(W.perf)
  W.ind <- W.ind/sum(W.ind)
  
  return(data.frame(W,W.perf,W.ind))
}

load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_CMIP6_hist_vs_fut.rda")

names(sets.reans)
names(sets.hists)
names(sets.futures)
#########################################################################################
# CMIP 6
#########################################################################################
respective.grid <-"CMIP6Amon_CanESM5_r1i1p1f1" 
respective.grid <-"interim_10d_akima_cubic"  
respective.grid <-"ncep_10d" 
respective.grid <- "CMIP6Amon_EC.Earth3.Veg_r1i1p1f1" 
respective.grid <- "CMIP6Amon_CNRM.CM6.1_r1i1p1f2"
respective.grid <- "CMIP6Amon_MIROC6_r1i1p1f1" 
respective.grid <- "CMIP6Amon_IPSL.CM6A.LR_r1i1p1f1"

sets.hists.CMIP6 <- sets.hists[grep("CMIP6",names(sets.hists))]
sets.futures.CMIP6 <- sets.futures[grep("CMIP6",names(sets.futures))]
prop.hist.6.ind <- namescubics6%in% names(sets.hists.CMIP6)
sets.hists.CMIP6<- sets.hists.CMIP6[namescubics6[prop.hist.6.ind]]
names(sets.hists.CMIP6)

prop_dif_pos_and_neg.CMIP6 <- prop_dif_pos_and_neg[grep("CMIP6",names(prop_dif_pos_and_neg))]
if(respective.grid%in%names(sets.hists.CMIP6)){
  resp.ind<-  which(names(sets.hists.CMIP6)==respective.grid)
  prop.hist.6.sets.ind <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.CMIP6)))%in%names(sets.hists.CMIP6)[-resp.ind]
} else {prop.hist.6.sets.ind <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.CMIP6)))%in%names(sets.hists.CMIP6)}
prop_dif_pos_and_neg.CMIP6 <- prop_dif_pos_and_neg.CMIP6[prop.hist.6.sets.ind]
names(prop_dif_pos_and_neg.CMIP6)

prop_dif_pos_and_neg.fut.CMIP6 <- prop_dif_pos_and_neg[grep("CMIP6",names(prop_dif_pos_and_neg))]
if(respective.grid%in%names(sets.futures.CMIP6)){
  resp.ind<-  which(names(sets.futures.CMIP6)==respective.grid)
  prop.futures.6.sets.ind <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.fut.CMIP6)))%in%names(sets.futures.CMIP6)[-resp.ind]
} else {prop.futures.6.sets.ind <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.fut.CMIP6)))%in%names(sets.futures.CMIP6)}
prop_dif_pos_and_neg.fut.CMIP6 <- prop_dif_pos_and_neg.fut.CMIP6[prop.futures.6.sets.ind]
names(prop_dif_pos_and_neg.fut.CMIP6)

helcoef <- hellinger_coefficients
rownames(helcoef) <- gsub("_historical","", rownames(helcoef))
colnames(helcoef) <- gsub("_historical","", colnames(helcoef))
names(sets.hists.CMIP6) %in% rownames(helcoef)
if(respective.grid%in%names(sets.hists.CMIP6)){
resp.ind<-  which(names(sets.hists.CMIP6)==respective.grid)
b <- -log(helcoef)[c(names(sets.hists.CMIP6)[-resp.ind],respective.grid),c(names(sets.hists.CMIP6)[-resp.ind],respective.grid)]
} else {b <- -log(helcoef)[c(names(sets.hists.CMIP6),respective.grid),c(names(sets.hists.CMIP6),respective.grid)]}

# -log(helcoef)[c(names(sets.hists.CMIP6),"interim_10d_akima_cubic"),"interim_10d_akima_cubic"][order(-log(helcoef)[c(names(sets.hists.CMIP6),"interim_10d_akima_cubic"),"interim_10d_akima_cubic"])]
# modelweights[order(modelweights)]
# # CMIP 5 and 6
# rownames(helcoef) <- gsub("_historical","", rownames(helcoef))
# colnames(helcoef) <- gsub("_historical","", colnames(helcoef))
# names(clim.hists) %in% rownames(helcoef)
# b <- -log(helcoef)[c(names(clim.hists),"interim_10d_akima_cubic"),c(names(clim.hists),"interim_10d_akima_cubic")]
# 

modelweights.split <- IPweight.split(b,respective.grid,scaleDS = TRUE,sD = 0.5,sS = 0.5, meth = "exponential")
modelweights.split
modelweights <- modelweights.split[,"W.ind"]
names(modelweights)<-rownames(modelweights.split)
-log(helcoef)[c(names(sets.hists.CMIP6),respective.grid),respective.grid][order(-log(helcoef)[c(names(sets.hists.CMIP6),respective.grid),respective.grid])]

hist(modelweights)
sum(modelweights)
###########################################################################
# CMIP6 hist propagation
###########################################################################
prop.hist.CMIP6.weighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg.CMIP6, y = modelweights,SIMPLIFY = FALSE)
prop.hist.CMIP6.weighted <- Reduce('+', prop.hist.CMIP6.weighted.grids)

lapply(prop_dif_pos_and_neg.CMIP6,function(x) x$Data)

prop.hist.CMIP6.grids.Data <- lapply(prop_dif_pos_and_neg.CMIP6,function(x) subsetGrid(x,drop = TRUE)$Data)
dim(prop.hist.CMIP6.grids.Data$prop_CMIP6Amon_CanESM5_r1i1p1f1_hc_1700_1800i_V81_equal2)

newgrid <- prop_dif_pos_and_neg.CMIP6[[1]]
newgrid$Data <- prop.hist.CMIP6.weighted
spatialPlot(newgrid,names.attr = names(modelweights),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg.CMIP6[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop.hist.CMIP6.unweighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg.CMIP6, y = 1/length(prop_dif_pos_and_neg.CMIP6),SIMPLIFY = FALSE)
prop.hist.CMIP6.unweighted <- Reduce('+', prop.hist.CMIP6.unweighted.grids)

oldgrid <- prop_dif_pos_and_neg.CMIP6[[1]]
oldgrid$Data <- prop.hist.CMIP6.unweighted
spatialPlot(oldgrid,names.attr =  names(modelweights),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg.CMIP6[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

compareoldnew <- bindGrid(oldgrid,newgrid, dimension = c("member"))
compareoldnew$Members <- c("unweighted","weighted")
spatialPlot(compareoldnew,names.attr = c("unweighted","weighted"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

all.equal(oldgrid,newgrid)

modelweights[order(modelweights)]
binds <- prop_dif_pos_and_neg.CMIP6[order(modelweights,decreasing = TRUE)]
length(binds)
compare.all <- bindGrid(binds, dimension = "member")

compare.all.int <- bindGrid(prop_dif_pos_and_neg[[respective.grid]],compare.all,dimension = "member")

spatialPlot(compare.all.int,names.attr = c(respective.grid,names(modelweights)[order(modelweights,decreasing = TRUE)]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

##################################################################################################################
# Compare sd weighted versus unweighted
##################################################################################################################
ensemble.unw.prop <- t(sapply(prop_dif_pos_and_neg.CMIP6,function(x) array3Dto2Dmat(x$Data)))
sd.box.unw <-apply(ensemble.unw.prop, MARGIN = 2, FUN = sd)
w <- modelweights
# sd.cors.fut.weighted <- sqrt(sum(modelweights*(keltoC(RMSE.cors.gcms.futures.CMIP5.int)-mean.cors.fut.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
# sd.cors.fut.weighted <- sqrt(sum(w*((x)-x)^2)/((length(w)-1)/(length(w)*sum(w))))

sd.box.w <- apply(ensemble.unw.prop, MARGIN = 2, FUN = function(x) sqrt(sum(w*(x-mean(x))^2)/((length(w)-1)/(length(w)*sum(w)))))# HIER weights integreren. 
clim.sd.box.w <- quantity2clim(sd.box.w, "sd per gridbox",ref.grid = cubicsets$interim_10d_akima_cubic)
clim.sd.box.unw <- quantity2clim(sd.box.unw, "sd per gridbox",ref.grid = cubicsets$interim_10d_akima_cubic)

comparesdoldnew <- bindGrid(clim.sd.box.unw,clim.sd.box.w, dimension = c("member"))
comparesdoldnew$Members <- c("unweighted","weighted")

spatialPlot(comparesdoldnew,backdrop.theme = "coastline", 
            lonCenter = 180,col = "Reds",rev.colors = FALSE,
            names.attr = c(paste0("unweighted sd = ",round(sum(sd.box.unw))),paste0("weighted sd = ",round(sum(sd.box.w))))
            # at = seq(-0.85,0.85,0.01),region = TRUE, 
            # col.regions= col,
            #set.min = -0.85,set.max = 0.85, 
            #colorkey = list(width = 0.6, lables = list(cex = 0.5))
)

spatialPlot(clim.sd.box.w,backdrop.theme = "coastline", lonCenter = 180,rev.colors = TRUE,
            # at = seq(-0.85,0.85,0.01),region = TRUE, 
            # col.regions= col,
            #set.min = -0.85,set.max = 0.85, 
            #colorkey = list(width = 0.6, lables = list(cex = 0.5))
)
spatialPlot(clim.sd.box.unw,backdrop.theme = "coastline", lonCenter = 180 
            # at = seq(-0.85,0.85,0.01),region = TRUE, 
            # col.regions= col,
            #set.min = -0.85,set.max = 0.85, 
            #colorkey = list(width = 0.6, lables = list(cex = 0.5))
)





# spatialPlot(prop_dif_pos_and_neg$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2,names.attr = c("interim",names(modelweights)[order(modelweights,decreasing = TRUE)]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
#             col.regions= col,
#             set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# 
# bindGrid(prop.hist.CMIP6.weighted.grids,dimension = "member")
# spatialPlot(compare.all.int,names.attr = c("interim",names(modelweights)[order(modelweights,decreasing = TRUE)]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
#             col.regions= col,
#             set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

##########################################################################
# CMIP6 fut con hist modelweights (To do: with future pesos. )
##########################################################################
prop.fut.CMIP6.weighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg.fut.CMIP6, y = modelweights,SIMPLIFY = FALSE)
prop.fut.CMIP6.weighted <- Reduce('+', prop.fut.CMIP6.weighted.grids)

newgrid <- prop_dif_pos_and_neg.fut.CMIP6[[1]]
newgrid$Data <- prop.fut.CMIP6.weighted
spatialPlot(newgrid,names.attr = names(modelweights),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg.fut.CMIP6[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop.fut.CMIP6.unweighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg.fut.CMIP6, y = 1/length(prop_dif_pos_and_neg.fut.CMIP6),SIMPLIFY = FALSE)
prop.fut.CMIP6.unweighted <- Reduce('+', prop.fut.CMIP6.unweighted.grids)

oldgrid <- prop_dif_pos_and_neg.fut.CMIP6[[1]]
oldgrid$Data <- prop.fut.CMIP6.unweighted
spatialPlot(oldgrid,names.attr =  names(modelweights),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg.fut.CMIP6[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

compareoldnew <- bindGrid(oldgrid,newgrid, dimension = c("member"))
compareoldnew$Members <- c("unweighted","weighted")
spatialPlot(compareoldnew,names.attr = c("unweighted","weighted"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

modelweights[order(modelweights)]
binds <- prop_dif_pos_and_neg.fut.CMIP6[order(modelweights,decreasing = TRUE)]
length(binds)
compare.all <- bindGrid(binds, dimension = "member")


compare.all.int <- bindGrid(prop_dif_pos_and_neg[[paste0("prop_",respective.grid,"_hc_1700_1800i_V81_equal2")]],compare.all,dimension = "member",skip.temporal.check = TRUE)

spatialPlot(compare.all.int,names.attr = c(respective.grid,names(modelweights)[order(modelweights,decreasing = TRUE)]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

###########################################################################
# CMIP 5
###########################################################################
sets.hists.CMIP5 <-sets.hists[grep("CMIP5",names(sets.hists))]
sets.futures.CMIP5 <-sets.futures[grep("CMIP5",names(sets.futures))]
helcoef <- hellinger_coefficients
b <- -log(helcoef)[c(names(sets.hists.CMIP5),respective.grid),c(names(sets.hists.CMIP5),respective.grid)]

respective.grid <- "JRA55_10d_akima_cubic"
respective.grid <-  "interim_10d_akima_cubic"  

prop_dif_pos_and_neg.CMIP5 <- prop_dif_pos_and_neg[grep("CMIP5",names(prop_dif_pos_and_neg))]
if(respective.grid%in%names(sets.hists.CMIP5)){
  resp.ind <- which(names(sets.hists.CMIP5)==respective.grid)
  prop_dif_pos_and_neg.CMIP5 <- prop_dif_pos_and_neg.CMIP5[-resp.ind]
  # a <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.CMIP5)))
  prop.fut.5.ind <- grep("rcp85",names(prop_dif_pos_and_neg.CMIP5))
} else {prop.fut.5.ind <- grep("rcp85",names(prop_dif_pos_and_neg.CMIP5))}
prop_dif_pos_and_neg.CMIP5 <- prop_dif_pos_and_neg.CMIP5[-prop.fut.5.ind]
names(prop_dif_pos_and_neg.CMIP5)
b <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.CMIP5)))%in%shortnames5
names(prop_dif_pos_and_neg.CMIP5)[b] <- paste0("prop_",namescubics5,"_hc_1700_1800i_V81_equal2")
c <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.CMIP5)))%in%names(sets.hists.CMIP5)
prop_dif_pos_and_neg.CMIP5<- prop_dif_pos_and_neg.CMIP5[c]
names(prop_dif_pos_and_neg.CMIP5)

helcoef <- hellinger_coefficients
rownames(helcoef) <- gsub("_historical","", rownames(helcoef))
colnames(helcoef) <- gsub("_historical","", colnames(helcoef))
names(sets.hists.CMIP5) %in% rownames(helcoef)
if(respective.grid%in%names(sets.hists.CMIP5)){
  resp.ind<-  which(names(sets.hists.CMIP5)==respective.grid)
  b <- -log(helcoef)[c(names(sets.hists.CMIP5)[-resp.ind],respective.grid),c(names(sets.hists.CMIP5)[-resp.ind],respective.grid)]
} else {b <- -log(helcoef)[c(names(sets.hists.CMIP5),respective.grid),c(names(sets.hists.CMIP5),respective.grid)]}


modelweights.split <- IPweight.split(b,respective.grid,scaleDS = TRUE,sD = 0.5,sS = 1, meth = "extreme")
modelweights.split
modelweights <- modelweights.split[,"W.ind"]
names(modelweights)<-rownames(modelweights.split)

# To DO: 
prop_dif_pos_and_neg.fut.CMIP5 <- prop_dif_pos_and_neg[grep("CMIP5",names(prop_dif_pos_and_neg))]
if(respective.grid%in%names(sets.futures.CMIP5)){
  resp.ind<-  which(names(sets.futures.CMIP5)==respective.grid)
  prop.futures.5.sets.ind <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.fut.CMIP5)))%in%names(sets.futures.CMIP5)[-resp.ind]
} else {prop.futures.5.sets.ind <- gsub("_hc_1700_1800i_V81_equal2","",gsub("prop_","",names(prop_dif_pos_and_neg.fut.CMIP5)))%in%names(sets.futures.CMIP5)}
prop_dif_pos_and_neg.fut.CMIP5 <- prop_dif_pos_and_neg.fut.CMIP5[prop.futures.5.sets.ind]
names(prop_dif_pos_and_neg.fut.CMIP5)




###########################################################################
# hist propagation CMIP5
###########################################################################
length(modelweights)
length(prop_dif_pos_and_neg.CMIP5)
prop.hist.CMIP5.weighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg.CMIP5, y = modelweights,SIMPLIFY = FALSE)
prop.hist.CMIP5.weighted <- Reduce('+', prop.hist.CMIP5.weighted.grids)

newgrid <- prop_dif_pos_and_neg.CMIP5[[1]]
newgrid$Data <- prop.hist.CMIP5.weighted
spatialPlot(newgrid,names.attr = names(modelweights),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg.CMIP5[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop.hist.CMIP5.unweighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg.CMIP5, y = 1/length(prop_dif_pos_and_neg.CMIP5),SIMPLIFY = FALSE)
prop.hist.CMIP5.unweighted <- Reduce('+', prop.hist.CMIP5.unweighted.grids)

oldgrid <- prop_dif_pos_and_neg.CMIP5[[1]]
oldgrid$Data <- prop.hist.CMIP5.unweighted
spatialPlot(oldgrid,names.attr =  names(modelweights),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg.CMIP5[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

compareoldnew <- bindGrid(oldgrid,newgrid, dimension = c("member"))
compareoldnew$Members <- c("unweighted","weighted")
spatialPlot(compareoldnew,names.attr = c("unweighted","weighted"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

all.equal(oldgrid,newgrid)

modelweights[order(modelweights)]
binds <- prop_dif_pos_and_neg.CMIP5[order(modelweights,decreasing = TRUE)]
length(binds)
compare.all <- bindGrid(binds, dimension = "member")


compare.all.int <- bindGrid(prop_dif_pos_and_neg[[respective.grid]],compare.all,dimension = "member")

spatialPlot(compare.all.int,names.attr = c(respective.grid,names(modelweights)[order(modelweights,decreasing = TRUE)]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))



# spatialPlot(prop_dif_pos_and_neg$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2,names.attr = c("interim"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
#             col.regions= col,
#             set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# 
# bindGrid(prop.hist.CMIP5.weighted.grids,dimension = "member")
# spatialPlot(compare.all.int,names.attr = c("interim",names(modelweights)[order(modelweights,decreasing = TRUE)]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
#             col.regions= col,
#             set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

##################################################################################################################
# Compare sd weighted versus unweighted CMIP5
##################################################################################################################
ensemble.unw.prop <- t(sapply(prop_dif_pos_and_neg.CMIP5,function(x) array3Dto2Dmat(x$Data)))
sd.box.unw <-apply(ensemble.unw.prop, MARGIN = 2, FUN = sd)
w <- modelweights
# sd.cors.fut.weighted <- sqrt(sum(modelweights*(keltoC(RMSE.cors.gcms.futures.CMIP5.int)-mean.cors.fut.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
# sd.cors.fut.weighted <- sqrt(sum(w*((x)-x)^2)/((length(w)-1)/(length(w)*sum(w))))

sd.box.w <- apply(ensemble.unw.prop, MARGIN = 2, FUN = function(x) sqrt(sum(w*(x-mean(x))^2)/((length(w)-1)/(length(w)*sum(w)))))# HIER weights integreren. 
clim.sd.box.w <- quantity2clim(sd.box.w, "sd per gridbox",ref.grid = cubicsets$interim_10d_akima_cubic)
clim.sd.box.unw <- quantity2clim(sd.box.unw, "sd per gridbox",ref.grid = cubicsets$interim_10d_akima_cubic)

comparesdoldnew <- bindGrid(clim.sd.box.unw,clim.sd.box.w, dimension = c("member"))
comparesdoldnew$Members <- c("unweighted","weighted")

spatialPlot(comparesdoldnew,backdrop.theme = "coastline", 
            lonCenter = 180,col = "Reds",rev.colors = FALSE,
            names.attr = c(paste0("unweighted sd = ",round(sum(sd.box.unw))),paste0("weighted sd = ",round(sum(sd.box.w))))
            # at = seq(-0.85,0.85,0.01),region = TRUE, 
            # col.regions= col,
            #set.min = -0.85,set.max = 0.85, 
            #colorkey = list(width = 0.6, lables = list(cex = 0.5))
)

##########################################################################
# TO DO: fut con hist modelweights (do with future pesos. )
##########################################################################
prop.fut.CMIP5.weighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg.fut.CMIP5, y = modelweights,SIMPLIFY = FALSE)
prop.fut.CMIP6.weighted <- Reduce('+', prop.fut.CMIP6.weighted.grids)

newgrid <- prop_dif_pos_and_neg.fut.CMIP6[[1]]
newgrid$Data <- prop.fut.CMIP6.weighted
spatialPlot(newgrid,names.attr = names(modelweights),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg.fut.CMIP6[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop.fut.CMIP6.unweighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg.fut.CMIP6, y = 1/length(prop_dif_pos_and_neg.fut.CMIP6),SIMPLIFY = FALSE)
prop.fut.CMIP6.unweighted <- Reduce('+', prop.fut.CMIP6.unweighted.grids)

oldgrid <- prop_dif_pos_and_neg.fut.CMIP6[[1]]
oldgrid$Data <- prop.fut.CMIP6.unweighted
spatialPlot(oldgrid,names.attr =  names(modelweights),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg.fut.CMIP6[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

compareoldnew <- bindGrid(oldgrid,newgrid, dimension = c("member"))
compareoldnew$Members <- c("unweighted","weighted")
spatialPlot(compareoldnew,names.attr = c("unweighted","weighted"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

modelweights[order(modelweights)]
binds <- prop_dif_pos_and_neg.fut.CMIP6[order(modelweights,decreasing = TRUE)]
length(binds)
compare.all <- bindGrid(binds, dimension = "member")


compare.all.int <- bindGrid(prop_dif_pos_and_neg[[paste0("prop_",respective.grid,"_hc_1700_1800i_V81_equal2")]],compare.all,dimension = "member",skip.temporal.check = TRUE)

spatialPlot(compare.all.int,names.attr = c(respective.grid,names(modelweights)[order(modelweights,decreasing = TRUE)]),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

#####################################################################################################
# Plot propagation fields with respect to clustering
#####################################################################################################
load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_CMIP6_hist_vs_fut.rda")
selecteddists <- hellinger_coefficients
mat<- -log(selecteddists)
rownames(mat)
names(sets.hists.CMIP5)%in%rownames(mat)
names(sets.reans)%in%rownames(mat)
rownames<- c(names(sets.hists.CMIP5),names(sets.reans))
submat.hist5 <- mat[rownames,rownames]
longdata <- reshape2::melt(submat.hist5)

mat[which(mat == Inf, arr.ind = TRUE)]<- 100
matclust <-hclust(dist(submat.hist5, method = "minkowski", p = 2),method = "complete")

matclust$labels[matclust$order]

binds <- prop_dif_pos_and_neg[matclust$labels[matclust$order]]
compare.all <- bindGrid(binds, dimension = "member")
compare.all.int <- bindGrid(prop_dif_pos_and_neg[matclust$labels[matclust$order]],compare.all,dimension = "member",skip.temporal.check = TRUE)


plotname <-  paste0("figs/evidence_propagation/CMIP5_cluster_propagation/",whichiterations,"_plus",whichnode,"pos.pdf")
pdf(plotname,width= 3, height= 33)


spatialPlot(compare.all,names.attr = matclust$labels[matclust$order],
            as.table = FALSE,
            layout = c(1,33),
           cex = 0.01,
            backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()


plotname <-  paste0("figs/evidence_propagation/CMIP5_cluster_propagation/dendogram",whichiterations,"_plus",whichnode,"pos.pdf")
pdf(plotname,width= 20, height= 7)

plot(matclust)

dev.off()


###########################################################################
# ex.subset propagation
###########################################################################
ex.subset.prop <- c("prop_CMIP5_BNU.ESM_r1i1p1_hc_1700_1800i_V81_equal2",
                    "prop_CMIP5_ACCESS1.0_r1i1p1_hc_1700_1800i_V81_equal2",
                    "prop_CMIP5_ACCESS1.3_r1i1p1_hc_1700_1800i_V81_equal2",
                    "prop_CMIP5_CMCC.CMS_r1i1p1_hc_1700_1800i_V81_equal2")

prop_dif_pos_and_neg[ex.subset.prop]

prop.ex.subset.weighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg[ex.subset.prop], y = modelweights.ex,SIMPLIFY = FALSE)
prop.ex.subset.weighted <- Reduce('+', prop.ex.subset.weighted.grids)

newgrid <- prop_dif_pos_and_neg[[1]]
newgrid$Data <- prop.ex.subset.weighted
spatialPlot(newgrid,names.attr = sort(combinames)[-c(ind6,indearth)],as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop.ex.subset.unweighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg[ex.subset.prop], y = 1/length(ex.subset.prop),SIMPLIFY = FALSE)
prop.ex.subset.unweighted <- Reduce('+', prop.ex.subset.unweighted.grids)

oldgrid <- prop_dif_pos_and_neg[[1]]
oldgrid$Data <- prop.ex.subset.unweighted

spatialPlot(oldgrid,names.attr = sort(combinames)[-c(ind6,indearth)],as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

compareoldnew <- bindGrid(oldgrid,newgrid, dimension = c("member"))
spatialPlot(compareoldnew,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

#################################################################################
# RMSE propagation fields;
# ERA.interim hist vs CMIPs. 
# hist CMIPS vs fut CMIPs. 
#################################################################################
prop_dif_pos_and_neg$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2
props.cmip5s.reans<- prop_dif_pos_and_neg[-c(ind6)]
props.cmip5s <- props.cmip5s.reans[1:(length(props.cmip5s.reans)-3)]


x <-props.cmip5s.reans$prop_CMIP5_ACCESS1.0_r1i1p1_hc_1700_1800i_V81_equal2
x <- redim(y,drop = TRUE)

gmt <- function(x){
  
  mean.pro.lat<- rowSums(x$Data)/ncol(x$Data)
  v <- as.vector(cos(x$xyCoords$y/(180)*pi))
  weight.pro.lat <- v/sum(v)
  weighted.lat.mean <- sum(mean.pro.lat*weight.pro.lat)
  return(keltoC(weighted.lat.mean))
}

keltoC(273)
keltoC <-function(x) {if (x>=0){
  x <- (x+273.15)%%273.15
} else {x <- x}
return(x)}

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

weighted.props.cmip5s <- lapply(props.cmip5s,gtimes.grid)
weighted.prop.ncep <- gtimes.grid(props.cmip5s.reans$prop_ncep_10d_hc_1700_1800i_V81_equal2)
weighted.prop.interim <- gtimes.grid(props.cmip5s.reans$prop_interim_10d_akima_cubic_hc_1700_1800i_V81_equal2)
weighted.prop.selrean <- weighted.prop.ncep
RMSE.int.w <- numeric(length(props.cmip5s))
RMSE.ncep.w <- numeric(length(props.cmip5s))
for (i in 1:length(weighted.props.cmip5s)){
x1 <- redim(weighted.props.cmip5s[[i]],drop = TRUE)$Data
x2 <- redim(weighted.prop.interim,drop = TRUE)$Data
x3 <- redim(weighted.prop.ncep,drop = TRUE)$Data
RMSE.int.w[i] <- sqrt(sum((x1-x2)^2)/length(x1))
RMSE.ncep.w[i] <- sqrt(sum((x1-x3)^2)/length(x1))
}

names(RMSE.int.w)<- names(cubicsets)[-ind6][1:(length(names(cubicsets)[-ind6])-3)]
save(RMSE.int.w,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/RMSE.int.w.rda")

names(RMSE.ncep.w)<- names(cubicsets)[-ind6][1:(length(names(cubicsets)[-ind6])-3)]
save(RMSE.ncep.w,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/RMSE.ncep.w.rda")

