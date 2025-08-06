#############################################################################
# Model Evaluation CMIP5/CMIP6
#############################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(sparsebn)
library(sparsebnUtils)
library(RColorBrewer)
library(gaussDiff)
########################################################################################################
# Model Evaluation CMIP5/CMIP6
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
load("../Data/tas_ncep_10d.rda")
load("data/tas_JRA55_10d_akima_cubic.rda")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")
load("data/tas_historical_cmip6_left_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda")
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
cubicsets6 <- tas_historical_cmip6_10d_akima_cubic_corrected
cubicsets6left <- tas_historical_cmip6_left_10d_akima_cubic_corrected
cubicsetsCSIRO <- tas_historical_cmip5_CSIRO_10d_akima_cubic

namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
namescubics6left <- names(cubicsets6left)
namescubicsearth <- names(cubicsetsearth)
namescubicsCSIRO <- gsub("-",".",gsub(".ncml","",gsub("https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip5/historical/historical_day_","",sapply(cubicsetsCSIRO,function(x)attr(x,"dataset")))))
names(cubicsetsCSIRO) <- namescubicsCSIRO
# With CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,cubicsets6,cubicsets6left,listinterim,listjra55,listncep)
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
# Without earth, without CMIP6, without ncep, with CSIROS
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsCSIRO,listinterim,listjra55)
names(cubicsets)

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
hc_gcmsCSIRO <- lapply(paste0("CMIP5_CSIRO/",namescubicsCSIRO), loadIterations, permused = permk, it = it)
hc_gcms5 <- lapply(paste0(namescubics5), loadIterations, permused = permk, it = it)
hc_gcms5extra <- lapply(paste0("CMIP5_extra/",namescubics5extra), loadIterations, permused = permk, it = it)
hc_gcms5left<- lapply(paste0("CMIP5_left/",namescubics5left), loadIterations, permused = permk, it = it)
hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk, it = it)
hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk, it = it)
hc_gcms6left <- lapply(paste0("CMIP6_left/",namescubics6left), loadIterations, permused = permk, it = it)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_jra55 <- lapply(c("JRA55_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_ncep <- lapply(c("../Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 
# With cmip6
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_gcms6,hc_gcms6left,hc_interim,hc_jra55,hc_ncep)
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
# Without earth, without CMIP6, without ncep, with CSIROS
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsCSIRO,hc_interim,hc_jra55)
names(hc_gcms)
# Choose between 'own optimums' or constant magnitude

# modelsize <- own_optimums
modelsize <- 18
selection_hc_gcms <-  mapply(function(x,y) x[[grep((y),x)]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
selection_hc_gcms
# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
#names(selection_fits_gcms_scaled) <- names(hc_gcms)
##########################################################################
# Global dataset evaluation from BN 
##########################################################################
whichBNFits <- selection_fits_gcms_scaled
names(selection_fits_gcms_scaled)
whichData <- data_gcms_anom_scaled
whichSet <- "interim_10d_akima_cubic"
whichSet <- "ncep_10d"
whichSet <- "JRA55_10d_akima_cubic"
#########################################################################
# pickedEOFsgrids is prepared under EOFs logliks preparation. 
#####################################################################
EOFS_preparation <- function(whichPCs){
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

varpercentages <- c(0.25,0.50,0.75,0.90,0.95,0.99)
# varpercentages <- c(0.50,0.75,0.825,0.90,0.95,0.99)
pickedEOFs <- lapply(varpercentages, FUN = function(x) which(expvars > x)[[1]])
pickedEOFsdata <- lapply(pickedEOFs, FUN = function(x) EOFS_preparation(1:x))

templatedata <- gridforEOF[[whichSet]]
templatesdata <- rep(list(templatedata),length(varpercentages))
pickedEOFsgrids <- mapply(function(x,y) {z<- x; z$Data[1,,,] <- y; return(z)},templatesdata, pickedEOFsdata, SIMPLIFY  = FALSE)
pickedEOFsnames <- paste0(whichSet,"_v.exp_",as.character(varpercentages))
names(pickedEOFsgrids) <- names(pickedEOFsdata) <- pickedEOFsnames

data_pickedEOFs <- lapply(pickedEOFsdata, function(x) as.data.frame(x))
whichData <- c(data_gcms_anom_scaled,data_pickedEOFs)

ordersets <- TRUE

loglik_selection_datasets <- matrix(data = NA, nrow = length(whichBNFits), ncol = length(whichData), dimnames = list(names(whichBNFits),names(whichData)))
for(i in 1:length(whichBNFits)){
  lo <- sapply(X = whichData, logLik, object = whichBNFits[[i]])
  loglik_selection_datasets[i,] <- lo
}

institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor")
institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","ACCESS","bcc","BNU","CCSM4","CESM1","CMCC","CSIRO","inmcm4","MRI")
institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","ACCESS","BCC","BNU","CCSM4","CESM2","CMCC","CSIRO","inmcm4","MRI","NESM3","UKESM1","FGOALS","CAMS")

# with cmip6
institutions <- rep(institutions1,each =2)
CMIPS <- rep(c("CMIP5","CMIP6"),length(institutions1))
# without cmip6
institutions <- rep(institutions1,each =1)
CMIPS <- rep(c("CMIP5"),length(institutions1))
# without cmip5
institutions <- rep(institutions1,each =1)
CMIPS <- rep(c("CMIP6"),length(institutions1))

combinations <- cbind(CMIPS,institutions)
indearth <- grep("EC.EARTH_r1i|EC.EARTH_r2i",rownames(loglik_selection_datasets))

namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(loglik_selection_datasets[-indearth,-indearth]))))
namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(loglik_selection_datasets))))

shortnames[-indearth][namesort]
# For variance percentage interim:
namesort1 <- c(namesort,36)
namesort2 <- c(namesort,36,37,38,39,40,41,42)
# For variance percentage interim plus CMIP6 and ALL:
namesort1 <- c(namesort,49,50,51)
namesort2 <- c(namesort,49,50,51,52,53,54,55,56,57)
# for including interim/ncep/jrr55:
namesort1 <- c(namesort,36,37,38)
namesort2 <- c(namesort,36,37,38)
# for including cmip5 left:
namesort1 <- c(namesort,49,50,51)
namesort2 <- c(namesort,49,50,51)
# for excluding cmip6 and including cmip5 left:
namesort1 <- c(namesort,35,36,37)
namesort2 <- c(namesort,35,36,37)
# for excluding cmip6 and including cmip5 left + eofs
namesort1 <- c(namesort,35,36,37)
namesort2 <- c(namesort,35,36,37,38,39,40,41,42,43)
# for excluding cmip6 and excluding earth and including cmip5 left + eofs
namesort1 <- c(namesort,32,33,34)
namesort2 <- c(namesort,32,33,34,35,36,37,38,39,40)
# for excluding cmip6 and excluding earth and excluding eofs and including cmip5 left 
namesort1 <- c(namesort,32,33,34)
namesort2 <- c(namesort,32,33,34)
# for excluding cmip6 and excluding earth and excluding eofs and including cmip5 left 
namesort1 <- c(namesort,30,31,32)
namesort2 <- c(namesort,30,31,32)
# for excluding cmip5  
namesort1 <- c(namesort,25,26,27)
namesort2 <- c(namesort,25,26,27)

orderby <- "interim_10d_akima_cubic"
orderby <- "CMIP5_HadGEM2.ES_r1i1p1"
orderby <- "ncep_10d"
orderby <- "ncep_10d_v.exp_0.99"
orderby <- "ncep_10d_v.exp_0.9"
orderby <- "interim_10d_akima_cubic_v.exp_0.95"
sort(colSums(loglik_selection_datasets[-indearth,-indearth]))
sort(rowSums(loglik_selection_datasets[-indearth,-indearth]))
sort(colSums(loglik_selection_datasets[-indearth,-indearth]) + rowSums(loglik_selection_datasets[-indearth,-indearth]))

logsort <- sort(loglik_selection_datasets[-indearth,-indearth][,orderby],index.return = TRUE)
logclust <-hclust(dist(loglik_selection_datasets[-indearth,-indearth], method = "minkowski", p = 2),method = "complete")
logclust <-hclust(dist(t(loglik_selection_datasets[-indearth,-indearth]), method = "minkowski", p = 2),method = "complete")

loglik_selection_datasets_ordered <- loglik_selection_datasets[namesort1,namesort2]
loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][rev(logclust$order),rev(logclust$order)]
loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][namesort1,namesort2]
loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][logsort$ix,logsort$ix]
loglik_selection_datasets_ordered <- loglik_selection_datasets[logsort$ix,namesort2]
loglik_selection_datasets_ordered <- loglik_selection_datasets[logsort$ix,c(pickedEOFsnames,whichSet)]
# with hellinger distance ordering:
loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][mclust$labels[mclust$order],mclust$labels[mclust$order]]
# with hellinger distance and only eofs
loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][mclust$labels[mclust$order],c(pickedEOFsnames,whichSet)]

if (ordersets == TRUE){ longdata <- melt(loglik_selection_datasets_ordered )
} else {longdata <- melt(loglik_selection_datasets)}

########################################################
# Visualize logliks global set
########################################################

# Plot normal value 
bottom <- -300000
top <- 0
top <- max(longdata$value)

#min(longdata$value)
b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  geom_text(aes(label = round(-value/10^4,0))) +
  scale_fill_gradient2(name = "log P(data|model)\nx -10^4",
                       # low="white", 
                       # high="red",
                       limits = c(bottom,top), labels = function(x)-x/10^4) +
  labs(x="data", y="models", title=paste0(modelsize*100," model-dataset evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=15, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=15),
        plot.title=element_text(size=12)) + 
  scale_x_discrete(position = "top")
  

b             

dev.off()

########################################################
# Visualize logliks example set 
########################################################
#CMCC.CM_|
rownames(loglik_selection_datasets)
#subr <-grep("ACCESS|BNU|interim",rownames(loglik_selection_datasets),fixed = FALSE)
#subc <- grep("ACCESS|BNU|interim",colnames(loglik_selection_datasets),fixed = FALSE)

ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")
ex.var.subset <- c("interim_10d_akima_cubic_v.exp_0.25","interim_10d_akima_cubic_v.exp_0.5",
  "interim_10d_akima_cubic_v.exp_0.75", "interim_10d_akima_cubic_v.exp_0.9", 
  "interim_10d_akima_cubic_v.exp_0.95", "interim_10d_akima_cubic_v.exp_0.99")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/exampleclustering/logliks_ex_ACCESS|BNU|CMCC.CM|interim_vs_variance_2times.png")
png(plotname,units = "in", width = 10, height = 6,  res =180)
par(mfrow = (c(1,2)))
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/exampleclustering/logliks_ex_ACCESS|BNU|CMCC.CM|interim_vs_variance_2times.pdf")
pdf(plotname, width = 10, height = 6)


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/exampleclustering/logliks_ex_ACCESS|BNU|CMCC.CMS|interim_vs_variance.png")
png(plotname,units = "in", width = 11, height = 7,  res =180)
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/exampleclustering/logliks_ex_ACCESS|BNU|CMCC.CMS|interim_vs_variance.pdf")
pdf(plotname, width = 11, height = 7)
n <-8

ex.var.logliks <- loglik_selection_datasets[ex.subset,c(ex.subset,ex.var.subset)]
colnames(ex.var.logliks)
ex.longdata <- melt(ex.var.logliks)
ex.longdata$d.value <- quant.discreter(n,ex.longdata$value)
br.lims.var <- quantile(ex.longdata$value, seq(0,1,1/n))

b <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value),interpolate = FALSE) + 
  geom_text(aes(label = round(-value/10^4,0)), size = 7) +
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Reds")), guide = "legend",name = "log P(data|model)\nx -10^4",n.breaks =n,labels = round(-br.lims.var[2:(n+1)]/10^4,0))+
  labs(x="data", y="models", title=paste0(modelsize*100," model-dataset evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=15, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=15),
        plot.title=element_text(size=12)) + 
  scale_x_discrete(position = "top")

b
dev.off()
########################################################
# Visualize logliks example set 
########################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/exampleclustering/logliks_ex_ACCESS|BNU|interim.png")
png(plotname,units = "in", width = 9, height = 6,  res =180)
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/exampleclustering/logliks_ex_ACCESS|BNU|interim.pdf")
pdf(plotname, width = 9, height = 6)
n <-7
ex.longdata <- melt(loglik_selection_datasets[ex.subset,ex.subset])
ex.longdata$d.value <- quant.discreter(n,ex.longdata$value)
# br.lims <- quantile(ex.longdata$value, seq(0,1,1/n))

b <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value),interpolate = FALSE) + 
  geom_text(aes(label = round(-value/10^4,0)), size = 7) +
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Reds")), guide = "legend",name = "log P(data|model)\nx -10^4",n.breaks =n,labels = round(-br.lims.var[2:(n+1)]/10^4,0))+
  labs(x="data", y="models", title=paste0(modelsize*100," model-dataset evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=15, angle=90, vjust=1, hjust = 1),
        axis.text.y=element_text(size=15),
         plot.title=element_text(size=12))# + 
  # scale_x_discrete(position = "top")

b
dev.off()
########################################################################################
# Calculate hellinger distance bayesian networks
# KL distances give NaN (because there are x such that P(x) = 0 does not imply Q(x) = 0)
########################################################################################
fits_NEL <- lapply(selection_fits_gcms_scaled,bnlearn:::as.graphNEL)
edgelistssparse <- lapply(fits_NEL,edgeList)
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


rownames(hellinger_coefficients) 
save(hellinger_coefficients, file ="results/hellinger_coefficient/hellinger_coefficients.rda")
testkl <- normdiff(mu1 = fits_Means[[1]],mu2 = fits_Means[[2]],sigma1 = fits_COVmats[[1]], sigma2 = fits_COVmats[[2]],method = "Mahalanobis", inv = FALSE)

hellinger_coefficients_CMIP6 <-hellinger_coefficients 
save(hellinger_coefficients_CMIP6, file ="results/hellinger_coefficient/hellinger_coefficients_CMIP6.rda")
testkl <- normdiff(mu1 = fits_Means[[1]],mu2 = fits_Means[[2]],sigma1 = fits_COVmats[[1]], sigma2 = fits_COVmats[[2]],method = "Mahalanobis", inv = FALSE)

grep("CSIRO",names(fits_Means))
hellinger_coefficients<- matrix(data = NA, nrow = length(fits_COVmats[grep("CSIRO",names(fits_Means))]), ncol = length(fits_COVmats), dimnames = list(names(fits_Means[grep("CSIRO",names(fits_Means))]),names(fits_Means)))

fits_Means_CSIRO <- fits_Means[grep("CSIRO",names(fits_Means))]
fits_COVmats_CSIRO <- fits_COVmats[grep("CSIRO",names(fits_Means))]
gc()
for(i in 1:length(fits_Means_CSIRO)){
  kl <- mapply(function(b,d) normdiff(mu1 =fits_Means_CSIRO[[i]], mu2 = b, sigma1 = fits_COVmats_CSIRO[[i]], sigma2 = d, method = "Hellinger"),b = fits_Means, d = fits_COVmats)
  hellinger_coefficients[i,] <- kl
}

hellinger_coefficients_CSIRO_CMIP5 <-hellinger_coefficients 
save(hellinger_coefficients_CSIRO_CMIP5, file ="results/hellinger_coefficient/hellinger_coefficients_CSIRO_CMIP5.rda")
#testkl <- normdiff(mu1 = fits_Means[[1]],mu2 = fits_Means[[2]],sigma1 = fits_COVmats[[1]], sigma2 = fits_COVmats[[2]],method = "Mahalanobis", inv = FALSE)

#########################################################################
# Analysis hellinger distances
#########################################################################
# Partitional clustering
library("cluster")
library("factoextra")
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
load(file ="results/hellinger_coefficient/hellinger_coefficients_CMIP6.rda")

h1 <- hellinger_coefficients[grep("CMIP5",rownames(hellinger_coefficients)),grep("CMIP5",rownames(hellinger_coefficients))]
which(h1 == 0, arr.ind = TRUE)
row.names(h1[-32,-32])[order(rowSums(-log(h1[-32,-32])),decreasing = TRUE)]
row.names(h1[-32,-32])[order(apply(-log(h1[-32,-32]),MARGIN = 1, FUN =sd),decreasing = TRUE)]
par(mar = c(15,4,4,2)+0.1)
boxplot(-log(h1[-32,-32]),xaxt = "n")
axis(1,at = 1:nrow(h1[-32,-32]),labels = row.names(h1[-32,-32]),las = 2)
lines(x = c(0,40), y = c(42,42),col = 'red')

earthgrep <- unlist(apply(matrix(c("CMIP5_EC.EARTH_r1i1p1_historical","CMIP5_EC.EARTH_r2i1p1_historical","CMIP5_EC.EARTH_r1i1p1","CMIP5_EC.EARTH_r2i1p1","CMIP5_EC.EARTH_r12i1p1_historical"),nrow = 5, ncol = 1),MARGIN = 1,function(x) grep(x,rownames(hellinger_coefficients))))
hel_coef_without_earth <- hellinger_coefficients[-earthgrep,-earthgrep]
namesorthc <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(hel_coef_without_earth))))
reangrep <- unlist(apply(matrix(c("interim","ncep","JRA55"),nrow = 3, ncol = 1),MARGIN = 1,function(x) grep(x,rownames(hel_coef_without_earth))))

selecteddists <- hellinger_coefficients_CMIP6
selecteddists <- hellinger_coefficients_CMIP6[namesort1,namesort2]
selecteddists <- hel_coef_without_earth[c(namesorthc,reangrep),c(namesorthc,reangrep)]
selecteddists <- hel_coef_without_earth[namesorthc,namesorthc]
nrow(selecteddists)

klsort <- sort(selecteddists,index.return = TRUE)
# select optimal amount of clusters
fviz_nbclust(-log(selecteddists), 
             FUNcluster = cluster::pam, 
             method = "silhouette", 
             k.max =25,medoids = NULL)
fviz_nbclust(-log(selecteddists),FUNcluster = cluster::pam,method = "gap_stat",k.max = 10)
# Use PAM to perform k-medoids
pamcluster <- pam(-log(selecteddists),k =2, 
                  medoids = c(30,31), 
                  do.swap = FALSE)
# info
names(selecteddists)
pamcluster$medoids
pamcluster$id.med
# visualize k-medoids
par(mar = c(5,10,3,4)+0.1)
plot(silhouette(pamcluster),max.strlen = 200, nmax.lab =40,cex = 0.7)

fviz_cluster(pamcluster, repel = TRUE, axes = c(1,2),show.clust.cent = TRUE)
###############################################################
# Hierarchical clustering batachari
###############################################################
library(proxy)
mclust <-hclust(dist(-log(selecteddists),method = "minkowski", p = 2),method = "complete")
mclust <-hclust(dist(-log(selecteddists),method = "euclidean"))
plot(mclust, hang = -1)
rect.hclust(mclust,k=11)
dev.off()

orderby <- "ncep_10d"
orderby <- "interim_10d_akima_cubic"
ordersort <- sort(-log(selecteddists)[,orderby],index.return = TRUE,decreasing = TRUE)

ordersets <- TRUE
if (ordersets == TRUE){ longdata <- melt(-log(selecteddists)[ordersort$ix,ordersort$ix])
} else {longdata <- melt(-log(selecteddists))}

# order by clustering.

get_upper_tri <- function(data2){
  data2[upper.tri(data2)]<- NA
  return(data2)
}
orderedmat<- -log(selecteddists)[mclust$order,mclust$order]
upper_tri <- get_upper_tri(orderedmat)
longdata <- melt(upper_tri, na.rm = TRUE)
###############################################################
# Hierarchical clustering hellinger distance
###############################################################
library(proxy)
mclust <-hclust(dist(sqrt(1-(selecteddists)),method = "minkowski", p = 2),method = "complete")
mclust <-hclust(dist(-log(selecteddists),method = "euclidean"))
plot(mclust, hang = -1)
rect.hclust(mclust,k=11)
dev.off()

names(data_gcms_anom_scaled[-indearth])
phy <-as.phylo(mclust)
phy$tip.label
datamat<-do.call(cbind, lapply(lapply(data_gcms_anom_scaled[-indearth],rowSums),mean))
treeSeg(datamat,phy,checkOrder = FALSE, fam = "gauss", alpha = 0.05)
row
orderby <- "ncep_10d"
orderby <- "interim_10d_akima_cubic"
ordersort <- sort(-log(selecteddists)[,orderby],index.return = TRUE,decreasing = TRUE)

ordersets <- TRUE
if (ordersets == TRUE){ longdata <- melt(-log(selecteddists)[ordersort$ix,ordersort$ix])
} else {longdata <- melt(-log(selecteddists))}

# order by clustering.

get_upper_tri <- function(data2){
  data2[upper.tri(data2)]<- NA
  return(data2)
}
orderedmat<- -log(selecteddists)[mclust$order,mclust$order]
upper_tri <- get_upper_tri(orderedmat)
longdata <- melt(upper_tri, na.rm = TRUE)

apply(orderedmat,1,sd)[order(apply(orderedmat,1,sd))]
hist(orderedmat)
colMeans(orderedmat)[order(colMeans(orderedmat))]
# ########################################################
# # Visualize hellinger distance example set 
# ########################################################
# ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")
# library(RColorBrewer)
# plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/Hellinger_CMIP5_ex_ACCESS|BNU|CMCC.CMS|interim.png")
# png(plotname,units = "in", width = 9, height = 7,  res =180)
# plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/Hellinger_CMIP5_ex_ACCESS|BNU|CMCC.CMS|interim.pdf")
# pdf(plotname, width = 9, height =7)
# # n <-7
# n <-7
# br.lims <- c(0,30,40,50,60,70,80)
# n <-9
# br.lims <- c(0,33,37,41,44,47,50,57,89)
#   #quantile(ex.longdata$value, seq(0,1,1/n))
# ex.longdata <- melt(orderedmat[ex.subset,ex.subset])
# 
# quant <- c(0,30,40,50,60,70,80)
# quant <- c(0,30,40,50,60,70,80)
# quant  <- c(0,33,37,41,44,47,50,57,89)
# 
# quant.discreter <- function(n,data.vec,quant = "quant"){
#   if (quant == "quant"){br.lims <- quantile(data.vec, seq(0,1,1/n))} else {br.lims <- quant}
#   d.vec <- numeric(length = length(data.vec))
#   # j <- 1
#   for(j in 1:length(data.vec)){
#     y <- data.vec[j]
#     i <- 1
#     while(i <= (length(br.lims)-1)){
#       if(y <= br.lims[1]){
#         d.vec[j] <- 1
#         break
#       } else if (y >= br.lims[i] & y <= br.lims[i+1]){
#         # if(y >= br.lims[i] & y <= br.lims[i+1]){
#         d.vec[j] <- i
#         break
#       } else {i <- i +1}
#     } 
#   }
#   return(d.vec)
# }
# 
# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value,  quant = c(0,33,37,41,44,47,50,57,89))
# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value,  quant = "quant")
# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value, quant = c(0,30,40,50,60,70,80))
# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value)
# 
# b <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
#   geom_raster(aes(fill=d.value)) + 
#   scale_fill_gradientn(aesthetics = c("fill"), breaks = NULL,colours = rev(brewer.pal(length(quant),"Blues")), guide = "legend",name = "Hellinger distance",labels = c(0,33,37,41,44,47,50,57,89))+
#   #scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[2:(n+1)],0))+
#   geom_text(aes(label = round(value,0)),size = 7) +
#   labs(x="models", y="models", title=paste0("hellinger distance ",modelsize*100," CMIP models")) +
#   theme_bw() + 
#   theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
#         axis.text.y=element_text(size=17),
#         plot.title=element_text(size=17)) + 
#   scale_y_discrete(position = "left") + 
#   scale_x_discrete(position = "top")
# 
# b
# 
# dev.off()
########################################################
# Visualize  Hellinger Distance CSIRO subset
########################################################

load(file ="results/hellinger_coefficient/hellinger_coefficients_CSIRO_CMIP5.rda")
ver <- colnames(hellinger_coefficients_CSIRO_CMIP5)
ind_namescluster_csiro <- which(ver%in%rownames(hellinger_coefficients_CSIRO_CMIP5)[-1])[c(2:10,1)]
ind_namescluster_cmippaper <- which(ver%in%namescluster_cmippaper)

ver2 <- ver[c(ind_namescluster_cmippaper,ind_namescluster_csiro)]
out <- which(ver2 %in% c("CMIP5_bcc.csm1.1.m_r1i1p1","CMIP5_CSIRO.Mk3.6.0_r1i1p1"))
ver3 <- ver2[-out]
ver3 <- c(sort(ver3[1:27]),ver3[28:length(ver3)])
batmat<- -log(hellinger_coefficients_CSIRO_CMIP5[-1,c(42,43,1:41)])
batmat <- -log(hellinger_coefficients_CSIRO_CMIP5[rownames(hellinger_coefficients_CSIRO_CMIP5)[-1][c(2:10,1)],ver3])
ex.longdata <- melt(t(batmat))

quant.discreter <- function(n,data.vec,quant = NULL){
  if (is.null(quant)){br.lims <- quantile(data.vec, seq(0,1,1/n))} else {br.lims <- quant}
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y < br.lims[1]){
        d.vec[j] <- 1
        break
      } else if (y >= br.lims[i] & y < br.lims[i+1]){
        # if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}



# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value, quant = c(0,30,40,50,60,70,80))
# ex.longdata$d.value <- quant.discreter(7,ex.longdata$value)
n <- 8

ex.longdata$d.value <- quant.discreter(n,round(ex.longdata$value), quant = c(0,33,37,41,44,47,50,57,89))


# br.lims<- quantile(ex.longdata$value, seq(0,1,1/n))

b <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  geom_text(aes(label = round(value,0)), size =7) +
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",breaks = 1:(n),labels = round(br.lims[2:(n+1)],0))+
labs(x="models", y="models", title=paste0("Bhattacharya distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")


plotname <- paste0("figs/hellinger_CMIP5_ordering/Hellinger_CSIRO_CMIP5.png")
png(plotname,units = "in", width = 10, height = 20,  res =180)
plotname <- paste0("figs/hellinger_CMIP5_ordering/Hellinger_CSIRO_CMIP5.pdf")
pdf(plotname, width = 10, height = 20)
b

dev.off()

########################################################
# Visualize  Hellinger Distance global set diagonal
########################################################
# Plot normal value 
bottom <- 32
top <- 0
top <- max(longdata$value)
library(RColorBrewer)
display.brewer.all()
# min(longdata$value)
n <- 4
br.lims <- quantile(longdata$value, seq(0,1,1/n))
data.vec <- longdata$value

quant.discreter <- function(n,data.vec){
  br.lims <- quantile(data.vec, seq(0,1,1/n))
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}

longdata$d.value <- quant.discreter(n,longdata$value)

b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  #discrete_scale(scale.name = d.value,palette = brewer.pal(5,"Reds"))
  #guides(fill = guide_colorsteps())+
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[2:(n+1)],0))+
  geom_text(aes(label = round(value,0))) +
 # scale_fill_brewer(aesthetics = longdata$d.value,
                 #scale_name = "fill_gradient2",
                 #name = "value",
                 #na.value = "white",
  #               palette = "Reds" )+
                # labels = function(x)x) +
                 #limits = c(bottom,top))+
  labs(x="data", y="models", title=paste0("hellinger distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=12, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=12),
        plot.title=element_text(size=12)) +
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_wrt_quantileclustering_withoutearth_diagonal.pdf"
pdf(plotname,width = 10, height = 10)
b             
dev.off()

##############################################################################
# Visualize  Hellinger Distance global set full matrix
##############################################################################
code <- FALSE
if(code == TRUE) {r.i. <- sapply(rownames(orderedmat),function(x)which(x == codes$rn.orderedmat))
toChange <- as.logical(sapply(r.i., length))
rownames(orderedmat)[toChange]<- codes$abrev.y.atm[unlist(r.i.[toChange])]
colnames(orderedmat)[toChange]<- codes$abrev.y.atm[unlist(r.i.[toChange])]}

rownames(orderedmat) <- gsub("CMIP6Amon_","",gsub("historical","",rownames(orderedmat)))
colnames(orderedmat) <- gsub("CMIP6Amon_","",gsub("historical","",colnames(orderedmat)))
longdata <- melt(orderedmat)
n <- 7
br.lims <- quantile(longdata$value, seq(0,1,1/n))
data.vec <- longdata$value

quant.discreter <- function(n,data.vec){
  br.lims <- quantile(data.vec, seq(0,1,1/n))
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}

longdata$d.value <- quant.discreter(n,longdata$value)

b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",
                       name = "value",labels = round(br.lims[2:(n+1)],0))+
  geom_text(aes(label = round(value,0))) +
  labs(x="data", y="models", title=paste0("hellinger distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=12, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=12),
        plot.title=element_text(size=12)) +
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")


if (code == TRUE){
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_wrt_quantileclustering_withoutearth_full_code.pdf"
pdf(plotname,width = 10, height = 10)} else {
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_wrt_quantileclustering_withoutearth_full.pdf"
pdf(plotname,width = 10, height = 10)}
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_wrt_quantileclustering_withoutearth_full.png"
png(plotname, width = 10, height = 10, units = "in", res = 180)
b             
dev.off()

# CMIP6
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP6_ordering/hellinger_CMIP6_wrt_quantileclustering_withoutearth_full.pdf"
pdf(plotname,width = 10, height = 10)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP6_ordering/hellinger_CMIP6_wrt_quantileclustering_withoutearth_full.png"
png(plotname, width = 10, height = 10, units = "in", res = 180)

b             
dev.off()
###################################################################
# CODE
###################################################################
codes <- read.csv(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIPs Metadata - Blad3.csv", header = TRUE, sep = ",",row.names = 1)
# eq.codes <- gsub(codes$Model.Name,pattern = "-",replacement = ".")
eq.codes <- gsub(rownames(codes),pattern = "-",replacement = ".")
eq.codes1 <- sapply(eq.codes, function(x)grep(x,abrev,ignore.case = TRUE),simplify = TRUE)
eq.codes.abrev <- sapply(eq.codes1, function(x) if(length(x)==1){return(x)} else {x <- which(min(abrev[x])==abrev); return(x)})
codes$eq.name <- abrev[eq.codes.abrev]

eq.codes1 <- sapply(eq.codes, function(x)grep(x,rownames(orderedmat),ignore.case = TRUE),simplify = TRUE)
eq.codes.orderedmat <- sapply(eq.codes1, function(x) if(length(x)==1){return(x)} else {x <- x[which(nchar(rownames(orderedmat))[x]==min(nchar(rownames(orderedmat))[x]))]; return(x)})
x <- eq.codes1$CMCC.CM
x[which(nchar(rownames(orderedmat))[x]==min(nchar(rownames(orderedmat))[x]))]

codes$rn.orderedmat <- rownames(orderedmat)[eq.codes.orderedmat]

codes$abrev.y.atm <- apply(codes, MARGIN = 1, function(x)paste0(x["eq.name"]," (",x["Atmosphere.CODE"],")"))
########################################
# CMIP5 CMIP6 only models for comparison 
########################################
library("cluster")
library("factoextra")
load(file ="results/hellinger_coefficient/hellinger_coefficients_CMIP6.rda")

institutions.overlap <- c("CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","CESM","MRI")



institutions <- rep(institutions.overlap,each =1)
CMIPS.6 <- rep(c("CMIP6"),length(institutions.overlap))
combinations.6 <- cbind(CMIPS.6,institutions)
namesort.6 <- unlist(apply(combinations.6,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(hellinger_coefficients_CMIP6))))
selecteddists <- hellinger_coefficients_CMIP6
orderedmat<- -log(selecteddists)[namesort.6,namesort.6]
longdata <- melt(orderedmat)
rownames(orderedmat) <- gsub("CMIP6Amon_","",gsub("historical","",rownames(orderedmat)))
colnames(orderedmat) <- gsub("CMIP6Amon_","",gsub("historical","",colnames(orderedmat)))

longdata <- melt(orderedmat)

quant <- c(0,30,40,50,60,70)


quant.discreter <- function(n,data.vec,quant = NULL){
  if (is.null(quant)){
    br.lims <- quantile(data.vec, seq(0,1,1/n))
  } else {br.lims <- quant}
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y <= br.lims[1]){
        d.vec[j] <- 1
        break
      } else if (y >= br.lims[i] & y <= br.lims[i+1]){
        # if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}

longdata$d.value <- quant.discreter(n,longdata$value, quant = c(0,30,40,50,60,70))

b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(length(quant),"Blues")), guide = "legend",name = "Hellinger distance",labels = c(0,30,40,50,60))+
  # scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[2:(n+1)],0))+
  geom_text(aes(label = round(value,0)),size = 7) +
  labs(x="models", y="models", title=paste0("hellinger distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

b

# cmip6
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP6_ordering/hellinger_CMIP6_historical_overlap.pdf"
pdf(plotname,width = 12, height = 10)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP6_ordering/hellinger_CMIP6_historical_overlap.png"
png(plotname, width = 12, height = 10, units = "in", res = 180)

b             
dev.off()


########################################
# CMIP5 CMIP6 only models for comparison: CMIP5
########################################
library("cluster")
library("factoextra")
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")


institutions.overlap <- c("CNRM","Can","EC.EARTH_r12i1p1","GFDL","IPSL","MIROC","MPI","Nor","CESM","MRI")
#earthgrep <- unlist(apply(matrix(c("CMIP5_EC.EARTH_r1i1p1_historical","CMIP5_EC.EARTH_r2i1p1_historical","CMIP5_EC.EARTH_r1i1p1","CMIP5_EC.EARTH_r2i1p1","CMIP5_EC.EARTH_r12i1p1_historical"),nrow = 5, ncol = 1),MARGIN = 1,function(x) grep(x,rownames(hellinger_coefficients))))
#hel_coef_without_earth <- hellinger_coefficients[-earthgrep,-earthgrep]
#namesorthc <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(hel_coef_without_earth))))
#reangrep <- unlist(apply(matrix(c("interim","ncep","JRA55"),nrow = 3, ncol = 1),MARGIN = 1,function(x) grep(x,rownames(hel_coef_without_earth))))


# cmip5
institutions <- rep(institutions.overlap,each =1)
CMIPS.5 <- rep(c("CMIP5"),length(institutions.overlap))
combinations.5 <- cbind(CMIPS.5,institutions)
namesort.5 <- unlist(apply(combinations.5,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(hellinger_coefficients))))
selecteddists <- hellinger_coefficients
orderedmat<- -log(selecteddists)[namesort.5,namesort.5]
longdata <- melt(orderedmat)
rownames(orderedmat) <- gsub("CMIP5_","",gsub("historical","",rownames(orderedmat)))
colnames(orderedmat) <- gsub("CMIP5_","",gsub("historical","",colnames(orderedmat)))

longdata <- melt(orderedmat)

quant <- c(0,30,40,50,60,70)


quant.discreter <- function(n,data.vec,quant = NULL){
  if (is.null(quant)){
    br.lims <- quantile(data.vec, seq(0,1,1/n))
  } else {br.lims <- quant}
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y <= br.lims[1]){
        d.vec[j] <- 1
        break
      } else if (y >= br.lims[i] & y <= br.lims[i+1]){
        # if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}

longdata$d.value <- quant.discreter(n,longdata$value, quant = c(0,30,40,50,60,70))

b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(length(quant),"Blues")), guide = "legend",name = "Hellinger distance",labels = c(0,30,40,50,60))+
  # scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[2:(n+1)],0))+
  geom_text(aes(label = round(value,0)),size = 7) +
  labs(x="models", y="models", title=paste0("hellinger distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

b

# cmip5
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_CMIP5_historical_overlap.pdf"
pdf(plotname,width = 12, height = 10)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_CMIP5_historical_overlap.png"
png(plotname, width = 12, height = 10, units = "in", res = 180)

b             
dev.off()

#########################################################################
# Analysis hellinger distances w.r.t Interim /NCEP /JRA55 (supp material table)
#########################################################################
# Partitional clustering
library("cluster")
library("factoextra")
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
load(file ="results/hellinger_coefficient/hellinger_coefficients_CMIP6.rda")

names(sort(-log(hellinger_coefficients)["interim_10d_akima_cubic", ]))
names(sort(-log(hellinger_coefficients)["JRA55_10d_akima_cubic", ]))
names(sort(-log(hellinger_coefficients)["ncep_10d", ]))

names(sort(-log(hellinger_coefficients_CMIP6)["interim_10d_akima_cubic", ]))
names(sort(-log(hellinger_coefficients_CMIP6)["JRA55_10d_akima_cubic", ]))
names(sort(-log(hellinger_coefficients_CMIP6)["ncep_10d", ]))

###########################################################################
# give weights by hellinger distance
###########################################################################
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")

bexp <- -log(hellinger_coefficients[ex.subset,ex.subset])[-2,-2]
ref.perf <- "interim_10d_akima_cubic"
refID <- which(rownames(bdists) == ref.perf)
D <- bdists[ref.perf,-refID]
S <- bdists[-refID,-refID]
namei <- "CMIP5_BNU.ESM_r1i1p1"
sS <- 30
sD <- 25

IPweight.old <- function(name.i,D,S,sD,sS){
  nameID.S <- which(rownames(S) == name.i)
  nameID.D <- which(names(D) == name.i)
  Di <- D[nameID.D]
  Sij <- S[nameID.S,][-nameID.S]
  sumS <- sum(exp(-Sij^2/sS^2))
  exp(-Di^2/sD^2)/(1+sumS)
}

IPweight <- function(bdists,ref.perf,sD,sS){
  refID <- which(rownames(bdists) == ref.perf)
  D <- bdists[ref.perf,-refID]
  S <- bdists[-refID,-refID]
  W <- numeric(length=length(D))
  for(i in 1:length(names(D))){
    name.i <-names(D)[i]
    nameID.S <- which(rownames(S) == name.i)
    nameID.D <- which(names(D) == name.i)
    Di <- D[nameID.D]
    Sij <- S[nameID.S,][-nameID.S]
    sumS <- sum(exp(-Sij^2/sS^2))
    W[i] <- exp(-Di^2/sD^2)/(1+sumS)
  }
  W <- W/sum(W)
  names(W) <- names(D)
  return(W)
}

# IPweight <- function(bdists,ref.perf,sD,sS){
#   refID <- which(rownames(bdists) == ref.perf)
#   D <- bdists[ref.perf,-refID]
#   S <- bdists[ref.perf,-refID]
#   W <- numeric(length=length(D))
#   for(i in 1:length(names(D))){
#     name.i <-names(D)[i]
#     nameID.S <- which(names(S) == name.i)
#     nameID.D <- which(names(D) == name.i)
#     Di <- D[nameID.D]
#     Sij <- min(abs(S[nameID.S]-S[-nameID.S]))
#     sumS <- sum(exp(-Sij^2/sS^2))
#     W[i] <- exp(-Di^2/sD^2)/(1+sumS)
#   }
#   W <- W/sum(W)
#   names(W) <- names(D)
#   return(W)
# }

##########################################
# Model weights cmip5 set
# 40 25 used
##########################################
matri <--log(hellinger_coefficients)[1:(nrow(hellinger_coefficients)-2),1:(ncol(hellinger_coefficients)-2)]
modelweights <- IPweight(matri,"interim_10d_akima_cubic",100,10)

par(mar = c(8, 4, 0, 2))
plot(modelweights,xaxt =  "n",xlab = "")
xlabels=gsub("r1i1p1","",gsub("historical","",names(modelweights)))
axis(1, at=1:length(modelweights),labels = xlabels, col.axis="red", las=2)
abline(h = 1/length(modelweights), lty = 2, col = "red")
###########################################
# Model weights subset example 
# 40 25 used so far
###########################################
modelweights.ex <- IPweight(-log(hellinger_coefficients[ex.subset,ex.subset]),
                            "interim_10d_akima_cubic",sD = 10, sS = 100)
bdists <- hellinger_coefficients
par(mar = c(5,4,4,2)+0.1)
plot(modelweights.ex,xaxt =  "n", ylim = c(0,0.6), ylab = expression("modelweights w"[i]),xlab = "", pch = 16, col = "blue")

abline(h = 0.25, lty = 2, col = "red")
xlabels=gsub("r1i1p1","",gsub("historical","",names(modelweights.ex)))
axis(1, at=1:length(modelweights.ex),labels = xlabels, col.axis="red", las=2)

###############################################################
# Hierarchical clustering hellinger distance JJRA NCEP Interim apart. 
###############################################################
r.ids <-sapply(c("interim_10d_akima_cubic","ncep_10d","JRA55_10d_akima_cubic"),function(x)which(rownames(selecteddists)==x))
selecteddists.cmips <-selecteddists[-r.ids,-r.ids]
selecteddists.reans <- selecteddists[r.ids,r.ids]
library(proxy)
mclust.cmips <-hclust(dist(-log(selecteddists.cmips),method = "euclidean"))
mclust.cmips <-hclust(as.dist(-log(selecteddists.cmips)))
plot(-log(selecteddists.cmips))
as.matrix(dist(-log(selecteddists.cmips),method = "euclidean"))
as.matrix(as.dist(-log(selecteddists.cmips)))
as.matrix(-log(selecteddists.cmips))
plot(mclust.cmips, hang = -1)
rect.hclust(mclust,k=11)
dev.off()
##############################################################################
# Visualize  Hellinger Distance global set full matrix
##############################################################################

cmips.ids <-sapply(rownames(selecteddists.cmips)[mclust.cmips$order],function(x)which(rownames(selecteddists)==x))
r.ids
orderedmat.cmips<- -log(selecteddists)[c(cmips.ids,r.ids),c(cmips.ids,r.ids)]
namescluster_cmippaper <- rownames(orderedmat.cmips)

longdata <- melt(orderedmat.cmips)
n <- 8
br.lims <- quantile(longdata$value, seq(0,1,1/n))
data.vec <- longdata$value

quant.discreter <- function(n,data.vec){
  br.lims <- quantile(data.vec, seq(0,1,1/n))
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}

longdata$d.value <- quant.discreter(n,longdata$value)

b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",breaks = 1:(n),labels = round(br.lims[2:(n+1)],0))+
  # scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),breaks = 1:(n),guide = "legend",name = "value",labels = round(br.lims[1:(n+1)],0))+
  geom_text(aes(label = round(value,0))) +
  labs(x="",y ="",title=paste0("Bhattacharyi distance",modelsize*100," CMIP models")) +
  # theme_bw() + 
  theme(axis.text.x=element_text(size=12, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=12),
        plot.title=element_text(size=12)) +
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")


#if (code == TRUE){
 # plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_wrt_quantileclustering_withoutearth_full_code.pdf"
  #pdf(plotname,width = 10, height = 10)} else {
    plotname <- "figs/hellinger_CMIP5_ordering/bhattacharyi_CMIP5clustering.pdf"
    pdf(plotname,width = 10, height = 10)
    #}
plotname <- "figs/hellinger_CMIP5_ordering/bhattacharyi_CMIP5clustering.png"
png(plotname, width = 10, height = 10, units = "in", res = 180)
b             
dev.off()
########################################################################
#
########################################################################

  ########################################################
# Visualize hellinger distance example set 
########################################################
ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")
library(RColorBrewer)
plotname <- paste0("figs/hellinger_CMIP5_ordering/Hellinger_CMIP5_ex_ACCESS|BNU|CMCC.CMS|interim_bottom33.png")
png(plotname,units = "in", width = 9, height = 7,  res =180)
plotname <- paste0("figs/hellinger_CMIP5_ordering/Hellinger_CMIP5_ex_ACCESS|BNU|CMCC.CMS|interim_bottom33.pdf")
pdf(plotname, width = 9, height =7)
# # n <-7
# n <-7
# br.lims <- c(0,30,40,50,60,70,80)
# #quantile(ex.longdata$value, seq(0,1,1/n))

ex.longdata <- melt(-log(selecteddists)[ex.subset,ex.subset])
ex.longdata <- melt(orderedmat[ex.subset,ex.subset])

# quant = c(0,30,40,50,60,70,80)

# quant.discreter <- function(n,data.vec,quant = NULL){
#   if (is.null(quant)){br.lims <- quantile(data.vec, seq(0,1,1/n))} else {br.lims <- quant}
#   d.vec <- numeric(length = length(data.vec))
#   # j <- 1
#   for(j in 1:length(data.vec)){
#     y <- data.vec[j]
#     i <- 1
#     while(i <= (length(br.lims)-1)){
#       if(y <= br.lims[1]){
#         d.vec[j] <- 1
#         break
#       } else if (y >= br.lims[i] & y <= br.lims[i+1]){
#         # if(y >= br.lims[i] & y <= br.lims[i+1]){
#         d.vec[j] <- i
#         break
#       } else {i <- i +1}
#     } 
#   }
#   return(d.vec)
# }
quant.discreter <- function(n,data.vec,quant = NULL){
  if (is.null(quant)){br.lims <- quantile(data.vec, seq(0,1,1/n))} else {br.lims <- quant}
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y < br.lims[1]){
        d.vec[j] <- 1
        break
      } else if (y >= br.lims[i] & y < br.lims[i+1]){
        # if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}

ex.longdata$d.value <- quant.discreter(n,round(ex.longdata$value,0), quant = c(0,33,37,41,44,47,50,57,89))
# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value)

b <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=round(d.value,0))) + 
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",breaks = 1:(n),labels = round(br.lims[2:(n+1)],0))+
  
  # scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(length(quant),"Blues")), guide = "legend",name = "Hellinger distance",labels = c(0,30,40,50,60,70))+
  #scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[2:(n+1)],0))+
  geom_text(aes(label = round(value,0)),size = 7) +
  labs(x="models", y="models", title=paste0("hellinger distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

b

dev.off()
