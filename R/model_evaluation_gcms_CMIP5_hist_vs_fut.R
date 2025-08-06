#############################################################################
# Model Evaluation CMIP5 hist vs fut
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


# cubicsets6 <- tas_historical_cmip6_10d_akima_cubic_corrected
# cubicsets6left <- tas_historical_cmip6_left_10d_akima_cubic_corrected
# # With CMIP6
# cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,cubicsets6,listinterim,listjra55,listncep)
# names(cubicsets)


cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsets5.rcp85,cubicsets5left.rcp85,listinterim,listjra55,listncep)
names(cubicsets)


namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics5.rcp85 <- names(cubicsets5.rcp85)
namescubics5left.rcp85 <- names(cubicsets5left.rcp85)

# namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
# namescubics6left <- names(cubicsets6left)
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


# without cmip6
institutions <- rep(institutions1,each =1)
CMIPS <- rep(c("CMIP5"),length(institutions1))

combinations <- cbind(CMIPS,institutions)
indearth <- grep("EC.EARTH_r1i|EC.EARTH_r2i",rownames(loglik_selection_datasets))

namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(loglik_selection_datasets[-indearth,-indearth]))))
namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(loglik_selection_datasets))))

shortnames[-indearth][namesort]

# for excluding cmip6 and including cmip5 left + eofs
namesort1 <- c(namesort,63,64,65)
namesort2 <- c(namesort,63,64,65,66,67,68,69,70,71)


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

hellinger_coefficients_CMIP5_hist_vs_fut <-hellinger_coefficients 
save(hellinger_coefficients_CMIP5_hist_vs_fut, file ="results/hellinger_coefficient/hellinger_coefficients_CMIP5_hist_vs_fut.rda")





####################################
# Visualize
####################################
load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_hist_vs_fut.rda")

selecteddists <- hellinger_coefficients_CMIP5_hist_vs_fut
orderedmat<- -log(selecteddists)
longdata <- melt(orderedmat)
rownames(orderedmat) <- gsub("CMIP5_","",gsub("historical","",rownames(orderedmat)))
colnames(orderedmat) <- gsub("CMIP5_","",gsub("historical","",colnames(orderedmat)))

longdata <- melt(orderedmat)
n = 7
br.lims <- quantile(longdata$value, seq(0,1,1/n))
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

longdata$d.value <- quant.discreter(n,longdata$value)

b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  # scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(length(quant),"Blues")), guide = "legend",name = "Hellinger distance",labels = c(0,30,40,50,60))+
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[2:(n+1)],0))+
  geom_text(aes(label = round(value,0)),size = 7) +
  labs(x="models", y="models", title=paste0("hellinger distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

b

#############################################################################################
# Emergent constraints: (1) Hellinger distance vs global mean temperature increase
#############################################################################################

futures <- grep("rcp85",names(cubicsets))
names(cubicsets)[futures]
reans <- futures[length(futures)]+ 1:3
hists <- 1:(futures[1]-1)
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
names(sets.hists)
names(sets.futures)

clim.sd.hists <- lapply(sets.hists, climatology,clim.fun = list(FUN = sd, na.rm = TRUE))
clim.sd.hists <- lapply(clim.sd.hists, redim, drop = TRUE)
clim.hists <- lapply(sets.hists, climatology)
clim.hists <- lapply(clim.hists,redim,drop = TRUE)

clim.sd.futures <- lapply(sets.futures, climatology,clim.fun = list(FUN = sd, na.rm = TRUE))
clim.sd.futures <- lapply(clim.sd.futures, redim, drop = TRUE)
clim.futures <- lapply(sets.futures,climatology)
clim.futures <- lapply(clim.futures,redim,drop = TRUE)

clim.sd.reans <- lapply(sets.reans, climatology,clim.fun = list(FUN = sd, na.rm = TRUE))
clim.sd.reans <- lapply(clim.sd.reans, redim, drop = TRUE)
clim.reans <- lapply(sets.reans, climatology)
clim.reans <- lapply(clim.reans,redim,drop = TRUE)


gmt <- function(x){
  mean.pro.lat<- rowSums(x$Data)/ncol(x$Data)
  v <- as.vector(cos(x$xyCoords$y/(180)*pi))
  weight.pro.lat <- v/sum(v)
  weighted.lat.mean <- sum(mean.pro.lat*weight.pro.lat)
  return(weighted.lat.mean)
}

gtimes <- function(x){
# x <- sets.hists$CMIP5_ACCESS1.0_r1i1p1
v <- as.vector(cos(x$xyCoords$y/(180)*pi))
m <- x$Data
for(i in 1:dim(x$Data)[2]){
m[,i,] <- x$Data[,i,]*v[i]}
weighted.gtimeseries <- apply(m, MARGIN= c(1),FUN =sum)/sum(v*dim(x$Data)[3])
return(weighted.gtimeseries)
}


# gmt.sd <- function(x){
#   data.pro.lat <- apply(x$Data,MARGIN = c(1,3), function(y)1/dim(y)[2]*sum(y))
#   nlat <- dim(x$Data)[2]
#   sd.pro.lat <- numeric(nlat)
#   i <- 1
#   for(i in 1:nlat){
#     data.pro.lat <- x$Data[,i,]
#     times.pro.lat <- rowSums(data.pro.lat)/ncol(x$Data[,i,])
#   sd.pro.lat[i] <-  sum(cov(data.pro.lat))*(1/nlat)^2
#   }
#   
#   v <- as.vector(cos(x$xyCoords$y/(180)*pi))
#   weight.pro.lat <- v/sum(v)
#   weighted.lat.sd <- sum(weight.pro.lat%o% weight.pro.lat*cov(sd.pro.lat))
#   return(weighted.lat.sd)
# }

gt.hists.cmips <- lapply(sets.hists,gtimes)
gmt.hists.cmips2 <- sapply(gt.hists.cmips,mean)
gsdt.hists.cmips <- sapply(gt.hists.cmips,sd) # sd with respect to timeseries, in 30 years.

keltoC <-function(x)(x+273)%%273
gt.futures.cmips <- lapply(sets.futures,gtimes)
gt.futures.cmips <- lapply(gt.futures.cmips,keltoC)
gmt.futures.cmips2 <- sapply(gt.futures.cmips,mean)
gsdt.futures.cmips <- sapply(gt.futures.cmips,sd)


gt.reans <- lapply(sets.reans,gtimes)
gt.reans <- lapply(gt.reans,keltoC)
gmt.reans2 <- sapply(gt.reans,mean)
gsdt.reans <- sapply(gt.reans,sd)

gmt.diffs.cmips2 <- gmt.futures.cmips2 - gmt.hists.cmips2
sd(gmt.diffs.cmips2) # sd of Delta gmt with respect to cmips unweigthed
mean(gmt.diffs.cmips2) # mean of Delta gmt with respect to cmips unweighted

#gmt.sd.hists.cmips <- sapply(clim.sd.hists,gmt)
gmt.hists.cmips <- sapply(clim.hists,gmt)
all.equal(gmt.hists.cmips,gmt.hists.cmips2)
#gmt.sd.futures.cmips <- sapply(clim.sd.futures,gmt)
gmt.futures.cmips <- sapply(clim.futures,gmt)
gmt.futures.cmips <- keltoC(gmt.futures.cmips)
all.equal(gmt.futures.cmips,gmt.futures.cmips2)
#gmt.sd.reans <- sapply(clim.sd.reans,gmt)
gmt.reans <- sapply(clim.reans,gmt)
gmt.reans <- keltoC(gmt.reans)
all.equal(gmt.reans,gmt.reans2)

# gmt.sd.difs.cmips.vs.int <- gmt.sd.hists.cmips - gmt.sd.reans["interim_10d_akima_cubic"]
gmt.mean.difs.cmips.vs.int <- gmt.hists.cmips - gmt.reans["interim_10d_akima_cubic"]
gmt.mean.difs.cmips.vs.ncep <- gmt.hists.cmips  - gmt.reans["ncep_10d"]
gmt.mean.difs.cmips.vs.jra55 <- gmt.hists.cmips  - gmt.reans["JRA55_10d_akima_cubic"]
gmt.mean.difs.cmips.vs.cmips <- gmt.futures.cmips - gmt.hists.cmips

load("results/hellinger_coefficient/hellinger_coefficients.rda")
load("results/hellinger_coefficient/hellinger_coefficients.rcp85.rda")
load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_hist_vs_fut.rda")

1-(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.hists),names(clim.futures)]
bd.cmips.fut.vs.hist <- diag(-log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.hists),names(clim.futures)])
bd.cmips.fut.vs.hist.wrt.int <- ((-log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.futures),"interim_10d_akima_cubic"])-(-log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.hists),"interim_10d_akima_cubic"]))
bd.cmips.fut.vs.hist.wrt.ncep<- (-log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.futures),"ncep_10d"])-(-log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.hists),"ncep_10d"])
bd.cmips.fut.vs.hist.wrt.jra55 <- (-log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.futures),"JRA55_10d_akima_cubic"])-(-log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.hists),"JRA55_10d_akima_cubic"])

bd.cmips.fut.vs.ncep <- -log(hellinger_coefficients_CMIP5_hist_vs_fut)[names(clim.futures),"ncep_10d"]
bd.cmips.hist.vs.int <- -log(hellinger_coefficients)[names(clim.hists),"interim_10d_akima_cubic"]
bd.cmips.hist.vs.ncep <- -log(hellinger_coefficients)[names(clim.hists),"ncep_10d"]
bd.cmips.hist.vs.JRA55 <- -log(hellinger_coefficients)[names(clim.hists),"JRA55_10d_akima_cubic"]
x <- "CMIP5_inmcm4_r1i1p1"
bd.cmips.hist.vs.x <- -log(hellinger_coefficients)[names(clim.hists),x ]


plot(bd.cmips.hist.vs.JRA55,gmt.mean.difs.cmips.vs.cmips)
cor(bd.cmips.hist.vs.int,gmt.mean.difs.cmips.vs.cmips, method = "pearson")

# change in Global Mean Temp versus change in spatial dependency structure CMIP5 BNs
plotname<- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/deltaGMT_vs_Bdist_CMIP5.pdf"
pdf(plotname,height = 5, width = 10)
plot(bd.cmips.fut.vs.hist,gmt.mean.difs.cmips.vs.cmips, 
     ylab = bquote(Delta ~ "Global Mean Temp (fut-hist)"),
     xlab = expression("d"[B]*"(BN"[fut]*",BN"[hist]*")"),
     main = "change in Global Mean Temp versus change in spatial dependency structure CMIP5 BNs \n(= Bhattacharya distance between future and historical BN)"
     )
text(40,3, paste("cor = ",round(cor(bd.cmips.fut.vs.hist,gmt.mean.difs.cmips.vs.cmips, method = "spearman"),2)))
lm1 <- lm(gmt.mean.difs.cmips.vs.cmips~bd.cmips.fut.vs.hist,data.frame(bd.cmips.fut.vs.hist,gmt.mean.difs.cmips.vs.cmips))
lines(x=c(18,36) , y = c(lm1$coefficients[2]*18+lm1$coefficients[1],lm1$coefficients[2]*36+lm1$coefficients[1]),col = 'red')
text(40,2.9, paste("R^2 = ",round(summary(lm1)$r.squared,3)))
dev.off()
###############################################################################################
# change in dependency pattern versus historical difference between cmip model and reanalysis
# Igual plot (F1 score vs delta precipitation in Causal networks for climate model evaluation..
# Nowack et al. )
###############################################################################################
plot(bd.cmips.fut.vs.hist,bd.cmips.hist.vs.int)
plot(gmt.mean.difs.cmips.vs.int,bd.cmips.hist.vs.int)
plot(gmt.mean.difs.cmips.vs.int,bd.cmips.fut.vs.hist)
plot(gmt.mean.difs.cmips.vs.cmips,bd.cmips.hist.vs.int)

df1 <- data.frame("gmt.change" = gmt.mean.difs.cmips.vs.cmips,"Bscore" = bd.cmips.hist.vs.int)
gmt.change2 <- gmt.mean.difs.cmips.vs.cmips^2
quadratic.model <-lm(Bscore ~ gmt.change +gmt.change2,df1)
summary(quadratic.model)
degreevalues <- seq(2, 5, 0.1)
predicteddegrees <- predict(quadratic.model,list(gmt.change=degreevalues, gmt.change2= degreevalues^2))

plot(gmt.mean.difs.cmips.vs.cmips,bd.cmips.hist.vs.int, pch=16, xlab = "GMT change °C", ylab = expression("Bscore: d"[B]*"(BN"[CMIP]*",BN"[ERAint]*")"), cex.lab = 1.3, col = "blue")
lines(degreevalues, predicteddegrees, col = "darkgreen", lwd = 3)


lin.model <-lm(Bscore ~ gmt.change,df1)
summary(lin.model)
degreevalues <- seq(2, 5, 0.1)
predicteddegrees <- predict(lin.model,list(gmt.change=degreevalues))

plot(gmt.mean.difs.cmips.vs.cmips,bd.cmips.hist.vs.int, pch=16, xlab = "GMT change °C", ylab = expression("Bscore: d"[B]*"(BN"[CMIP]*",BN"[ERAint]*")"), cex.lab = 1.3, col = "blue")
lines(degreevalues, predicteddegrees, col = "darkgreen", lwd = 3)
text(gmt.mean.difs.cmips.vs.cmips,bd.cmips.hist.vs.int, labels = names(gmt.mean.difs.cmips.vs.cmips))
######################################################
# Check temperature increase, weighted model ensemble
######################################################
sD <- 50
sS <- 50
modelweights <- IPweight(-log(hellinger_coefficients)[c(names(clim.hists),"interim_10d_akima_cubic"),c(names(clim.hists),"interim_10d_akima_cubic")],"interim_10d_akima_cubic",sD = sD,sS =sS)

plotname<- paste0("figs/weights/weights_CMIP5_sD",sD,"_sS",sS,".pdf")
pdf(plotname,height = 5, width = 10)
par(mar = c(9, 4, 1, 2))
plot(modelweights,xaxt =  "n",xlab = "", pch = 16)
xlabels=gsub("_r1i1p1","",gsub("CMIP5_","",gsub("historical","",names(modelweights))))
axis(1, at=1:length(modelweights),labels = xlabels, col.axis="red", las=2)
abline(h = 1/length(modelweights), lty = 2, col = "grey")
text(25,0.15, bquote(sigma~"D = "));text(27.5,0.15,paste0(sD))
text(25,0.11, bquote(sigma~"S = "));text(27.5,0.11,paste0(sS))
dev.off()

mean.hist.weighted <- sum(modelweights*gmt.hists.cmips)
mean.hist.unweighted <- sum(gmt.hists.cmips)/length(gmt.hists.cmips)
mean.fut.weighted <- sum(modelweights*keltoC(gmt.futures.cmips))
mean.fut.unweighted <- sum(keltoC(gmt.futures.cmips))/length(gmt.futures.cmips)
mean.diff.weighted <- sum(modelweights*gmt.diffs.cmips2)
mean.diff.unweighted <- sum(gmt.diffs.cmips2)/length(gmt.diffs.cmips2)



sd.hist.weighted <- sqrt(sum(modelweights*(gmt.hists.cmips-mean.hist.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
sd.hist.unweighted <- sqrt(sum((gmt.hists.cmips-mean.hist.unweighted)^2)/((length(gmt.hists.cmips)-1)))
sd.fut.weighted <- sqrt(sum(modelweights*(keltoC(gmt.futures.cmips)-mean.fut.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
sd.fut.unweighted <- sqrt(sum((keltoC(gmt.futures.cmips)-mean.fut.unweighted)^2)/((length(gmt.futures.cmips)-1)))
sd.diff.weighted <- sqrt(sum(modelweights*(keltoC(gmt.diffs.cmips2)-mean.diff.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
sd.diff.unweighted <- sqrt(sum((gmt.diffs.cmips2-mean.diff.unweighted)^2)/((length(gmt.diffs.cmips2)-1)))

# p.fh <- cor(gmt.hists.cmips,gmt.futures.cmips)
# sqrt(modelweights^2*(sd.fut.unweighted^2+sd.hist.unweighted^2+2*p.fh))

# sd.diff.unweighted <- sqrt(sd.fut.unweighted^2+sd.hist.unweighted^2+2*p.fh)
# sd.diff.weighted <- sqrt(sd.fut.weighted^2+sd.hist.weighted^2)

# sqrt(sum(modelweights^2*(sd.fut.unweighted^2+sd.hist.unweighted^2)))


mean.hist.weighted
mean.hist.unweighted
mean.fut.weighted
mean.fut.unweighted

# change in Global Mean Temp versus change in spatial dependency structure CMIP5 BNs
percs.w <- qnorm(c(0.05,0.33,0.5,0.66,0.9),mean = mean.diff.weighted,sd = sd.diff.weighted)
percs.uw <- qnorm(c(0.05,0.33,0.5,0.66,0.9),mean = mean.diff.unweighted,sd = sd.diff.unweighted)

# plotname<- paste0("figs/weights/deltaGMT_weighted_vs_unweighted_CMIP5_sD",sD,"_sS",sS,".pdf")
# pdf(plotname,height = 5, width = 10)
# par(mar = c(9, 4, 1, 2))
# plot(gmt.mean.difs.cmips.vs.cmips,xaxt =  "n",xlab = "",pch = 16,
#      ylab = bquote(Delta ~ "Global Mean Temp (fut-hist) °C"), main = "Weighted versus unweighted projection")
# xlabels=gsub("_r1i1p1","",gsub("CMIP5_","",gsub("historical","",names(modelweights))))
# axis(1, at=1:length(modelweights),labels = xlabels, col.axis="red", las=2)
# abline(h = mean.diff.weighted, col = "blue",lty = 1)
# abline(h = percs.w[2], col = "blue",lty = 2)
# abline(h = percs.w[4], col = "blue",lty = 2)
# abline(h = mean.diff.unweighted, col = "grey",lty = 1)
# abline(h = percs.uw[2], col = "grey",lty = 2)
# abline(h = percs.uw[4], col = "grey",lty = 2)
# 
# text(1,4.7, bquote(sigma~"D = "));text(5,4.7,paste0(sD))
# text(1,4.5, bquote(sigma~"S = "));text(5,4.5,paste0(sS))
# 
# legend("bottomleft",legend = c("weighted","unweighted"), col = c("blue","grey"),lty = 2)
# 
# dev.off()

gmt.reans <- sapply(clim.reans,gmt)

#########################################################################
# Plot weights and temperature rise in one figure
#########################################################################

plotname<- paste0("figs/weights/deltaGMT_and_weights_CMIP5_sD",sD,"_sS",sS,".pdf")
pdf(plotname,height = 5, width = 10)
par(mar = c(9, 4, 2, 2))
plot(gmt.mean.difs.cmips.vs.cmips,xaxt =  "n",xlab = "",pch = 16,
     ylab = bquote(Delta ~ "Global Mean Temp (fut-hist) °C"), main = "Weighted versus unweighted projection")
xlabels=gsub("_r1i1p1","",gsub("CMIP5_","",gsub("historical","",names(modelweights))))
axis(1, at=1:length(modelweights),labels = xlabels, col.axis="red", las=2)
abline(h = mean.diff.weighted, col = "blue",lty = 1)
abline(h = percs.w[2], col = "blue",lty = 2)
abline(h = percs.w[4], col = "blue",lty = 2)
abline(h = mean.diff.unweighted, col = "grey",lty = 1)
abline(h = percs.uw[2], col = "grey",lty = 2)
abline(h = percs.uw[4], col = "grey",lty = 2)

text(1,4.7, bquote(sigma~"D = "));text(2,4.7,paste0(sD))
text(1,4.5, bquote(sigma~"S = "));text(2,4.5,paste0(sS))

par(new = T)
plot(modelweights,
     axes = F,
     xlab = NA,
     ylab = NA, col = "brown" ,pch = 16)

abline(h = 1/length(modelweights), lty = 2, col = "brown")
axis(side = 4, col.axis = "brown", col = "brown")
mtext(side = 4, line = 3, 'Weights',col = "brown")

legend("bottomleft",legend = c("weighted","unweighted"), col = c("blue","grey"),lty = 2)

dev.off()


######################################################################################
# Check if Bhatachary distance is diagnosis for 
# climatology and standard deviation
######################################################################################
clim.sd.futures
clim.sd.reans
clim.sd.hists

x <- clim.hists$CMIP5_ACCESS1.0_r1i1p1
m.grid <- function(x){
  # x <- sets.hists$CMIP5_ACCESS1.0_r1i1p1
  v <- as.vector(cos(x$xyCoords$y/(180)*pi))
  m <- x$Data
  for(i in 1:dim(x$Data)[1]){
    m[i,] <- x$Data[i,]*v[i]}
  norm.m  <- m/(sum(v)*dim(x$Data)[2])
  x$Data <- norm.m
  return(x)
}


RMSE.sd.grid.int <- numeric(length(clim.sd.hists))
for (i in 1:length(clim.sd.hists)){
   x1 <- m.grid(clim.sd.hists[[i]])$Data
   x2 <- m.grid(clim.sd.reans$interim_10d_akima_cubic)$Data
   # x1 <- clim.sd.hists[[i]]$Data
   # x2 <- clim.sd.reans$interim_10d_akima_cubic$Data
   RMSE.sd.grid.int[i] <- sqrt(sum((x1-x2)^2)/length(x1))
}
names(RMSE.sd.grid.int) <- names(clim.sd.hists)

RMSE.clim.int <- numeric(length(clim.hists))
i <- 1

clim.interim.c <- clim.reans$interim_10d_akima_cubic
clim.interim.c$Data <- clim.reans$interim_10d_akima_cubic$Data-273.15

for (i in 1:length(clim.hists)){
  x1 <- m.grid(clim.hists[[i]])$Data
  x2 <- m.grid(clim.interim.c)$Data
  # x1 <- clim.sd.hists[[i]]$Data
  # x2 <- clim.sd.reans$interim_10d_akima_cubic$Data
  RMSE.clim.int[i] <- sqrt(sum((x1-x2)^2)/length(x1))
}

names(RMSE.clim.int)<- names(clim.hists)

# save(RMSE.int.w,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/RMSE.int.w.rda")
bd.cmips.hist.vs.int
plot(bd.cmips.hist.vs.int,RMSE.clim.int)
cor(bd.cmips.hist.vs.int,RMSE.clim.int)
text(bd.cmips.hist.vs.int,RMSE.clim.int,labels = names(RMSE.clim.int))
plot(bd.cmips.hist.vs.int,RMSE.sd.grid.int)
cor(bd.cmips.hist.vs.int,RMSE.sd.grid.int)
text(bd.cmips.hist.vs.int,RMSE.sd.grid.int,labels = names(RMSE.sd.grid.int))

########################################
# Library TreeSeg
########################################
names(sets.hists)
rownames(selecteddists)
colnames(selecteddists)

hist.dists <-selecteddists[names(sets.hists),names(sets.hists)]
fut.dists <- selecteddists[names(sets.futures),names(sets.futures)]

# For RMSE clim int
obs <- RMSE.clim.int
names(obs)

# For Observations are differences future-hist mean temperatures
obs <- gmt.diffs.cmips2
length(obs)
names(obs) <- rownames(hist.dists)

# For observations are sd hist mean tempeartures
obs <- gsdt.hists.cmips
length(obs)
names(obs)<- rownames(hist.dists)

# For Observations are historical mean temperatures
obs <- gmt.hists.cmips
names(obs)<- rownames(hist.dists)
all.equal(names(gmt.hists.cmips),rownames(hist.dists))

# For Observations are future mean temperatures
obs <-gmt.futures.cmips[names(sets.futures)]
names(obs)<- rownames(hist.dists)
names(obs)<- rownames(fut.dists)

###############################################################
# Hierarchical clustering batachari
###############################################################
library(proxy)
library(devtools)
library(treeSeg)
library(phytools)
library("vegan")
library(ape)
library('TreeTools')
##############################################################
# Making dendogram and visualize with hclust
##############################################################
# For dendogram made of historical BNs
mclust <-hclust(dist(-log(hist.dists),method = "minkowski", p = 2),method = "complete")
mclust <-hclust(dist(-log(hist.dists),method = "minkowski", p = 2),method = "single")
mclust <-hclust(dist(-log(hist.dists),method = "minkowski", p = 2),method = "ward.D")

# For dendogram made of Future BNs
mclust <-hclust(dist(-log(fut.dists),method = "minkowski", p = 2),method = "complete")
mclust <-hclust(dist(-log(fut.dists),method = "minkowski", p = 2),method = "single")
mclust <-hclust(dist(-log(fut.dists),method = "minkowski", p = 2),method = "ward.D")

# Visualize

plot(mclust, hang = -1)

######################################################
# Converting to phylo object
######################################################
phy <-as.phylo(mclust)
# attr(phy,"order")
#attributes(phy)
#phy$tip.label
#phy$edge
#mclust$merge
#mclust$order
# Get orders tips in plot
phy$tip.label[mclust$order]
#################################
# Treeseg: Combine observations tips with treestructure
#################################
# y <- obs
# tree <- phy
# tipOrder <- "cladewise"
# checkOrder <- TRUE

# #phy.re <-reorder.phylo(phy, order = "postorder",index.only =TRUE)
# #phy.o <-phy.re[(length(phy.re)-29):length(phy.re)]
# new.phy <-phy
# new.phy$tip.label<-phy$tip.label[mclust$order]
# new.phy$tip.label

# nodepath(phy)

plot.phylo(phy, direction = 'downwards')
nodelabels(cex =1, frame = "none")
tiplabels(cex=1, frame ="none")

#new.tS <-treeSeg(obs,new.phy,checkOrder = FALSE, fam = "gauss", alpha = 0.05)
tS<-treeSeg(obs,phy,checkOrder = FALSE, fam = "gauss", alpha = 0.1)

node.cols <-rainbow(length(tS$mlAN))
node.vec <-rep("black",1,length(phy$tip.label)+phy$Nnode)
tip.vec <- rep("black",1,length(phy$tip.label))

for(i in 1:length(tS$mlAN)){
  desc <-getDescendants(phy,tS$mlAN[i])
  node.vec[desc]<-node.cols[i]
  tip.vec[which(desc <= length(phy$tip.label))]<-node.cols[i]
}

tip.vec
par(mar = c(0,0,0,0))
plot.phylo(phy,node.color = node.vec, direction = 'downwards',no.margin =TRUE, cex = 0.7)
nodelabels(cex =0.8, adj = c(0,0),frame = "none")
tiplabels(cex=0.8, adj = c(0,0),frame ="none")

###############################
#
##############################
tiporder <-phy$tip.label[mclust$order]
par(mar = c(1,4, 13, 2))

plot(obs[tiporder],xlab="", ylab = "Mean Temp in 2071-2100 °C", xaxt ="n", col =node.vec[1:length(phy$tip.label)][mclust$order],pch = 16)
xlabels=tiporder
axis(3, at=1:length(tiporder),labels = xlabels, col.axis= "red", las=2)



