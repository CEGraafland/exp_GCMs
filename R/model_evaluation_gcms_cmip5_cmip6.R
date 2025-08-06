#############################################################################
# Model Evaluation CMIP6
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
########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("data/tas_JRA55_10d_akima_cubic.rda")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")
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
# cubicsets <- c(cubicsets5,cubicsets5extra,cubicsetsearth,cubicsets6,listinterim,listjra55)
# With CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,cubicsets6,listinterim,listjra55,listncep)
names(cubicsets)
# Without CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,listinterim,listjra55,listncep)
names(cubicsets)


namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
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
loadIterations <- function(pattern,permused, hctype = NULL, ncep = FALSE){
  if(is.null(hctype)){hctype <- ""}
  if(ncep == FALSE){fullpattern <- paste0("data/hciterations/",pattern)} else {fullpattern <- pattern}
  hc_interim_list <- list.files(paste0(fullpattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0(fullpattern,"/perm",permused,hctype))
  hc_interim_names <- gsub(".rda", "", hc_interim_names)
  
  hc_interim_networks <- list()
  
  hc_interim_networks <- lapply(hc_interim_list, function(x){get(load(x))})
  names(hc_interim_networks) <- hc_interim_names
  interimsizes <- sapply(hc_interim_networks,narcs)
  hc_interims <- hc_interim_networks[order(interimsizes)]
  return(hc_interims)
}

##################################################################
# Make al datasets standardized
##################################################################
permk <- 3
# Anomalies and standarization over seasons (networks are learned with this set)
data_gcms_out <- lapply(cubicsets,function(x) TimeCoordsAnom_from_Grid_rms(x, rms = TRUE))
data_gcms_out_df <- lapply(data_gcms_out, function(x) as.data.frame(x))
data_gcms_anom_out <- lapply(data_gcms_out, function(x) mat2Dto3Darray(x,attr(x,"Xcoords"), attr(x,"Ycoords")))
grid_gcms_anom <- mapply(function(x,y) {x$Data <- y 
return(x)}, 
x = cubicsets, 
y = data_gcms_anom_out,SIMPLIFY = FALSE)
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
# data_gcms_anom_scaled_3 <- lapply(data_gcms_anom_scaled, function(x) x[,permutations[[3]]])
assign(paste0("data_gcms_anom_scaled_",permk),lapply(data_gcms_anom_scaled, function(x) x[,permutations[[permk]]]))



hc_gcms5 <- lapply(paste0(namescubics5), loadIterations, permused = permk)
hc_gcms5extra <- lapply(paste0("CMIP5_extra/",namescubics5extra), loadIterations, permused = permk)
hc_gcms5left<- lapply(paste0("CMIP5_left/",namescubics5left), loadIterations, permused = permk)
hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk)
hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk) 
hc_jra55 <- lapply(c("JRA55_10d_akima_cubic"), loadIterations, permused = permk) 
hc_ncep <- lapply(c("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE) 
# With cmip6
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_gcms6,hc_interim,hc_jra55,hc_ncep)
# Without cmip6
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_interim,hc_jra55,hc_ncep)
length(hc_gcms)
# names(hc_gcms) <- shortnames
names(hc_gcms) <- shortnames2

# Choose between 'own optimums' or constant magnitude
# modelsize <- own_optimums
modelsize <- 18
selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
# selection_fits_gcms_scaled_3 <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled_3, SIMPLIFY = FALSE)
# all.equal(selection_fits_gcms_scaled,selection_fits_gcms_scaled_3)
# selection_fits_gcms_scaled_3$CMIP5_CanESM2$V550
hc_jra55
hc_interim
##########################################################################
# Global dataset evaluation from BN 
##########################################################################
whichBNFits <-selection_fits_gcms_scaled
whichData <- data_gcms_anom_scaled

#########################################################################
# pickedEOFsgrids is prepared under EOFs logliks preparation. 
#####################################################################
data_pickedEOFs <- lapply(pickedEOFsgrids,
                                function(x) { x <- redim(x,drop = TRUE)
                                xdata <- array3Dto2Dmat(x$Data) 
                                return(xdata)}) 
data_pickedEOFs <- lapply(data_pickedEOFs, function(x) as.data.frame(x))
whichData <- c(data_gcms_anom_scaled,data_pickedEOFs)

ordersets <- TRUE

loglik_selection_datasets <- matrix(data = NA, nrow = length(whichBNFits), ncol = length(whichData), dimnames = list(names(whichBNFits),names(whichData)))
for(i in 1:length(whichBNFits)){
  lo <- sapply(X = whichData, logLik, object = whichBNFits[[i]])
  loglik_selection_datasets[i,] <- lo
}

loglik_selection_datasets[,"interim_10d_akima_cubic"]

institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor")
institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","ACCESS","bcc","BNU","CCSM4","CESM1","CMCC","CSIRO","inmcm4","MRI")

# with cmip6
institutions <- rep(institutions1,each =2)
CMIPS <- rep(c("CMIP5","CMIP6"),length(institutions1))
# without cmip6
institutions <- rep(institutions1,each =1)
CMIPS <- rep(c("CMIP5"),length(institutions1))
combinations <- cbind(CMIPS,institutions)
nrow(loglik_selection_datasets)
namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(loglik_selection_datasets))))
length(namesort)
# For variance percentage interim:
namesort1 <- c(namesort,36)
namesort2 <- c(namesort,36,37,38,39,40,41,42)
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

orderby <- "interim_10d_akima_cubic"
orderby <- "CMIP5_HadGEM2.ES_r1i1p1"
orderby <- "ncep_10d"
logsort <- sort(loglik_selection_datasets[,orderby],index.return = TRUE)
logsort
rep_cnrm <- c(1,2,7,3,4,5,6,8,9,10)


loglik_selection_datasets_ordered <- loglik_selection_datasets[namesort1,namesort2]
loglik_selection_datasets_ordered <- loglik_selection_datasets[logsort$ix,logsort$ix]

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
  geom_text(aes(label = round(value/10^4,0))) +
  scale_fill_gradient2(name = "value x 10^4",
                      # low="white", 
                      # high="red",
                      limits = c(bottom,top), labels = function(x)x/10^4) +
  labs(x="data", y="models", title=paste0(modelsize*100," optimum model evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

b             

dev.off()
##########################################################################
# determinant valuation
##########################################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# BiocManager::install("graph")
fits_NEL <- lapply(selection_fits_gcms_scaled,bnlearn:::as.graphNEL)
# fits_NEL$CMIP5_CanESM2@edgeData
# fits_NEL$CMIP5_CanESM2@nodes
# fits_NEL$CMIP5_CanESM2@edgeL

# selection_fits_gcms_scaled$CMIP5_CanESM2$'V619'

# bn_NEL <- lapply(fits_NEL,as.bn)
# bnNEL_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = bn_NEL, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
# all.equal(selection_fits_gcms_scaled,bnNEL_fits_gcms_scaled)
edgelistssparse <- lapply(fits_NEL,edgeList)

# print(edgelistssparse$CMIP5_CanESM2[[1]])
# selection_fits_gcms_scaled$CMIP5_CanESM$V303
# (1:648)[permutations[[3]]][edgelistssparse$CMIP5_CanESM2[[4]]]

# as.matrix(get.adjacency.matrix(edgelistssparse$interim_10d_akima_cubic))
# blub <- lapply(data_gcms_anom_scaled, function(x) x[permutations[[3]]])
data_gcms_anom_scaled_3 <- lapply(data_gcms_anom_scaled, function(x) x[,permutations[[3]]])
# all.equal(blub, data_gcms_anom_scaled_3)
# data_gcms_anom_scaled_3
data_gcms_anom_scaled_sparse <- lapply(data_gcms_anom_scaled_3,sparsebnData, type = "continuous")



fits_COVS <- mapply(get.covariance, x = edgelistssparse, data = data_gcms_anom_scaled_sparse)

# fits_PRECS<- mapply(get.precision, x = edgelistssparse, data = data_gcms_anom_scaled_sparse)
# fits_COVS$CMIP5_CanESM2
fits_COVmats <- lapply(fits_COVS, as.matrix)
fits_invmats <- lapply(fits_COVmats, solve)
det_COVmats <- lapply(fits_COVmats, determinant, logarithm = TRUE)
unlist(det_COVmats)

sapply(det_COVmats, function(x)x)
# det_COVmats2 <- sapply(fits_COVmats,det)
# log(det_COVmats2)
# cor_COVmats <- lapply(fits_COVmats,cov2cor)
# det_CORmats <- sapply(cor_COVmats,det)
# log(det_CORmats)
# all.equal(cor_COVmats,fits_COVmats)

logLik_const <- sapply(det_COVmats,function(y) {-0.5*(y$modulus+648*log(2*pi))})
logLik_const*360
# lapply(fits_COVmats,function(x) {-0.5*(log(det(x))+648*log(2*pi))})

fits_Means <- lapply(data_gcms_anom_scaled, function(x) colMeans(x[permutations[[3]]]) )

logliksSP <- matrix(nrow = length(data_gcms_anom_scaled_3), ncol =  length(data_gcms_anom_scaled_3), dimnames = list(names(data_gcms_anom_scaled_3),names(data_gcms_anom_scaled_3)))
nedgesSP<- c()
library(mvtnorm)
for (j in 1:length(data_gcms_anom_scaled)){
  for (i in 1:length(fits_NEL)){
    littles <- dmvnorm(data_gcms_anom_scaled[[j]][permutations[[3]]], mean = fits_Means[[i]], sigma = fits_COVmats[[i]], log = TRUE)
    logliksSP[j,i] <- sum(littles)
  }
}


# SINGLE CASE 
# # for (j in 1:length(data_gcms_anom_scaled)){
# # for (i in 1:length(fits_NEL)){
#     testvec <- data_gcms_anom_scaled[[10]][permutations[[3]]][1,]-fits_Means[[10]]
#     prod1 <- fits_invmats[[10]]%*%t(testvec)
#     dev1 <- -(1/2)*(as.matrix(testvec) %*% as.matrix(prod1))
#     dim(as.matrix(testvec))
#     dim(as.matrix(prod1))
# 
#     logLik_const[[10]]

gcms_event_devs <- array(data = NA, 
                         dim = c(nrow(data_gcms_anom_scaled[[length(fits_invmats)]]),length(fits_invmats)),
                         dimnames = list(NULL,names(data_gcms_anom_scaled)))
for (j in 1:ncol(gcms_event_devs)){
  testvec2 <- data_gcms_anom_scaled[[length(fits_invmats)]][permutations[[3]]]-fits_Means[[j]]
  events_devs <- c()
  for (i in 1:nrow(testvec2)){
    prod2 <- fits_invmats[[j]]%*%t(testvec2[i,])
    dev2 <- -(1/2)*(as.matrix(testvec2[i,]) %*% as.matrix(prod2))
    events_devs[i]<- dev2
  }
  gcms_event_devs[,j] <- events_devs
}

# CHECK 
colSums(gcms_event_devs)+360*logLik_const
logliksSP[36,]

sum(events_devs)+360*logLik_const[[ncol(gcms_event_devs)]]

colnames(gcms_event_devs) <- abrev
# all.equal(testvec,testvec2[1,])
# all.equal(prod1,as.matrix(prod2[,1]))
# dim(prod2[,1])

# logliks_event <- array(data = NA, dim = c(nrow(respectiveDataSet),length(whichBNFits)))
# 
# for (i in 1:length(whichBNFits)){
#   logliks_event[,i] <- sapply(single_dfs, FUN = logLik, object = whichBNFits[[i]])
# }

# Total data:
gcms_event_devs
logLik_const

cmat <- matrix(rep(logLik_const, each = nrow(gcms_event_devs)), ncol = length(logLik_const),nrow = nrow(gcms_event_devs))
gcms_event_method2 <- gcms_event_devs +cmat

################################################################
# globaldata evaluation from SIGMA method
###############################################################
longSPdata <- melt(logliksSP)
# Plot normal value 
bottom <- -300000
top <- max(longdata$value)

#min(longdata$value)
a <- ggplot(longSPdata, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  geom_text(aes(label = round(value/10^4,0))) +
  scale_fill_gradient(name = "value x 10^4",low="white", high="red",
                      limits = c(bottom,top), 
                      labels = function(x)x/10^4) +
  labs(x="data", y="models", title=paste0(modelsize*100," model evaluation via Sigma representation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

a
################################################################
# globaldata evaulation from SIGMA method Deviation part
################################################################
longDEVdata <- melt(gcms_event_devs)
# Plot event value 
# bottom <- -5000
bottom <- min(longDEVdata$value)/10
top <- max(longDEVdata$value)

#min(longdata$value)
ggplot(longDEVdata, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  # geom_text(aes(label = round(value/10^4,0))) +
  scale_fill_gradient(name = "Deviation",low="white", high="red",
                      limits = c(bottom,top)
                      # labels = function(x)x/10^4
  ) +
  labs(x="data", y="models", title=paste0(modelsize*100," optimum model evaluation")) +
  # scale_y_discrete(name="models", labels = shortnames) #+
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))


plotdata <- gcms_event_devs
brewer.pal.info
n <- 36
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

par(mar=c(6.1, 4.1, 4.1, 9.1), xpd=NA)
cols <- col_vector
cols <- brewer.pal(9,name = "Set1")
maxl <- max(plotdata)
maxl <- 100
minl <- min(plotdata)
minl <- -2000
minl <- minl / 5
plot(plotdata[,36],pch = 16,col = cols[1],
     ylim = c(minl,maxl), 
     xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = paste0("Loglik constants and deviations data realizations (wrt interim) METHOD2"))
for(i in 1:35){
  points(plotdata[,i], col = cols[i+1],pch = 16)
  abline(h = logLik_const[i],col = cols[i+1] )
}

legend(370,300, c(abrev[36],abrev[1:35]),bty = "n", pch = 16, col =cols[1:36], cex = 0.6)

##################################################################
# Event logliks preparation
##################################################################
whichBNFits <- selection_fits_gcms_scaled
whichData <- data_gcms_anom_scaled
data_gcms_anom_scaled$interim_10d_akima_cubic
names(cubicsets)
whichSet <- "interim_10d_akima_cubic"
whichSet <- "ncep_10d"
# whichSet <- "CMIP5_EC.EARTH_r12i1p1"

# # EOFS preparation
# whichPC  <- "PC1"
# prinCompSet <- prinComp(grid_gcms_anom_scaled[[whichSet]], n.eofs = 360)
# EOFs <- prinCompSet[[1]][[1]]$EOFs
# PCs <- prinCompSet[[1]][[1]]$PCs
# selected_EOF <- EOFs[,whichPC]
# selected_PC <- PCs[,whichPC]
# projection_EOF <- matrix(nrow = length(selected_PC),ncol = length(selected_EOF))
# for(i in 1:length(selected_PC)){
#   projection_EOF[i,] <- selected_EOF*selected_PC[i]
# }

# EOFS preparation

# whichPCs  <- c("PC1","PC2","PC3","PC4","PC5")
# whichPCs <- 1:10
grid_gcms_anom_scaled$CMIP5_CanESM2_r1i1p1$Data
gridforEOF <- grid_gcms_anom_scaled
prinCompSet <- prinComp(gridforEOF[[whichSet]], n.eofs = 360)
expvars <- attr(prinCompSet[[1]][[1]],"explained_variance")
EOFs <- prinCompSet[[1]][[1]]$EOFs
PCs <- prinCompSet[[1]][[1]]$PCs

varpercentages <- c(0.25,0.50,0.75,0.90,0.95,0.99)
pickedEOFs <- lapply(varpercentages, FUN = function(x) which(expvars > x)[[1]])
pickedEOFsdata <- lapply(pickedEOFs, FUN = function(x) EOFS_preparation(1:x))

templatedata <- gridforEOF[[whichSet]]
templatesdata <- rep(list(templatedata),6)
templatesdata[[3]]
templatedata$Data[1,,,] <- pickedEOFsdata[[1]]00
pickedEOFsgrids <- mapply(function(x,y) {z<- x; z$Data[1,,,] <- y; return(z)},templatesdata, pickedEOFsdata, SIMPLIFY  = FALSE)
templatenames <- paste0("interim_10d_akima_cubic_v.exp_",as.character(varpercentages))
names(pickedEOFsgrids) <- names(pickedEOFsdata) <- templatenames
pickedEOFsgrids$interim_10d_akima_cubic_v.exp_0.25

# now equalize first new longitude with last already existing longitude
newdata.historical3[,,,,1]<-data.historical3[,,,,144]
# give same dimensions
attr(newdata.historical3,"dimensions") <-  attr(historical3$Data,"dimensions") 
# replace data 
historical3$Data <- newdata.historical3


prinCompset_byvar <- prinComp(gridforEOF[[whichSet]], n.eofs = varpercentages)
prinCompSet_byvar <- apply(varpercentages, MARGIN = 1, FUN = function(x)prinComp(gridforEOF[[whichSet]],v.exp = x))
prinCompSet_byvar

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

eventLoglikCalc <- function(respectiveDataSet, manipulateWith){
  
  single_dfs <- apply(respectiveDataSet, MARGIN = 1, FUN = function(x) as.data.frame(x))
  single_dfs  <- lapply(single_dfs , t)
  single_dfs  <- lapply(single_dfs , as.data.frame)
  
  # Manipulate data if necesarry
  if(!is.null(manipulateWith)){
    for(i in 1:length(single_dfs)){
      single_dfs[[i]] <- single_dfs[[i]]*manipulateWith[[i]]
    }
  }
  
  logliks_event <- array(data = NA, dim = c(nrow(respectiveDataSet),length(whichBNFits)))
  
  for (i in 1:length(whichBNFits)){
    logliks_event[,i] <- sapply(single_dfs, FUN = logLik, object = whichBNFits[[i]])
  }
  
  logliks_event_scaled <- scale(logliks_event, center = TRUE, scale = TRUE)
  colnames(logliks_event_scaled) <- abrev
  colnames(logliks_event) <- abrev
  
  return(list("logliks_event" = logliks_event,"logliks_event_scaled" = logliks_event_scaled))
}

#################################################################################
# Visualize models vs events
#################################################################################
melteddfn <- melt(gcms_event_devs)
# melteddfn <- melt(logliks_event)
meltedfn <- melt(Logliks_event_normal$logliks_event_scaled)
Logliks_event_normal$logliks_event_scaled
titlename <- ""
limits = NA

plot_loglik_event_raster <- function(meltedfn,titlename, limits = c(min(meltedfn$value),max(meltedfn$value))){
  slde <- ggplot(meltedfn, aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient2(high="blue", low="red", midpoint =0, 
                         limits = limits, 
                         na.value = "black") +
    labs(x= paste0("events"), y="models", title = titlename) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
          axis.text.y=element_text(size=9),
          plot.title=element_text(size=11))
  return(slde)
}

min(meltedfn$value)
max(meltedfn$value)
dev.off()
###############################################################################
# Visualize correlation matrix loglikelhood data events models
###############################################################################
# loglikscor <- cor(logliks_event)
# loglikscor <- cor(logliks_event_scaled)
# meltloglikscor <- melt(loglikscor)

plot_loglikscor <- function(meltloglikscor,titlename, limits = c(0,1), decimals = 3){
  if (is.null(titlename)){titlename <- paste0("GCMs ", modelsize*100," loglikelihood correlation w.r.t ",whichSet)}
  ggplot(meltloglikscor, aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill=value)) + 
    geom_text(aes(label = round(value,decimals))) +
    scale_fill_gradient(high="blue", low="white", 
                        #midpoint =0, 
                        limits = limits, 
                        na.value = "black") +
    labs(x= paste0("models"), y="models", title=titlename) +
    theme_bw() + 
    theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
          axis.text.y=element_text(size=9),
          plot.title=element_text(size=11))
}

plot_loglikscor(meltloglikscor,titlename = "hola")
################################################################################
# Visualize Events one by one
################################################################################


plot_loglik_event_points <- function(logliks_event_matrix, loglik_constants = NULL, minl = NULL, maxl = NULL, titlename = NULL){
  if(is.null(titlename)){titlename <-  paste0("Loglik data realizations given models")}
  if(is.null(maxl)){maxl <- max(logliks_event_matrix)}
  if(is.null(minl)){minl <- min(logliks_event_matrix)}
  
  totnum <- ncol(logliks_event_matrix)
  
  plot(logliks_event_matrix[,totnum],pch = 16,col = rainbow(totnum)[1],
       ylim = c(minl,maxl), 
       xlab = "D_i", ylab = bquote("log P(D_i|"~Mu~")"), main = titlename)
  if (!is.null(loglik_constants)){abline(h = loglik_constants[totnum], col = rainbow(totnum)[1])}
  for(i in 1:(totnum-1)){
    points(logliks_event_matrix[,i], col = rainbow(totnum)[i+1],pch = 16)
    if (!is.null(loglik_constants)){abline(h = loglik_constants[i], col = rainbow(totnum)[i+1])}
  }
  legend("bottomleft", c(abrev[totnum],abrev[(1:(totnum-1))]),bty = "n", pch = 16, col =rainbow(totnum)[1:totnum], cex = 0.6) 
}


Logliks_event_normal <- eventLoglikCalc(respectiveDataSet = whichData[[whichSet]], manipulateWith =NULL)
Logliks_event_0.01 <- eventLoglikCalc(respectiveDataSet = whichData[[whichSet]], manipulateWith = rep(0.01,360))
Logliks_event_PC1 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(c(1)), manipulateWith = NULL)
Logliks_event_PC1 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(c(1)), manipulateWith = NULL)
Logliks_event_PC1_3 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:3), manipulateWith = NULL)
Logliks_event_PC1_10 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:10), manipulateWith = NULL)
Logliks_event_PC1_360 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:360), manipulateWith = NULL)

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMS/figs/logliks_CMIP_5_6_intercomparison/Loglik_models_variouseventsEOF.pdf"
pdf(plotname, width = 7, height = 7)
par(mfrow = c(3,2))


plot_loglik_event_points(Logliks_event_normal$logliks_event,loglik_constants = NULL, titlename = "Loglik data realizations", minl = -2000, maxl = 100)
plot_loglik_event_points(Logliks_event_0.01$logliks_event,loglik_constants = NULL, titlename = "Loglik data realizations factor 0.01")
plot_loglik_event_points(Logliks_event_PC1$logliks_event,loglik_constants = NULL, titlename = "Loglik EOF1 part of data ")
plot_loglik_event_points(Logliks_event_PC1_3$logliks_event,loglik_constants = NULL, titlename = "Loglik EOF1 + ... + EOF3 part of data ")
plot_loglik_event_points(Logliks_event_PC1_10$logliks_event,loglik_constants = NULL, titlename = "Loglik EOF1 + ... + EOF10 part of data ")
plot_loglik_event_points(Logliks_event_PC1_360$logliks_event,loglik_constants = NULL, titlename = "Loglik EOF1 + ... + EOF360 part of data ", minl = -2000, maxl = 100)
dev.off()

b1 <- plot_loglik_event_raster(melt(Logliks_event_normal$logliks_event_scaled), titlename = "Normal data")
b2 <- plot_loglik_event_raster(melt(Logliks_event_0.01$logliks_event_scaled), titlename = "Scaled 0.01")
b3 <- plot_loglik_event_raster(melt(Logliks_event_PC1$logliks_event_scaled), titlename = "PC1")
b4 <- plot_loglik_event_raster(melt(Logliks_event_PC1_3$logliks_event_scaled), titlename = "PC 1:3")
b5 <- plot_loglik_event_raster(melt(Logliks_event_PC1_10$logliks_event_scaled), titlename = "PC 1:10")
b6 <- plot_loglik_event_raster(melt(Logliks_event_PC1_360$logliks_event_scaled), titlename = "PC 1:360")
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/logliks_CMIP_5_6_intercomparison/Loglik_scaled_models_variouseventsEOF.pdf"

pdf(plotname, width = 20, height = 20)
joep1<- arrangeGrob(b1,b2,b3,b4,b5,b6,nrow = 3)
grid.arrange(joep1, newpage = FALSE)
dev.off()

if (ordersets == TRUE){o <- namesort} else{o <- 1:ncol(Logliks_event_normal$logliks_event)}
c1 <- plot_loglikscor(melt(cor(Logliks_event_normal$logliks_event)[o,o]),titlename = NULL, decimals = 2)
c2 <- plot_loglikscor(melt(cor(Logliks_event_0.01$logliks_event)[o,o]),titlename = "Normal scaled 0.01")
c3 <- plot_loglikscor(melt(cor(Logliks_event_PC1$logliks_event)[o,o]),titlename = "PC1")
c4 <- plot_loglikscor(melt(cor(Logliks_event_PC1_3$logliks_event)[o,o]),titlename = "PC 1:3", limits = c(0.75,1))
c5 <- plot_loglikscor(melt(cor(Logliks_event_PC1_10$logliks_event)[o,o]),titlename = "PC 1:10", limits = c(0.75,1))
c6 <- plot_loglikscor(melt(cor(Logliks_event_PC1_360$logliks_event)[o,o]),titlename = "PC 1:360")
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/logliks_CMIP_5_6_intercomparison/correlationLoglik_models_variouseventsEOF_cmip5_cmip6.pdf"
pdf(plotname, width = 30, height = 20)
joep2 <- arrangeGrob(c1,c4,c5,c6,nrow = 2)
grid.arrange(joep2, newpage = FALSE)
dev.off()

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/logliks_CMIP_5_6_intercomparison/correlationLoglik_models_interimEOF1_3_cmip5_cmip6.pdf"
pdf(plotname, width = 17, height = 12)
c4
dev.off()


####################################################################################
# Submatrices
####################################################################################
subm <-grep("CESM|CCSM|Nor",colnames(Logliks_event_normal$logliks_event))
subm <-grep("MPI|CMCC",colnames(Logliks_event_normal$logliks_event))
subm <-grep("CNRM|EC.EARTH",colnames(Logliks_event_normal$logliks_event))
subm <-grep("CCSM|CESM|bcc",colnames(Logliks_event_normal$logliks_event))
plot_loglikscor(melt(cor(Logliks_event_normal$logliks_event)[,subm]),titlename = NULL, decimals = 2)

#######################################################################################
# investigate rare modeled events
#######################################################################################
library('RColorBrewer')
library(gridExtra)
col.blue <- rev(brewer.pal(8,"Blues"))
col.blue
col.red <- brewer.pal(8,"Reds")
col.div <- c(col.blue, col.red)
col.div[8:9] <- "white"
col.b <- c(col.blue,col.red)

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/backpermutations.rda")
Logliks_event_PC1_3[,"CMIP6Amon_GFDL.ESM4f1_historical"]
Logliks_event_PC1_3$logliks_event_scaled
lowtohighGFDL1 <- order(Logliks_event_PC1_3$logliks_event_scaled[,"CMIP6Amon_GFDL.ESM4f1_historical"])
lowtohighGFDL2 <- order(Logliks_event_PC1_3$logliks_event_scaled[,"CMIP6Amon_EC.Earth3f1_historical"])
lowtohighINT <- order(Logliks_event_PC1_3$logliks_event_scaled[,"interim_10d"])
lowtohighGFDL1
lowtohighGFDL2
lowtohighINT
datasetpc110 <- EOFS_preparation(1:3)
lowev<- datasetpc110[lowtohighINT,][358,]
lowev <- datasetpc110[232,]
climlow <- quantity2clim(lowev,"lowevent",ref.grid = cubicsets$CMIP6Amon_GFDL.ESM4_r1i1p1f1_historical,backperm = NULL)


spatialPlot(climlow,backdrop.theme = "coastline",lonCenter = 180, col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,0.5))
#######################################################################################
# clustering of loglike events? euc plus ward
#######################################################################################
earthhist <- grep("historical",colnames(Logliks_event_normal$logliks_event_scaled))
reanalysis <- grep("_",colnames(Logliks_event_normal$logliks_event_scaled))

distmethod <- "euclidian"
exclude <- -c(earthhist)
d <- dist(cor(Logliks_event_normal$logliks_event_scaled[exclude, exclude]), method = distmethod)
clustmethod = "ward"
fit <- hclust(d, method = clustmethod)

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/modelclustering/",distmethod,"_",clustmethod,".pdf")
pdf(plotname)
plot(fit)
# rect.hclust(fit, k=8) 
dev.off()
