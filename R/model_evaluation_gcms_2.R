#############################################################################
# Model Evaluation CMIP5
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
citation("sparsebn")
########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"
cubicsets <- c(tas_historical_10d_akima_cubic_corrected,listinterim)

namescubics <- names(cubicsets)
shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
                   pattern = "_r12i1p1", replacement = "")
abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")

############################################################################
# Function to load hciterations of models in a list 
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
##############################################################################
# load results loglikelihood train and tests and 
# combine loglikelihood results of every single model
##############################################################################
loglik_traintest_gcms <- lapply(shortnames, function(x){get(load(paste0("results/loglik_traintest/loglik_traintest_",x,".rda")))})
names(loglik_traintest_gcms) <- shortnames
# then select the statistical optimums of the models
own_optimums <-sapply(loglik_traintest_gcms, function(x) which(x[,"testtrain"] == max(x[,"testtrain"])))
rm(loglik_traintest_gcms)
##################################################################
# Make al datasets standardized
##################################################################
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
data_gcms_anom_scaled_3 <- lapply(data_gcms_anom_scaled, function(x) x[,permutations[[3]]])



hc_gcms <- lapply(namescubics, loadIterations, permused = 3)
names(hc_gcms) <- shortnames
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
##########################################################################
# determinant valuation
##########################################################################
fits_NEL <- lapply(selection_fits_gcms_scaled,as.graphNEL)
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
# data_gcms_anom_scaled_sparse <- lapply(data_gcms_anom_scaled,sparsebnData, type = "continuous")
# data_gcms_anom_out_sparse$CMIP5_CanESM2_r1i1p1$data
# data_gcms_out_sparse <- lapply(data_gcms_out,sparsebnData, type = "continuous")
# fit_COV_int <- get.covariance(x = edgelistssparse$interim_10d_akima_cubic,data = data_gcms_anom_scaled_sparse$interim_10d_akima_cubic)
# all.equal(fit_COV_int,fits_COVS$interim_10d_akima_cubic)


fits_COVS <- mapply(get.covariance, x = edgelistssparse, data = data_gcms_anom_scaled_sparse)

# fits_PRECS<- mapply(get.precision, x = edgelistssparse, data = data_gcms_anom_scaled_sparse)
# fits_COVS$CMIP5_CanESM2
fits_COVmats <- lapply(fits_COVS, as.matrix)
fits_invmats <- lapply(fits_COVmats, solve)
det_COVmats <- lapply(fits_COVmats, determinant, logarithm = TRUE)
det_COVmats

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

logliksSP <- matrix(nrow = 10, ncol = 10)
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

gcms_event_devs <- array(data = NA, dim = c(nrow(data_gcms_anom_scaled[[10]]),length(fits_invmats)))
for (j in 1:ncol(gcms_event_devs)){
  testvec2 <- data_gcms_anom_scaled[[10]][permutations[[3]]]-fits_Means[[j]]
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
logliksSP[10,]

sum(events_devs)+360*logLik_const[[10]]

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
# globaldata evaulation from SIGMA method
###############################################################
longSPdata <- melt(logliksSP)
# Plot normal value 
bottom <- -500000
top <- max(longdata$value)

#min(longdata$value)
a <- ggplot(longSPdata, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  geom_text(aes(label = round(value/10^4,0))) +
  scale_fill_gradient(name = "value x 10^4",low="white", high="red",
                      #limits = c(bottom,top), 
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
bottom <- -500000
top <- max(longDEVdata$value)

#min(longdata$value)
ggplot(longDEVdata, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  # geom_text(aes(label = round(value/10^4,0))) +
  scale_fill_gradient(name = "Deviation",low="white", high="red"
                      # ,
                      #limits = c(bottom,top), 
                      # labels = function(x)x/10^4
                      ) +
  labs(x="data", y="models", title=paste0(modelsize*100," optimum model evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))


plotdata <- gcms_event_devs
maxl <- max(plotdata)
maxl <- 100
minl <- min(plotdata)
# minl <- -2000
 minl <- minl / 5
plot(plotdata[,10],pch = 16,col = rainbow(10)[1],
     ylim = c(minl,maxl), 
     xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = paste0("Loglik constants and deviations data realizations (wrt interim) METHOD2"))
for(i in 1:9){
  points(plotdata[,i], col = rainbow(10)[i+1],pch = 16)
  abline(h = logLik_const[i],col = rainbow(10)[i+1] )
}

legend("bottomleft", c(abrev[10],abrev[1:9]),bty = "n", pch = 16, col =rainbow(10)[1:10], cex = 0.6)


##########################################################################
# Global dataset evaluation from BN 
##########################################################################
whichBNFits <-selection_fits_gcms_scaled
whichData <- data_gcms_anom_scaled
ordersets <- FALSE

loglik_selection_datasets <- matrix(data = NA, nrow = length(whichBNFits), ncol = length(whichData), dimnames = list(names(whichBNFits),names(whichData)))
for(i in 1:length(whichBNFits)){
  lo <- sapply(X = whichData, logLik, object = whichBNFits[[i]])
  loglik_selection_datasets[i,] <- lo
}

orderby <- "interim_10d_akima_cubic"
logsort <- sort(loglik_selection_datasets[,orderby],index.return = TRUE)
logsort
rep_cnrm <- c(1,2,7,3,4,5,6,8,9,10)

loglik_selection_datasets_ordered <- loglik_selection_datasets[logsort$ix,logsort$ix][rep_cnrm,rep_cnrm]

if (ordersets == TRUE){ longdata <- melt(loglik_selection_datasets_ordered )
} else {longdata <- melt(loglik_selection_datasets)}

########################################################
# Visualize logliks global set
########################################################
# Plot normal value 
bottom <- -300000
top <- max(longdata$value)

#min(longdata$value)
b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  geom_text(aes(label = round(value/10^4,0))) +
  scale_fill_gradient(name = "value x 10^4",low="white", high="red",limits = c(bottom,top), labels = function(x)x/10^4) +
  labs(x="data", y="models", title=paste0(modelsize*100," optimum model evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

b

##################################################################
# Event logliks preparation
##################################################################
whichBNFits <-selection_fits_gcms_scaled
whichData <- data_gcms_anom_scaled
names(cubicsets)
whichSet <- "interim_10d_akima_cubic"
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

gridforEOF <- grid_gcms_anom_scaled
prinCompSet <- prinComp(gridforEOF[[whichSet]], n.eofs = 360)
EOFs <- prinCompSet[[1]][[1]]$EOFs
PCs <- prinCompSet[[1]][[1]]$PCs

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


# randomSet preparation
randomSet <- whichData[[whichSet]]
for(i in 1:ncol(randomSet)){
  randomSet[[i]]<- sample(randomSet[[i]],replace = FALSE)
}


# scalars <- c(0.01,0.1,1,10,100)
# scalar <- 0.01
# 
# respectiveDataSet <- whichData[[whichSet]]
# respectiveDataSet <- projection_EOF
# respectiveDataSet <- projection_EOFs
# respectiveDataSet <- randomSet
# 
# manipulateWith <- rep(scalar,360)
# manipulateWith <- selected_PC
# manipulateWith <- NULL

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

# ev <- eventLoglikCalc(projection_EOFs, NULL)
#################################################################################
# Visualize models vs events
#################################################################################
melteddfn <- melt(gcms_event_devs)
# melteddfn <- melt(logliks_event)
# melteddfn <- melt(logliks_event_scaled)


plot_loglik_event_raster <- function(meltedfn,titlename, limits = NA){
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


###############################################################################
# Visualize correlation matrix loglikelhood data events models
###############################################################################
# loglikscor <- cor(logliks_event)
# loglikscor <- cor(logliks_event_scaled)
# meltloglikscor <- melt(loglikscor)

plot_loglikscor <- function(meltloglikscor,titlename, limits = c(0,1)){
  if (is.null(titlename)){titlename <- paste0("GCMs ", modelsize*100," loglikelihood correlation w.r.t ",whichSet)}
  ggplot(meltloglikscor, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  geom_text(aes(label = round(value,3))) +
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
  
  plot(logliks_event_matrix[,10],pch = 16,col = rainbow(10)[1],
     ylim = c(minl,maxl), 
     xlab = "D_i", ylab = bquote("log P(D_i|"~Mu~")"), main = titlename)
if (!is.null(loglik_constants)){abline(h = loglik_constants[10], col = rainbow(10)[1])}
for(i in 1:9){
  points(logliks_event_matrix[,i], col = rainbow(10)[i+1],pch = 16)
if (!is.null(loglik_constants)){abline(h = loglik_constants[i], col = rainbow(10)[i+1])}
}
  legend("bottomleft", c(abrev[10],abrev[1:9]),bty = "n", pch = 16, col =rainbow(10)[1:10], cex = 0.6) 
}



###############################################################################
# Make Figures
###############################################################################

Logliks_event_normal <- eventLoglikCalc(respectiveDataSet = whichData[[whichSet]], manipulateWith =NULL)
Logliks_event_0.01 <- eventLoglikCalc(respectiveDataSet = whichData[[whichSet]], manipulateWith = rep(0.01,360))
Logliks_event_PC1 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(c(1)), manipulateWith = NULL)
Logliks_event_PC1 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(c(1)), manipulateWith = NULL)
Logliks_event_PC1_3 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:3), manipulateWith = NULL)
Logliks_event_PC1_10 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:10), manipulateWith = NULL)
Logliks_event_PC1_360 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:360), manipulateWith = NULL)

Logliks_event_PC1_2 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:2), manipulateWith = NULL)
Logliks_event_PC1_4 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:4), manipulateWith = NULL)
Logliks_event_PC1_5 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:5), manipulateWith = NULL)
Logliks_event_PC1_6 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:6), manipulateWith = NULL)
Logliks_event_PC1_7 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:7), manipulateWith = NULL)
Logliks_event_PC1_8 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:8), manipulateWith = NULL)
Logliks_event_PC1_9 <- eventLoglikCalc(respectiveDataSet = EOFS_preparation(1:9), manipulateWith = NULL)

Logliks_event_PC_list <- list(Logliks_event_PC1, Logliks_event_PC1_2,Logliks_event_PC1_3,Logliks_event_PC1_4,
     Logliks_event_PC1_5,Logliks_event_PC1_6,Logliks_event_PC1_7,
     Logliks_event_PC1_8,Logliks_event_PC1_9,Logliks_event_PC1_10)

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/Resumen8/figures/Loglik_models_variouseventsEOF.pdf"
pdf(plotname, width = 7, height = 7)
par(mfrow = c(3,2))

plot_loglik_event_points(Logliks_event_normal$logliks_event,loglik_constants = logLik_const, titlename = "Loglik data realizations", minl = -2000, maxl = 100)
plot_loglik_event_points(Logliks_event_0.01$logliks_event,loglik_constants = logLik_const, titlename = "Loglik data realizations factor 0.01")
plot_loglik_event_points(Logliks_event_PC1$logliks_event,loglik_constants = logLik_const, titlename = "Loglik EOF1 part of data ")
plot_loglik_event_points(Logliks_event_PC1_3$logliks_event,loglik_constants = logLik_const, titlename = "Loglik EOF1 + ... + EOF3 part of data ")
plot_loglik_event_points(Logliks_event_PC1_10$logliks_event,loglik_constants = logLik_const, titlename = "Loglik EOF1 + ... + EOF10 part of data ")
plot_loglik_event_points(Logliks_event_PC1_360$logliks_event,loglik_constants = logLik_const, titlename = "Loglik EOF1 + ... + EOF360 part of data ", minl = -2000, maxl = 100)
dev.off()



b1 <- plot_loglik_event_raster(melt(Logliks_event_normal$logliks_event_scaled), titlename = "Normal data")
b2 <- plot_loglik_event_raster(melt(Logliks_event_0.01$logliks_event_scaled), titlename = "Scaled 0.01")
b3 <- plot_loglik_event_raster(melt(Logliks_event_PC1$logliks_event_scaled), titlename = "PC1")
b4 <- plot_loglik_event_raster(melt(Logliks_event_PC1_3$logliks_event_scaled), titlename = "PC 1:3")
b5 <- plot_loglik_event_raster(melt(Logliks_event_PC1_10$logliks_event_scaled), titlename = "PC 1:10")
b6 <- plot_loglik_event_raster(melt(Logliks_event_PC1_360$logliks_event_scaled), titlename = "PC 1:360")
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/Resumen8/figures/Loglik_scaled_models_variouseventsEOF.pdf"
pdf(plotname, width = 7, height = 7)
joep1<- arrangeGrob(b1,b2,b3,b4,b5,b6,nrow = 3)
grid.arrange(joep1, newpage = FALSE)
dev.off()


c1 <- plot_loglikscor(melt(cor(Logliks_event_normal$logliks_event)),titlename = NULL)
c2 <- plot_loglikscor(melt(cor(Logliks_event_0.01$logliks_event)),titlename = "Normal scaled 0.01")
c3 <- plot_loglikscor(melt(cor(Logliks_event_PC1$logliks_event)),titlename = "PC1")
c4 <- plot_loglikscor(melt(cor(Logliks_event_PC1_3$logliks_event)),titlename = "PC 1:3", limits = c(0.75,1))
c5 <- plot_loglikscor(melt(cor(Logliks_event_PC1_10$logliks_event)),titlename = "PC 1:10", limits = c(0.75,1))
c6 <- plot_loglikscor(melt(cor(Logliks_event_PC1_360$logliks_event)),titlename = "PC 1:360")
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/Resumen8/figures/correlationLoglik_models_variouseventsEOF.pdf"
pdf(plotname, width = 13, height = 15)
joep2 <- arrangeGrob(c1,c2,c3,c4,c5,c6,nrow = 3)
grid.arrange(joep2, newpage = FALSE)
dev.off()

d1 <- plot_loglikscor(melt(cor(Logliks_event_PC1_2$logliks_event)),titlename = "PC 1:2", limits = c(0.75,1))
d2 <- plot_loglikscor(melt(cor(Logliks_event_PC1_3$logliks_event)),titlename = "PC 1:3", limits = c(0.75,1))
d3 <- plot_loglikscor(melt(cor(Logliks_event_PC1_4$logliks_event)),titlename = "PC 1:4", limits = c(0.75,1))
d4 <- plot_loglikscor(melt(cor(Logliks_event_PC1_5$logliks_event)),titlename = "PC 1:5", limits = c(0.75,1))
d5 <- plot_loglikscor(melt(cor(Logliks_event_PC1_6$logliks_event)),titlename = "PC 1:6", limits = c(0.75,1))
d6 <- plot_loglikscor(melt(cor(Logliks_event_PC1_7$logliks_event)),titlename = "PC 1:7", limits = c(0.75,1))
d7 <- plot_loglikscor(melt(cor(Logliks_event_PC1_8$logliks_event)),titlename = "PC 1:8", limits = c(0.75,1))
d8 <- plot_loglikscor(melt(cor(Logliks_event_PC1_9$logliks_event)),titlename = "PC 1:9", limits = c(0.75,1))
d9 <- plot_loglikscor(melt(cor(Logliks_event_PC1_10$logliks_event)),titlename = "PC 1:10", limits = c(0.75,1))

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/Resumen8/figures/correlationLoglik_models_caidaEOF.pdf"
pdf(plotname, width = 20, height = 13)
joep3 <- arrangeGrob(d1,d2,d3,d4,d5,d6,d7,d8,d9,nrow = 3)
grid.arrange(joep3, newpage = FALSE)
dev.off()


plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/Resumen8/figures/Loglik_models_2reps.pdf"
pdf(plotname, width = 13, height = 5)
joep <- arrangeGrob(b,a,nrow = 1)
grid.arrange(joep, newpage = FALSE)
dev.off()

##################################################################################
# Differences PC 1:10
##################################################################################
Logliks_event_PC_list[[1]]$logliks_event
Logliks_event_PC_list_cors <- lapply(Logliks_event_PC_list, function(x) cor(x$logliks_event))
diff_cors <- list()
for (i in 1:(length(Logliks_event_PC_list_cors)-1)){
  diff_cors[[i]] <- Logliks_event_PC_list_cors[[i+1]]-Logliks_event_PC_list_cors[[i]]
}
dev.off()
diff_cors_interim <- sapply(1:9,function(x) diff_cors[[x]][,"interim_10d"])

cols <- brewer.pal(9,name = "Set1")

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/logliks_scaleddata/correlationfallgcms_PC1_10_interim.pdf"
pdf(plotname)
plot(diff_cors_interim[10,],pch = 16, col = "red",
     ylim = c(min(diff_cors_interim),max(diff_cors_interim)),
     main = "Correlation Fall w.r.t Interim timeseries for every PC added",xaxt = "n",
     xlab = "PC added",
     ylab = paste0("Cor( log P(Dpc|M), log P(Dpc|",abrev[10],"))"))
for(i in 1:9){lines(diff_cors_interim[i,], col = cols[i])}
legend("bottomleft", c(abrev[10],abrev[1:9]),bty = "n", pch = 16, col = c("red",cols), cex = 0.6) 
axis(side = 1, at = 1:9, labels = as.character(2:10))
dev.off()
