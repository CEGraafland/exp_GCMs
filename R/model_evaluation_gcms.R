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
shortnames
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


###################################################################################################
# Select for each Model its optimum with respect to train and test (random permutation)
###################################################################################################
# load indices of the random Permuation with which train and test are learned
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
samplesize <- dim(cubicsets[[1]]$Data)[1]
indTEST1 <- (1:samplesize)[-indTRAIN1]
# load train and test hcmodels for each (gc)m in list
hc_train_gcms <- lapply(namescubics, loadIterations, permused = 3, hctype = "train")
hc_test_gcms <- lapply(namescubics, loadIterations, permused = 3, hctype = "test")
names(hc_train_gcms) <- shortnames
names(hc_test_gcms) <- shortnames
# create data
data_train_gcms <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE, subind = indTRAIN1)))
data_test_gcms <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE, subind = indTEST1)))
############################################################################
# make data fits one by one (to much memory)
# and calculate loglikelihood values one by one (to much memory)
# and save loglikelihood 
for(k in 1:length(shortnames)){
i <- shortnames[k]
j <- namescubics[k]
hctrains <- hc_train_gcms[[i]]
hctests <- hc_test_gcms[[i]]
#loglikresults <- matrix(data = NA, nrow = length(hctrains), ncol = 3, dimnames = list(names(hctrains),c("traintrain","testtrain","testtest")))
loglikresults <- matrix(data = NA, nrow = length(hctrains), ncol = 2, dimnames = list(names(hctrains),c("traintrain","testtrain")))
fitstrain <- lapply(hctrains, bn.fit, data_train_gcms[[j]])
fitstest <- lapply(hctests, bn.fit, data_test_gcms[[j]])

loglikstrain <- sapply(X = fitstrain, logLik, data = data_train_gcms[[j]])
loglikstest <- sapply(X = fitstrain, logLik, data = data_test_gcms[[j]])
#loglikstesttest <- sapply(X = fitstest, logLik, data = data_test_gcms[[j]])

loglikresults[,"traintrain"] <- loglikstrain
loglikresults[,"testtrain"] <- loglikstest
#loglikresults[,"testtest"] <- loglikstesttest

assign(paste0("loglik_traintest_",i), loglikresults)
save(list = paste0("loglik_traintest_",i), 
     file = paste0("results/loglik_traintest/loglik_traintest_",i,".rda"))

loglikstrain <- loglikstest <- fitstrain <- fitstest <- hctrains <- hctests <- NULL
}
##################################################################################
# Skipped the above: 
# load results loglikelihood train and tests and 
# combine loglikelihood results of every single model
loglik_traintest_gcms <- lapply(shortnames, function(x){get(load(paste0("results/loglik_traintest/loglik_traintest_",x,".rda")))})
names(loglik_traintest_gcms) <- shortnames
# then select the statistical optimums of the models
own_optimums <-sapply(loglik_traintest_gcms, function(x) which(x[,"testtrain"] == max(x[,"testtrain"])))

################################################################################
# Compare the models on their generalization capacities
# How do the hc models explain data of other models ? 
# First figures (under, equal, optimum fits)
# EVALUATED ON TRAIN DATASET????
################################################################################
hc_gcms <- lapply(namescubics, loadIterations, permused = 3)
names(hc_gcms) <- shortnames
data_gcms <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE)))

under_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = 6, SIMPLIFY = FALSE)
equal_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = 18, SIMPLIFY = FALSE)
optimum_hc_gcms <- mapply(function(x,y) x[[y]], x = hc_gcms, y = own_optimums, SIMPLIFY = FALSE)

under_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = under_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)
equal_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = equal_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)
optimum_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = optimum_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)

loglik_under_datasets <- matrix(data = NA, nrow = length(under_fits_gcms), ncol = length(data_gcms), dimnames = list(names(under_fits_gcms),names(data_gcms)))
loglik_optimum_datasets <- matrix(data = NA, nrow = length(optimum_fits_gcms), ncol = length(data_gcms), dimnames = list(names(optimum_fits_gcms),names(data_gcms)))
loglik_equal_datasets <- matrix(data = NA, nrow = length(equal_fits_gcms), ncol = length(data_gcms), dimnames = list(names(equal_fits_gcms),names(data_gcms)))



for(i in 1:length(under_fits_gcms)){
  lo <- sapply(X = data_train_gcms, logLik, object = under_fits_gcms[[i]])
  loglik_under_datasets[i,] <- lo
}
loglik_under_datasets

for(i in 1:length(optimum_fits_gcms)){
lo <- sapply(X = data_train_gcms, logLik, object = optimum_fits_gcms[[i]])
loglik_optimum_datasets[i,] <- lo
}
loglik_optimum_datasets

for(i in 1:length(equal_fits_gcms)){
  lo <- sapply(X = data_train_gcms, logLik, object = equal_fits_gcms[[i]])
  loglik_equal_datasets[i,] <- lo
}
loglik_equal_datasets



par(mar=c(5, 4, 4, 2) + 0.1)
par(xpd=T, mar=par()$mar+c(2,0,0,0))
plot(loglik_optimum_datasets[10,], xaxt = "n", xlab = c(""), ylab = c("Loglik(x|model)"), col = rainbow(10)[10],ylim = c(-600000,max(loglik_optimum_datasets)),pch = 16)
for(i in 1:9){
  points(loglik_optimum_datasets[i,], col = rainbow(10)[i],pch = 16)
}
legend("bottomright",legend = gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = ""), fill = rainbow(10),cex = 0.7)
axis(1,seq(1,10,1),gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = ""),las = 2)


par(mfrow = c(1,2))
plot(loglik_equal_datasets[10,], xaxt = "n", xlab = c(""), ylab = c("Loglik(x|model)"), col = rainbow(10)[10],ylim = c(-200000,max(loglik_equal_datasets)),pch = 16)
for(i in 1:9){
  points(loglik_equal_datasets[i,], col = rainbow(10)[i],pch = 16)
}
legend("bottomright",legend = gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = ""), fill = rainbow(10),cex = 0.7)
axis(1,seq(1,10,1),gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = ""),las = 2)

plot(loglik_under_datasets[10,], xaxt = "n", xlab = c(""), ylab = c("Loglik(x|model)"), col = rainbow(10)[10],ylim = c(-200000,max(loglik_equal_datasets)),pch = 16)
for(i in 1:9){
  points(loglik_under_datasets[i,], col = rainbow(10)[i],pch = 16)
}
legend("bottomright",legend = gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = ""), fill = rainbow(10),cex = 0.7)
axis(1,seq(1,10,1),gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = ""),las = 2)

abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")

# Short cluster analísis; corresponds to general results later
fivemeans <- kmeans(loglik_optimum_datasets,5)
fivemeans$cluster
# heat plot 
longdata <- melt(loglik_optimum_datasets) 


ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red",limits = c(-100000,max(longdata$value))) +
  labs(x="data", y="models", title="optimum model evaluation") +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))


################################################################################
# How do the models of various gcms explain the data of various gcms?
# Generalization
# Note that selection_fits_gcms uses unscaled (no Field) data
################################################################################
hc_gcms <- lapply(namescubics, loadIterations, permused = 3)
names(hc_gcms) <- shortnames
data_gcms <- lapply(cubicsets,function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE)))

modelsize <- 18

selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms, SIMPLIFY = FALSE)
loglik_selection_datasets <- matrix(data = NA, nrow = length(selection_fits_gcms), ncol = length(data_gcms), dimnames = list(names(selection_fits_gcms),names(data_gcms)))
# hc_gcms <- NULL

for(i in 1:length(selection_fits_gcms)){
  lo <- sapply(X = data_gcms, logLik, object = selection_fits_gcms[[i]])
  loglik_selection_datasets[i,] <- lo
}

abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")

orderby <- "interim_10d_akima_cubic"
logsort <- sort(loglik_selection_datasets[,orderby],index.return = TRUE)
logsort
rep_cnrm <- c(1,2,7,3,4,5,6,8,9,10)

loglik_selection_datasets_ordered <- loglik_selection_datasets[logsort$ix,logsort$ix][rep_cnrm,rep_cnrm]
longdata <- melt(loglik_selection_datasets_ordered )
longdata$logvalue <- log10(-longdata$value)
longdata

# bottom <- log10(-min(longdata$value))
bottom <- log10(200000)
top <- log10(-max(longdata$value))

# Plot logvalue (not necesarilly -> second plot)
a <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=logvalue)) + 
  scale_fill_gradient(high="grey90", low="red",limits = c(top,bottom)) +
  labs(x="data", y="models", title=paste0(modelsize*100," optimum model evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

assign(paste0("model",modelsize*100),a)
model1800
a

# Plot normal value 
bottom <- -300000
top <- max(longdata$value)

#min(longdata$value)
a <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  geom_text(aes(label = round(value/10^4,0))) +
  scale_fill_gradient(name = "value x 10^4",low="white", high="red",limits = c(bottom,top), labels = function(x)x/10^4) +
  labs(x="data", y="models", title=paste0(modelsize*100," optimum model evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))
assign(paste0("model",modelsize*100),a)
a


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/model_data_matrix_gcms_",modelsize*100,"_",orderby,".pdf")
pdf(plotname, height = 6, width = 7)

model1800
dev.off()

######################################################################
# Figure with loglikelihoods of under, equal and over estimated models
# If generated below models with different sizes
######################################################################
grid.arrange(model600,model1800,model7000)

######################################################################
# plot loglikelihood different method than heatmap. 
#####################################################################
plot(loglik_selection_datasets[10,], xaxt = "n", xlab = c(""), ylab = c("Loglik(x|model)"), col = rainbow(10)[10],ylim = c(-600000,max(loglik_selection_datasets)),pch = 16)
for(i in 1:9){
  points(loglik_selection_datasets[i,], col = rainbow(10)[i],pch = 16)
}
legend("bottomright",legend = gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = ""), fill = rainbow(10),cex = 0.7)
axis(1,seq(1,10,1),gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = ""),las = 2)


##################################################################
# Event logliks data
##################################################################
##################################################################
# Make al datasets standardized
##################################################################
# Anomalies and standarization over seasons (networks are learned with this set)
data_gcms_out <- lapply(cubicsets,function(x) TimeCoordsAnom_from_Grid_rms(x, rms = TRUE))
data_gcms_anom_out <- lapply(data_gcms_out, function(x)mat2Dto3Darray(x,attr(x,"Xcoords"), attr(x,"Ycoords")))
grid_gcms_anom <- mapply(function(x,y) {x$Data <- y 
                                        return(x)}, x = cubicsets, y = data_gcms_anom_out,SIMPLIFY = FALSE)
# Global standarazation over the above: FIELD
grid_gcms_anom_scaled <- lapply(grid_gcms_anom, function(x) scaleGrid(x,spatial.frame ="field",type = "standardize"))

# Extract the FIELD standarized data
data_gcms_anom_scaled <- lapply(grid_gcms_anom_scaled,function(x) {x <- redim(x,drop = TRUE)
                                                                    xdata <- array3Dto2Dmat(x$Data) 
                                                                    return(xdata)}) 
data_gcms_anom_scaled <- lapply(data_gcms_anom_scaled, function(x) as.data.frame(x))

# Make fits
# unscaled
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out, SIMPLIFY = FALSE)
# scaled
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)

# all.equal(selection_fits_gcms,selection_fits_gcms_scaled)
# all.equal(selection_fits_gcms_scaled[[1]],selection_fits_gcms[[1]])
# selection_fits_gcms_scaled[[10]]$V215$sd
# selection_fits_gcms[[10]]$V215$sd

# m1<-lapply(data_gcms_anom_scaled,colMeans)
# lapply(m1,mean)

# whichData <- cubicsets$interim_10d_akima_cubic
whichBNFits <-selection_fits_gcms_scaled

scalars <- c(0.01,0.1,1,10,100)
j <- 2
loglik_gcms_event_vs_int_scalars <- array(data = NA, dim = c(nrow(dataRMS),length(whichBNFits),length(scalars)), dimnames = list(NULL,names(whichBNFits),as.character(scalars)))
for (j in 1:length(scalars)){
  scalar <- scalars[j]
  dataRMS <- data_gcms_anom_scaled$interim_10d_akima_cubic
  # dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(whichData, rms = TRUE))
  dataRMSframes <- apply(dataRMS, MARGIN = 1, FUN = function(x) as.data.frame(x))
  dataRMSframes <- lapply(dataRMSframes, t)
  dataRMSframes <- lapply(dataRMSframes, as.data.frame)

  # Manipulate data with scalars
  for(i in 1:length(dataRMSframes)){
    dataRMSframes[[i]] <- dataRMSframes[[i]]*scalar
  }
  # # Manipulate data with PCs
  # for(i in 1:length(dataRMSframes)){
  #    dataRMSframes[[i]] <- dataRMSframes[[i]]*PCs[i]


  # dataRMSframes[[2]] <- dataRMSframes[[1]]*0.8
  logliks_gcms_vs_int_event <- array(data = NA, dim = c(nrow(dataRMS),length(whichBNFits)))

  for (i in 1:length(whichBNFits)){
    logliks_gcms_vs_int_event[,i] <- sapply(dataRMSframes, FUN = logLik, object = whichBNFits[[i]])
  }
  loglik_gcms_event_vs_int_scalars[,,j] <-  logliks_gcms_vs_int_event
}

all.equal(loglik_gcms_event_vs_int_scalars[,,4],loglik_gcms_event_vs_int_scalars[,,5])
############################################################################
# DIfferent scalar plots
############################################################################
# plot(logliks_gcms_vs_int_event[1,],logliks_gcms_vs_int_event[2,])
plotname <- paste0("figs/logliks_scaleddata/Loglik_datarealizations_gcms_vs_int_databy",scalar,".pdf")
pdf(plotname)

maxl <- max(logliks_gcms_vs_int_event)
# maxl <- 100
minl <- min(logliks_gcms_vs_int_event)
# minl <- -2000
# minl <- minl / 5
plot(logliks_gcms_vs_int_event[,10],pch = 16,col = rainbow(10)[1],
     ylim = c(minl,maxl), 
     xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = paste0("Loglik data realizations given models scaled by factor",scalar))
for(i in 1:9){
  points(logliks_gcms_vs_int_event[,i], col = rainbow(10)[i+1],pch = 16)
}

legend("bottomleft", c(abrev[10],abrev[1:9]),bty = "n", pch = 16, col =rainbow(10)[1:10], cex = 0.6)

dev.off()



#########################################################################
# Scaled logliks from data events
#########################################################################

loglik_gcms_event_vs_int_scalars_scaled <- apply(loglik_gcms_event_vs_int_scalars, MARGIN = c(2), FUN = function(x){x <- scale(x,center = TRUE,scale = TRUE)
return(x)}) 
apply(loglik_gcms_event_vs_int_scalars[,,1], MARGIN = c(2), FUN = function(x){x <- scale(x,center = TRUE,scale = TRUE)
return(x)})                                                                     
dim(loglik_gcms_event_vs_int_scalars_scaled)
all.equal(loglik_gcms_vs_int_event_scalar_scaled[2,,],loglik_gcms_vs_int_event_scalar_scaled[3,,])


logliks_gcms_vs_int_event_scaled <- scale(logliks_gcms_vs_int_event, center = TRUE, scale = TRUE)
loglik_gcms_event_vs_int_scalars_scaled <- array(data = NA, dim = c(nrow(dataRMS),length(whichBNFits),length(scalars)), dimnames = list(NULL,names(whichBNFits),as.character(scalars)))
for (j in 1:length(scalars)){
loglik_gcms_event_vs_int_scalars_scaled[,,j] <- scale(loglik_gcms_event_vs_int_scalars[,,j], center = TRUE, scale = TRUE)
}

all.equal(loglik_gcms_event_vs_int_scalars_scaled[,,4],loglik_gcms_event_vs_int_scalars_scaled[,,3])
colnames(logliks_gcms_vs_int_event_scaled) <- abrev
melteddfn <- melt(logliks_gcms_vs_int_event_scaled)
colnames(loglik_gcms_event_vs_int_scalars_scaled) <- abrev
melteddfn <- melt(loglik_gcms_event_vs_int_scalars_scaled[,,5])
# plot(logliks_gcms_vs_int_event_scaled[,8])


slde <- ggplot(melteddfn, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(high="blue", low="red", midpoint =0, 
                       #limits = c(-200,100), 
                       na.value = "black") +
  labs(x= paste0("Interim events"), y="models", title=paste0("GCMs ", modelsize*100,"evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

plotname <- paste0("figs/logliks_scaleddata/Scaled_Loglik_data_gcms_vs_int_databy",scalar,".pdf")
pdf(plotname)
slde
dev.off()
#######################################################################
# correlation matrix from scaled logliks normal events. 
#######################################################################
plotname <- paste0("figs/logliks_scaleddata/Scaled_Loglik_correlation_gcms_vs_int_databy",scalar,".pdf")
pdf(plotname)
loglikscor <- cor(logliks_gcms_vs_int_event_scaled)
meltloglikscor <- melt(loglikscor)
ggplot(meltloglikscor, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(high="blue", low="white", 
                       #midpoint =0, 
                       #limits = c(-200,100), 
                       na.value = "black") +
  labs(x= paste0("Interim events"), y="models", title=paste0("GCMs ", modelsize*100,"evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))
dev.off()
}
###########################################################################
# Normal plot (scalar = 1) 
###########################################################################
plotname <- paste0("figs/Loglik_datarealizations_gcms_vs_int_color.pdf")
pdf(plotname)

maxl <- max(logliks_gcms_vs_int_event)
maxl
minl <- min(logliks_gcms_vs_int_event)
minl
minl <- minl / 5
plot(logliks_gcms_vs_int_event[,10],pch = 16,col = rainbow(10)[1],
     ylim = c(minl,maxl), 
     xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = "Loglik data realizations given models")
for(i in 1:9){
  points(logliks_gcms_vs_int_event[,i], col = rainbow(10)[i+1],pch = 16)
}

legend("bottomleft", c(abrev[10],abrev[1:9]),bty = "n", pch = 16, col =rainbow(10)[1:10], cex = 0.6)

dev.off()

######################################################################
# Old plot blue, gray. 
plotname <- paste0("figs/Loglik_datarealizations_gcms_vs_int.pdf")
pdf(plotname)

bluebn<- colorRampPalette(c("light blue","navy"), alpha = TRUE)(9)
blackcn <- colorRampPalette(c("grey","black"), alpha = TRUE)(4)

plot(logliks_gcms_vs_int_event[,10],col = blackcn[3],ylim = c(-1000,500), xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = "Loglik data realizations given models")
for(i in 1:9){
  points(logliks_gcms_vs_int_event[,i], col = bluebn[i],pch = i)
}

legend("topleft", c(abrev[10],abrev[1:9]),bty = "n", pch = c(1,1:9), col =c(blackcn[3],bluebn), cex = 0.6)

dev.off()



# #############################################################################
# # EOFs logliks without global standardization. 
# #############################################################################
# whichSet <- data_gcms_out$interim_10d_akima_cubic
# data_Int_anom_out <- mat2Dto3Darray(whichset, attr(whichset,"Xcoords"), attr(whichset,"Ycoords"))
# grid_Int_anom <- cubicsets$interim_10d_akima_cubic
# grid_Int_anom$Data <- data_Int_anom_out
# 
# prinCompInt <- prinComp(grid_Int_anom, n.eofs = 100, scaling = "gridbox")
# prinCompInt
# 
# EOFs_Int <- prinCompInt$`2T`[[1]]$EOFs
# PCs_Int <- prinCompInt$`2T`[[1]]$PCs
# thirdEOFint<- EOFs_Int[,"PC10"]
# thirdPCint <- PCs_Int[,"PC10"]
# dim(thirdPCint) <- c(length(thirdPCint),1)
# dim(thirdEOFint) <- c(1,length(thirdEOFint))
# 
# third_EOFs_int <- t(sapply(thirdPCint, function(x) thirdEOFint*x))
# plotEOF(prinCompInt,backdrop.theme = "coastline")
# dataRMSframes
# 
# dataEOFframes <- apply(third_EOFs_int, MARGIN = 1, FUN = function(x) as.data.frame(x))
# dataEOFframes <- lapply(dataEOFframes, t)
# dataEOFframes <- lapply(dataEOFframes, as.data.frame)
# 
# 
# 
# logliks_gcms_vs_int_EOFs <- array(data = NA, dim = c(nrow(third_EOFs_int),length(selection_fits_gcms)))
# 
# 
# for (i in 1:length(selection_fits_gcms)){
#   logliks_gcms_vs_int_EOFs[,i] <- sapply(dataEOFframes, FUN = logLik, object = selection_fits_gcms[[i]])
# }
# 
# bluebn<- colorRampPalette(c("light blue","navy"), alpha = TRUE)(9)
# blackcn <- colorRampPalette(c("grey","black"), alpha = TRUE)(4)
# 
# plotname <- paste0("figs/Loglik_datarealizations_gcms_vs_int.pdf")
# pdf(plotname)
# 
# plot(logliks_gcms_vs_int_EOFs[,10],col = blackcn[3],ylim = c(-200,100), xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = "Loglik data realizations given models")
# for(i in 1:9){
#   points(logliks_gcms_vs_int_EOFs[,i], col = bluebn[i],pch = i)
# }
# 
# legend("topleft", c(abrev[10],abrev[1:9]),bty = "n", pch = c(1,1:9), col =c(blackcn[3],bluebn), cex = 0.6)
# 
# dev.off()


#############################################################################
# EOFs logliks with global standarization: Field
#############################################################################

# whichset <- grid_gcms_anom_scaled[[10]]
# data_Int_anom_out <- mat2Dto3Darray(whichset, attr(whichset,"Xcoords"), attr(whichset,"Ycoords"))
# grid_Int_anom <- cubicsets$interim_10d_akima_cubic
# grid_Int_anom$Data <- data_Int_anom_out
# grid_Int_anom <- scaleGrid(grid_Int_anom,spatial.frame = "field",type = "standardize")
# prinCompInt <- prinComp(grid_Int_anom, n.eofs = 360, scaling = "gridbox")
# prinCompInt
# all.equal(grid_Int_anom,grid_gcms_anom_scaled[[10]])

prinCompsScaled <- lapply(grid_gcms_anom_scaled, function(x) prinComp(x,n.eofs = 10))
###################################################################################
# Visuzlize EOFS
##################################################################################
for(i in 1:length(prinCompsScaled)){
  plotname <- paste0("figs/EOFs/EOFs_10_",names(prinCompsScaled)[i],".pdf")
  pdf(plotname)
  plotEOF(prinCompsScaled[[i]], backdrop.theme = "coastline",n.eofs = 10, main = names(prinCompsScaled)[i], at = seq(-0.18,0.18,0.02) )
  dev.off()
}
  
##########################################################################
#
##########################################################################
whichSet <- "interim_10d_akima_cubic"
whichPC  <- "PC1"
whichBNFits <- selection_fits_gcms_scaled
prinCompSet <- prinComp(grid_gcms_anom_scaled[[whichSet]], n.eofs = 360)
EOFs <- prinCompSet[[1]][[1]]$EOFs
PCs <- prinCompSet[[1]][[1]]$PCs

selected_EOF <- EOFs[,whichPC]
selected_PC <- PCs[,whichPC]

projection_EOF <- matrix(nrow = length(selected_PC),ncol = length(selected_EOF))
for(i in 1:length(selected_PC)){
  projection_EOF[i,] <- selected_EOF*selected_PC[i]
}

dim(selected_PC) <- c(length(selected_PC),1)
dim(selected_EOF) <- c(1,length(selected_EOF))

projection_EOF2 <- t(sapply(selected_PC, function(x) selected_EOF*x))
all.equal(projection_EOF,projection_EOF2)

# plot(projection_EOF[3,])
same <- array(data = selected_EOF, dim = c(length(selected_EOF),360))
same <- t(same)
all.equal(same[1,],same[2,])

dataEOFframes <- apply(same, MARGIN = 1, FUN = function(x) as.data.frame(x))
dataEOFframes <- lapply(dataEOFframes, t)
dataEOFframes <- lapply(dataEOFframes, as.data.frame)
same[2,]
dataEOFframes[[2]]

dataEOFframes[[1]]
projection_EOF[]

for(i in 1:length(dataEOFframes)){
  dataEOFframes[[i]] <- dataEOFframes[[i]]*selected_PC[i]
}

all.equal(dataEOFframes[[3]],projection_EOF2[3,])

logliks_gcms_vs_int_EOFs <- array(data = NA, dim = c(nrow(projection_EOF),length(whichBNFits)))
# dataEOFframes[[1]]<- dataEOFframes[[2]]*0.5

for (i in 1:length(whichBNFits)){
  logliks_gcms_vs_int_EOFs[,i] <- sapply(dataEOFframes, FUN = logLik, object = whichBNFits[[i]])
}
abrev
dev.off()
plot(logliks_gcms_vs_int_EOFs[6,],logliks_gcms_vs_int_EOFs[7,])
logLik(dataEOFframes[[3]], object = selection_fits_gcms_scaled[[1]])
logLik(dataEOFframes[[3]], object = selection_fits_gcms_scaled[[2]])

logliks_gcms_vs_int_EOFs[1:3,1:2]
cor(logliks_gcms_vs_int_EOFs[1:3,1:2])

logLik(dataRMSframes[[1]],  object = selection_fits_gcms_scaled[[1]])
logLik(dataRMSframes[[1]],  object = selection_fits_gcms_scaled[[2]])

logliks_gcms_vs_int_event_scaled[,2]

cor(logliks_gcms_vs_int_event_scaled, method = "pearson")
cor(logliks_gcms_vs_int_EOFs, method = "spearman")
cov2cor(v)


logliks_gcms_vs_int_EOFs_scaled <-scale(logliks_gcms_vs_int_EOFs, center = TRUE, scale = TRUE)
#plotname <- paste0("figs/Loglik_datarealizations_gcms_vs_int.pdf")
#pdf(plotname)

bluebn<- colorRampPalette(c("light blue","navy"), alpha = TRUE)(9)
blackcn <- colorRampPalette(c("grey","black"), alpha = TRUE)(4)

plot(logliks_gcms_vs_int_EOFs[,10],col = blackcn[3],ylim = c(-200,100), xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = "Loglik data realizations given models")
for(i in 1:9){
  points(logliks_gcms_vs_int_EOFs[,i], col = bluebn[i],pch = i)
}

legend("topleft", c(abrev[10],abrev[1:9]), bty = "n", pch = c(1,1:9), col =c(blackcn[3],bluebn), cex = 0.6)

dev.off()



colnames(logliks_gcms_vs_int_EOFs_scaled) <- abrev
melteddf <- melt(logliks_gcms_vs_int_EOFs_scaled)

ggplot(melteddf, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(high="blue", low="red", midpoint =0, 
                       #limits = c(-50,50), 
                       na.value = "black") +
  labs(x= paste0("Interim projected EOF ",whichPC, " events"), y="models", title=paste0("GCMs ", modelsize*100," EOF ",whichPC," ",whichSet," evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

plot(thirdPCint)
#############################################################################
# Correlation between loglikelihood values models. 
cor()


#############################################################################
# random Event logliks globally scaled
#############################################################################

dataRMS<- data_gcms_anom_scaled$interim_10d_akima_cubic
dataRMSrandom <- dataRMS

for(i in 1:ncol(dataRMS)){
  dataRMSrandom[[i]]<- sample(dataRMS[[i]],replace = FALSE)
}
all.equal(dataRMSrandom,dataRMS)
all.equal(colSums(dataRMSrandom),colSums(dataRMS))


# dataRMS.2 <- as.data.frame(TimeCoordsAnom_from_Grid_rms(cubicsets$interim_10d_akima_cubic, rms = TRUE))
dataRMSframes <- apply(dataRMSrandom, MARGIN = 1, FUN = function(x) as.data.frame(x))
dataRMSframes <- lapply(dataRMSframes, t)
dataRMSframes <- lapply(dataRMSframes, as.data.frame)

logliks_gcms_vs_int_event_r <- array(data = NA, dim = c(nrow(dataRMS),length(selection_fits_gcms_scaled)))

for (i in 1:length(selection_fits_gcms)){
  logliks_gcms_vs_int_event_r[,i] <- sapply(dataRMSframes, FUN = logLik, object = selection_fits_gcms_scaled[[i]])
}

bluebn<- colorRampPalette(c("light blue","navy"), alpha = TRUE)(9)
blackcn <- colorRampPalette(c("grey","black"), alpha = TRUE)(4)

plotname <- paste0("figs/Loglik_random_datarealizations_gcms_vs_int.pdf")
pdf(plotname)

maxl <- max(logliks_gcms_vs_int_event_r)
maxl
minl <- min(logliks_gcms_vs_int_event_r)
minl
minl <- minl / 5
plot(logliks_gcms_vs_int_event_r[,10],pch = 16,col = rainbow(10)[1],
     ylim = c(-500000,maxl), 
     xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = "Loglik random data realizations given models")
for(i in 1:9){
  points(logliks_gcms_vs_int_event_r[,i], col = rainbow(10)[i+1],pch = 16)
}

legend("bottomleft", c(abrev[10],abrev[1:9]),bty = "n", pch = 16, col =rainbow(10)[1:10], cex = 0.6)

dev.off()


########################################################################
# Mean and SD original grid
########################################################################

ens <- bindGrid(cubicsets, dimension = "member")
ens$Members <- abrev
ens2 <- bindGrid(cubicsets2, dimension = "member")
ens2$Members <- abrev

cubicsets2 <- cubicsets
cubicsets2$interim_10d_akima_cubic$Data <- cubicsets2$interim_10d_akima_cubic$Data-273.15
cubicsets$CMIP5_CanESM2_r1i1p1$Data

spatialPlot(climatology(ens2, clim.fun = list(FUN = function(x){mean(x)})), 
            backdrop.theme = "coastline", as.table = TRUE,
            rev.colors = TRUE,
            lonCenter = 180)

spatialPlot(climatology(ens, clim.fun = list(FUN = sd)), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,
            rev.colors = TRUE)

spatialPlot(climatology(ens2, clim.fun = list(FUN = sd)), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,
            rev.colors = TRUE)


######################################################################
# mean and sd after TimeCoordsAnom function. 
######################################################################
# All are equal:
sdgcms <- lapply(data_gcms, function(x) apply(X = x, MARGIN = 2, FUN = sd))
# Means differ: TimeCoordsAnom substracts 
meangcms <- lapply(data_gcms, function(x) colMeans(x))
clims <- mapply(quantity2clim, quantity = meangcms, ref.grid = cubicsets, MoreArgs = list(backperm = backpermutations[[1]], what = "mean"), SIMPLIFY = FALSE)
ensclims <- bindGrid(clims, dimension = c("member"))
ensclims$Members <- abrev
spatialPlot(ensclims,backdrop.theme = "coastline", as.table = TRUE, lonCenter = 180, rev.colors = TRUE)

######################################################################
# Hier verder gaan met clustering bijv?
######################################################################
fivemeans <- kmeans(loglik_optimum_datasets,5)
fivemeans$cluster
loglik_optimum_datasets
# fits_gcms <- mapply(hc_gcms, function(sel) lapply(FUN = bn.fit, x = sel, data = data_gcms)
# fits: order:
aver1 <- lapply(hc_gcms[[1]], function(x) lapply(data_gcms, function(y) bn.fit(x = x, data = y)))
rm(aver1)
# all fits:
fits_gcms <- lapply(hc_gcms, function(z) lapply(z, function(x) lapply(data_gcms, function(y) bn.fit(x = x, data = y))))
