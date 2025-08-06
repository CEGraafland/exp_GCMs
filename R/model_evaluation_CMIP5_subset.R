#######################################################################
# Model evaluation ex subset 5d
#######################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(sparsebnUtils)
library(gridExtra)
library(RColorBrewer)
library(gaussDiff)

########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_cmip5_exsubset_5d_akima_cubic.rda")
load("data/hciterations/permutations/perm1.rda")
load("data/hciterations/permutations/5d/indTRAIN1.rda")
load("data/hciterations/permutations/5d/indTEST1.rda")
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
cubicsets <- get(load("data/tas_historical_cmip5_exsubset_5d_akima_cubic.rda"))

namescubics <- names(cubicsets)
# shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
#                    pattern = "_r12i1p1", replacement = "")
# shortnames
############################################################################
# Function to load hciterations of models in a list 
############################################################################
loadIterations <- function(pattern,permused, hctype = NULL,it = NULL){
  if(is.null(hctype)){hctype <- ""}
  hc_interim_list <- list.files(paste0("data/hciterations/CMIP5_5d/",pattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0("data/hciterations/CMIP5_5d/",pattern,"/perm",permused,hctype))
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


###################################################################################################
# Select for each Model its optimum with respect to train and test (random permutation)
###################################################################################################
# load indices of the random Permuation with which train and test are learned

samplesize <- dim(cubicsets[[1]]$Data)[1]
indTEST1 <- (1:samplesize)[-indTRAIN1]
# load train and test hcmodels for each (gc)m in list
hc_train_gcms <- lapply(namescubics[5], loadIterations, permused = 1, hctype = "train")
hc_test_gcms <- lapply(namescubics[5], loadIterations, permused = 1, hctype = "test")
names(hc_train_gcms) <- namescubics[5]
names(hc_test_gcms) <- namescubics[5]
# create data

data_train_gcms <- lapply(list(cubicsets[[5]]),function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE, subind = indTRAIN1)))
data_test_gcms <- lapply(list(cubicsets[[5]]),function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE, subind = indTEST1)))
####################################################################
#
####################################################################
k <- 5
#for(k in 1:length(shortnames)){
  i <- namescubics[k]
  j <- namescubics[k]
  hctrains <- hc_train_gcms[[i]]
  hctests <- hc_test_gcms[[i]]
  #loglikresults <- matrix(data = NA, nrow = length(hctrains), ncol = 3, dimnames = list(names(hctrains),c("traintrain","testtrain","testtest")))
  loglikresults <- matrix(data = NA, nrow = length(hctrains), ncol = 2, dimnames = list(names(hctrains),c("traintrain","testtrain")))
  fitstrain <- lapply(hctrains, bn.fit, data_train_gcms[[1]])
  fitstest <- lapply(hctests, bn.fit, data_test_gcms[[1]])
  
  loglikstrain <- sapply(X = fitstrain, logLik, data = data_train_gcms[[1]])
  loglikstest <- sapply(X = fitstrain, logLik, data = data_test_gcms[[1]])
  narcstrain <- sapply(X = hc_train_gcms[[i]], narcs)
  #loglikstesttest <- sapply(X = fitstest, logLik, data = data_test_gcms[[j]])
  
  loglikresults[,"traintrain"] <- loglikstrain
  loglikresults[,"testtrain"] <- loglikstest
  #loglikresults[,"testtest"] <- loglikstesttest
  
  assign(paste0("loglik_traintest_",i), loglikresults)
  # save(list = paste0("loglik_traintest_",i), 
       # file = paste0("results/loglik_traintest/loglik_traintest_",i,".rda"))
  
  loglikstrain <- loglikstest <- fitstrain <- fitstest <- hctrains <- hctests <- NULL
# }
gc()
loglik_traintest_tas_2T_interim
plot(narcstrain,loglik_traintest_tas_2T_Interim[,1], 
     #xaxt = "n", 
     xlab = c("|E|"), ylab = c("Loglik(x|model)"), 
     col = rainbow(10)[10],
     #ylim = c(-600000,max(loglik_optimum_datasets)),
     pch = 16)

points(narcstrain,loglik_traintest_tas_2T_Interim[,2], col = rainbow(10)[9],pch = 16)
##############################################################
#
##############################################################
##################################################################
# Make al datasets standardized
##################################################################
permk <- perm1
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
assign(paste0("data_gcms_anom_scaled_1"),lapply(data_gcms_anom_scaled, function(x) x[,permk]))
length(data_gcms_anom_scaled)

it <- "2900_3000"
hc_gcms5sub <- lapply(paste0(namescubics), loadIterations, permused = permk, it = it)
modelsize <- 30
selection_hc_gcms <-  mapply(function(x,y) x[[grep(y,x)]], x = hc_gcms5sub, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
#selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
########################################################################################
# Calculate hellinger distance bayesian networks
# KL distances give NaN (because there are x such that P(x) = 0 does not imply Q(x) = 0)
########################################################################################
fits_NEL <- lapply(selection_fits_gcms_scaled,bnlearn:::as.graphNEL)
edgelistssparse <- lapply(X =fits_NEL,FUN = edgeList)
data_gcms_anom_scaled_sparse <- lapply(data_gcms_anom_scaled_1,sparsebnData, type = "continuous")
fits_COVS <- mapply(get.covariance, x = edgelistssparse, data = data_gcms_anom_scaled_sparse)
fits_COVmats <- lapply(fits_COVS, as.matrix)
fits_Means <- lapply(data_gcms_anom_scaled, function(x) colMeans(x[perm1]) )

hellinger_coefficients <- matrix(data = NA, nrow = length(fits_COVmats), ncol = length(fits_COVmats), dimnames = list(names(fits_COVmats),names(fits_COVmats)))
i <- 1
for(i in 1:length(fits_COVmats)){
  kl <- mapply(function(b,d) normdiff(mu1 =fits_Means[[i]], mu2 = b, sigma1 = fits_COVmats[[i]], sigma2 = d, method = "Hellinger"),b = fits_Means, d = fits_COVmats)
  hellinger_coefficients[i,] <- kl
}
str(fits_Means)

rownames(hellinger_coefficients) <- colnames(hellinger_coefficients) <-names(fits_Means)
hellinger_coefficients_sub5d_3000<- hellinger_coefficients
save(hellinger_coefficients_sub5d_3000, file ="results/hellinger_coefficient/hellinger_coefficients_sub5d_3000.rda")
####################################################################################
#
####################################################################################
library("gmp")
load(file ="results/hellinger_coefficient/hellinger_coefficients_sub5d_3000.rda")
log(hellinger_coefficients_sub5d_3000)
num <- as.bigz(hellinger_coefficients_sub5d_3000[3,2])
num1 <- as.bigz(1)
num1-num 
########################################################
# Visualize hellinger distance example set 
########################################################
#ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")

plotname <- paste0("figs/hellinger_CMIP5_ordering/Bhat_CMIP5_ex_ACCESS|BNU|CMCC.CMS|interim_5d_3000.png")
png(plotname,units = "in", width = 9, height = 7,  res =180)
# plotname <- paste0("figs/hellinger_CMIP5_ordering/Hellinger_CMIP5_ex_ACCESS|BNU|CMCC.CMS|interim.pdf")
# pdf(plotname, width = 9, height =7)
# n <-7
# n <-7
# br.lims <- c(0,30,40,50,60,70,80)
#quantile(ex.longdata$value, seq(0,1,1/n))
batmat<- -log(hellinger_coefficients_sub5d_3000)
ex.longdata <- melt(batmat)

# quant = c(0,30,40,50,60,70,80)

quant.discreter <- function(n,data.vec,quant = "quant"){
  if (quant == "quant"){br.lims <- quantile(data.vec, seq(0,1,1/n))} else {br.lims <- quant}
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

# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value, quant = c(0,30,40,50,60,70,80))
ex.longdata$d.value <- quant.discreter(6,ex.longdata$value)
n <- 6
br.lims<- quantile(ex.longdata$value, seq(0,1,1/n))

b <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
 #scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(length(quant),"Blues")), guide = "legend",name = "Bhat distance",labels = c(0,30,40,50,60,70,80))+
   # scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[1:(n+1)],0))+
  geom_text(aes(label = round(value,0)),size = 7) +
  labs(x="models", y="models", title=paste0("Bhattacharya distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

b

dev.off()

#######################################################################
# Model evaluation ex subset 20d
#######################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(sparsebnUtils)
library(gridExtra)
library(RColorBrewer)
library(gaussDiff)

########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_cmip5_exsubset_20d_akima_cubic.rda")
# load("data/hciterations/permutations/perm1.rda")
# load("data/hciterations/permutations/5d/indTRAIN1.rda")
# load("data/hciterations/permutations/5d/indTEST1.rda")
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
cubicsets <- get(load("data/tas_historical_cmip5_exsubset_20d_akima_cubic.rda"))

namescubics <- names(cubicsets)
# shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
#                    pattern = "_r12i1p1", replacement = "")
# shortnames
############################################################################
# Function to load hciterations of models in a list 
############################################################################
loadIterations <- function(pattern,permused, hctype = NULL,it = NULL){
  if(is.null(hctype)){hctype <- ""}
  hc_interim_list <- list.files(paste0("data/hciterations/CMIP5_20d/",pattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0("data/hciterations/CMIP5_20d/",pattern,"/perm",permused,hctype))
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


###################################################################################################
# Select for each Model its optimum with respect to train and test (random permutation)
###################################################################################################
# load indices of the random Permuation with which train and test are learned
ntime <- dim(tas_historical_cmip5_exsubset_20d_akima_cubic$tas_2T_Interim$Data)[1]
nlat <- dim(tas_historical_cmip5_exsubset_20d_akima_cubic$tas_2T_Interim$Data)[2]
nlon <- dim(tas_historical_cmip5_exsubset_20d_akima_cubic$tas_2T_Interim$Data)[3]
set.seed(5)
indTRAIN1 <- sample(1:ntime, size =ntime/2)
# samplesize is amount of months (t in grid) nvars is 18 * 36
samplesize <- dim(cubicsets[[1]]$Data)[1]
nvars <- dim(cubicsets[[1]]$Data)[2]*dim(cubicsets[[1]]$Data)[3]
# indices for test is what is left
indTEST1 <- (1:samplesize)[-indTRAIN1]
# load train and test hcmodels for each (gc)m in list
hc_train_gcms <- lapply(namescubics[5], loadIterations, permused = 1, hctype = "train")
hc_test_gcms <- lapply(namescubics[5], loadIterations, permused = 1, hctype = "test")
names(hc_train_gcms) <- namescubics[5]
names(hc_test_gcms) <- namescubics[5]
# create data

data_train_gcms <- lapply(list(cubicsets[[5]]),function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE, subind = indTRAIN1)))
data_test_gcms <- lapply(list(cubicsets[[5]]),function(x) as.data.frame(TimeCoordsAnom_from_Grid_rms(x, rms = TRUE, subind = indTEST1)))
####################################################################
#
####################################################################
k <- 5
#for(k in 1:length(shortnames)){
i <- namescubics[k]
j <- namescubics[k]
hctrains <- hc_train_gcms[[i]]
hctests <- hc_test_gcms[[i]]
#loglikresults <- matrix(data = NA, nrow = length(hctrains), ncol = 3, dimnames = list(names(hctrains),c("traintrain","testtrain","testtest")))
loglikresults <- matrix(data = NA, nrow = length(hctrains), ncol = 2, dimnames = list(names(hctrains),c("traintrain","testtrain")))
fitstrain <- lapply(hctrains, bn.fit, data_train_gcms[[1]])
fitstest <- lapply(hctests, bn.fit, data_test_gcms[[1]])

loglikstrain <- sapply(X = fitstrain, logLik, data = data_train_gcms[[1]])
loglikstest <- sapply(X = fitstrain, logLik, data = data_test_gcms[[1]])
narcstrain <- sapply(X = hc_train_gcms[[i]], narcs)
#loglikstesttest <- sapply(X = fitstest, logLik, data = data_test_gcms[[j]])

loglikresults[,"traintrain"] <- loglikstrain
loglikresults[,"testtrain"] <- loglikstest
#loglikresults[,"testtest"] <- loglikstesttest

assign(paste0("loglik_traintest_",i), loglikresults)
# save(list = paste0("loglik_traintest_",i), 
# file = paste0("results/loglik_traintest/loglik_traintest_",i,".rda"))

loglikstrain <- loglikstest <- fitstrain <- fitstest <- hctrains <- hctests <- NULL
# }
gc()
loglik_traintest_tas_2T_interim
plot(narcstrain,loglik_traintest_tas_2T_Interim[,1], 
     #xaxt = "n", 
     xlab = c("|E|"), ylab = c("Loglik(x|model)"), 
     col = rainbow(10)[10],
     #ylim = c(-600000,max(loglik_optimum_datasets)),
     pch = 16)


points(narcstrain,loglik_traintest_tas_2T_Interim[,2], col = rainbow(10)[9],pch = 16)
##############################################################
#
##############################################################
##################################################################
# Make al datasets standardized
##################################################################
permk <- perm1
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
assign(paste0("data_gcms_anom_scaled_1"),lapply(data_gcms_anom_scaled, function(x) x[,permk]))
length(data_gcms_anom_scaled)

it <- "400_500"
hc_gcms5sub <- lapply(paste0(namescubics), loadIterations, permused = permk, it = it)
modelsize <- 5
selection_hc_gcms <-  mapply(function(x,y) x[[grep(y,x)]], x = hc_gcms5sub, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
#selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
########################################################################################
# Calculate hellinger distance bayesian networks
# KL distances give NaN (because there are x such that P(x) = 0 does not imply Q(x) = 0)
########################################################################################
fits_NEL <- lapply(selection_fits_gcms_scaled,bnlearn:::as.graphNEL)
edgelistssparse <- lapply(X =fits_NEL,FUN = edgeList)
data_gcms_anom_scaled_sparse <- lapply(data_gcms_anom_scaled_1,sparsebnData, type = "continuous")
fits_COVS <- mapply(get.covariance, x = edgelistssparse, data = data_gcms_anom_scaled_sparse)
fits_COVmats <- lapply(fits_COVS, as.matrix)
fits_Means <- lapply(data_gcms_anom_scaled, function(x) colMeans(x[perm1]) )

hellinger_coefficients <- matrix(data = NA, nrow = length(fits_COVmats), ncol = length(fits_COVmats), dimnames = list(names(fits_COVmats),names(fits_COVmats)))
i <- 1
for(i in 1:length(fits_COVmats)){
  kl <- mapply(function(b,d) normdiff(mu1 =fits_Means[[i]], mu2 = b, sigma1 = fits_COVmats[[i]], sigma2 = d, method = "Hellinger"),b = fits_Means, d = fits_COVmats)
  hellinger_coefficients[i,] <- kl
}
str(fits_Means)

rownames(hellinger_coefficients) <- colnames(hellinger_coefficients) <-names(fits_Means)
hellinger_coefficients_sub20d_500<- hellinger_coefficients
save(hellinger_coefficients_sub20d_500, file ="results/hellinger_coefficient/hellinger_coefficients_sub20d_500.rda")
####################################################################################
#
####################################################################################
load(file ="results/hellinger_coefficient/hellinger_coefficients_sub20d_500.rda")


########################################################
# Visualize hellinger distance example set 
########################################################
#ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")

plotname <- paste0("figs/hellinger_CMIP5_ordering/Bhat_CMIP5_ex_ACCESS|BNU|CMCC.CMS|interim_20d_500.png")
png(plotname,units = "in", width = 9, height = 7,  res =180)
# plotname <- paste0("figs/hellinger_CMIP5_ordering/Hellinger_CMIP5_ex_ACCESS|BNU|CMCC.CMS|interim.pdf")
# pdf(plotname, width = 9, height =7)
# n <-7
# n <-7
# br.lims <- c(0,30,40,50,60,70,80)
#quantile(ex.longdata$value, seq(0,1,1/n))
batmat<- -log(hellinger_coefficients_sub20d_500)
ex.longdata <- melt(batmat)

# quant = c(0,30,40,50,60,70,80)

quant.discreter <- function(n,data.vec,quant = "quant"){
  if (quant == "quant"){br.lims <- quantile(data.vec, seq(0,1,1/n))} else {br.lims <- quant}
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

# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value, quant = c(0,30,40,50,60,70,80))
ex.longdata$d.value <- quant.discreter(6,ex.longdata$value)
n <- 6
br.lims<- quantile(ex.longdata$value, seq(0,1,1/n))

b <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  #scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(length(quant),"Blues")), guide = "legend",name = "Bhat distance",labels = c(0,30,40,50,60,70,80))+
  # scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[1:(n+1)],0))+
  geom_text(aes(label = round(value,2)),size = 7) +
  labs(x="models", y="models", title=paste0("Bhattacharya distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

b

dev.off()


#######################################################################
# Model evaluation CSIRO different runs same permutation
#######################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(sparsebnUtils)
library(gridExtra)
library(RColorBrewer)
library(gaussDiff)
########################################################################################################
# Model Evaluation CSIRO
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
load("../Data/Struct_learn/permutations.rda")
load("data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
cubicsets <- get(load("data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda"))
namescubics <- gsub("-",".",gsub(".ncml","",gsub("https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip5/historical/historical_day_","",sapply(cubicsets,function(x)attr(x,"dataset")))))

# shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
#                    pattern = "_r12i1p1", replacement = "")
# shortnames
############################################################################
# Function to load hciterations of models in a list 
############################################################################
loadIterations <- function(pattern,permused, hctype = NULL,it = NULL){
  if(is.null(hctype)){hctype <- ""}
  hc_interim_list <- list.files(paste0("data/hciterations/CMIP5_CSIRO/",pattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0("data/hciterations/CMIP5_CSIRO/",pattern,"/perm",permused,hctype))
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
assign(paste0("data_gcms_anom_scaled_perm",permk),lapply(data_gcms_anom_scaled, function(x) x[,permutations[[permk]]]))
length(data_gcms_anom_scaled)

it <- "1700_1800"
hc_gcmsCSIRO <- lapply(paste0(namescubics), loadIterations, permused = permk, it = it)
modelsize <- 18
selection_hc_gcms <-  mapply(function(x,y) x[[grep(y,x)]], x = hc_gcmsCSIRO, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
#selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
########################################################################################
# Calculate hellinger distance bayesian networks
# KL distances give NaN (because there are x such that P(x) = 0 does not imply Q(x) = 0)
########################################################################################
fits_NEL <- lapply(selection_fits_gcms_scaled,bnlearn:::as.graphNEL)
edgelistssparse <- lapply(X =fits_NEL,FUN = edgeList)
data_gcms_anom_scaled_sparse <- lapply(data_gcms_anom_scaled_perm3,sparsebnData, type = "continuous")
fits_COVS <- mapply(get.covariance, x = edgelistssparse, data = data_gcms_anom_scaled_sparse)
fits_COVmats <- lapply(fits_COVS, as.matrix)
fits_Means <- lapply(data_gcms_anom_scaled, function(x) colMeans(x[permutations[[permk]]]) )

hellinger_coefficients <- matrix(data = NA, nrow = length(fits_COVmats), ncol = length(fits_COVmats), dimnames = list(names(fits_COVmats),names(fits_COVmats)))
i <- 1
for(i in 1:length(fits_COVmats)){
  kl <- mapply(function(b,d) normdiff(mu1 =fits_Means[[i]], mu2 = b, sigma1 = fits_COVmats[[i]], sigma2 = d, method = "Hellinger"),b = fits_Means, d = fits_COVmats)
  hellinger_coefficients[i,] <- kl
}
str(fits_Means)

rownames(hellinger_coefficients) <- colnames(hellinger_coefficients) <-namescubics
hellinger_coefficients_CSIRO<- hellinger_coefficients
save(hellinger_coefficients_CSIRO, file ="results/hellinger_coefficient/hellinger_coefficients_CSIRO.rda")


#######################################################################
# Model evaluation CSIRO different runs different permutations NEED BACKPERMUTATIONS

#######################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(sparsebnUtils)
library(gridExtra)
library(RColorBrewer)
library(gaussDiff)
########################################################################################################
# Model Evaluation CSIRO
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
load("../Data/Struct_learn/permutations.rda")
load("data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
cubicsets <- get(load("data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda"))
namescubics <- gsub("-",".",gsub(".ncml","",gsub("https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip5/historical/historical_day_","",sapply(cubicsets,function(x)attr(x,"dataset")))))

# shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
#                    pattern = "_r12i1p1", replacement = "")
# shortnames
############################################################################
# Function to load hciterations of models in a list 
############################################################################
loadIterations <- function(pattern,permused, hctype = NULL,it = NULL){
  if(is.null(hctype)){hctype <- ""}
  hc_interim_list <- list.files(paste0("data/hciterations/CMIP5_CSIRO/",pattern,"/perm",permused,hctype), full.names = T)
  hc_interim_names <- list.files(paste0("data/hciterations/CMIP5_CSIRO/",pattern,"/perm",permused,hctype))
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
##################################################################
# Make al datasets standardized
##################################################################
for(permk in 1:5){
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
assign(paste0("data_gcms_anom_scaled_perm",permk),lapply(data_gcms_anom_scaled, function(x) x[,permutations[[permk]]]))
length(data_gcms_anom_scaled)
}

it <- "1700_1800"
hc_gcmsCSIRO <- lapply(1:5,function(x)lapply(paste0(namescubics), loadIterations, permused = x, it = it))
modelsize <- 18
selection_hc_gcms <- lapply(
  hc_gcmsCSIRO, function(x) mapply(function(z,y) x[[z]][[grep(y,x[[z]])]], z = 1:5, y = modelsize, SIMPLIFY = FALSE))
# Make fits
# unscaled fit (uses data_gcms_out)
#selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- lapply(selection_hc_gcms,function(x) mapply(function (z,y) bn.fit(x = x[[z]], data = y), z = 1:5, y = data_gcms_anom_scaled, SIMPLIFY = FALSE))
########################################################################################
# Calculate hellinger distance bayesian networks
# KL distances give NaN (because there are x such that P(x) = 0 does not imply Q(x) = 0)
########################################################################################
selection_fits_gcms_scaled[[1]][[10]]
fits_NEL <- lapply(1:5, function(x)lapply(selection_fits_gcms_scaled[[x]],bnlearn:::as.graphNEL))
edgelistssparse <- lapply(1:5,function(x)lapply(X =fits_NEL[[x]],FUN = edgeList))
data_gcms_anom_scaled_sparse <- lapply(1:5,function(x) lapply(eval(get(paste0("data_gcms_anom_scaled_perm",x))),sparsebnData, type = "continuous"))
fits_COVS <- lapply(1:5,function(x) mapply(get.covariance, x = edgelistssparse[[x]], data = data_gcms_anom_scaled_sparse[[x]]))
fits_COVmats <- lapply(1:5,function(x)lapply(fits_COVS[[x]], as.matrix))
fits_Means <- lapply(1:5,function(z)lapply(data_gcms_anom_scaled, function(x) colMeans(x[permutations[[z]]]) ))

run <- 1

for (run in 1:10){
hellinger_coefficients <- matrix(data = NA, nrow = length(fits_COVmats), ncol = length(fits_COVmats))
i <- 1
for(i in 1:length(fits_COVmats)){
  kl <- mapply(function(b,d) normdiff(mu1 =fits_Means[[i]][[run]], mu2 = b, sigma1 = fits_COVmats[[i]][[run]], sigma2 = d, method = "Hellinger"),b = lapply(fits_Means,function(x)x[[run]]), d = lapply(fits_COVmats,function(x)x[[run]]))
  hellinger_coefficients[i,] <- kl
}
rownames(hellinger_coefficients) <- colnames(hellinger_coefficients) <-assign(paste("perm",1:5))

assign(paste0("hellinger_coefficients_run",run), value = hellinger_coefficients)
}
-log(hellinger_coefficients)
# Need backpermutations. 

save(hellinger_coefficients_CSIRO, file ="results/hellinger_coefficient/hellinger_coefficients_CSIRO.rda")

hellinger_coefficients_run1
