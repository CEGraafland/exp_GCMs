#############################################################################
# Model Evaluation CMIP5 rcp 85
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
library(RColorBrewer)
library(gaussDiff)
########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
load("data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda")
load("data/tas_rcp85_cmip5_left_10d_akima_cubic.rda")
load("data/tas_ssp585_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_ssp585_cmip6_left_10d_akima_cubic_corrected.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models
###########################################################################################
cubicsets5.rcp85 <- tas_rcp85_cmip5_10d_akima_cubic_corrected
# Leave out bcc "CMIP5_bcc.csm1.1.m_r1i1p1_rcp85" because does not reach 2100 (one year atrasado)
xout <- which(names(tas_rcp85_cmip5_left_10d_akima_cubic) == "CMIP5_bcc.csm1.1.m_r1i1p1_rcp85")
cubicsets5left.rcp85 <- tas_rcp85_cmip5_left_10d_akima_cubic[-xout]

cubicsets6.ssp585 <- tas_ssp585_cmip6_10d_akima_cubic_corrected
cubicsets6left.ssp585 <- tas_ssp585_cmip6_left_10d_akima_cubic_corrected

cubicsets.future <- c(cubicsets6.ssp585,cubicsets6left.ssp585)
cubicsets.future <- c(cubicsets5.rcp85,cubicsets5left.rcp85)
names(cubicsets.future)

namescubics5.rcp85 <- names(cubicsets5.rcp85)
namescubics5left.rcp85 <- names(cubicsets5left.rcp85)
namescubics6.ssp585 <- names(cubicsets6.ssp585)
namescubics6left.ssp585 <- names(cubicsets6left.ssp585)

shortnames.future <- gsub("CMIP5_",gsub(names(cubicsets.future), pattern = "_r1i1p1", replacement = ""),replacement = "")

#,
#pattern = "_r12i1p1", replacement = "")
# shortnames2 <- names(cubicsets)
abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames.future,replacement = ""),replacement = "")

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
data.future_gcms_out <- lapply(cubicsets.future,function(x) TimeCoordsAnom_from_Grid_rms(x, rms = TRUE))
data.future_gcms_out_df <- lapply(data.future_gcms_out, function(x) as.data.frame(x))
data.future_gcms_anom_out <- lapply(data.future_gcms_out, function(x) mat2Dto3Darray(x,attr(x,"Xcoords"), attr(x,"Ycoords")))
grid.future_gcms_anom <- mapply(function(x,y) {x$Data <- y ;return(x)}, x = cubicsets.future, y = data.future_gcms_anom_out,SIMPLIFY = FALSE)
# Global standarazation over the above: FIELD
grid.future_gcms_anom_scaled <- lapply(grid.future_gcms_anom, 
                                function(x) scaleGrid(x,
                                                      spatial.frame ="field",
                                                      type = "standardize"))

# Extract the FIELD standarized data
data.future_gcms_anom_scaled <- lapply(grid.future_gcms_anom_scaled,
                                function(x) { x <- redim(x,drop = TRUE)
                                xdata <- array3Dto2Dmat(x$Data) 
                                return(xdata)}) 
data.future_gcms_anom_scaled <- lapply(data.future_gcms_anom_scaled, function(x) as.data.frame(x))
# adapt to k-th permutation
assign(paste0("data.future_gcms_anom_scaled_",permk),lapply(data.future_gcms_anom_scaled, function(x) x[,permutations[[permk]]]))

it <- "1700_1800"
hc_gcms5.rcp85 <- lapply(paste0("FUTURE_",namescubics5.rcp85), loadIterations, permused = permk, it = it)
hc_gcms5left.rcp85<- lapply(paste0("FUTURE_CMIP5_left/FUTURE_",namescubics5left.rcp85), loadIterations, permused = permk, it = it)
hc_gcms6.ssp585 <- lapply(paste0("FUTURE_CMIP6/FUTURE_",namescubics6.ssp585), loadIterations, permused = permk, it = it)
hc_gcms6left.ssp585 <- lapply(paste0("FUTURE_CMIP6_left/FUTURE_",namescubics6left.ssp585), loadIterations, permused = permk, it = it)
# hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk, it = it)
# hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk, it = it)

hc_gcms.future <- c(hc_gcms5.rcp85,hc_gcms5left.rcp85)
# hc_gcms.future <- c(hc_gcms6.ssp585,hc_gcms6left.ssp585)

names(hc_gcms.future) <- names(cubicsets.future)
# names(hc_gcms.future) <- shortnames.future


# Choose between 'own optimums' or constant magnitude
# modelsize <- own_optimums
modelsize <- 18
selection_hc_gcms.future <-  mapply(function(x,y) x[[grep((y),x)]], x = hc_gcms.future, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms.future <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms.future, y = data.future_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms.future_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms.future, y = data.future_gcms_anom_scaled, SIMPLIFY = FALSE)
##########################################################################
# Global dataset evaluation from BN 
##########################################################################
whichBNFits <- selection_fits_gcms.future_scaled
whichData <- data.future_gcms_anom_scaled
# whichSet <- "interim_10d_akima_cubic"
# whichSet <- "ncep_10d"
# whichSet <- "JRA55_10d_akima_cubic"
#########################################################################
# pickedEOFsgrids is prepared under EOFs logliks preparation. 
#####################################################################
# EOFS_preparation <- function(whichPCs){
#   selected_EOFs <- EOFs[,whichPCs]
#   selected_PCs <- PCs[,whichPCs]
#   if(is.null(dim(selected_PCs))){
#     dim(selected_PCs) <- c(length(selected_PCs),1)
#     dim(selected_EOFs) <- c(length(selected_EOFs),1)
#   }
#   
#   projection_EOFs <- matrix(data = 0,nrow = nrow(selected_PCs),ncol = nrow(selected_EOFs))
#   for(i in 1:nrow(selected_PCs)){
#     for(j in 1:ncol(selected_PCs)){
#       projection_EOFs[i,] <- projection_EOFs[i,] + selected_EOFs[,j]*selected_PCs[i,j]
#     }
#   }
#   return(projection_EOFs)
# }
# 
# gridforEOF <- grid_gcms_anom_scaled
# prinCompSet <- prinComp(gridforEOF[[whichSet]], n.eofs = 360)
# expvars <- attr(prinCompSet[[1]][[1]],"explained_variance")
# EOFs <- prinCompSet[[1]][[1]]$EOFs
# PCs <- prinCompSet[[1]][[1]]$PCs
# 
# varpercentages <- c(0.25,0.50,0.75,0.90,0.95,0.99)
# # varpercentages <- c(0.50,0.75,0.825,0.90,0.95,0.99)
# pickedEOFs <- lapply(varpercentages, FUN = function(x) which(expvars > x)[[1]])
# pickedEOFsdata <- lapply(pickedEOFs, FUN = function(x) EOFS_preparation(1:x))
# 
# templatedata <- gridforEOF[[whichSet]]
# templatesdata <- rep(list(templatedata),length(varpercentages))
# pickedEOFsgrids <- mapply(function(x,y) {z<- x; z$Data[1,,,] <- y; return(z)},templatesdata, pickedEOFsdata, SIMPLIFY  = FALSE)
# pickedEOFsnames <- paste0(whichSet,"_v.exp_",as.character(varpercentages))
# names(pickedEOFsgrids) <- names(pickedEOFsdata) <- pickedEOFsnames
# 
# data_pickedEOFs <- lapply(pickedEOFsdata, function(x) as.data.frame(x))
# whichData <- c(data_gcms_anom_scaled,data_pickedEOFs)

ordersets <- TRUE

loglik_selection_datasets.future <- matrix(data = NA, nrow = length(whichBNFits), ncol = length(whichData), dimnames = list(names(whichBNFits),names(whichData)))
for(i in 1:length(whichBNFits)){
  lo <- sapply(X = whichData, logLik, object = whichBNFits[[i]])
  loglik_selection_datasets.future[i,] <- lo
}

institutions.future <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","ACCESS","BCC","BNU","CCSM4","CESM2","CMCC","CSIRO","inmcm4","MRI","UKESM1","NESM3")
institutions.future <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","ACCESS","bcc","BNU","CCSM4","CESM1","CMCC","CSIRO","inmcm4","MRI")
dim(institutions.future)<- c(length(institutions.future),1)


# institutions <- rep(institutions.rcp85,each =1)
# CMIPS.rcp85 <- rep(c("CMIP5"),length(institutions.rcp85))

# combinations.rcp85 <- cbind(CMIPS.rcp85,institutions.rcp85)
indearth.rcp85 <- grep("EC.EARTH_r1i|EC.EARTH_r2i",rownames(loglik_selection_datasets.future))

# namesort.rcp85 <- unlist(apply(combinations.rcp85,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(loglik_selection_datasets.rcp85[-indearth.rcp85,-indearth.rcp85]))))
namesort.future <- unlist(apply(institutions.future,MARGIN = 1,function(x) grep(x,rownames(loglik_selection_datasets.future[-indearth.rcp85,-indearth.rcp85]))))
namesort.future <- unlist(apply(institutions.future,MARGIN = 1,function(x) grep(x,rownames(loglik_selection_datasets.future))))
shortnames.future[-indearth.rcp85][namesort.future]


# orderby <- "interim_10d_akima_cubic"
# orderby <- "CMIP5_HadGEM2.ES_r1i1p1"
# orderby <- "ncep_10d"
# orderby <- "ncep_10d_v.exp_0.99"
# orderby <- "ncep_10d_v.exp_0.9"
# orderby <- "interim_10d_akima_cubic_v.exp_0.95"
# sort(colSums(loglik_selection_datasets[-indearth,-indearth]))
# sort(rowSums(loglik_selection_datasets[-indearth,-indearth]))
# sort(colSums(loglik_selection_datasets[-indearth,-indearth]) + rowSums(loglik_selection_datasets[-indearth,-indearth]))
# 
# logsort <- sort(loglik_selection_datasets[-indearth,-indearth][,orderby],index.return = TRUE)
# logclust <-hclust(dist(loglik_selection_datasets[-indearth,-indearth], method = "minkowski", p = 2),method = "complete")
# logclust <-hclust(dist(t(loglik_selection_datasets[-indearth,-indearth]), method = "minkowski", p = 2),method = "complete")


# loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][rev(logclust$order),rev(logclust$order)]
loglik_selection_datasets.future_ordered <- loglik_selection_datasets.future[-indearth.rcp85,-indearth.rcp85][namesort.future,namesort.future]
loglik_selection_datasets.future_ordered <- loglik_selection_datasets.future[namesort.future,namesort.future]
# loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][logsort$ix,logsort$ix]
# loglik_selection_datasets_ordered <- loglik_selection_datasets[logsort$ix,namesort2]
# loglik_selection_datasets_ordered <- loglik_selection_datasets[logsort$ix,c(pickedEOFsnames,whichSet)]
# with hellinger distance ordering:
# loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][mclust$labels[mclust$order],mclust$labels[mclust$order]]
# with hellinger distance and only eofs
# loglik_selection_datasets_ordered <- loglik_selection_datasets[-indearth,-indearth][mclust$labels[mclust$order],c(pickedEOFsnames,whichSet)]
ordersets <- TRUE
if (ordersets == TRUE){ longdata <- melt(loglik_selection_datasets.future_ordered)
} else {longdata <- melt(loglik_selection_datasets)}

########################################################
# Visualize logliks future set
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
########################################################################################
fits_NEL.future <- lapply(selection_fits_gcms.future_scaled,bnlearn:::as.graphNEL)
edgelistssparse.future <- lapply(fits_NEL.future,edgeList)
data.future_gcms_anom_scaled_sparse <- lapply(data.future_gcms_anom_scaled,sparsebnData, type = "continuous")
fits_COVS.future <- mapply(get.covariance, x = edgelistssparse.future, data = data.future_gcms_anom_scaled_sparse)
fits_COVmats.future <- lapply(fits_COVS.future, as.matrix)
fits_Means.future <- lapply(data.future_gcms_anom_scaled, function(x) colMeans(x[permutations[[3]]]) )

hellinger_coefficients.future <- matrix(data = NA, nrow = length(fits_COVmats.future), ncol = length(fits_COVmats.future), dimnames = list(names(fits_COVmats.future),names(fits_COVmats.future)))
i <- 1
for(i in 1:length(fits_COVmats.future)){
  kl <- mapply(function(b,d) normdiff(mu1 =fits_Means.future[[i]], mu2 = b, sigma1 = fits_COVmats.future[[i]], sigma2 = d, method = "Hellinger"),b = fits_Means.future, d = fits_COVmats.future)
  hellinger_coefficients.future[i,] <- kl
}
hellinger_coefficients.future

hellinger_coefficients.rcp85 <- hellinger_coefficients.future
save(hellinger_coefficients.rcp85, file ="results/hellinger_coefficient/hellinger_coefficients.rcp85.rda")
hellinger_coefficients.ssp585 <-hellinger_coefficients.future
save(hellinger_coefficients.ssp585, file ="results/hellinger_coefficient/hellinger_coefficients.ssp585.rda")
#########################################################################
# Analysis hellinger distances
#########################################################################
# Partitional clustering
library("cluster")
library("factoextra")
load(file ="results/hellinger_coefficient/hellinger_coefficients.rcp85.rda")
load(file ="results/hellinger_coefficient/hellinger_coefficients.ssp585.rda")

hellinger_coefficients.rcp85
earthgrep <- unlist(apply(matrix(c("EC.EARTH_rcp85","EC.EARTH_r2i1p1_rcp85","bcc.csm1.1.m_rcp85"),nrow =3, ncol = 1),MARGIN = 1,function(x) grep(x,rownames(hellinger_coefficients.rcp85))))
hel_coef_without_earth.rcp85 <- hellinger_coefficients.rcp85[-earthgrep,-earthgrep]
namesorthc.rcp85 <- unlist(apply(institutions.future,MARGIN = 1,function(x) grep(x,rownames(hel_coef_without_earth.rcp85))))
selecteddists.rcp85 <- hel_coef_without_earth.rcp85[namesorthc.rcp85,namesorthc.rcp85]
nrow(selecteddists.rcp85)

namesorthc.ssp585 <- unlist(apply(institutions.future,MARGIN = 1,function(x) grep(x,rownames(hellinger_coefficients.ssp585))))
selecteddists.ssp585 <- hellinger_coefficients.ssp585[namesorthc.ssp585,namesorthc.ssp585]

selecteddists.future <- selecteddists.ssp585
selecteddists.future <- selecteddists.rcp85

klsort <- sort(selecteddists.future,index.return = TRUE)
# select optimal amount of clusters
fviz_nbclust(-log(selecteddists.future), 
             FUNcluster = cluster::pam, 
             method = "silhouette", 
             k.max =25,medoids = NULL)
fviz_nbclust(-log(selecteddists.future),FUNcluster = cluster::pam,method = "gap_stat",k.max = 10)
# Use PAM to perform k-medoids
pamcluster <- pam(-log(selecteddists.future),k =2, 
                  medoids = c(1,2), 
                  do.swap = FALSE)
# info
names(selecteddists.future)
pamcluster$medoids
pamcluster$id.med
# visualize k-medoids
par(mar = c(5,10,3,4)+0.1)
plot(silhouette(pamcluster),max.strlen = 200, nmax.lab =40,cex = 0.7)

fviz_cluster(pamcluster, repel = TRUE, axes = c(1,2),show.clust.cent = TRUE)
###############################################################
# Hierarchical clustering 
library(proxy)
mclust.future <-hclust(dist(-log(selecteddists.future),method = "minkowski", p = 2),method = "complete")
plot(mclust.future, hang = -1)
rect.hclust(mclust.future,k=11)
dev.off()

# ordersets <- TRUE
# if (ordersets == TRUE){ longdata <- melt(-log(selecteddists)[ordersort$ix,ordersort$ix])
# } else {longdata <- melt(-log(selecteddists))}

orderedmat.future<- -log(selecteddists.future)[mclust.future$order,mclust.future$order]


##############################################################################
# Visualize  Hellinger Distance global set full matrix
##############################################################################
code <- FALSE
if(code == TRUE) {r.i. <- sapply(rownames(orderedmat),function(x)which(x == codes$rn.orderedmat))
toChange <- as.logical(sapply(r.i., length))
rownames(orderedmat)[toChange]<- codes$abrev.y.atm[unlist(r.i.[toChange])]
colnames(orderedmat)[toChange]<- codes$abrev.y.atm[unlist(r.i.[toChange])]}


longdata.future <- melt(orderedmat.future)
n <- 7
br.lims.future <- quantile(longdata.future$value, seq(0,1,1/n))
data.future.vec <- longdata.future$value

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

longdata.future$d.value <- quant.discreter(n,longdata.future$value)

b <- ggplot(longdata.future, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims.future[2:(n+1)],0))+
  geom_text(aes(label = round(value,0))) +
  labs(x="data", y="models", title=paste0("hellinger distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=12, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=12),
        plot.title=element_text(size=12)) +
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_wrt_quantileclustering_withoutearth_rcp85.pdf"
pdf(plotname,width = 10, height = 10)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_wrt_quantileclustering_withoutearth_rcp85.png"
png(plotname, width = 10, height = 10, units = "in", res = 180)
b
dev.off()

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP6_ordering/hellinger_wrt_quantileclustering_ssp585.pdf"
pdf(plotname,width = 10, height = 10)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP6_ordering/hellinger_wrt_quantileclustering_ssp585.png"
png(plotname, width = 10, height = 10, units = "in", res = 180)
b
dev.off()
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
rownames(loglik_selection_datasets.future)
colnames(loglik_selection_datasets.future)
#subr <-grep("ACCESS|BNU|interim",rownames(loglik_selection_datasets),fixed = FALSE)
#subc <- grep("ACCESS|BNU|interim",colnames(loglik_selection_datasets),fixed = FALSE)
r.ex.rcp85.subset <- c("BNU.ESM_rcp85","ACCESS1.0_rcp85","ACCESS1.3_rcp85","CMCC.CMS_rcp85")
c.ex.rcp85.subset <- c("CMIP5_BNU.ESM_r1i1p1_rcp85","CMIP5_ACCESS1.0_r1i1p1_rcp85","CMIP5_ACCESS1.3_r1i1p1_rcp85","CMIP5_CMCC.CMS_r1i1p1_rcp85")
n <-8

ex.rcp85.logliks <- loglik_selection_datasets.future[r.ex.rcp85.subset,c.ex.rcp85.subset]
ex.rcp85.longdata <- melt(ex.rcp85.logliks)
ex.rcp85.longdata$d.value <- quant.discreter(n,ex.rcp85.longdata$value)
br.lims.rcp85 <- quantile(ex.rcp85.longdata$value, seq(0,1,1/n))

b <- ggplot(ex.rcp85.longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value),interpolate = FALSE) + 
  geom_text(aes(label = round(-value/10^4,0)), size = 7) +
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Reds")), guide = "legend",name = "log P(data|model)\nx -10^4",n.breaks =n,labels = round(-br.lims.rcp85[2:(n+1)]/10^4,0))+
  labs(x="data", y="models", title=paste0(modelsize*100," model-dataset evaluation")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=15, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=15),
        plot.title=element_text(size=12)) + 
  scale_x_discrete(position = "top")

b
dev.off()




########################################
# CMIP5 CMIP6 only models for comparison 
########################################
library("cluster")
library("factoextra")
load(file ="results/hellinger_coefficient/hellinger_coefficients.ssp585.rda")

out <- which(rownames(hellinger_coefficients.ssp585) == "CMIP6Amon_CanESM5_ssp585_r1i1p2f1")
hellinger_coefficients.ssp585 <-hellinger_coefficients.ssp585[-out,-out]
institutions.overlap <- c("IPSL","CNRM","CESM","EC","CanESM5","MRI","MIROC","MPI","GFDL","Nor")
institutions <- rep(institutions.overlap,each =1)
CMIPS.6 <- rep(c("CMIP6"),length(institutions.overlap))
combinations.6 <- cbind(CMIPS.6,institutions)
namesort.6 <- unlist(apply(combinations.6,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(hellinger_coefficients.ssp585))))
selecteddists <- hellinger_coefficients.ssp585
orderedmat<- -log(selecteddists)[namesort.6,namesort.6]

# namesorthc.ssp585 <- unlist(apply(institutions.overlap,MARGIN = 1,function(x) grep(x,rownames(hellinger_coefficients.ssp585))))
# selecteddists.ssp585 <- hellinger_coefficients.ssp585[namesorthc.ssp585,namesorthc.ssp585]

orderedmat<- -log(selecteddists)[namesort.6,namesort.6]
longdata <- melt(orderedmat)
rownames(orderedmat) <- gsub("CMIP6Amon_","",gsub("f1","",rownames(orderedmat)))
colnames(orderedmat) <- gsub("CMIP6Amon_","",gsub("f1","",colnames(orderedmat)))

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

# cmip6
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP6_ordering/hellinger_CMIP6_ssp585_overlap.pdf"
pdf(plotname,width = 12, height = 10)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP6_ordering/hellinger_CMIP6_ssp585_overlap.png"
png(plotname, width = 12, height = 10, units = "in", res = 180)

b             
dev.off()


########################################
# CMIP5 CMIP6 only models for comparison: CMIP5
########################################
library("cluster")
library("factoextra")
load(file ="results/hellinger_coefficient/hellinger_coefficients.rcp85.rda")


institutions.overlap <- c("IPSL.CM5A","CESM1","EC.EARTH_r12i1p1","Can","MRI","MIROC.ES","MIROC5","CNRM","MPI","GFDL.ESM","Nor")
hellinger_coefficients.rcp85


# cmip5
institutions <- rep(institutions.overlap,each =1)
CMIPS.5 <- rep(c("CMIP5"),length(institutions.overlap))
combinations.5 <- cbind(CMIPS.5,institutions)
namesort.5 <- unlist(apply(combinations.5,MARGIN = 1,function(x) grep(paste0(x[2]),rownames(hellinger_coefficients.rcp85))))
selecteddists <- hellinger_coefficients.rcp85
orderedmat<- -log(selecteddists)[namesort.5,namesort.5]
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

# cmip5
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_CMIP5_rcp85_overlap.pdf"
pdf(plotname,width = 12, height = 10)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_CMIP5_ordering/hellinger_CMIP5_rcp85_overlap.png"
png(plotname, width = 12, height = 10, units = "in", res = 180)

b             
dev.off()



