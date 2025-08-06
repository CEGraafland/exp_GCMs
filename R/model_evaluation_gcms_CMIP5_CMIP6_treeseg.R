#############################################################################
# Model Evaluation CMIP5 CMIP6 Treeseg
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
library(sparsebnUtils)
library(RColorBrewer)
library(gaussDiff)
library(visualizeR)
########################################################################################################
# Model Evaluation CMIP5/CMIP6
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
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


cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsets5.rcp85,cubicsets5left.rcp85,cubicsets6all, cubicsets6all.ssp585,listinterim,listjra55,listncep)
names(cubicsets)<- gsub(names(cubicsets), pattern = "_historical",replacement ="") 


namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics5.rcp85 <- names(cubicsets5.rcp85)
namescubics5left.rcp85 <- names(cubicsets5left.rcp85)

#namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
#namescubics6left <- gsub(names(cubicsets6left), pattern = "_historical",replacement ="") 
namescubics6all <- gsub(names(cubicsets6all), pattern = "_historical",replacement ="")
#namescubics6.ssp585 <- names(cubicsets6.ssp585)
#namescubics6left.ssp585 <- names(cubicsets6left.ssp585)
namescubics6all.ssp585 <- names(cubicsets6all.ssp585)


namescubicsearth <- names(cubicsetsearth)
namescubics <- names(cubicsets)
#shortnames <- #gsub(
#  gsub(names(cubicsets), pattern = "_r1i1p1", replacement = "")
#,
#pattern = "_r12i1p1", replacement = "")
#shortnames2 <- names(cubicsets)
#abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")




#############################################################################################
# Emergent constraints: (1) Hellinger distance vs global mean temperature increase
#############################################################################################
names(cubicsets)


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

keltoC <-function(x)(x)%%273

gt.hists.cmips <- lapply(sets.hists,gtimes)
gt.hists.cmips <- lapply(gt.hists.cmips,keltoC)
gmt.hists.cmips2 <- sapply(gt.hists.cmips,mean)
gsdt.hists.cmips <- sapply(gt.hists.cmips,sd) # sd with respect to timeseries, in 30 years.



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
gmt.hists.cmips <- sapply(clim.hists,gmt)%>%keltoC
all.equal(gmt.hists.cmips,gmt.hists.cmips2)
#gmt.sd.futures.cmips <- sapply(clim.sd.futures,gmt)
gmt.futures.cmips <- sapply(clim.futures,gmt)
gmt.futures.cmips <- keltoC(gmt.futures.cmips)
all.equal(gmt.futures.cmips,gmt.futures.cmips2)
#gmt.sd.reans <- sapply(clim.sd.reans,gmt)
gmt.reans <- sapply(clim.reans,gmt)
gmt.reans <- keltoC(gmt.reans)
all.equal(gmt.reans,gmt.reans2)





########################################
# Library TreeSeg
########################################
load("results/hellinger_coefficient/hellinger_coefficients.rda")
load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_CMIP6_hist_vs_fut.rda")
selecteddists <- hellinger_coefficients
#orderedmat<- -log(selecteddists)
#longdata <- melt(orderedmat)
rownames(selecteddists) <- gsub("_historical","",rownames(selecteddists))
colnames(selecteddists) <- gsub("_historical","",colnames(selecteddists))

names(sets.hists)
rownames(selecteddists)
colnames(selecteddists)

names(sets.hists)[!names(sets.hists)%in%rownames(selecteddists)]
names(sets.futures)[!names(sets.futures)%in%rownames(selecteddists)]

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

#names(obs)<- rownames(fut.dists)

# For observations are sd fut mean tempeartures
obs <- gsdt.futures.cmips
length(obs)
names(obs)<- rownames(hist.dists)


###############################################################
# Hierarchical clustering batachari
###############################################################
library(proxy)
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

as.dist(-log(hist.dists))
mclust <- hclust(as.dist(-log(hist.dists)),method ="complete")
mclust <- hclust(as.dist(-log(hist.dists)),method ="single")
mclust <- hclust(as.dist(-log(hist.dists)),method ="ward.D")

#For dendogram made of Future BNs
mclust <-hclust(dist(-log(fut.dists),method = "minkowski", p = 2),method = "complete")
#mclust <-hclust(dist(-log(fut.dists),method = "minkowski", p = 2),method = "single")
#mclust <-hclust(dist(-log(fut.dists),method = "minkowski", p = 2),method = "ward.D")


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
tS<-treeSeg(obs,phy,checkOrder = FALSE, fam = "gauss", alpha = 0.05)

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

names(obs) <- row.names(fut.dists)
plot(obs[tiporder],xlab="", ylab = "Mean Temp in 2071-2100 °C", xaxt ="n", col =node.vec[1:length(phy$tip.label)][mclust$order],pch = 16)
xlabels=tiporder
axis(3, at=1:length(tiporder),labels = xlabels, col.axis= "red", las=2)

