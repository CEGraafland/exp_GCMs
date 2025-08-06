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
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics5.rcp85 <- names(cubicsets5.rcp85)
namescubics5left.rcp85 <- names(cubicsets5left.rcp85)

namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
#namescubics6 <- names(cubicsets6)
namescubics6left <- names(cubicsets6left)
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
hc_gcms6left <- lapply(paste0("CMIP6_left/",namescubics6left), loadIterations, permused = permk, it = it)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_jra55 <- lapply(c("JRA55_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_ncep <- lapply(c("../Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 

hc_gcms5.rcp85 <- lapply(paste0("FUTURE_",namescubics5.rcp85), loadIterations, permused = permk, it = it)
hc_gcms5left.rcp85<- lapply(paste0("FUTURE_CMIP5_left/FUTURE_",namescubics5left.rcp85), loadIterations, permused = permk, it = it)
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

# institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor")
# institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","ACCESS","bcc","BNU","CCSM4","CESM1","CMCC","CSIRO","inmcm4","MRI")
institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","ACCESS","BCC","bcc","BNU","CCSM4","CESM1","CESM2","CMCC","CSIRO","inmcm4","MRI","NESM3","UKESM1","FGOALS","CAMS")


# without cmip6
institutions <- rep(institutions1,each =2)
rep("CMIP5",length(institutions1))

CMIPS <- rep(c("CMIP5","CMIP6"),length(institutions1))
# combinations_hf <-rbind(cbind(institutions1,"historical"),cbind(institutions1,"rcp85"),cbind(institutions1,"ssp585"))
combinations <- cbind(CMIPS,institutions)
indearth <- grep("EC.EARTH_r1i|EC.EARTH_r2i",rownames(loglik_selection_datasets))

namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(loglik_selection_datasets[-indearth,-indearth]))))
namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(loglik_selection_datasets))))

unique(rownames(loglik_selection_datasets)[namesort])

# for including cmip6 and including cmip5 left + eofs
namesort1 <- c(namesort,length(namesort)+1:3)
namesort2 <- c(namesort,length(namesort)+1:9)


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


####################################
# Visualize Bhattacharya distance
####################################
load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_CMIP6_hist_vs_fut.rda")

selecteddists <- hellinger_coefficients

mat<- -log(selecteddists)

namesort_hel <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),
                                                                  rownames(mat))))
orderedmat <- mat[c(namesort_hel,length(namesort_hel)+1:3),c(namesort_hel,length(namesort_hel)+1:3)]
namesort_f <- unlist(sapply(c("rcp85","ssp585"),function(x) grep(x,rownames(orderedmat))))
namesort_h <- (1:nrow(orderedmat))[-namesort_f]
namesort_hf <- c(namesort_h,namesort_f)
orderedmat <- orderedmat[c(namesort_hf,length(namesort_hel)+1:3),c(namesort_hf,length(namesort_hel)+1:3)]
longdata <- melt(orderedmat)

mat[which(mat == Inf, arr.ind = TRUE)]<- 100
matclust <-hclust(dist(mat[], method = "minkowski", p = 2),method = "complete")
longdata <- melt(mat[matclust$order,matclust$order])

# rownames(orderedmat) <- gsub("CMIP5_","",gsub("historical","",rownames(orderedmat)))
# colnames(orderedmat) <- gsub("CMIP5_","",gsub("historical","",colnames(orderedmat)))

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
  labs(x="models", y="models", title=paste0("Bhattacharya distance ",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

plotname <-  "figs/Bhattacharya_CMIP5_CMIP6/Bhattacharya_CMIP5_CMIP6_hist_vs_fut.pdf"
pdf(plotname,width= 45, height= 45)
b
dev.off()

plotname <-  "figs/Bhattacharya_CMIP5_CMIP6/Bhattacharya_CMIP5_CMIP6_hist_vs_fut_clustcomplete.pdf"
pdf(plotname,width= 45, height= 45)
b
dev.off()

#############################################################################################
# Emergent constraints: (1) Hellinger distance vs global mean temperature increase
############################################################################################
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


# gmt.sd.difs.cmips.vs.int <- gmt.sd.hists.cmips - gmt.sd.reans["interim_10d_akima_cubic"]
gmt.mean.difs.cmips.vs.int <- gmt.hists.cmips - gmt.reans["interim_10d_akima_cubic"]
gmt.mean.difs.cmips.vs.ncep <- gmt.hists.cmips  - gmt.reans["ncep_10d"]
gmt.mean.difs.cmips.vs.jra55 <- gmt.hists.cmips  - gmt.reans["JRA55_10d_akima_cubic"]
gmt.mean.difs.cmips.vs.cmips <- gmt.futures.cmips - gmt.hists.cmips

###################################################################################
# hellinger coefficients: weighted vs unweighted.
###################################################################################
# load("results/hellinger_coefficient/hellinger_coefficients.rda")
# load("results/hellinger_coefficient/hellinger_coefficients.rcp85.rda")
# load("results/hellinger_coefficient/hellinger_coefficients_CMIP5_hist_vs_fut.rda")
# 
helcoef <- hellinger_coefficients
rownames(helcoef) <-gsub("_historical","",rownames(hellinger_coefficients)) 
colnames(helcoef)<-gsub("_historical","",colnames(hellinger_coefficients))

bd.cmips.fut.vs.hist <- diag(-log(helcoef)[names(clim.hists),names(clim.futures)])
bd.cmips.fut.vs.hist.wrt.int <- ((-log(helcoef)[names(clim.futures),"interim_10d_akima_cubic"])-(-log(helcoef)[names(clim.hists),"interim_10d_akima_cubic"]))
bd.cmips.fut.vs.hist.wrt.ncep<- (-log(helcoef)[names(clim.futures),"ncep_10d"])-(-log(helcoef)[names(clim.hists),"ncep_10d"])
bd.cmips.fut.vs.hist.wrt.jra55 <- (-log(helcoef)[names(clim.futures),"JRA55_10d_akima_cubic"])-(-log(helcoef)[names(clim.hists),"JRA55_10d_akima_cubic"])

bd.cmips.fut.vs.ncep <- -log(helcoef)[names(clim.futures),"ncep_10d"]
bd.cmips.hist.vs.int <- -log(helcoef)[names(clim.hists),"interim_10d_akima_cubic"]
bd.cmips.fut.vs.int <- -log(helcoef)[names(clim.futures),"interim_10d_akima_cubic"]
bd.cmips.hist.vs.ncep <- -log(helcoef)[names(clim.hists),"ncep_10d"]
bd.cmips.hist.vs.JRA55 <- -log(helcoef)[names(clim.hists),"JRA55_10d_akima_cubic"]
# x <- "CMIP5_inmcm4_r1i1p1"
# bd.cmips.hist.vs.x <- -log(hellinger_coefficients)[names(clim.hists),x ]
# 
# 
plot(bd.cmips.hist.vs.JRA55,gmt.mean.difs.cmips.vs.cmips)
cor(bd.cmips.hist.vs.int,gmt.mean.difs.cmips.vs.cmips, method = "pearson")

plot(bd.cmips.fut.vs.hist,gmt.mean.difs.cmips.vs.cmips)
cor(bd.cmips.fut.vs.hist,gmt.mean.difs.cmips.vs.cmips)

plot(bd.cmips.hist.vs.JRA55,gmt.mean.difs.cmips.vs.cmips)
cor(bd.cmips.hist.vs.int,gmt.mean.difs.cmips.vs.cmips, method = "pearson")

plot(bd.cmips.hist.vs.int,bd.cmips.fut.vs.int)
cor(bd.cmips.hist.vs.int,bd.cmips.fut.vs.int, method = "pearson")

# 
# # change in Global Mean Temp versus change in spatial dependency structure CMIP5 BNs
plotname<- "figs/Bhattacharya_CMIP5_CMIP6/deltaGMT_vs_Bdist_CMIP5_CMIP6.pdf"
pdf(plotname,height = 5, width = 10)
plot(bd.cmips.fut.vs.hist,gmt.mean.difs.cmips.vs.cmips, 
ylab = bquote(Delta ~ "Global Mean Temp (fut-hist)"),
      xlab = expression("d"[B]*"(BN"[fut]*",BN"[hist]*")"), pch = 16,
     main = "change in Global Mean Temp versus change in spatial dependency structure CMIP5 BNs \n(= Bhattacharya distance between future and historical BN)"
 )
text(40,3, paste("cor = ",round(cor(bd.cmips.fut.vs.hist,gmt.mean.difs.cmips.vs.cmips, method = "spearman"),2)))
lm1 <- lm(gmt.mean.difs.cmips.vs.cmips~bd.cmips.fut.vs.hist,data.frame(bd.cmips.fut.vs.hist,gmt.mean.difs.cmips.vs.cmips))
lines(x=c(18,40) , y = c(lm1$coefficients[2]*18+lm1$coefficients[1],lm1$coefficients[2]*40+lm1$coefficients[1]),col = 'red')
text(40,2.9, paste("R^2 = ",round(summary(lm1)$r.squared,3)))
dev.off()
# ###############################################################################################
# # change in dependency pattern versus historical difference between cmip model and reanalysis
# # Igual plot (F1 score vs delta precipitation in Causal networks for climate model evaluation..
# # Nowack et al. )
# ###############################################################################################
plot(bd.cmips.fut.vs.hist,bd.cmips.hist.vs.int)
plot(gmt.mean.difs.cmips.vs.int,bd.cmips.hist.vs.int)
plot(gmt.mean.difs.cmips.vs.int,bd.cmips.fut.vs.hist)
plot(gmt.mean.difs.cmips.vs.cmips,bd.cmips.hist.vs.int)
# 
# df1 <- data.frame("gmt.change" = gmt.mean.difs.cmips.vs.cmips,"Bscore" = bd.cmips.hist.vs.int)
# gmt.change2 <- gmt.mean.difs.cmips.vs.cmips^2
# quadratic.model <-lm(Bscore ~ gmt.change +gmt.change2,df1)
# summary(quadratic.model)
# degreevalues <- seq(2, 5, 0.1)
# predicteddegrees <- predict(quadratic.model,list(gmt.change=degreevalues, gmt.change2= degreevalues^2))
# 
plot(gmt.mean.difs.cmips.vs.cmips,bd.cmips.hist.vs.int, pch=16, xlab = "GMT change °C", ylab = expression("Bscore: d"[B]*"(BN"[CMIP]*",BN"[ERAint]*")"), cex.lab = 1.3, col = "blue")
 lines(degreevalues, predicteddegrees, col = "darkgreen", lwd = 3)
# 
# 
# lin.model <-lm(Bscore ~ gmt.change,df1)
# summary(lin.model)
# degreevalues <- seq(2, 5, 0.1)
# predicteddegrees <- predict(lin.model,list(gmt.change=degreevalues))
# 
# plot(gmt.mean.difs.cmips.vs.cmips,bd.cmips.hist.vs.int, pch=16, xlab = "GMT change °C", ylab = expression("Bscore: d"[B]*"(BN"[CMIP]*",BN"[ERAint]*")"), cex.lab = 1.3, col = "blue")
# lines(degreevalues, predicteddegrees, col = "darkgreen", lwd = 3)
# text(gmt.mean.difs.cmips.vs.cmips,bd.cmips.hist.vs.int, labels = names(gmt.mean.difs.cmips.vs.cmips))
 ###############################################################################################
 # variance explained 
 ###############################################################################################
 rownames(loglik_selection_datasets)
 colnames(loglik_selection_datasets)
 rownames(loglik_selection_datasets)<- gsub("historical_","",rownames(loglik_selection_datasets))
 colnames(loglik_selection_datasets)<- gsub("historical_","",colnames(loglik_selection_datasets))
 
 bd.cmips.fut.vs.int
 
 which()
 
 names(sets.hists)[40]
 
 plot(bd.cmips.hist.vs.int, -loglik_selection_datasets[names(sets.hists),"interim_10d_akima_cubic_v.exp_0.95"])
 cor(bd.cmips.hist.vs.int, -loglik_selection_datasets[names(sets.hists),"interim_10d_akima_cubic_v.exp_0.95"],method = "spearman")
 
names(bd.cmips.hist.vs.int)
 
 
######################################################
# Check temperature increase, weighted model ensemble
######################################################
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

##################################################################################
# Exponential weighting
##################################################################################
ref.perf <- "interim_10d_akima_cubic"

clim.hists.CMIP5 <- clim.hists[grep("CMIP5",names(clim.hists))]
helcoef <- hellinger_coefficients
bdists <- -log(helcoef)[c(names(clim.hists.CMIP5),"interim_10d_akima_cubic"),c(names(clim.hists.CMIP5),"interim_10d_akima_cubic")]
refID <- which(rownames(bdists) == ref.perf)
D <- bdists[ref.perf,-refID]
D <- D/(max(D))

d1 <- D/max(D)
hist(exp(-d1))
hist(d1)

S <- bdists[-refID,-refID]
S <- S/max(S)

sums <- numeric(length=nrow(S))
for( i in 1:nrow(S)){
sums[i] <- 1+sum(exp(-(S[i,][-i])))
}


sS <- 0.175
sD <- 1

W.ind <- numeric(length=length(D))
W.perf <- numeric(length=length(D))
names(W.ind)<-names(W.perf)<- names(D)

i <- 1
for(i in 1:length(names(D))){
  name.i <-names(D)[i]
  nameID.S <- which(rownames(S) == name.i)
  nameID.D <- which(names(D) == name.i)
  Di <- D[nameID.D]
  Sij <- S[nameID.S,][-nameID.S]
  
  W.perf[i] <- exp(-Di^2/sD^2)
  sumS <- sum(exp(-Sij^2/sS^2))
  W.ind[i] <- 1/(1+sumS)
}

W.perf
W <- W.perf*W.ind
W <- W/sum(W)
W.perf<- W.perf/sum(W.perf)
W.ind <- W.ind/sum(W.ind)


sd.W.perf <- sd(W.perf)
sd.W.ind <- sd(W.ind)

mean(W.perf)
mean(W.ind)
mean(W)

plotname<- paste0("figs/weights/weights_CMIP5_exponntial_sD_",sD,"_sS_",sS,".pdf")
pdf(plotname,height = 7, width = 10)

par(mar = c(15, 4, 1, 4))
plot(W,main = "exponential weighting",xaxt =  "n",xlab = "", pch = 16,cex = 2, col = "red",ylim = c(min(W,W.perf,W.ind),max(W,W.perf,W.ind)))
xlabels=names(W)
axis(1, at=1:length(W),labels = xlabels, col.axis="black", las=2)
# abline(h = 1/length(modelweights), lty = 2, col = "grey")

# text(1,min(W)+0.95*(max(W)-min(W)), bquote(sigma~"D = "));text(4,min(W)+0.95*(max(W)-min(W)),paste0(sD))
# text(1,min(W)+0.90*(max(W)-min(W)), bquote(sigma~"S = "));text(4,min(W)+0.90*(max(W)-min(W)),paste0(sS))

# par(new = T)
# plot(W.ind,axes = F,
     # xlab = NA,
     # ylab = NA,col = "dark green",pch = 16)
# axis(side = 4, col.axis = "dark green", col = "dark green")
# mtext(side = 4, line = 3, 'W.ind',col = "dark green")
points(W.ind, col = "dark green", pch = 16)

# par(new = T)
# plot(W.perf,axes = F,
     # xlab = NA,
     # ylab = NA,col = "blue",pch = 16)
# axis(side = 4, col.axis = "blue", col = "blue", line = 2)
# mtext(side = 4, line = 3, 'W.ind',col = "blue")
points(W.perf, col = "blue", pch = 16)
abline(h = 1/length(W), lty = 2, col = "grey")

legend("bottomleft",legend = c(paste0("W.perf sD  = ",sD," stand.dev = ", round(sd.W.perf,3)),paste0("W.ind sS = ",sS," stand.dev = ", round(sd.W.ind,3))), col = c("blue","dark green"),pch = 16)
                
dev.off()

###################################################################################
# First estandarizar W.perf and W.ind, than multiplicate (W), than estandarizar deneuvo
###################################################################################
W2 <- W.perf*W.ind
W2 <- W2/sum(W2)

W2-W

plotname<- paste0("figs/weights/weights_CMIP5_exponntial_sD_",sD,"_sS_",sS,"_2xnorm.pdf")
pdf(plotname,height = 7, width = 10)

par(mar = c(15, 4, 1, 4))
plot(W2,main = "exponential weighting",xaxt =  "n",xlab = "", pch = 16,cex = 2, col = "red",ylim = c(min(W2,W.perf,W.ind),max(W2,W.perf,W.ind)))
xlabels=names(W2)
axis(1, at=1:length(W),labels = xlabels, col.axis="black", las=2)
points(W.ind, col = "dark green", pch = 16)
points(W.perf, col = "blue", pch = 16)
abline(h = 1/length(W2), lty = 2, col = "grey")

legend("bottomleft",legend = c(paste0("W.perf sD  = ",sD," stand.dev = ", round(sd.W.perf,3)),paste0("W.ind sS = ",sS," stand.dev = ", round(sd.W.ind,3))), col = c("blue","dark green"),pch = 16)

dev.off()


##################################################################################
# Hyperbolic weighting
##################################################################################
ref.perf <- "interim_10d_akima_cubic"

clim.hists.CMIP5 <- clim.hists[grep("CMIP5",names(clim.hists))]

bdists <- -log(helcoef)[c(names(clim.hists.CMIP5),"interim_10d_akima_cubic"),c(names(clim.hists.CMIP5),"interim_10d_akima_cubic")]
refID <- which(rownames(bdists) == ref.perf)
D <- bdists[ref.perf,-refID]
D <- D/(max(D))

S <- bdists[-refID,-refID]
S <- S/max(S)

sS <- 1
sD <- 1

W.ind <- numeric(length=length(D))
W.perf <- numeric(length=length(D))
names(W.ind)<-names(W.perf)<- names(D)

i <- 1
for(i in 1:length(names(D))){
  name.i <-names(D)[i]
  nameID.S <- which(rownames(S) == name.i)
  nameID.D <- which(names(D) == name.i)
  Di <- D[nameID.D]
  Sij <- S[nameID.S,][-nameID.S]
  
  W.perf[i] <- 1/Di
  sumS <- sum(1/Sij)
  W.ind[i] <- 1/(1+sumS)
}

W.perf
W <- W.perf*W.ind
W <- W/sum(W)
W.perf<- W.perf/sum(W.perf)
W.ind <- W.ind/sum(W.ind)

plotname<- paste0("figs/weights/weights_CMIP5_hyperbolic_sD_",sD,"_sS_",sS,".pdf")
pdf(plotname,height = 7, width = 10)

par(mar = c(15, 4, 1, 4))
plot(W,xaxt =  "n",xlab = "", main = "hyperbolic weighting", pch = 16,cex = 2, col = "red",ylim = c(min(W,W.perf,W.ind),max(W,W.perf,W.ind)))
xlabels=names(W)
axis(1, at=1:length(W),labels = xlabels, col.axis="black", las=2)

points(W.ind, col = "dark green", pch = 16)
points(W.perf, col = "blue", pch = 16)
abline(h = 1/length(W), lty = 2, col = "grey")

legend("bottomleft",legend = c(paste0("W.perf sD  = ",sD),paste0("W.ind sS = ",sS)), col = c("blue","dark green"),pch = 16)
dev.off()

##################################################################################
# exponential performance weighting + step function independence weighting
##################################################################################
ref.perf <- "interim_10d_akima_cubic"

clim.hists.CMIP5 <- clim.hists[grep("CMIP5",names(clim.hists))]

bdists <- -log(helcoef)[c(names(clim.hists.CMIP5),"interim_10d_akima_cubic"),c(names(clim.hists.CMIP5),"interim_10d_akima_cubic")]
refID <- which(rownames(bdists) == ref.perf)
D <- bdists[ref.perf,-refID]
D <- D/(max(D))

S <- bdists[-refID,-refID]
S <- S/max(S)

sS <- 0.4
sD <- 1

W.ind <- numeric(length=length(D))
W.perf <- numeric(length=length(D))
names(W.ind)<-names(W.perf)<- names(D)

i <- 1
for(i in 1:length(names(D))){
  name.i <-names(D)[i]
  nameID.S <- which(rownames(S) == name.i)
  nameID.D <- which(names(D) == name.i)
  Di <- D[nameID.D]
  Sij <- S[nameID.S,][-nameID.S]
  Step.Sij <-Sij
  Step.Sij[Sij>sS]<- 0
  Step.Sij[Sij<=sS]<- 1
  
  W.perf[i] <- exp(-Di^2/sD^2)
  sumS <- sum(Step.Sij)
  W.ind[i] <- 1/(1+sumS)
}

W.perf
W <- W.perf*W.ind
W <- W/sum(W)
W.perf<- W.perf/sum(W.perf)
W.ind <- W.ind/sum(W.ind)


plotname<- paste0("figs/weights/weights_CMIP5_expstep_sD_",sD,"_sS_",sS,".pdf")
pdf(plotname,height = 7, width = 10)

par(mar = c(15, 4, 1, 4))
plot(W,xaxt =  "n",xlab = "", main = "exponential performance weighting + step function independence weighting", pch = 16,cex = 2, col = "red",ylim = c(min(W,W.perf,W.ind),max(W,W.perf,W.ind)))
xlabels=names(W)
axis(1, at=1:length(W),labels = xlabels, col.axis="black", las=2)

points(W.ind, col = "dark green", pch = 16)
points(W.perf, col = "blue", pch = 16)
abline(h = 1/length(W), lty = 2, col = "grey")

legend("bottomleft",legend = c(paste0("W.perf sD  = ",sD),paste0("W.ind sS = ",sS)), col = c("blue","dark green"),pch = 16)

dev.off()
##################################################################################
# One for one
##################################################################################
sD <- 0.2
sS <- 0.8
helcoef <- hellinger_coefficients
rownames(helcoef) <-gsub("_historical","",rownames(hellinger_coefficients)) 
colnames(helcoef)<-gsub("_historical","",colnames(hellinger_coefficients))


modelweights <- IPweight(-log(helcoef)[c(names(clim.hists),"interim_10d_akima_cubic"),c(names(clim.hists),"interim_10d_akima_cubic")],"interim_10d_akima_cubic",sD = sD,sS =sS)

center <- FALSE
modelweights <- IPweight(scale(-log(helcoef)[c(names(clim.hists),"interim_10d_akima_cubic"),c(names(clim.hists),"interim_10d_akima_cubic")],center = center),"interim_10d_akima_cubic",sD = sD,sS =sS)

scale(-log(helcoef)[c(names(clim.hists),"interim_10d_akima_cubic"),c(names(clim.hists),"interim_10d_akima_cubic")],center = center)

if (center == TRUE) {plotname<- paste0("figs/weights/weights_CMIP5_CMIP6_sD",sD,"_sS",sS,".pdf")} else {plotname<- paste0("figs/weights/weights_CMIP5_CMIP6_sD",sD,"_sS",sS,"_nocenter.pdf")}
pdf(plotname,height = 5, width = 10)
par(mar = c(9, 4, 1, 2))
plot(modelweights,xaxt =  "n",xlab = "", pch = 16)
xlabels=gsub("_r1i1p1","",gsub("CMIP5_","",gsub("historical","",names(modelweights))))
axis(1, at=1:length(modelweights),labels = xlabels, col.axis="red", las=2)
abline(h = 1/length(modelweights), lty = 2, col = "grey")
text(25,min(modelweights) + 0.10*(max(modelweights)-min(modelweights)), bquote(sigma~"D = "));text(27.5,min(modelweights) + 0.10*(max(modelweights)-min(modelweights)),paste0(sD))
text(25,min(modelweights) + 0.15*(max(modelweights)-min(modelweights)), bquote(sigma~"S = "));text(27.5,min(modelweights) + 0.14*(max(modelweights)-min(modelweights)),paste0(sS))
dev.off()


##################################################################################
# Generate modelweights figures
##################################################################################
sDs <- seq(2,1,0.1)
sSs <- rev(seq(2,1,0.1))

helcoef <- hellinger_coefficients
rownames(helcoef) <-gsub("_historical","",rownames(hellinger_coefficients)) 
colnames(helcoef)<-gsub("_historical","",colnames(hellinger_coefficients))

center <- FALSE

if (center == TRUE) {plotname<- paste0("figs/weights/weights_CMIP5_CMIP6_sDs_",sDs[1],"_",sDs[length(sDs)],"_sSs",sSs[1],"_",sSs[length(sSs)],".pdf")} else {plotname<- paste0("figs/weights/weights_CMIP5_CMIP6_sDs_",sDs[1],"_",sDs[length(sDs)],"_sSs",sSs[1],"_",sSs[length(sSs)],"nocenter.pdf")}
pdf(plotname,height = 20, width = 25)
par(mfrow = c(3,3))
for(i in 1:length(sDs)){
modelweights <- IPweight(scale(-log(helcoef)[c(names(clim.hists),"interim_10d_akima_cubic"),c(names(clim.hists),"interim_10d_akima_cubic")],center = center),"interim_10d_akima_cubic",sD = sDs[i],sS =sSs[i])
par(mar = c(9, 4, 1, 2))
plot(modelweights,xaxt =  "n",xlab = "", pch = 16,ylim = c(0,0.2))
xlabels=gsub("_r1i1p1","",gsub("CMIP5_","",gsub("historical","",names(modelweights))))
axis(1, at=1:length(modelweights),labels = xlabels, col.axis="red", las=2)
abline(h = 1/length(modelweights), lty = 2, col = "grey")
text(1,0.80*0.2, bquote(sigma~"D = "));text(4,0.80*0.2,paste0(sDs[i]))
text(1,0.90*0.2, bquote(sigma~"S = "));text(4,0.90*0.2,paste0(sSs[i]))
}
dev.off()

##############################################################################################################
# 
##############################################################################################################
IPweight.split <- function(bdists,ref.perf,sD,sS,meth = c("exponential","hyperbolic","expstep")){
  refID <- which(rownames(bdists) == ref.perf)
  D <- bdists[ref.perf,-refID]
  D <- D/(max(D))
  S <- bdists[-refID,-refID]
  S <- S/max(S)
  
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
    W.ind[i] <- 1/(1+sumS)}
  }
  
  W <- W.perf*W.ind
  W <- W/sum(W)
  W.perf<- W.perf/sum(W.perf)
  W.ind <- W.ind/sum(W.ind)
  
  return(data.frame(W,W.perf,W.ind))
}


# CMIP 5
clim.hists.CMIP5 <- clim.hists[grep("CMIP5",names(clim.hists))]
helcoef <- hellinger_coefficients
b <- -log(helcoef)[c(names(clim.hists.CMIP5),"interim_10d_akima_cubic"),c(names(clim.hists.CMIP5),"interim_10d_akima_cubic")]


# CMIP 6
clim.hists.CMIP6 <- clim.hists[grep("CMIP6",names(clim.hists))]
helcoef <- hellinger_coefficients
rownames(helcoef) <- gsub("_historical","", rownames(helcoef))
colnames(helcoef) <- gsub("_historical","", colnames(helcoef))
names(clim.hists.CMIP6) %in% rownames(helcoef)
b <- -log(helcoef)[c(names(clim.hists.CMIP6),"interim_10d_akima_cubic"),c(names(clim.hists.CMIP6),"interim_10d_akima_cubic")]

# CMIP 5 and 6
rownames(helcoef) <- gsub("_historical","", rownames(helcoef))
colnames(helcoef) <- gsub("_historical","", colnames(helcoef))
names(clim.hists) %in% rownames(helcoef)
b <- -log(helcoef)[c(names(clim.hists),"interim_10d_akima_cubic"),c(names(clim.hists),"interim_10d_akima_cubic")]


modelweights.split <- IPweight.split(b,"interim_10d_akima_cubic",sD = 1,sS = 0.175, meth = "exponential")
modelweights <- modelweights.split[,"W.perf"]
names(modelweights)<-rownames(modelweights.split)
##############################################################################################################
# Mean and standarddeviation CMIP5
##############################################################################################################
ind <- 1:length(gmt.hists.cmips) # UITKIJKEN ALS JE EERST DE LIJN HIERONDER DOET EN DAN WEER DEZE DAN IS IS HIJ AL OVERSCHREVEN
ind <-grep("CMIP5",names(gmt.hists.cmips))
gmt.hists.cmips <- gmt.hists.cmips[ind]
gmt.futures.cmips <- gmt.futures.cmips[ind]
gmt.diffs.cmips2 <- gmt.diffs.cmips2[ind]

mean.hist.weighted <- sum(modelweights*gmt.hists.cmips)
mean.hist.unweighted <- sum(gmt.hists.cmips)/length(gmt.hists.cmips)
mean.fut.weighted <- sum(modelweights*gmt.futures.cmips)
mean.fut.unweighted <- sum(gmt.futures.cmips)/length(gmt.futures.cmips)
mean.diff.weighted <- sum(modelweights*gmt.diffs.cmips2)
mean.diff.unweighted <- sum(gmt.diffs.cmips2)/length(gmt.diffs.cmips2)

sd.hist.weighted <- sqrt(sum(modelweights*(gmt.hists.cmips-mean.hist.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
sd.hist.unweighted <- sqrt(sum((gmt.hists.cmips-mean.hist.unweighted)^2)/((length(gmt.hists.cmips)-1)))
sd.fut.weighted <- sqrt(sum(modelweights*(keltoC(gmt.futures.cmips)-mean.fut.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
sd.fut.unweighted <- sqrt(sum((keltoC(gmt.futures.cmips)-mean.fut.unweighted)^2)/((length(gmt.futures.cmips)-1)))
sd.diff.weighted <- sqrt(sum(modelweights*(keltoC(gmt.diffs.cmips2)-mean.diff.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
sd.diff.weighted2 <- sqrt(sd.hist.weighted^2 + sd.fut.weighted^2-2*cov(modelweights*(keltoC(gmt.futures.cmips)),modelweights*(gmt.hists.cmips)))
sd.diff.unweighted <- sqrt(sum((gmt.diffs.cmips2-mean.diff.unweighted)^2)/((length(gmt.diffs.cmips2)-1)))
sd.diff.unweighted2 <- sqrt(sd.hist.unweighted^2 + sd.fut.unweighted^2-2*cov((keltoC(gmt.futures.cmips)),(gmt.hists.cmips)))

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
plot(bd.cmips.fut.vs.hist,sd.diff)
#########################################################################
# Plot weights and temperature rise in one figure
#########################################################################

# plotname<- paste0("figs/weights/deltaGMT_and_weights_CMIP5_CMIP6_sD",sD,"_sS",sS,".pdf")
# pdf(plotname,height = 5, width = 12)
par(mar = c(11, 4, 2, 2))
plot(gmt.mean.difs.cmips.vs.cmips[ind],xaxt =  "n",xlab = "",pch = 16,
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

legend("topleft",legend = c("weighted","unweighted"), col = c("blue","grey"),lty = 2)

dev.off()
#################################################################################
# Calculate correlation between climatologies
#################################################################################
clim.hists$CMIP5_ACCESS1.0_r1i1p1



#################################################################################
# Calculate spatial correlation between gridpoints each dataset.
#################################################################################
cors_gcms <- lapply(data_gcms_anom_scaled,cor, method = "spearman")

cors.cors.gcms <- matrix(nrow = length(cors_gcms),ncol = length(cors_gcms),dimnames = list(names(cors_gcms),names(cors_gcms)))
RMSE.cors.gcms <- matrix(nrow = length(cors_gcms),ncol = length(cors_gcms),dimnames = list(names(cors_gcms),names(cors_gcms)))
for (j in 1:length(cors_gcms)){
  cors.cors.gcms.i<- numeric(length(cors_gcms))
  RMSE.cors.gcms.i <- numeric(length(cors_gcms))
for (i in 1:length(cors_gcms)){
  x1 <- cors_gcms[[i]]
  x2 <- cors_gcms[[j]]
  # x2 <- cors_gcms[["ncep_10d"]]
 cors.cors.gcms.i[i]<- cor(as.vector(x1),as.vector(x2))
  RMSE.cors.gcms.i[i] <- sqrt(sum((x1-x2)^2)/length(x1))
}
RMSE.cors.gcms[,j] <- RMSE.cors.gcms.i
cors.cors.gcms[,j] <- cors.cors.gcms.i
}


rownames(RMSE.cors.gcms)
RMSE.cors.gcms.int <- RMSE.cors.gcms["interim_10d_akima_cubic",]


spatialPlot(quantity2clim(cors_gcms$CMIP6Amon_UKESM1.0.LL_r1i1p1f2[81,],"cor v81",sets.hists$CMIP6Amon_UKESM1.0.LL_r1i1p1f2),backdrop.theme = "coastline", at = seq(-1,1,0.1))
# ind <- 1:length(gmt.hists.cmips) # UITKIJKEN ALS JE EERST DE LIJN HIERONDER DOET EN DAN WEER DEZE DAN IS IS HIJ AL OVERSCHREVEN

RMSE.cors.gcms.hists.int <- RMSE.cors.gcms.int[names(sets.hists)]
RMSE.cors.gcms.futures.int <- RMSE.cors.gcms.int[names(sets.futures)]
RMSE.cors.gcms.reans.int <- RMSE.cors.gcms.int[names(sets.reans)]

RMSE.cors.gcms.d <- RMSE.cors.gcms
bhat.gcms.d <- -log(helcoef)
diag(bhat.gcms.d) <- diag(RMSE.cors.gcms.d)<-NA
bhat.gcms.d[which(bhat.gcms.d== Inf, arr.ind = TRUE)]<- NA

cor(as.vector(bhat.gcms.d),as.vector(RMSE.cors.gcms.d),use = "pairwise.complete.obs")

plot(RMSE.cors.gcms[names(sets.futures),names(sets.futures)],-log(helcoef)[names(sets.futures),names(sets.futures)])
cor(as.vector(RMSE.cors.gcms[names(sets.futures),names(sets.futures)]),as.vector(-log(helcoef)[names(sets.futures),names(sets.futures)]),method = "spearman")
points(RMSE.cors.gcms[names(sets.hists),names(sets.hists)],-log(helcoef)[names(sets.hists),names(sets.hists)], col = "blue")
cor(use = "complete.obs",as.vector(RMSE.cors.gcms[names(sets.hists),names(sets.hists)]),as.vector(bhat.gcms.d[names(sets.hists),names(sets.hists)]), method = "spearman")

points(RMSE.cors.gcms[names(sets.hists),names(sets.futures)],-log(helcoef)[names(sets.hists),names(sets.futures)], col = "purple")
cor(use = "complete.obs",as.vector(RMSE.cors.gcms[names(sets.hists),names(sets.futures)]),as.vector(bhat.gcms.d[names(sets.hists),names(sets.futures)]))

points(RMSE.cors.gcms[names(sets.reans),names(sets.reans)],-log(helcoef)[names(sets.reans),names(sets.reans)], col = "green")

points(RMSE.cors.gcms[names(sets.reans),names(sets.futures)],-log(helcoef)[names(sets.reans),names(sets.futures)], col = "red")
cor(as.vector(RMSE.cors.gcms[names(sets.reans),names(sets.futures)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.futures)]),method = "spearman")
points(RMSE.cors.gcms[names(sets.reans),names(sets.hists)],-log(helcoef)[names(sets.reans),names(sets.hists)], col = "orange")
cor(as.vector(RMSE.cors.gcms[names(sets.reans),names(sets.hists)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.hists)]),method = "spearman")
    
    




plot(RMSE.cors.gcms.hists.int,bd.cmips.hist.vs.int)

text(RMSE.cors.gcms.hists.int,bd.cmips.hist.vs.int,names(RMSE.cors.gcms.hists.int))

ind.CMIP5.hist <- 1:length(RMSE.cors.gcms.hists.int)
# ind.CMIP5.hist <-grep("CMIP5",names(RMSE.cors.gcms.hists.int))
RMSE.cors.gcms.hists.CMIP5.int <- RMSE.cors.gcms.hists.int[ind.CMIP5.hist]
ind.CMIP5.fut <-  1:length(RMSE.cors.gcms.hists.int)
# ind.CMIP5.fut <-grep("CMIP5",names(RMSE.cors.gcms.futures.int))
RMSE.cors.gcms.futures.CMIP5.int <- RMSE.cors.gcms.futures.int[ind.CMIP5.fut]


mean.cors.hist.weighted <- sum(modelweights*RMSE.cors.gcms.hists.CMIP5.int)
mean.cors.hist.unweighted <- sum(RMSE.cors.gcms.hists.CMIP5.int)/length(RMSE.cors.gcms.hists.CMIP5.int)
mean.cors.fut.weighted <- sum(modelweights*RMSE.cors.gcms.futures.CMIP5.int)
mean.cors.fut.unweighted <- sum(RMSE.cors.gcms.futures.CMIP5.int)/length(RMSE.cors.gcms.futures.CMIP5.int)
# mean.diff.weighted <- sum(modelweights*gmt.diffs.cmips2)
# mean.diff.unweighted <- sum(gmt.diffs.cmips2)/length(gmt.diffs.cmips2)

sd.cors.hist.weighted <- sqrt(
  sum(modelweights*(
    RMSE.cors.gcms.hists.CMIP5.int-mean.cors.hist.weighted
    )^2)
  /
    ((length(modelweights)-1)/(
      length(modelweights)*sum(modelweights))
     )
  )
sd.cors.hist.unweighted <- sqrt(sum((RMSE.cors.gcms.hists.CMIP5.int-mean.cors.hist.unweighted)^2)/((length(RMSE.cors.gcms.hists.CMIP5.int)-1)))
sd.cors.fut.weighted <- sqrt(sum(modelweights*(keltoC(RMSE.cors.gcms.futures.CMIP5.int)-mean.cors.fut.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
sd.cors.fut.unweighted <- sqrt(sum((RMSE.cors.gcms.futures.CMIP5.int-mean.cors.fut.unweighted)^2)/((length(RMSE.cors.gcms.futures.CMIP5.int)-1)))
# sd.diff.weighted <- sqrt(sum(modelweights*(keltoC(gmt.diffs.cmips2)-mean.diff.weighted)^2)/((length(modelweights)-1)/(length(modelweights)*sum(modelweights))))
# sd.diff.unweighted <- sqrt(sum((gmt.diffs.cmips2-mean.diff.unweighted)^2)/((length(gmt.diffs.cmips2)-1)))

#######################################################################################
# cors cors analyse
#######################################################################################
plotname<- paste0("figs/weights/Bhattacharya_vs_spatcorrelation_CMIP5_CMIP6.pdf")
pdf(plotname,height = 7, width = 7) 

cor.method <- "pearson"
plot(cors.cors.gcms[names(sets.futures),names(sets.futures)],-log(helcoef)[names(sets.futures),names(sets.futures)],
     xlab = "Correlation spatial correlation field [Model A, Model B]",ylab = "Bhattacharya distance BNs Model A Model B",pch = 16)
corfuts <- cor(as.vector(cors.cors.gcms[names(sets.futures),names(sets.futures)]),as.vector(-log(helcoef)[names(sets.futures),names(sets.futures)]),method = cor.method)%>% round(digits = 2)

points(cors.cors.gcms[names(sets.hists),names(sets.hists)],-log(helcoef)[names(sets.hists),names(sets.hists)], col = "blue",pch = 16)
corhists <- cor(use = "complete.obs",as.vector(cors.cors.gcms[names(sets.hists),names(sets.hists)]),as.vector(bhat.gcms.d[names(sets.hists),names(sets.hists)]), method = cor.method)%>% round(digits = 2)

points(cors.cors.gcms[names(sets.hists),names(sets.futures)],-log(helcoef)[names(sets.hists),names(sets.futures)], col = "purple",pch = 16)
corhistfut <- cor(use = "complete.obs",as.vector(cors.cors.gcms[names(sets.hists),names(sets.futures)]),as.vector(bhat.gcms.d[names(sets.hists),names(sets.futures)]),method = cor.method)%>% round(digits = 2)

points(cors.cors.gcms[names(sets.reans),names(sets.reans)],-log(helcoef)[names(sets.reans),names(sets.reans)], col = "green",pch = 16)
correans<-cor(as.vector(cors.cors.gcms[names(sets.reans),names(sets.reans)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.reans)]),method = cor.method)%>% round(digits = 2)

points(cors.cors.gcms[names(sets.reans),names(sets.futures)],-log(helcoef)[names(sets.reans),names(sets.futures)], col = "red",pch = 16)
correansfut <- cor(as.vector(cors.cors.gcms[names(sets.reans),names(sets.futures)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.futures)]),method = cor.method)%>% round(digits = 2)

points(cors.cors.gcms[names(sets.reans),names(sets.hists)],-log(helcoef)[names(sets.reans),names(sets.hists)], col = "orange",pch = 16)
correanshist <- cor(as.vector(cors.cors.gcms[names(sets.reans),names(sets.hists)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.hists)]),method = cor.method)%>% round(digits = 2)


legend("topright",legend = c(paste0("A = fut, B = fut cor = ",corfuts),
                             paste0("A = hist, B = hist cor = ",corhists),
                             paste0("A = hist, B = fut cor = ",corhistfut),
                             paste0("A = reans, B = fut cor = ",correansfut),
                             paste0("A = reans, B = hist cor = ",correanshist),
                             paste0("A = reans, B = reans cor = ",correans)),
                             col = c("black","blue","purple","red","orange","green"),pch = 16)

dev.off()


which(names(sets.hists)=="CMIP6Amon_FGOALS.f3.L_r1i1p1f1")
plot(as.vector(cors.cors.gcms["interim_10d_akima_cubic",names(sets.hists)[-40]]),as.vector(-log(helcoef)["interim_10d_akima_cubic",names(sets.hists)][-40]),col = "brown")
cor(as.vector(cors.cors.gcms["interim_10d_akima_cubic",names(sets.hists)[-40]]),as.vector(-log(helcoef)["interim_10d_akima_cubic",names(sets.hists)[-40]]),method = "pearson")
cors.cors.gcms["interim_10d_akima_cubic",names(sets.hists)[-40]][order(cors.cors.gcms["interim_10d_akima_cubic",names(sets.hists)[-40]])]
names(sets.reans)

dev.off()
#################################################################################
# w.r.t EOFS
#################################################################################
cors_EOFs <- lapply(data_pickedEOFs,cor, method = "spearman")

RMSE.cors.EOFs <- matrix(nrow = length(cors_gcms),ncol = length(cors_EOFs),dimnames = list(names(cors_gcms),names(cors_EOFs)))
cors.cors.EOFs <- matrix(nrow = length(cors_gcms),ncol = length(cors_EOFs),dimnames = list(names(cors_gcms),names(cors_EOFs)))

for (j in 1:length(cors_EOFs)){
  

  RMSE.cors.EOFs.i <- numeric(length(cors_gcms))
  cors.cors.EOFs.i <- numeric(length(cors_gcms))
  for (i in 1:length(cors_gcms)){
    x1 <- cors_gcms[[i]]
    x2 <- cors_EOFs[[j]]
    # x2 <- cors_gcms[["ncep_10d"]]
    
    RMSE.cors.EOFs.i[i] <- sqrt(sum((x1-x2)^2)/length(x1))
    cors.cors.EOFs.i[i] <- cor(as.vector(x1),as.vector(x2))
  }
  RMSE.cors.EOFs[,j] <- RMSE.cors.EOFs.i
  cors.cors.EOFs[,j] <- cors.cors.EOFs.i
}


rownames(RMSE.cors.gcms)
RMSE.cors.gcms.int <- RMSE.cors.gcms["interim_10d_akima_cubic",]


spatialPlot(
  quantity2clim(cors_gcms$CMIP5_IPSL.CM5A.LR_r1i1p1[81,],"cor v81",sets.hists$CMIP5_IPSL.CM5A.LR_r1i1p1),
  backdrop.theme = "coastline", at = seq(-1,1,0.1), rev.colors = TRUE)
# ind <- 1:length(gmt.hists.cmips) # UITKIJKEN ALS JE EERST DE LIJN HIERONDER DOET EN DAN WEER DEZE DAN IS IS HIJ AL OVERSCHREVEN

# RMSE.cors.gcms.hists.int <- RMSE.cors.gcms.int[names(sets.hists)]
# RMSE.cors.gcms.futures.int <- RMSE.cors.gcms.int[names(sets.futures)]
# RMSE.cors.gcms.reans.int <- RMSE.cors.gcms.int[names(sets.reans)]


plot(RMSE.cors.EOFs[names(clim.hists.CMIP5),"interim_10d_akima_cubic_v.exp_0.75"],-log(helcoef)[names(clim.hists.CMIP5),"interim_10d_akima_cubic"])
cor(RMSE.cors.EOFs[names(clim.hists.CMIP5),"interim_10d_akima_cubic_v.exp_0.75"],-log(helcoef)[names(clim.hists.CMIP5),"interim_10d_akima_cubic"],method = "spearman")
text(RMSE.cors.EOFs[names(clim.hists.CMIP5),"interim_10d_akima_cubic_v.exp_0.75"],-log(helcoef)[names(clim.hists.CMIP5),"interim_10d_akima_cubic"],names(clim.hists.CMIP5))

plot(RMSE.cors.EOFs[names(clim.hists),"interim_10d_akima_cubic_v.exp_0.9"],loglik_selection_datasets[names(clim.hists),"interim_10d_akima_cubic_v.exp_0.9"])
cor(RMSE.cors.EOFs[names(clim.hists),"interim_10d_akima_cubic_v.exp_0.9"],loglik_selection_datasets[names(clim.hists),"interim_10d_akima_cubic_v.exp_0.9"],method = "kendall")
text(RMSE.cors.EOFs[names(clim.hists),"interim_10d_akima_cubic_v.exp_0.9"],loglik_selection_datasets[names(clim.hists),"interim_10d_akima_cubic_v.exp_0.9"],names(clim.hists))
     

plot(RMSE.cors.EOFs[names(sets.futures),"interim_10d_akima_cubic_v.exp_0.99"],-log(helcoef)[names(sets.futures),"interim_10d_akima_cubic"])
cor(RMSE.cors.EOFs[names(sets.futures),"interim_10d_akima_cubic_v.exp_0.99"],-log(helcoef)[names(sets.futures),"interim_10d_akima_cubic"])

points(RMSE.cors.gcms[names(sets.hists),names(sets.hists)],-log(helcoef)[names(sets.hists),names(sets.hists)], col = "blue")
cor(as.vector(RMSE.cors.gcms[names(sets.hists),names(sets.hists)]),as.vector(-log(helcoef)[names(sets.hists),names(sets.hists)]))

#############################################################################################################################
# Plot correlation spatial correlation CMIP5 dataset and EOFs vs loglikelihood EOFs under CMIP5 BN
################################################################################################################################
whichEOF <- "interim_10d_akima_cubic_v.exp_0.95"
whichOUT <- c(40)

plot(cors.cors.EOFs[names(clim.hists.CMIP5),whichEOF],-log(helcoef)[names(clim.hists.CMIP5),"interim_10d_akima_cubic"])
cor(cors.cors.EOFs[names(clim.hists.CMIP5),whichEOF],-log(helcoef)[names(clim.hists.CMIP5),"interim_10d_akima_cubic"],method = "spearman")
text(cors.cors.EOFs[names(clim.hists.CMIP5),whichEOF],-log(helcoef)[names(clim.hists.CMIP5),"interim_10d_akima_cubic"],names(clim.hists.CMIP5))


plot(cors.cors.EOFs[names(clim.hists.CMIP5),whichEOF],-log(helcoef)[names(clim.hists.CMIP5),"interim_10d_akima_cubic"], xlim = c(0.5,0.7))
cors.EOF.BHAT.CMIP5[1] <- cor(cors.cors.EOFs[names(clim.hists.CMIP5),whichEOF],-log(helcoef)[names(clim.hists.CMIP5),"interim_10d_akima_cubic"],method = "spearman")

for (i in 2:length(colnames(cors.cors.EOFs))){
  whichEOF <- colnames(cors.cors.EOFs)[i]
  points(cors.cors.EOFs[names(clim.hists.CMIP5)[-whichOUT],whichEOF],
         -log(helcoef)[names(clim.hists.CMIP5)[-whichOUT],"interim_10d_akima_cubic"],
         pch = 16, col = rainbow(i))
  cors.EOF.BHAT.CMIP5[i]<- cor(cors.cors.EOFs[names(clim.hists.CMIP5)[-whichOUT],whichEOF],-log(helcoef)[names(clim.hists.CMIP5)[-whichOUT],"interim_10d_akima_cubic"],method = "spearman")
}


mapply(function(x,y) paste0(x," cor = ", y), x= colnames(cors.cors.EOFs),y = round(cors.EOF.BHAT.CMIP5,digits = 2))



#######################################################################################################################################3
whichEOF <- "interim_10d_akima_cubic_v.exp_0.75"
whichOUT <- c(40)
plot(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],loglik_selection_datasets[names(clim.hists)[-whichOUT],whichEOF])
cor(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],loglik_selection_datasets[names(clim.hists)[-whichOUT],whichEOF],method = "spearman")
text(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],loglik_selection_datasets[names(clim.hists)[-whichOUT],whichEOF],names(clim.hists)[-40])

##################################################################################################
# Plot correlation spatial correlation CMIP dataset and EOFs vs loglikelihood EOFs under CMIP BN
##################################################################################################
cors.EOF.BHAT<- c()
whichEOF<-colnames(cors.cors.EOFs)[1]
plot(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],-log(helcoef)[names(clim.hists)[-whichOUT],"interim_10d_akima_cubic"],col = rainbow(1),pch = 16)
cors.EOF.BHAT[1]<- cor(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],-log(helcoef)[names(clim.hists)[-whichOUT],"interim_10d_akima_cubic"],method = "spearman")

for (i in 2:length(colnames(cors.cors.EOFs))){
  whichEOF <- colnames(cors.cors.EOFs)[i]
  points(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],
         -log(helcoef)[names(clim.hists)[-whichOUT],"interim_10d_akima_cubic"],
         pch = 16, col = rainbow(i))
  cors.EOF.BHAT[i]<- cor(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],-log(helcoef)[names(clim.hists)[-whichOUT],"interim_10d_akima_cubic"],method = "spearman")

}

mapply(function(x,y) paste0(x," cor = ", y), x= colnames(cors.cors.EOFs),y = round(cors.EOF.BHAT,digits = 2))

plot(-cors.EOF.BHAT, xaxt = "n", xlab = "", ylab = "correlation", main = "correlation Bhattacharya and spatial correlation w.r.t. EOFs" )
xlabels = colnames(cors.cors.EOFs)
axis(1, at=1:length(cors.EOF.BHAT),labels = xlabels, col.axis="black", las=2)


# legend("topright",legend = c(paste0(whichEOFs,cors.EOF.BHAT[i]),
#                              paste0("A = hist, B = hist cor = ",cors.EOF.BHAT[i]),
#                              paste0("A = hist, B = fut cor = ",cors.EOF.BHAT[i]),
#                              paste0("A = reans, B = fut cor = ",cors.EOF.BHAT[i]),
#                              paste0("A = reans, B = hist cor = ",cors.EOF.BHAT[i]),
#                              paste0("A = reans, B = reans cor = ",cors.EOF.BHAT[i])),
#        col = c("black","blue","purple","red","orange","green"),pch = 16)

plot(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],-log(helcoef)[names(clim.hists)[-whichOUT],"interim_10d_akima_cubic"])
cor(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],-log(helcoef)[names(clim.hists)[-whichOUT],"interim_10d_akima_cubic"],method = "spearman")

#################################################################################################
# Plot correlation spatial correlation CMIP dataset and EOFs vs loglikelihood EOFs under CMIP BN
#################################################################################################
cors.EOF.LOGLIK <- c()
whichEOF <-colnames(cors.cors.EOFs)[6]
plot(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],loglik_selection_datasets[names(clim.hists)[-whichOUT],whichEOF],col = rainbow(1),pch = 16)
cors.EOF.LOGLIK[6] <- cor(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF], loglik_selection_datasets[names(clim.hists)[-whichOUT],whichEOF],method = "spearman")

for (i in 1:(length(colnames(cors.cors.EOFs))-1)){
  whichEOF <- colnames(cors.cors.EOFs)[i]
  points(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF],
         loglik_selection_datasets[names(clim.hists)[-whichOUT],whichEOF],
         pch = 16, col = rainbow(i))
  cors.EOF.LOGLIK[i] <- cor(cors.cors.EOFs[names(clim.hists)[-whichOUT],whichEOF], loglik_selection_datasets[names(clim.hists)[-whichOUT],whichEOF],method = "spearman")

}

mapply(function(x,y) paste0(x," cor = ", y), x= colnames(cors.cors.EOFs),y = round(cors.EOF.LOGLIK,digits = 2))
plot(cors.EOF.LOGLIK, xaxt = "n", xlab = "", ylab = "correlation", main = "correlation loglikelihood and spatial correlation w.r.t. EOFs" )
xlabels = colnames(cors.cors.EOFs)
axis(1, at=1:length(cors.EOF.LOGLIK),labels = xlabels, col.axis="black", las=2)

#####################################################################################
# Check if Bhattacharya distance has relationship with distance principal components: Method 1 EOF projection
#####################################################################################
Models_prinCompSet <- lapply(grid_gcms_anom_scaled, function(x) prinComp(x, n.eofs = 360))
Models_expvars <- lapply(Models_prinCompSet, function(x) attr(x[[1]][[1]],"explained_variance"))
Models_EOFs <- lapply(Models_prinCompSet, function(x) x[[1]][[1]]$EOFs)
Models_PCs <- lapply(Models_prinCompSet, function(x) x[[1]][[1]]$PCs)

Models_EOFS_preparation <- function(whichPCs, EOFs, PCs){
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

Models_EOFS_preparation(whichPCs = c(1:5), EOFs = Models_EOFs$CMIP5_CanESM2_r1i1p1, PCs = Models_PCs$CMIP5_CanESM2_r1i1p1)

firstsEOFs <- mapply(function(x,y) Models_EOFS_preparation(whichPCs = 1:30, EOFs =x, PCs = y), x = Models_EOFs, y = Models_PCs,SIMPLIFY = FALSE)

cors_firstsEOFs <- lapply(firstsEOFs,cor, method = "spearman")


cors.cors.firstsEOFs <- matrix(nrow = length(cors_firstsEOFs),ncol = length(cors_firstsEOFs),dimnames = list(names(cors_firstsEOFs),names(cors_firstsEOFs)))
RMSE.cors.firstsEOFs <- matrix(nrow = length(cors_firstsEOFs),ncol = length(cors_firstsEOFs),dimnames = list(names(cors_firstsEOFs),names(cors_firstsEOFs)))
for (j in 1:length(cors_firstsEOFs)){
  cors.cors.firstsEOFs.i<- numeric(length(cors_firstsEOFs))
  RMSE.cors.firstsEOFs.i <- numeric(length(cors_firstsEOFs))
  for (i in 1:length(cors_firstsEOFs)){
    x1 <- cors_firstsEOFs[[i]]
    x2 <- cors_firstsEOFs[[j]]
    # x2 <- cors_gcms[["ncep_10d"]]
    cors.cors.firstsEOFs.i[i]<- cor(as.vector(x1),as.vector(x2))
    RMSE.cors.firstsEOFs.i[i] <- sqrt(sum((x1-x2)^2)/length(x1))
  }
  RMSE.cors.firstsEOFs[,j] <- RMSE.cors.firstsEOFs.i
  cors.cors.firstsEOFs[,j] <- cors.cors.firstsEOFs.i
}


# ind <- 1:length(gmt.hists.cmips) # UITKIJKEN ALS JE EERST DE LIJN HIERONDER DOET EN DAN WEER DEZE DAN IS IS HIJ AL OVERSCHREVEN

RMSE.cors.gcms.hists.int <- RMSE.cors.gcms.int[names(sets.hists)]
RMSE.cors.gcms.futures.int <- RMSE.cors.gcms.int[names(sets.futures)]
RMSE.cors.gcms.reans.int <- RMSE.cors.gcms.int[names(sets.reans)]

cors.cors.firstsEOFs.d <- cors.cors.firstsEOFs
bhat.gcms.d <- -log(helcoef)
diag(bhat.gcms.d) <- diag(cors.cors.firstsEOFs.d)<-NA
bhat.gcms.d[which(bhat.gcms.d== Inf, arr.ind = TRUE)]<- NA

cor(as.vector(bhat.gcms.d),as.vector(cors.cors.firstsEOFs.d),use = "pairwise.complete.obs")

plot(cors.cors.firstsEOFs[names(sets.futures),names(sets.futures)],-log(helcoef)[names(sets.futures),names(sets.futures)])
cor(use = "complete.obs",as.vector(cors.cors.firstsEOFs[names(sets.futures),names(sets.futures)]),as.vector(-log(helcoef)[names(sets.futures),names(sets.futures)]),method = "spearman")
points(cors.cors.firstsEOFs[names(sets.hists),names(sets.hists)],-log(helcoef)[names(sets.hists),names(sets.hists)], col = "blue")
cor(use = "complete.obs",as.vector(cors.cors.firstsEOFs[names(sets.hists),names(sets.hists)]),as.vector(bhat.gcms.d[names(sets.hists),names(sets.hists)]), method = "spearman")

points(cors.cors.firstsEOFs[names(sets.hists),names(sets.futures)],-log(helcoef)[names(sets.hists),names(sets.futures)], col = "purple")
cor(use = "complete.obs",as.vector(cors.cors.firstsEOFs[names(sets.hists),names(sets.futures)]),as.vector(bhat.gcms.d[names(sets.hists),names(sets.futures)]))

points(cors.cors.firstsEOFs[names(sets.reans),names(sets.reans)],-log(helcoef)[names(sets.reans),names(sets.reans)], col = "green")

points(cors.cors.firstsEOFs[names(sets.reans),names(sets.futures)],-log(helcoef)[names(sets.reans),names(sets.futures)], col = "red")
cor(use= "complete.obs",as.vector(cors.cors.firstsEOFs[names(sets.reans),names(sets.futures)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.futures)]),method = "spearman")
points(cors.cors.firstsEOFs[names(sets.reans),names(sets.hists)],-log(helcoef)[names(sets.reans),names(sets.hists)], col = "orange")
cor(as.vector(cors.cors.firstsEOFs[names(sets.reans),names(sets.hists)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.hists)]),method = "spearman")


cor(as.vector(bhat.gcms.d),as.vector(RMSE.cors.firstsEOFs),use = "pairwise.complete.obs")
plot(as.vector(bhat.gcms.d[,"interim_10d_akima_cubic"]),as.vector(cors.cors.firstsEOFs[,"interim_10d_akima_cubic"]))
cor(as.vector(bhat.gcms.d[,"interim_10d_akima_cubic"]),as.vector(cors.cors.firstsEOFs[,"interim_10d_akima_cubic"]), use = "complete.obs", method = "spearman")


plot(as.vector(-log(helcoef)["interim_10d_akima_cubic",]),as.vector(cors.cors.gcms["interim_10d_akima_cubic",]))
cor(as.vector(cors.cors.gcms["interim_10d_akima_cubic",]),as.vector(-log(helcoef)["interim_10d_akima_cubic",]),method = "pearson")

lm2 <- lm(as.vector(cors.cors.gcms["interim_10d_akima_cubic",])~as.vector(-log(helcoef)["interim_10d_akima_cubic",]),data.frame(as.vector(cors.cors.gcms["interim_10d_akima_cubic",]),as.vector(-log(helcoef)["interim_10d_akima_cubic",])))
lm2$`as.vector(cors.cors.gcms["interim_10d_akima_cubic", ])`
summary(lm2)$r.squared
# Models_EOFS_preparation$
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

##############################################################################################################
# Check if Bhattacharya distance has relationship with distance principal components: Method 2 EOF climatology
##############################################################################################################
Models_EOFs_clim <- lapply(Models_prinCompSet, EOF2clim,ind.var = 1,n.eofs = 2)
xy_grid_EOFs_clim <- lapply(Models_EOFs_clim,redim, drop = TRUE)
xy_mat_EOFs_clim <- lapply(xy_grid_EOFs_clim,function(x)x$Data)

xy_mat_EOFs_clim$CMIP5_CanESM2_r1i1p1[1,,]


quantity2clim

cors.CLIMs.firstsEOFs <- matrix(nrow = length(xy_mat_EOFs_clim),ncol = length(xy_mat_EOFs_clim),dimnames = list(names(xy_mat_EOFs_clim),names(xy_mat_EOFs_clim)))
RMSE.CLIMs.firstsEOFs <- matrix(nrow = length(xy_mat_EOFs_clim),ncol = length(xy_mat_EOFs_clim),dimnames = list(names(xy_mat_EOFs_clim),names(xy_mat_EOFs_clim)))

i <- 1
j <- 2


for (j in 1:length(xy_mat_EOFs_clim)){
  cors.CLIMs.firstsEOFs.i<- numeric(length(xy_mat_EOFs_clim))
  RMSE.CLIMs.firstsEOFs.i <- numeric(length(xy_mat_EOFs_clim))
  for (i in 1:length(xy_mat_EOFs_clim)){
    x1 <- xy_mat_EOFs_clim[[i]]
    x2 <- xy_mat_EOFs_clim[[j]]
    
    cors.CLIMs.firstsEOFs.i.k.l  <- matrix(nrow = dim(xy_mat_EOFs_clim[[1]])[1],ncol = dim(xy_mat_EOFs_clim[[1]])[1])
    RMSE.CLIMs.firstsEOFs.i.k.l <- matrix(nrow = dim(xy_mat_EOFs_clim[[1]])[1],ncol = dim(xy_mat_EOFs_clim[[1]])[1])
    
  
    for(k in 1:dim(xy_mat_EOFs_clim[[1]])[1]){
      for(l in 1:dim(xy_mat_EOFs_clim[[1]])[1]){
      cors.CLIMs.firstsEOFs.i.k.l[k,l]<- cor(as.vector(x1[k,,]),as.vector(x2[l,,]))
      RMSE.CLIMs.firstsEOFs.i.k.l[k,l]<- sqrt(sum((x1[k,,]-x2[l,,])^2)/length(x1[k,,]))
      }
    }
      
    # x2 <- cors_gcms[["ncep_10d"]]
    cors.CLIMs.firstsEOFs.i[i]<- sum(apply(cors.CLIMs.firstsEOFs.i.k.l, MARGIN = 2, FUN = function(x)max(abs((x)))))
    RMSE.CLIMs.firstsEOFs.i[i] <- sum(apply(RMSE.CLIMs.firstsEOFs.i.k.l[k,l], MARGIN = 2, FUN = min))
  }
  RMSE.CLIMs.firstsEOFs[,j] <- RMSE.CLIMs.firstsEOFs.i
  cors.CLIMs.firstsEOFs[,j] <- cors.CLIMs.firstsEOFs.i
}

cors.CLIMs.firstsEOFs.d <- cors.CLIMs.firstsEOFs
bhat.gcms.d <- -log(helcoef)
diag(bhat.gcms.d) <- diag(cors.CLIMS.firstsEOFs.d)<-NA
bhat.gcms.d[which(bhat.gcms.d== Inf, arr.ind = TRUE)]<- NA

cor(as.vector(bhat.gcms.d),as.vector(cors.CLIMs.firstsEOFs.d),use = "pairwise.complete.obs")

plot(cors.CLIMs.firstsEOFs[names(sets.futures),names(sets.futures)],-log(helcoef)[names(sets.futures),names(sets.futures)])
cor(use = "complete.obs",as.vector(cors.CLIMs.firstsEOFs[names(sets.futures),names(sets.futures)]),as.vector(-log(helcoef)[names(sets.futures),names(sets.futures)]),method = "spearman")
points(cors.CLIMs.firstsEOFs[names(sets.hists),names(sets.hists)],-log(helcoef)[names(sets.hists),names(sets.hists)], col = "blue")
cor(use = "complete.obs",as.vector(cors.CLIMs.firstsEOFs[names(sets.hists),names(sets.hists)]),as.vector(bhat.gcms.d[names(sets.hists),names(sets.hists)]), method = "spearman")

points(cors.cors.firstsEOFs[names(sets.hists),names(sets.futures)],-log(helcoef)[names(sets.hists),names(sets.futures)], col = "purple")
cor(use = "complete.obs",as.vector(cors.cors.firstsEOFs[names(sets.hists),names(sets.futures)]),as.vector(bhat.gcms.d[names(sets.hists),names(sets.futures)]))

points(cors.cors.firstsEOFs[names(sets.reans),names(sets.reans)],-log(helcoef)[names(sets.reans),names(sets.reans)], col = "green")

points(cors.cors.firstsEOFs[names(sets.reans),names(sets.futures)],-log(helcoef)[names(sets.reans),names(sets.futures)], col = "red")
cor(use= "complete.obs",as.vector(cors.cors.firstsEOFs[names(sets.reans),names(sets.futures)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.futures)]),method = "spearman")
points(cors.cors.firstsEOFs[names(sets.reans),names(sets.hists)],-log(helcoef)[names(sets.reans),names(sets.hists)], col = "orange")
cor(as.vector(cors.cors.firstsEOFs[names(sets.reans),names(sets.hists)]),as.vector(-log(helcoef)[names(sets.reans),names(sets.hists)]),method = "spearman")


cor(as.vector(bhat.gcms.d),as.vector(RMSE.cors.firstsEOFs),use = "pairwise.complete.obs")
plot(as.vector(bhat.gcms.d[,"interim_10d_akima_cubic"]),as.vector(cors.cors.firstsEOFs[,"interim_10d_akima_cubic"]))
cor(as.vector(bhat.gcms.d[,"interim_10d_akima_cubic"]),as.vector(cors.cors.firstsEOFs[,"interim_10d_akima_cubic"]), use = "complete.obs", method = "spearman")


plot(as.vector(-log(helcoef)["interim_10d_akima_cubic",]),as.vector( cors.CLIMs.firstsEOFs["interim_10d_akima_cubic",]))
cor(as.vector( cors.CLIMs.firstsEOFs["interim_10d_akima_cubic",]),as.vector(-log(helcoef)["interim_10d_akima_cubic",]),method = "pearson")

lm2 <- lm(as.vector(cors.cors.gcms["interim_10d_akima_cubic",])~as.vector(-log(helcoef)["interim_10d_akima_cubic",]),data.frame(as.vector(cors.cors.gcms["interim_10d_akima_cubic",]),as.vector(-log(helcoef)["interim_10d_akima_cubic",])))
lm2$`as.vector(cors.cors.gcms["interim_10d_akima_cubic", ])`
summary(lm2)$r.squared

#####################################################
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
cor.clim.int <- numeric(length(clim.hists))
i <- 1

clim.interim.c <- clim.reans$interim_10d_akima_cubic
clim.interim.c$Data <- clim.reans$interim_10d_akima_cubic$Data-273.15

for (i in 1:length(clim.hists)){
  if(max(clim.hists[[i]]$Data)>60){
    clim.hists[[i]]$Data <- keltoC(clim.hists[[i]]$Data)}

  x1 <- m.grid(clim.hists[[i]])$Data
  
  x2 <- m.grid(clim.interim.c)$Data
  # x1 <- clim.sd.hists[[i]]$Data
  # x2 <- clim.sd.reans$interim_10d_akima_cubic$Data
  RMSE.clim.int[i] <- sqrt(sum((x1-x2)^2)/length(x1))
  cor.clim.int[i] <- cor(as.vector(x1),as.vector(x2))
}

names(RMSE.clim.int)<- names(clim.hists)
names(cor.clim.int)<- names(clim.hists)


# save(RMSE.int.w,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/RMSE.int.w.rda")
bd.cmips.hist.vs.int
plot(bd.cmips.hist.vs.int[-which(RMSE.clim.int>0.12,arr.ind=TRUE)],RMSE.clim.int[-which(RMSE.clim.int>0.12,arr.ind=TRUE)])
cor(bd.cmips.hist.vs.int[-which(RMSE.clim.int>0.12,arr.ind=TRUE)],RMSE.clim.int[-which(RMSE.clim.int>0.12,arr.ind=TRUE)])
text(bd.cmips.hist.vs.int,RMSE.clim.int,labels = names(RMSE.clim.int))
plot(bd.cmips.hist.vs.int,RMSE.sd.grid.int)
cor(bd.cmips.hist.vs.int,RMSE.sd.grid.int)
text(bd.cmips.hist.vs.int,RMSE.sd.grid.int,labels = names(RMSE.sd.grid.int))

cor.clim.int[-which(cor.clim.int<0,arr.ind=TRUE)]

plot(bd.cmips.hist.vs.int[-which(cor.clim.int<0,arr.ind=TRUE)],cor.clim.int[-which(cor.clim.int<0,arr.ind=TRUE)])
cor(bd.cmips.hist.vs.int[-which(cor.clim.int<0,arr.ind=TRUE)],cor.clim.int[-which(cor.clim.int<0,arr.ind=TRUE)])
text(bd.cmips.hist.vs.int,cor.clim.int,labels = names(cor.clim.int))
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

