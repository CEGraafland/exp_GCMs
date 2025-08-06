#############################################################################
# Propgation of evidence CMIP5 rcp85 and CMIP6 ssp585
#############################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list=ls())
library(visualizeR)
library(RColorBrewer)
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
########################################################################################################
# Model Evaluation CMIP5
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda")
load("data/tas_rcp85_cmip5_left_10d_akima_cubic.rda")
load("data/tas_ssp585_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_ssp585_cmip6_left_10d_akima_cubic_corrected.rda")
load("../Data/Struct_learn/permutations.rda")
###########################################################################################
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
# Functions to load hciterations of models in a list 
############################################################################
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
################################################################################
# How do the models of various gcms explain the data of various gcms?
################################################################################
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
hc_gcms.future <- c(hc_gcms6.ssp585,hc_gcms6left.ssp585)
names(hc_gcms.future) <- shortnames.future
names(hc_gcms.future) <- names(cubicsets.future)

# Choose between 'own optimums' or constant magnitude
# modelsize <- own_optimums
modelsize <- 18
selection_hc_gcms.future <-  mapply(function(x,y) x[[grep((y),x)]], x = hc_gcms.future, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms.future <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms.future, y = data.future_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms.future_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms.future, y = data.future_gcms_anom_scaled, SIMPLIFY = FALSE)

# whichBNFits <-selection_fits_gcms_scaled
# whichData <- data_gcms_anom_scaled
# ordersets <- TRUE
# 
# loglik_selection_datasets <- matrix(data = NA, nrow = length(whichBNFits), ncol = length(whichData), dimnames = list(names(whichBNFits),names(whichData)))
# for(i in 1:length(whichBNFits)){
#   lo <- sapply(X = whichData, logLik, object = whichBNFits[[i]])
#   loglik_selection_datasets[i,] <- lo
# }

#######################################################################################
# Automatic V81 Overlap CMIP5 CMIP6
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


institutions.overlap.5 <- c("IPSL.CM5A","CESM1","EC.EARTH_r12i1p1","Can","MRI","MIROC.ES","MIROC5","CNRM","MPI","GFDL.ESM","Nor")

# cmip5
institutions <- rep(institutions.overlap.5,each =1)
CMIPS.5 <- rep(c("CMIP5"),length(institutions.overlap.5))
combinations.5 <- cbind(CMIPS.5,institutions)
namesort.5 <- unlist(apply(combinations.5,MARGIN = 1,function(x) grep(paste0(x[2]),names(cubicsets.future))))
names.overlap.rcp85 <- names(cubicsets.future)[namesort.5]


# cmip6
out <- which(names(cubicsets.future) == "CMIP6Amon_CanESM5_ssp585_r1i1p2f1")
cubicsets.future <-cubicsets.future[-out]
institutions.overlap <- c("IPSL","CNRM","CESM","EC","CanESM5","MRI","MIROC","MPI","GFDL","Nor")
institutions <- rep(institutions.overlap,each =1)
CMIPS.6 <- rep(c("CMIP6"),length(institutions.overlap))
combinations.6 <- cbind(CMIPS.6,institutions)
namesort.6 <- unlist(apply(combinations.6,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),names(cubicsets.future))))
names.overlap.ssp585 <- names(cubicsets.future)[namesort.6]



files5.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP5_rcp85/prop_",names(cubicsets.future)[namesort.5],"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files5.rcp85,function(x)load(file = x))
for (i in 1:length(propfiles)){load(file = files5.rcp85[i])}

files6.ssp585 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP6_ssp585/prop_",names(cubicsets.future)[namesort.6],"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files6.ssp585,function(x)load(file = x))
for (i in 1:length(propfiles)){load(file = files6.ssp585[i])}



posposmix <- list()
combinames <- names.overlap.rcp85
combinames <- names.overlap.ssp585
posposmix <- lapply(combinames, function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
all.equal(posposmix[[1]],posposmix[[2]])
for(i in 1:length(combinames)){
  posposmix[[i]] <- eval(parse(text = paste0("prop_",combinames[i],"_hc_",whichiterations,"i_",whichnode,"_equal2")))
}
names(posposmix) <- propfiles

griddata.future<- cubicsets.future[namesort.5]
griddata.future<- cubicsets.future[namesort.6]

prop_dif_pos <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, griddata = griddata.future, SIMPLIFY = FALSE)
enspropdifsposclims <- bindGrid(prop_dif_pos, dimension = c("member"))
enspropdifsposclims$Members <- combinames
enspropdifsposclims$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsposclims, names.attr = combinames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


files5.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP5_rcp85/propneg_",names(cubicsets.future)[namesort.5],"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files5.rcp85,function(x)load(file = x))
for (i in 1:length(propfiles)){load(file = files5.rcp85[i])}

files6.ssp585 <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/propagation/perm3/pos",whichnode,"neg/CMIP6_ssp585/propneg_",names(cubicsets.future)[namesort.6],"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files6.ssp585,function(x)load(file = x))
for (i in 1:length(propfiles)){load(file = files6.ssp585[i])}




posnegmix <- lapply(combinames, function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posnegmix
names(posnegmix) <- propfiles


prop_dif_neg <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata = griddata.future, SIMPLIFY = FALSE)
enspropdifsnegclims <- bindGrid(prop_dif_neg, dimension = c("member"))
enspropdifsnegclims$Members <- combinames
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsnegclims,names.attr = combinames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_neg[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))



prop_dif_pos_and_neg <- mapply(function(proppos,propneg,griddata) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),proppos = posposmix,propneg = posnegmix, griddata = griddata.future, SIMPLIFY = FALSE)
enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg, dimension = c("member"))

enspropdifsposandnegclims$Members <- combinames
col.b

col <- c(rev(col.b),col.r)
cmip6membernames<- gsub("Amon",combinames,replacement = "")

#cmip5
spatialPlot(enspropdifsposandnegclims,names.attr =combinames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

#cmip6
spatialPlot(enspropdifsposandnegclims,names.attr =cmip6membernames,as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

############################################################################################################
# Weighted vs unweighted
############################################################################################################
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
# "CMIP5_MPI.ESM.LR_r1i1p1" no está en historical. 
IDfutinhel <- unlist(sapply(gsub("_rcp85","",names.overlap.rcp85), function(x) which(row.names(hellinger_coefficients) == x)))
IDhelinfut <- unlist(sapply(row.names(hellinger_coefficients), function(x) which(gsub("_rcp85","",names.overlap.rcp85) == x)))
IDinterim <- which(row.names(hellinger_coefficients) == "interim_10d_akima_cubic")


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

sD <- 20
sS <- 25

overlap.rcp85.weights <- IPweight(-log(hellinger_coefficients)[c(IDfutinhel,IDinterim),c(IDfutinhel,IDinterim)],"interim_10d_akima_cubic",sD,sS)
IDproptoorderweights <- sapply(names(overlap.rcp85.weights), function(x) grep(x,names(prop_dif_pos_and_neg[IDhelinfut])))
names(prop_dif_pos_and_neg[IDhelinfut][IDproptoorderweights])

# 
# 
# redims_prop <- lapply(prop_dif_pos_and_neg,redim,drop = TRUE)
# mapply(function(x,y)quantity2clim(quantity = x, what = attr(y$Data,"climatology:fun"), ref.grid = y), x = redims_prop, y= prop_dif_pos_and_neg ) 
# 
# x <- prop_dif_pos_and_neg[IDhelinfut]$prop_CMIP5_CanESM2_r1i1p1_rcp85_hc_1700_1800i_V81_equal2
# y <- overlap.rcp85.weights[[1]]
# 
# x$Data * y 
# redim_x <- redim(x,drop = TRUE)
# spatialPlot(quantity2clim(redim_x$Data*y, what = attr(x$Data,"climatology:fun"), ref.grid = x),backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
#             col.regions= col,
#             set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# 
# all.equal(prop_dif_pos_and_neg[IDhelinfut],prop.cmip5.fut.weighted.grids)
#
# names(prop_dif_pos_and_neg[IDhelinfut])


prop.cmip5.fut.weighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg[IDhelinfut][IDproptoorderweights], y = overlap.rcp85.weights ,SIMPLIFY = FALSE)
prop.cmip5.fut.weighted <- Reduce('+', prop.cmip5.fut.weighted.grids)

newgrid <- prop_dif_pos_and_neg[[1]]
newgrid$Data <- prop.cmip5.fut.weighted
# all.equal(prop_dif_pos_and_neg$prop_CMIP5_CNRM.CM5_r1i1p1_rcp85_hc_1700_1800i_V81_equal2$Data,newgrid$Data)
spatialPlot(newgrid,names.attr = attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop.cmip5.fut.unweighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg[IDhelinfut][IDproptoorderweights], y = 1/length(overlap.rcp85.weights),SIMPLIFY = FALSE)
prop.cmip5.fut.unweighted <- Reduce('+', prop.cmip5.fut.unweighted.grids)

oldgrid <- prop_dif_pos_and_neg[[1]]
oldgrid$Data <- prop.cmip5.fut.unweighted

spatialPlot(oldgrid,names.attr = attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

compareoldnew <- bindGrid(oldgrid,newgrid, dimension = c("member"))
compareoldnew$Members <- c("unweighted mean propagation","weighted mean propagation")
spatialPlot(compareoldnew,names.attr = c("unweighted mean propagation","weighted mean propagation"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

###############################################################################
# Plot modelweights
###############################################################################
par(mar = c(8, 4, 0, 0))
plot(overlap.rcp85.weights,xaxt =  "n",xlab = "",ylab = "weights",ylim = c(0,0.5))
xlabels<-gsub("_r12i1p1","",gsub("_r1i1p1","",gsub("CMIP5_","",names(overlap.rcp85.weights))))
axis(1, at=1:length(overlap.rcp85.weights),labels = xlabels, col.axis="red", las=2)
abline(h = 1/length(overlap.rcp85.weights), lty = 2, col = "red")

text(13,0.12, expression(paste(sigma,"D=")))
text(14.2,0.12,paste(sD))
text(13,0.10, bquote(sigma*"S ="))
text(14.2,0.10,paste(sS))





#######################################################################################
# Automatic V81 All CMIP5 
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


files5.all.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"pos/CMIP5_rcp85/prop_",names(cubicsets.future),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files5.all.rcp85,function(x)load(file = x))
for (i in 1:length(propfiles)){load(file = files5.all.rcp85[i])}

posposmix <- list()
posposmix <- lapply(names(cubicsets.future), function(x) eval(parse(text = paste0("prop_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))

# for(i in 1:length(combinames)){
#   posposmix[[i]] <- eval(parse(text = paste0("prop_",combinames[i],"_hc_",whichiterations,"i_",whichnode,"_equal2")))
# }
names(posposmix) <- propfiles

griddata.future<- cubicsets.future[namesort.5]
griddata.future<- cubicsets.future[namesort.6]

prop_dif_pos <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posposmix, griddata = cubicsets.future, SIMPLIFY = FALSE)
enspropdifsposclims <- bindGrid(prop_dif_pos, dimension = c("member"))
enspropdifsposclims$Members <-names(cubicsets.future)
enspropdifsposclims$Members
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
spatialPlot(enspropdifsposclims, names.attr = names(cubicsets.future),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


files5.all.rcp85 <- paste0("results/propagation/perm3/pos",whichnode,"neg/CMIP5_rcp85/propneg_",names(cubicsets.future),"_hc_",whichiterations,"i_",whichnode,"_equal2.rda")
propfiles <- lapply(files5.all.rcp85,function(x)load(file = x))
for (i in 1:length(propfiles)){load(file = files5.all.rcp85[i])}



posnegmix <- lapply(names(cubicsets.future), function(x) eval(parse(text = paste0("propneg_",x,"_hc_",whichiterations,"i_",whichnode,"_equal2"))))
posnegmix
names(posnegmix) <- propfiles


prop_dif_neg <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = posnegmix, griddata =  cubicsets.future, SIMPLIFY = FALSE)
enspropdifsnegclims <- bindGrid(prop_dif_neg, dimension = c("member"))
enspropdifsnegclims$Members <- names(cubicsets.future)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(100)
spatialPlot(enspropdifsnegclims,names.attr = names(cubicsets.future),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_neg[[1]]$Data,"climatology:fun"))), at = seq(0.03,0.85,0.01),region = TRUE, col.regions= col.b,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))



prop_dif_pos_and_neg <- mapply(function(proppos,propneg,griddata) quantity2clim((proppos$with - proppos$without)-(propneg$with-propneg$without), paste0(attr(proppos$with, "probability"),"-", attr(proppos$without, "probability")," & ",attr(propneg$with, "probability"),"-", attr(propneg$without, "probability")), griddata),proppos = posposmix,propneg = posnegmix, griddata = cubicsets.future, SIMPLIFY = FALSE)
prop_dif_pos_and_neg_all <- prop_dif_pos_and_neg
enspropdifsposandnegclims <- bindGrid(prop_dif_pos_and_neg, dimension = c("member"))

enspropdifsposandnegclims$Members <- combinames
col.b

col <- c(rev(col.b),col.r)


#cmip5
spatialPlot(enspropdifsposandnegclims,names.attr =names(cubicsets.future),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos[[1]]$Data,"climatology:fun")," & ",attr(prop_dif_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


############################################################################################################
# Weighted vs unweighted
############################################################################################################
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
# "CMIP5_MPI.ESM.LR_r1i1p1" no está en historical. 
IDfutallinhel <- unlist(sapply(gsub("_rcp85","",names(cubicsets.future)), function(x) which(row.names(hellinger_coefficients) == x)))
IDhelinfutall <- unlist(sapply(row.names(hellinger_coefficients), function(x) which(gsub("_rcp85","",names(cubicsets.future)) == x)))
IDinterim <- which(row.names(hellinger_coefficients) == "interim_10d_akima_cubic")
IDncep <- which(row.names(hellinger_coefficients) == "ncep_10d")
IDjra55 <- which(row.names(hellinger_coefficients) == "JRA55_10d_akima_cubic")

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

sD <- 20
sS <- 25

rcp85.weights <- IPweight(-log(hellinger_coefficients)[c(IDfutallinhel,IDinterim),c(IDfutallinhel,IDinterim)],"interim_10d_akima_cubic",sD,sS)
IDproptoorderweights <- sapply(names(rcp85.weights), function(x) grep(x,names(prop_dif_pos_and_neg[IDhelinfutall])))
names(prop_dif_pos_and_neg[IDhelinfutall][IDproptoorderweights])

# MDS
fit <- cmdscale(-log(hellinger_coefficients)[c(IDfutallinhel,IDinterim),c(IDfutallinhel,IDinterim)],eig = TRUE)
fit # view results
fit
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x,y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS",type = "n")
text(x, y, labels = names(x), cex=.7) 

 
# redims_prop <- lapply(prop_dif_pos_and_neg,redim,drop = TRUE)
# mapply(function(x,y)quantity2clim(quantity = x, what = attr(y$Data,"climatology:fun"), ref.grid = y), x = redims_prop, y= prop_dif_pos_and_neg ) 
# 
# x <- prop_dif_pos_and_neg[IDhelinfut]$prop_CMIP5_CanESM2_r1i1p1_rcp85_hc_1700_1800i_V81_equal2
# y <- overlap.rcp85.weights[[1]]
# 
# x$Data * y 
# redim_x <- redim(x,drop = TRUE)
# spatialPlot(quantity2clim(redim_x$Data*y, what = attr(x$Data,"climatology:fun"), ref.grid = x),backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
#             col.regions= col,
#             set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# 
# all.equal(prop_dif_pos_and_neg[IDhelinfut],prop.cmip5.fut.weighted.grids)
#
# names(prop_dif_pos_and_neg[IDhelinfut])


prop.cmip5.fut.weighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg[IDhelinfutall][IDproptoorderweights], y = rcp85.weights ,SIMPLIFY = FALSE)
prop.cmip5.fut.weighted <- Reduce('+', prop.cmip5.fut.weighted.grids)

newgrid <- prop_dif_pos_and_neg[[1]]
newgrid$Data <- prop.cmip5.fut.weighted
# all.equal(prop_dif_pos_and_neg$prop_CMIP5_CNRM.CM5_r1i1p1_rcp85_hc_1700_1800i_V81_equal2$Data,newgrid$Data)
spatialPlot(newgrid,names.attr = attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop.cmip5.fut.unweighted.grids <-mapply(function(x,y) x$Data*y,x = prop_dif_pos_and_neg[IDhelinfutall][IDproptoorderweights], y = 1/length(rcp85.weights),SIMPLIFY = FALSE)
prop.cmip5.fut.unweighted <- Reduce('+', prop.cmip5.fut.unweighted.grids)

oldgrid <- prop_dif_pos_and_neg[[1]]
oldgrid$Data <- prop.cmip5.fut.unweighted

spatialPlot(oldgrid,names.attr = attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

compareoldnew <- bindGrid(oldgrid,newgrid, dimension = c("member"))
compareoldnew$Members <- c("unweighted mean propagation","weighted mean propagation")
w.prop.plot <- spatialPlot(compareoldnew,names.attr = c("unweighted mean propagation","weighted mean propagation"),as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0(attr(prop_dif_pos_and_neg[[1]]$Data,"climatology:fun"))), at = seq(-0.85,0.85,0.01),region = TRUE, 
            col.regions= col,
            set.min = -0.85,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/evidence_propagation/CMIP5_gmt_weighted/CMIP5_1700_1800_BNgcms_gmt_weighted_plusV81pos_and_neg_perm3_sD_",sD,"_sS_",sS,".pdf")
pdf(plotname)
w.prop.plot
dev.off()
###############################################################################
# Plot modelweights cmip5 all future
###############################################################################
par(mar = c(8, 4, 0, 0))
plot(rcp85.weights,xaxt =  "n",xlab = "",ylab = "weights",ylim = c(0,0.14))
xlabels<-gsub("_r12i1p1","",gsub("_r1i1p1","",gsub("CMIP5_","",names(rcp85.weights))))
axis(1, at=1:length(rcp85.weights),labels = xlabels, col.axis="red", las=2)
abline(h = 1/length(rcp85.weights), lty = 2, col = "red")

text(13,0.12, expression(paste(sigma,"D=")))
text(14.2,0.12,paste(sD))
text(13,0.10, bquote(sigma*"S ="))
text(14.2,0.10,paste(sS))

# Load RMSE:int
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
load("results/RMSE.int.rda")
load("results/RMSE.int.w.rda")
-log(hellinger_coefficients)[c(IDfutallinhel,IDinterim),c(IDfutallinhel,IDinterim)]
p1<-names(hellinger_coefficients[IDfutallinhel,IDinterim])
a<- sapply(p1, function(x) grep(x,names(RMSE.int.w)),simplify = TRUE)


plot(RMSE.int.w[a],-log(hellinger_coefficients[IDfutallinhel,IDinterim]))
cor(RMSE.int.w[a],-log(hellinger_coefficients[IDfutallinhel,IDinterim]), method = "pearson")

plot(RMSE.int.w[a],-log(hellinger_coefficients[IDfutallinhel,IDncep]))
cor(RMSE.int.w[a],-log(hellinger_coefficients[IDfutallinhel,IDncep]), method = "pearson")

##################################################################################
# area weighted RMSE between historical propagation and future propagation cmip5s :
# Make new R file with historical and future in one. 
##################################################################################
props.cmip5s.reans<- prop_dif_pos_and_neg_all$
props.cmip5s <- props.cmip5s.reans[1:(length(props.cmip5s.reans)-3)]



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

weighted.rcp85.props.cmip5s <- lapply(prop_dif_pos_and_neg_all,gtimes.grid)

RMSE.delta.w <- numeric(length(weighted.rcp85.props.cmip5s))
for (i in 1:length(weighted.rcp85.props.cmip5s)){
  x1 <- redim(weighted.rcp85.props.cmip5s[[i]],drop = TRUE)$Data
  RMSE.delta.w[i] <- sqrt(sum((x1-x1)^2)/length(x1))
}

names(RMSE.int.w)<- names(cubicsets)[-ind6][1:(length(names(cubicsets)[-ind6])-3)]
save(RMSE.int.w,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/RMSE.int.w.rda")

