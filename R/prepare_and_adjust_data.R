library(transformeR)
library(visualizeR)
library(gridExtra)
library(loadeR)

source("/oceano/gmeteo/WORK/lisette/Trabajo/creds")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/historical.rda")
load("data/tas_historical_10d_akima_cubic.rda")
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
######################################
# Make multiGrid of tas_historical_10d_akima_cubic
# View multiGrid
######################################

tas_historical_10d_akima_cubic.ens <- bindGrid(tas_historical_10d_akima_cubic, dimension = "member")
modnames <- gsub(gsub(UDG.datasets("CMIP5.*historical")$name[-6], pattern = "_historical", replacement = ""),
                 pattern = "-", replacement = ".")
                 
tas_historical_10d_akima_cubic.ens$Members <- modnames

spatialPlot(climatology(tas_historical_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)

####¿qué pasa con el IPSL en el 180?
# 'same type of coordinates es others'
di <- dataInventory(UDG.datasets("CMIP5.*IPSL.*historical")$name)
di2 <- dataInventory(UDG.datasets("CMIP5.*CNRM.CM5.*historical")$name)
di$tas$Dimensions
di2$tas$Dimensions

# After second time interpolation the error persists: 
tas_historical_10d_akima_cubic_3 <- interpGrid(historical[[3]],
                                               new.coordinates = getGrid(tas_interim_10dnew), 
                                               method = 'bilinear', bilin.method = 'akima', linear = FALSE,
                                               extrap = TRUE)


all.equal(tas_historical_10d_akima_cubic_3,tas_historical_10d_akima_cubic[[3]]) # = true
# In the original grid no habia un error en este sitio:
spatialPlot(climatology(historical[[3]]), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)



#################################################################################################
# Add extra first data slice which is the same as last dataslice (180) to IPSL (historical[[3]])
# Then add extra x coordinate to xyCoords in historical3
#################################################################################################
historical3 <- historical[[3]]
data.historical3 <- historical3$Data
# Make newdata with one dimension in longitude more
newdata.historical3 <- array(data = NA, dim = dim(data.historical3) + c(0,0,0,0,1), dimnames = dimnames(data.historical3))
newdata.historical3[,,,,2:145]<- data.historical3
# now equalize first new longitude with last already existing longitude
newdata.historical3[,,,,1]<-data.historical3[,,,,144]
# give same dimensions
attr(newdata.historical3,"dimensions") <-  attr(historical3$Data,"dimensions") 
# replace data 
historical3$Data <- newdata.historical3
# Transmit change to coordinates (create extra x coordinate -180)
newxCoords.historical3 <- vector(mode = "numeric", length = length(historical3$xyCoords$x)+1)
newxCoords.historical3[2:145] <- historical3$xyCoords$x
newxCoords.historical3[1] <- -180
historical3$xyCoords$x <- newxCoords.historical3

# do interpolation with new historical 3
tas_historical_10d_akima_cubic_3new <- interpGrid(historical3,
                                               new.coordinates = getGrid(tas_interim_10dnew), 
                                               method = 'bilinear', bilin.method = 'akima', linear = FALSE,
                                               extrap = TRUE)


# save(tas_historical_10d_akima_cubic_3new, file = "data/tas_historical_10d_akima_cubic_3new.rda")
load("data/tas_historical_10d_akima_cubic_3new.rda")


spatialPlot(climatology(tas_historical_10d_akima_cubic[[3]]))
spatialPlot(climatology(tas_historical_10d_akima_cubic_3new))
# observe little change in other values: Method bilinear akima ??
#################################
# Make new tas_historical_10d_akima_cubic_corrected
#################################

tas_historical_10d_akima_cubic_corrected <- tas_historical_10d_akima_cubic
tas_historical_10d_akima_cubic_corrected[[3]] <- tas_historical_10d_akima_cubic_3new

tas_historical_10d_akima_cubic_corrected.ens <- bindGrid(tas_historical_10d_akima_cubic_corrected, dimension = "member")
modnames <- gsub(gsub(UDG.datasets("CMIP5.*historical")$name[-6], pattern = "_historical", replacement = ""),
                 pattern = "-", replacement = ".")

tas_historical_10d_akima_cubic_corrected.ens$Members <- modnames
# seems to work:
spatialPlot(climatology(tas_historical_10d_akima_cubic_corrected.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)


modnames <- gsub(gsub(UDG.datasets("CMIP5.*historical")$name[-6], pattern = "_historical", replacement = ""),
                 pattern = "-", replacement = ".")
names(tas_historical_10d_akima_cubic_corrected) <- modnames
save(tas_historical_10d_akima_cubic_corrected, file = "data/tas_historical_10d_akima_cubic_corrected.rda")

#################################################################################################
# CMIP5 extra
#################################################################################################
######################################
# Make multiGrid of tas_historical_10d_akima_cubic
# View multiGrid
# Seems ok 
######################################
load("data/historical_cmip5_extra.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")

tas_historical_cmip5_extra_10d_akima_cubic
tas_historical_cmip5_extra_10d_akima_cubic.ens <- bindGrid(tas_historical_cmip5_extra_10d_akima_cubic, dimension = "member")

tas_historical_cmip5_extra_10d_akima_cubic.ens$Members <- names(historical_cmip5_extra)

spatialPlot(climatology(tas_historical_cmip5_extra_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rev.colors = TRUE)
#################################################################################################
# CMIP5 earth
#################################################################################################
######################################
# Make multiGrid of tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic
# View multiGrid
# Seems ok 
######################################
load("data/historical_cmip5_earth_r1r2r12.rda")
load("data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda")

tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic
tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.ens <- bindGrid(tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic, dimension = "member")


tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.ens$Members <- names(historical_cmip5_earth_r1r2r12)

spatialPlot(climatology(tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180, rev.colors = TRUE)
str(climatology(tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.ens))
dev.off()
#################################################################################################
# CMIP5 left
#################################################################################################
######################################
# Make multiGrid of tas_historical_cmip5_left_10d_akima_cubic
# View multiGrid
# Seems ok 
######################################
# load("data/historical_cmip5_extra.rda"):
dataset_list <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left", full.names = T)
dataset_names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left")
dataset_names <- gsub(".rda", "", dataset_names)
historical_cmip5_left <- lapply(dataset_list, function(x){get(load(x))})
names(historical_cmip5_left) <- dataset_names

load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")

tas_historical_cmip5_left_10d_akima_cubic
tas_historical_cmip5_left_10d_akima_cubic.ens <- bindGrid(tas_historical_cmip5_left_10d_akima_cubic, dimension = "member")

tas_historical_cmip5_left_10d_akima_cubic.ens$Members <- names(historical_cmip5_left)

spatialPlot(climatology(tas_historical_cmip5_left_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rev.colors = TRUE)


#################################################################################################
# CMIP5 FUTURE
#################################################################################################
######################################
# Make multiGrid of tas_rcp85_cmip5_10d_akima_cubic
# View multiGrid
# Seems ok 
######################################
# load("data/rcp85_cmip5.rda"):
rcp85_cmip5_list <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_2071_2100", full.names = T)
rcp85_cmip5_names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_2071_2100")
rcp85_cmip5_names <- gsub(".rda", "", rcp85_cmip5_names)
rcp85_cmip5 <- lapply(rcp85_cmip5_list, function(x){get(load(x))})
names(rcp85_cmip5) <- rcp85_cmip5_names

load("data/tas_rcp85_cmip5_10d_akima_cubic.rda")

tas_rcp85_cmip5_10d_akima_cubic
tas_rcp85_cmip5_10d_akima_cubic.ens <- bindGrid(tas_rcp85_cmip5_10d_akima_cubic, dimension = "member")

tas_rcp85_cmip5_10d_akima_cubic.ens$Members <- names(rcp85_cmip5)

spatialPlot(climatology(tas_rcp85_cmip5_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rev.colors = TRUE)

# MEMBER 12 has problems

#################################################################################################
# Add extra first data slice which is the same as last dataslice (180) to IPSL (CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85)
# Then add extra x coordinate to xyCoords in newCMIP5_IPSL.CM5A.MR_r1i1p1_rcp85
#################################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_2071_2100/CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85.rda")

CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$Data
historical3 <- historical[[3]]
data.historical3 <- historical3$Data
# Make newdata with one dimension in longitude more
newdata.rcp853.IPSL.CM5A.MR <- array(data = NA, dim = dim(CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$Data) + c(0,0,1))
newdata.rcp853.IPSL.CM5A.MR[,,2:145]<- CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$Data
# now equalize first new longitude with last already existing longitude
newdata.rcp853.IPSL.CM5A.MR[,,1]<-CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$Data[,,144]
# give same dimensions
attr(newdata.rcp853.IPSL.CM5A.MR,"dimensions") <-  attr(CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$Data,"dimensions") 
# replace data 
newCMIP5_IPSL.CM5A.MR_r1i1p1_rcp85 <- CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85
newCMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$Data <- newdata.rcp853.IPSL.CM5A.MR
# Transmit change to coordinates (create extra x coordinate -180)
newxCoords.CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85 <- vector(mode = "numeric", length = length(CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$xyCoords$x)+1)
newxCoords.CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85[2:145] <- CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$xyCoords$x
newxCoords.CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85[1] <- -180
newCMIP5_IPSL.CM5A.MR_r1i1p1_rcp85$xyCoords$x <- newxCoords.CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85

# do interpolation with newCMIP5_IPSL.CM5A.MR_r1i1p1_rcp85
tas_rcp85_newCMIP5_IPSL.CM5A.MR_r1i1p1_10d_akima_cubic <- interpGrid(newCMIP5_IPSL.CM5A.MR_r1i1p1_rcp85,
                                                  new.coordinates = getGrid(tas_interim_10dnew), 
                                                  method = 'bilinear', bilin.method = 'akima', linear = FALSE,
                                                  extrap = TRUE)


# save(tas_rcp85_newCMIP5_IPSL.CM5A.MR_r1i1p1_10d_akima_cubic, file = "data/tas_rcp85_newCMIP5_IPSL.CM5A.MR_r1i1p1_10d_akima_cubic.rda")
load("data/tas_rcp85_newCMIP5_IPSL.CM5A.MR_r1i1p1_10d_akima_cubic.rda")


spatialPlot(climatology(tas_rcp85_cmip5_10d_akima_cubic[[12]]))
spatialPlot(climatology(tas_rcp85_newCMIP5_IPSL.CM5A.MR_r1i1p1_10d_akima_cubic))
# observe little change in other values: Method bilinear akima ??
######################################################
# Make new tas_rcp85_cmip5_10d_akima_cubic_corrected
######################################################

tas_rcp85_cmip5_10d_akima_cubic_corrected <- tas_rcp85_cmip5_10d_akima_cubic
tas_rcp85_cmip5_10d_akima_cubic_corrected[[12]] <- tas_rcp85_newCMIP5_IPSL.CM5A.MR_r1i1p1_10d_akima_cubic

tas_rcp85_cmip5_10d_akima_cubic_corrected.ens <- bindGrid(tas_rcp85_cmip5_10d_akima_cubic_corrected, dimension = "member")

# seems to work:
spatialPlot(climatology(tas_rcp85_cmip5_10d_akima_cubic_corrected.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)

save(tas_rcp85_cmip5_10d_akima_cubic_corrected, file = "data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda")
#################################################################################################
# CMIP5 left FUTURE
#################################################################################################
######################################
# Make multiGrid of tas_rcp85_cmip5_left_10d_akima_cubic
# View multiGrid
# Seems ok 
######################################
# load("data/rcp85_cmip5.rda"):
rcp85_cmip5_left_list <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_left_2071_2100", full.names = T)
rcp85_cmip5_left_names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_left_2071_2100")
rcp85_cmip5_left_names <- gsub(".rda", "", rcp85_cmip5_left_names)
rcp85_cmip5_left <- lapply(rcp85_cmip5_left_list, function(x){get(load(x))})
names(rcp85_cmip5_left) <- rcp85_cmip5_left_names

load("data/tas_rcp85_cmip5_left_10d_akima_cubic.rda")

# CMIP5_bcc.csm1.1.m_r1i1p1_rcp85 heeft eerder begin en eind jaar; kan niet in het ensemble
sapply(tas_rcp85_cmip5_left_10d_akima_cubic,function(x)x$Dates$end[length(x$Dates$end)])
tas_rcp85_cmip5_left_10d_akima_cubic.ens <- bindGrid(tas_rcp85_cmip5_left_10d_akima_cubic[-3], dimension = "member")
tas_rcp85_cmip5_left_10d_akima_cubic.ens$Members <- names(rcp85_cmip5_left)[-3]

spatialPlot(climatology(tas_rcp85_cmip5_left_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rev.colors = TRUE)
spatialPlot(climatology(tas_rcp85_cmip5_left_10d_akima_cubic[[3]]), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rev.colors = TRUE)
#################################################################################################
# CMIP6 historical
#################################################################################################
######################################
# Make multiGrid of tas_historical_10d_akima_cubic
# View multiGrid
# IPSL is wrong.
######################################
load("data/historical_cmip6.rda")
load("data/tas_historical_cmip6_10d_akima_cubic.rda")

tas_historical_cmip6_10d_akima_cubic
names(tas_historical_cmip6_10d_akima_cubic)<-names(historical_cmip6)

tas_historical_cmip6_10d_akima_cubic.ens <- bindGrid(tas_historical_cmip6_10d_akima_cubic, dimension = "member")
tas_historical_cmip6_10d_akima_cubic.ens$Members <- names(historical_cmip6)

spatialPlot(climatology(tas_historical_cmip6_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180, rev.colors = TRUE)
str(climatology(tas_historical_cmip6_10d_akima_cubic.ens))
dev.off()
####¿qué pasa con el IPSL en el 180?
# 'same type of coordinates es others'
di <- dataInventory(UDG.datasets("CMIP6.*IPSL.*historical"))
din <- UDG.datasets("CMIP6Amon.*CanESM5.*historical")
din
di2 <- dataInventory(din$CMIP6Amon[1])
di$tas$Dimensions
di2$tas$Dimensions


# In the original grid no habia un error en este sitio:
historical_cmip6$CMIP6Amon_IPSL.CM6A.LR_r1i1p1f1_historical
spatialPlot(climatology(historical_cmip6$CMIP6Amon_IPSL.CM6A.LR_r1i1p1f1_historical), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)

#################################################################################################
# Add extra first data slice which is the same as last dataslice (180) to IPSL CMIP6 (historical_cmip6)
# Then add extra x coordinate to xyCoords in historical3
#################################################################################################
historical_cmip6$CMIP6Amon_IPSL.CM6A.LR_r1i1p1f1_historical
historical_ipsl <- historical_cmip6$CMIP6Amon_IPSL.CM6A.LR_r1i1p1f1_historical
data.historical_ipsl <- historical_ipsl$Data
data.historical_ipsl[1,,]
historical_ipsl$xyCoords
str(data.historical_ipsl)
# Make newdata with one dimension in longitude more
newdata.historical_ipsl <- array(data = NA, dim = dim(data.historical_ipsl) + c(0,0,1), dimnames = dimnames(data.historical_ipsl))
newdata.historical_ipsl[,,2:145]<- data.historical_ipsl
# now equalize first new longitude with last already existing longitude
newdata.historical_ipsl[,,1] <-data.historical_ipsl[,,144]
# give same dimensions
attr(newdata.historical_ipsl,"dimensions") <-  attr(historical_ipsl$Data,"dimensions") 
# replace data 
historical_ipsl$Data <- newdata.historical_ipsl
# Transmit change to coordinates (create extra x coordinate -180)
newxCoords.historical_ipsl <- vector(mode = "numeric", length = length(historical_ipsl$xyCoords$x)+1)
newxCoords.historical_ipsl[2:145] <- historical_ipsl$xyCoords$x
newxCoords.historical_ipsl[1] <- -180
historical_ipsl$xyCoords$x <- newxCoords.historical_ipsl

# do interpolation with new historical 3
tas_historical_cmip6_10d_akima_cubic
tas_historical_cmip6_10d_akima_cubic_ipsl_new <- interpGrid(historical_ipsl,
                                                  new.coordinates = getGrid(tas_interim_10dnew), 
                                                  method = 'bilinear', bilin.method = 'akima', linear = FALSE,
                                                  extrap = TRUE)


save(tas_historical_cmip6_10d_akima_cubic_ipsl_new, file = "data/tas_historical_cmip6_10d_akima_cubic_ipsl_new.rda")


spatialPlot(climatology(tas_historical_cmip6_10d_akima_cubic$CMIP6Amon_IPSL.CM6A.LR_r1i1p1f1_historical), lonCenter = 180,backdrop.theme = "coastline")
spatialPlot(climatology(tas_historical_cmip6_10d_akima_cubic_ipsl_new), lonCenter = 180,backdrop.theme = "coastline")
# observe little change in other values: Method bilinear akima ??
###############################################################################
# Make new tas_historical_cmip6_10d_akima_cubic_corrected
###############################################################################

tas_historical_cmip6_10d_akima_cubic_corrected <- tas_historical_cmip6_10d_akima_cubic
tas_historical_cmip6_10d_akima_cubic_corrected$CMIP6Amon_IPSL.CM6A.LR_r1i1p1f1_historical <- tas_historical_cmip6_10d_akima_cubic_ipsl_new
tas_historical_cmip6_10d_akima_cubic_corrected.ens <- bindGrid(tas_historical_cmip6_10d_akima_cubic_corrected, dimension = "member")

names(historical_cmip6)
# modnames <- gsub(gsub(UDG.datasets("CMIP5.*historical")$name[-6], pattern = "_historical", replacement = ""),
                 # pattern = "-", replacement = ".")

tas_historical_cmip6_10d_akima_cubic_corrected.ens$Members <- names(historical_cmip6)
# seems to work:
spatialPlot(climatology(tas_historical_cmip6_10d_akima_cubic_corrected.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)


save(tas_historical_cmip6_10d_akima_cubic_corrected, file = "data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
#################################################################################################
# CMIP6 left historical
#################################################################################################
######################################
# Make multiGrid of tas_historical_cmip6_left_10d_akima_cubic
# View multiGrid
# Seems ok except from BCC.CSM2.MR
######################################
library(convertR)
load("data/tas_historical_cmip6_left_10d_akima_cubic.rda")
# transform to Celcius
tas_historical_cmip6_left_10d_akima_cubic_degC <- lapply(tas_historical_cmip6_left_10d_akima_cubic, udConvertGrid, new.units = "degC")
tas_historical_cmip6_left_10d_akima_cubic.ens <- bindGrid(tas_historical_cmip6_left_10d_akima_cubic_degC, dimension = "member")

require(RColorBrewer)
spatialPlot(climatology(tas_historical_cmip6_left_10d_akima_cubic.ens, by.member = TRUE), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rev.colors = TRUE,
            at = seq(-60,80,10))
# filipines out of range
# member 1 max 5, 6 min out of range.
lapply(tas_historical_cmip6_left_10d_akima_cubic, function(x)max(x$Data)) # member 1
lapply(tas_historical_cmip6_left_10d_akima_cubic, function(x)min(x$Data)) # member 4???


# fix filipines coordinate CMIP6Amon_BCC.CSM2.MR_historical_r1i1p1f1
dmat <-tas_historical_cmip6_left_10d_akima_cubic$CMIP6Amon_BCC.CSM2.MR_historical_r1i1p1f1$Data
fail <- which(dmat == max(dmat),arr.ind = TRUE)

dmat[fail+c(0,1,0)]#lat up
dmat[fail+c(0,1,-1)]#lat up west
dmat[fail+c(0,1,1)]# lat up east
dmat[fail+c(0,-1,-1)]# lat down west
dmat[fail+c(0,-1,1)]# lat down east
dmat[fail+c(0,-1,0)]# lat down
dmat[fail+c(0,0,-1)]# lon west
dmat[fail+c(0,0,1)]# lon east

dmat[fail] <- (dmat[fail+c(0,1,0)] + 
  dmat[fail+c(0,1,-1)] + #lat up
  dmat[fail+c(0,1,1)] + #lat up west
  dmat[fail+c(0,-1,-1)]+ # lat down west
  dmat[fail+c(0,-1,1)] + # lat down east
  dmat[fail+c(0,-1,0)] + # lat down
  dmat[fail+c(0,0,-1)] + # lon west
  dmat[fail+c(0,0,1)]# lon east
)/8

tas_historical_cmip6_left_10d_akima_cubic_corrected <- tas_historical_cmip6_left_10d_akima_cubic
tas_historical_cmip6_left_10d_akima_cubic_corrected$CMIP6Amon_BCC.CSM2.MR_historical_r1i1p1f1$Data <- dmat

# check again
# transform to Celcius
tas_historical_cmip6_left_10d_akima_cubic_corrected_degC <- lapply(tas_historical_cmip6_left_10d_akima_cubic_corrected, udConvertGrid, new.units = "degC")
tas_historical_cmip6_left_10d_akima_cubic_corrected.ens <- bindGrid(tas_historical_cmip6_left_10d_akima_cubic_corrected_degC, dimension = "member")

require(RColorBrewer)
spatialPlot(climatology(tas_historical_cmip6_left_10d_akima_cubic_corrected.ens, by.member = TRUE), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rev.colors = TRUE,
            at = seq(-60,60,10))



save(tas_historical_cmip6_left_10d_akima_cubic_corrected, file = "data/tas_historical_cmip6_left_10d_akima_cubic_corrected.rda")
#################################################################################################
# CMIP6 Future ssp585
# CHECK the same as above 
#################################################################################################
######################################
# Make multiGrid of tas_ssp585_cmip6_10d_akima_cubic
# View multiGrid
# IPSL is wrong.
######################################
load("data/ssp585_cmip6.rda")
load("data/tas_ssp585_cmip6_10d_akima_cubic.rda")

sapply(tas_ssp585_cmip6_10d_akima_cubic,function(x)x$Dates$end[length(x$Dates$end)])
tas_ssp585_cmip6_10d_akima_cubic.ens <- bindGrid(tas_ssp585_cmip6_10d_akima_cubic, dimension = "member")
tas_ssp585_cmip6_10d_akima_cubic.ens$Members <- names(tas_ssp585_cmip6_10d_akima_cubic)

spatialPlot(climatology(tas_ssp585_cmip6_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)

####¿qué pasa con el IPSL en el 180?
# 'same type of coordinates es others'
din1<- UDG.datasets("CMIP6Amon.*IPSL.*ssp585")
di1 <- dataInventory(din1[[1]][4])
din2 <- UDG.datasets("CMIP6Amon.*CanESM5.*ssp585")
din2
di2 <- dataInventory(din2[[1]][43])
di1$tas$Dimensions
di2$tas$Dimensions

# In the original grid no habia un error en este sitio:
ssp585_cmip6$CMIP6Amon_IPSL.CM6A.LR_ssp585_r1i1p1f1
spatialPlot(climatology(ssp585_cmip6$CMIP6Amon_IPSL.CM6A.LR_ssp585_r1i1p1f1), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)
#################################################################################################
# Add extra first data slice which is the same as last dataslice (180) to IPSL CMIP6 (ssp85_cmip6)
# Then add extra x coordinate to xyCoords in ssp585
#################################################################################################
ssp585_cmip6$CMIP6Amon_IPSL.CM6A.LR_ssp585_r1i1p1f1
ssp585_ipsl <- ssp585_cmip6$CMIP6Amon_IPSL.CM6A.LR_ssp585_r1i1p1f1
data.ssp585_ipsl <- ssp585_ipsl$Data
data.ssp585_ipsl[1,,]
historical_ipsl$xyCoords
str(data.ssp585_ipsl)
# Make newdata with one dimension in longitude more
newdata.ssp585_ipsl <- array(data = NA, dim = dim(data.ssp585_ipsl) + c(0,0,1), dimnames = dimnames(data.ssp585_ipsl))
newdata.ssp585_ipsl[,,2:145]<- data.ssp585_ipsl
# now equalize first new longitude with last already existing longitude
newdata.ssp585_ipsl[,,1] <-data.ssp585_ipsl[,,144]
# give same dimensions
attr(newdata.ssp585_ipsl,"dimensions") <-  attr(ssp585_ipsl$Data,"dimensions") 
# replace data 
ssp585_ipsl$Data <- newdata.ssp585_ipsl
# Transmit change to coordinates (create extra x coordinate -180)
newxCoords.ssp585_ipsl <- vector(mode = "numeric", length = length(ssp585_ipsl$xyCoords$x)+1)
newxCoords.ssp585_ipsl[2:145] <- ssp585_ipsl$xyCoords$x
newxCoords.ssp585_ipsl[1] <- -180
ssp585_ipsl$xyCoords$x <- newxCoords.ssp585_ipsl

# do interpolation with new ssp585
tas_ssp585_cmip6_10d_akima_cubic
tas_ssp585_cmip6_10d_akima_cubic_ipsl_new <- interpGrid(ssp585_ipsl,
                                                            new.coordinates = getGrid(tas_interim_10dnew), 
                                                            method = 'bilinear', bilin.method = 'akima', linear = FALSE,
                                                            extrap = TRUE)


save(tas_ssp585_cmip6_10d_akima_cubic_ipsl_new, file = "data/tas_ssp585_cmip6_10d_akima_cubic_ipsl_new.rda")


spatialPlot(climatology(tas_ssp585_cmip6_10d_akima_cubic$CMIP6Amon_IPSL.CM6A.LR_ssp585_r1i1p1f1), lonCenter = 180,backdrop.theme = "coastline")
spatialPlot(climatology(tas_ssp585_cmip6_10d_akima_cubic_ipsl_new), lonCenter = 180,backdrop.theme = "coastline")
# observe little change in other values: Method bilinear akima ??
#####################################################################################
# Make new tas_ssp585_cmip6_10d_akima_cubic_corrected
#####################################################################################

tas_ssp585_cmip6_10d_akima_cubic_corrected <- tas_ssp585_cmip6_10d_akima_cubic
tas_ssp585_cmip6_10d_akima_cubic_corrected$CMIP6Amon_IPSL.CM6A.LR_ssp585_r1i1p1f1 <- tas_ssp585_cmip6_10d_akima_cubic_ipsl_new
tas_ssp585_cmip6_10d_akima_cubic_corrected.ens <- bindGrid(tas_ssp585_cmip6_10d_akima_cubic_corrected, dimension = "member")

names(ssp585_cmip6)
# modnames <- gsub(gsub(UDG.datasets("CMIP5.*historical")$name[-6], pattern = "_historical", replacement = ""),
# pattern = "-", replacement = ".")

tas_ssp585_cmip6_10d_akima_cubic_corrected.ens$Members <- names(ssp585_cmip6)
# seems to work:
spatialPlot(climatology(tas_ssp585_cmip6_10d_akima_cubic_corrected.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180)


save(tas_ssp585_cmip6_10d_akima_cubic_corrected, file = "data/tas_ssp585_cmip6_10d_akima_cubic_corrected.rda")
#######################################################################################
# CMIP6 ssp585 left FUTURE
#######################################################################################
######################################
# Make multiGrid of tas_ssp585_cmip6_left_10d_akima_cubic
# View multiGrid
# Seems ok 
######################################
load("data/tas_ssp585_cmip6_left_10d_akima_cubic.rda")

# begin en eind jaar;
sapply(tas_ssp585_cmip6_left_10d_akima_cubic,function(x)x$Dates$end[length(x$Dates$end)])
# transform to Celcius
tas_ssp585_cmip6_left_10d_akima_cubic_degC <- lapply(tas_ssp585_cmip6_left_10d_akima_cubic, udConvertGrid, new.units = "degC")
tas_ssp585_cmip6_left_10d_akima_cubic.ens <- bindGrid(tas_ssp585_cmip6_left_10d_akima_cubic_degC, dimension = "member")

spatialPlot(climatology(tas_ssp585_cmip6_left_10d_akima_cubic.ens), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rlonCenter = 180,rev.colors = TRUE,
            at = seq(-80,60,10))

# member 1 max 5, 6 min out of range.
lapply(tas_ssp585_cmip6_left_10d_akima_cubic, function(x)max(x$Data)) # member 1
lapply(tas_ssp585_cmip6_left_10d_akima_cubic, function(x)min(x$Data)) # member 5 & 6

# Fix member 1 maximum
dmat1 <- tas_ssp585_cmip6_left_10d_akima_cubic[[1]]$Data
max(dmat1)
max(dmat1[dmat1!=max(dmat1)] ) # only entry of this coordinate is a problem
fail1 <- which(dmat1 == max(dmat1),arr.ind = TRUE)

dmat1[fail1+c(0,1,0)]#lat up
dmat1[fail1+c(0,1,-1)]#lat up west
dmat1[fail1+c(0,1,1)]# lat up east
dmat1[fail1+c(0,-1,-1)]# lat down west; does not exist
dmat1[fail1+c(0,-1,1)]# lat down east; does not exist
dmat1[fail1+c(0,-1,0)]# lat down; does not exist
dmat1[fail1+c(0,0,-1)]# lon west
dmat1[fail1+c(0,0,1)]# lon east

dmat1[fail1] <- (dmat1[fail1+c(0,1,0)] + 
                 dmat1[fail1+c(0,1,-1)] + #lat up
                 dmat1[fail1+c(0,1,1)] + #lat up west
                 dmat1[fail1+c(0,0,-1)] + # lon west
                 dmat1[fail1+c(0,0,1)]# lon east
)/5

# Fix member 5 minimum
dmat5 <- tas_ssp585_cmip6_left_10d_akima_cubic[[5]]$Data
min(dmat5)
min(dmat5[dmat5!=min(dmat5)] ) # only entry of this coordinate is a problem
fail5 <- which(dmat5 == min(dmat5),arr.ind = TRUE)

dmat5[fail5+c(0,1,0)]#lat up
dmat5[fail5+c(0,1,-1)]#lat up west
dmat5[fail5+c(0,1,1)]# lat up east
dmat5[fail5+c(0,-1,-1)]# lat down west; does not exist
dmat5[fail5+c(0,-1,1)]# lat down east; does not exist
dmat5[fail5+c(0,-1,0)]# lat down; does not exist
dmat5[fail5+c(0,0,-1)]# lon west
dmat5[fail5+c(0,0,1)]# lon east

dmat5[fail5] <- (dmat5[fail5+c(0,1,0)] + 
                   dmat5[fail5+c(0,1,-1)] + #lat up
                   dmat5[fail5+c(0,1,1)] + #lat up west
                   dmat5[fail5+c(0,0,-1)] + # lon west
                   dmat5[fail5+c(0,0,1)]# lon east
)/5

# Fix member 6 minimum
dmat6 <- tas_ssp585_cmip6_left_10d_akima_cubic[[6]]$Data
min(dmat6)
min(dmat6[dmat6!=min(dmat6)] ) # only entry of this coordinate is a problem
fail6 <- which(dmat6 == min(dmat6),arr.ind = TRUE)

dmat6[fail6+c(0,1,0)]#lat up
dmat6[fail6+c(0,1,-1)]#lat up west
dmat6[fail6+c(0,1,1)]# lat up east
dmat6[fail6+c(0,-1,-1)]# lat down west
dmat6[fail6+c(0,-1,1)]# lat down east
dmat6[fail6+c(0,-1,0)]# lat down
dmat6[fail6+c(0,0,-1)]# lon west
dmat6[fail6+c(0,0,1)]# lon east

dmat6[fail6] <- (dmat6[fail6+c(0,1,0)] + 
                   dmat6[fail6+c(0,1,-1)] + #lat up
                   dmat6[fail6+c(0,1,1)] + #lat up west
                   dmat6[fail6+c(0,-1,-1)] + # lat down west
                    dmat6[fail6+c(0,-1,1)] + # lat down east
                    dmat6[fail6+c(0,-1,0)] + # lat down
                   dmat6[fail6+c(0,0,-1)] + # lon west
                   dmat6[fail6+c(0,0,1)]# lon east
)/8


tas_ssp585_cmip6_left_10d_akima_cubic_corrected <- tas_ssp585_cmip6_left_10d_akima_cubic
tas_ssp585_cmip6_left_10d_akima_cubic_corrected[[1]]$Data <- dmat1
tas_ssp585_cmip6_left_10d_akima_cubic_corrected[[5]]$Data <- dmat5
tas_ssp585_cmip6_left_10d_akima_cubic_corrected[[6]]$Data <- dmat6

# check again
# transform to Celcius
tas_ssp585_cmip6_left_10d_akima_cubic_corrected_degC <- lapply(tas_ssp585_cmip6_left_10d_akima_cubic_corrected, udConvertGrid, new.units = "degC")
tas_ssp585_cmip6_left_10d_akima_cubic_corrected.ens <- bindGrid(tas_ssp585_cmip6_left_10d_akima_cubic_corrected_degC, dimension = "member")

require(RColorBrewer)
spatialPlot(climatology(tas_ssp585_cmip6_left_10d_akima_cubic_corrected.ens, by.member = TRUE), 
            backdrop.theme = "coastline", as.table = TRUE,
            lonCenter = 180,rev.colors = TRUE,
            at = seq(-60,80,10))



save(tas_ssp585_cmip6_left_10d_akima_cubic_corrected, file = "data/tas_ssp585_cmip6_left_10d_akima_cubic_corrected.rda")


