library(devtools)
# install_github("SantanderMetGroup/loadeR@devel")
# install_github("SantanderMetGroup/climate4R.UDG@mai-devel")

library(loadeR)
library(transformeR)
library(visualizeR)

source("/oceano/gmeteo/WORK/lisette/Trabajo/creds")
####################################################################################
# JRA55
####################################################################################
AvailableJRA55 <- UDG.datasets("JRA55",full.info = TRUE)
a <-dataInventory("http://spock.meteo.unican.es/tds5/dodsC/jra55/daily/JRA55_daily_Dataset.ncml")
JRA55 <- loadGridData("http://spock.meteo.unican.es/tds5/dodsC/jra55/daily/JRA55_daily_Dataset.ncml", "tas", years = 1981:2010, aggr.m = "mean")
save("JRA55", file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/JRA55.rda")

###############################################################################
# CMIP5 models 
###############################################################################
# TWEEDE KEER
# download.file("https://raw.githubusercontent.com/SantanderMetGroup/ATLAS/master/AtlasHub-inventory/Hub/CMIP5_day_Hub.csv", 
               # destfile = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP5_day_Hub.csv", method = "curl")
cmip5inv <- read.csv("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP5_day_Hub.csv")
# EERSTE KEER
# download.file("https://raw.githubusercontent.com/SantanderMetGroup/ATLAS/master/AtlasHub-inventory/Hub/CMIP5_Hub_20191212.csv", 
              # destfile = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP5_Hub_20191212.csv", method = "curl")
#cmip5inv <- read.csv("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP5_Hub_20191212.csv")

# requisitos
subtas <- which(cmip5inv$tas == TRUE,arr.ind = TRUE)
subr1 <- grep("r1i1",cmip5inv$X)
subr2 <- grep("r2i1",cmip5inv$X)
subhist <- grep("historical",cmip5inv$X)
submodels <- grep("Had|CNRM|Can|EARTH|GFDL|IPSL|MIROC|MPI|Nor",cmip5inv$X)
sub85 <- grep("rcp85",cmip5inv$X)

allmodels <- 1:length(cmip5inv$X)
submodels_left <-  setdiff(allmodels,submodels)
select85_left <- intersect(subr1,intersect(sub85,submodels_left))
selecthist_left <- intersect(subr1,intersect(subhist,submodels_left))
cmip5inv$X[selecthist_left]


UDG.datasets("CMIP5.*rcp85")$CMIP5_subset
# datasets with requistitos
select85 <- intersect(intersect(subtas,subr1),intersect(sub85,submodels))
selecthist <- intersect(intersect(subtas,subr1),intersect(subhist,submodels))
# selecthist <- intersect(intersect(subtas,subr1),subhist) para ver que modeles hay además. 
# Datasets that have successor in rcp 85:
namesa <- gsub(cmip5inv$X[selecthist], pattern = "_historical", replacement = "")
namesb <- gsub(cmip5inv$X[select85], pattern = "_rcp85", replacement = "")
diffel <- paste0(setdiff(namesa,namesb),"_historical") 
dim(diffel) <- c(1,length(diffel))
subdif <- apply(diffel,MARGIN = 2, FUN = function(x) grep(x,cmip5inv$X))
selecthisteq <- setdiff(selecthist,subdif)
cmip5inv$X[selecthisteq]
cmip5inv$X[select85]
# Manually add r12i1p1 and r2i1p1 for earth 
subr2 <- grep("r2i1",cmip5inv$X)
subr12 <- grep("r12i1",cmip5inv$X)
subearth <- grep("EARTH",cmip5inv$X)
earthr2h <- intersect(intersect(subtas,subhist),intersect(subr2,subearth))
earthr2rcp85 <- intersect(intersect(subtas,sub85),intersect(subr2,subearth))
earthr12h <- intersect(intersect(subtas,subhist),intersect(subr12,subearth))
earthr12rcp85 <- intersect(intersect(subtas,sub85),intersect(subr12,subearth))
select85plus <- c(select85,earthr2rcp85,earthr12rcp85)
selecthistplus <- c(selecthisteq,earthr2h,earthr12h)


# Sets used in first CMIP5 analisis:
prueba <- UDG.datasets("CMIP5.*historical")
prueba.85 <- UDG.datasets("CMIP5.*rcp85")
prueba$CMIP5_subset
prueba.85$CMIP5_subset
# comprobobar que all elements in first elements are in amplified datasets:
is.element(prueba$CMIP5_subset,cmip5inv$X[selecthistplus])
is.element(prueba.85$CMIP5_subset,cmip5inv$X[select85plus])
# extradatasets for intercomparison cmip5:
datasets_cmip5 <-cmip5inv$X[selecthistplus]
datasets.85_cmip5 <-cmip5inv$X[select85plus]
datasets_cmip5_extra <- setdiff(cmip5inv$X[selecthistplus],prueba$CMIP5_subset)
datasets.85_cmip5_extra <- setdiff(cmip5inv$X[select85plus],prueba.85$CMIP5_subset)
datasets_cmip5_earth <- datasets_cmip5[grep("EC-EARTH",datasets_cmip5)]
datasets.85_cmip5_earth <- datasets.85_cmip5[grep("EC-EARTH",datasets_cmip5)]

# Check enddates historical models:
# cmip5_histInv <- lapply(datasets_cmip5,function(x) dataInventory(x))
# names(cmip5_histInv) <- datasets_cmip5
# timeranges <- sapply(cmip5_histInv, function(x) x$tas$Dimensions$time$Date_range)
timeranges

# Both CMIP5_HadGEM2-CC_r1i1p1_historical 
# and CMIP5_HadGEM2-ES_r1i1p1_historical until 2005-11-30
# CMIP5_EC-EARTH_r1i1p1_historical until 2009-12-31 
# CMIP5_EC-EARTH_r2i1p1_historical until 2012-12-31
# CMIP5_EC-EARTH_r12i1p1_historical until 2012-12-31
# Effectivamente en los nmcl's se ve que van hasta estas fechas:
# UDG.datasets("CMIP5",full.info = TRUE)$CMIP5_AR5_1run


historical_cmip5_extra <- lapply(1:length(datasets_cmip5_extra), function(x) {
  pre <- loadGridData(datasets_cmip5_extra[x], "tas", years = 1981:2004, aggr.m = "mean")
  if(datasets_cmip5_extra[x] == "CMIP5_HadGEM2-CC_r1i1p1_historical") {
    pre21 <- loadGridData(datasets_cmip5_extra[x], "tas", years = 2005, aggr.m = "mean", season = 1:11)
    pre22 <- loadGridData(datasets.85_cmip5_extra[x], "tas", years = 2005, aggr.m = "mean", season = 12)
    pre2 <- bindGrid(pre21, pre22, dimension = "time")
  } else {
    pre2 <- loadGridData(datasets_cmip5_extra[x], "tas", years = 2005, aggr.m = "mean")
  }
  hist1 <- bindGrid(pre, pre2, dimension = "time")
  hist2 <- loadGridData(datasets.85_cmip5_extra[x], "tas", years = 2006:2010, aggr.m = "mean")
  bindGrid(hist1, hist2, dimension = "time")
})

setnames <- gsub(gsub(datasets_cmip5_extra, pattern = "_historical", replacement = ""),
                 pattern = "-", replacement = ".")
names(historical_cmip5_extra) <- setnames
save(historical_cmip5_extra, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_extra.rda")
###############################################################################
# CMIP5 models met nieuwe instituties.
###############################################################################
# ALLEEN "TWEEDE KEER"
# download.file("https://raw.githubusercontent.com/SantanderMetGroup/ATLAS/master/AtlasHub-inventory/Hub/CMIP5_day_Hub.csv", 
# destfile = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP5_day_Hub.csv", method = "curl")
cmip5inv <- read.csv("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP5_day_Hub.csv")

# requisitos
subtas <- which(cmip5inv$tas == TRUE,arr.ind = TRUE)
subr1 <- grep("r1i1",cmip5inv$X)
subr2 <- grep("r2i1",cmip5inv$X)
subhist <- grep("historical",cmip5inv$X)
sub85 <- grep("rcp85",cmip5inv$X)
submodels <- grep("Had|CNRM|Can|EARTH|GFDL|IPSL|MIROC|MPI|Nor",cmip5inv$X)
allmodels <- 1:length(cmip5inv$X)
submodels_left <-  setdiff(allmodels,submodels)

# allmodels <- 1:length(Availablecmip5)
# leftoutsubmodels <-  setdiff(allmodels,submodels)
# lefoutselect85 <- intersect(subr1,intersect(sub85,leftoutsubmodels))
# lefoutselecthist <- intersect(subr1,intersect(subhist,leftoutsubmodels))
# Availablecmip5[lefoutselect85]
# Availablecmip5[lefoutselecthist]

# datasets with requistitos
select85_left <- intersect(subr1,intersect(sub85,submodels_left))
selecthist_left <- intersect(subr1,intersect(subhist,submodels_left))
cmip5inv$X[select85_left]
# Datasets that have successor in rcp 85:
namesa_left <- gsub(cmip5inv$X[selecthist_left], pattern = "_historical", replacement = "")
namesb_left <- gsub(cmip5inv$X[select85_left], pattern = "_rcp85", replacement = "")
diffel_left <- paste0(setdiff(namesa_left,namesb_left),"_historical") 
dim(diffel_left) <- c(1,length(diffel_left))
subdif_left <- apply(diffel_left,MARGIN = 2, FUN = function(x) grep(x,cmip5inv$X))
selecthisteq_left <- setdiff(selecthist_left,subdif_left)
cmip5inv$X[selecthisteq_left]
cmip5inv$X[subdif_left]

# extradatasets for intercomparison cmip5 (that have succesor rcp.85):
datasets_cmip5_left <-cmip5inv$X[selecthisteq_left]
datasets.85_cmip5_left <-cmip5inv$X[select85_left]

# Check enddates historical models:
# cmip5_left_histInv <- lapply(datasets_cmip5_left, function(x) dataInventory(x))
# names(cmip5_left_histInv) <- datasets_cmip5_left
# timeranges_left <- sapply(cmip5_left_histInv, function(x) x$tas$Dimensions$time$Date_range)
# timeranges_left
cmip5_left_histInv$`CMIP5_FGOALS-g2_r1i1p1_historical`$tas
# All CMIP5_left historical until 2005-12-31 (most common) 
# Except 
# CMIP5_bcc-csm1-1-m_r1i1p1_historical until 2012-12-31
# CMIP5_bcc-csm1-1_r1i1p1_historical until 2012-12-30
# AS CMIP5_CESM1-FASTCHEM_r1i1p1_historical CMIP5_MRI-ESM1_r1i1p1_historical don t have succesor rcp85 and go until 2005_12_31. We do left out these two. 
# Effectivamente en los nmcl's se ve que van hasta estas fechas:
# UDG.datasets("CMIP5",full.info = TRUE)$CMIP5_AR5_1run

for (i in 1:length(datasets_cmip5_left)){
  hist1 <- loadGridData(datasets_cmip5_left[i], "tas", years = 1981:2005, aggr.m = "mean")
  hist2 <- loadGridData(datasets.85_cmip5_left[i], "tas", years = 2006:2010, aggr.m = "mean")
  assign(gsub(gsub(datasets_cmip5_left, pattern = "_historical", replacement = ""),pattern = "-", replacement = ".")[i],
         bindGrid(hist1, hist2, dimension = "time"))
  save(list = gsub(gsub(datasets_cmip5_left, pattern = "_historical", replacement = ""),pattern = "-", replacement = ".")[i], file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left/",gsub(gsub(datasets_cmip5_left, pattern = "_historical", replacement = ""),pattern = "-", replacement = ".")[i],".rda"))
  hist1 <- hist2 <- NULL
}

datasets_cmip5_left[12]

# historical_cmip5_left <- lapply(1:length(datasets_cmip5_left), function(x) {
#   hist1 <- loadGridData(datasets_cmip5_left[x], "tas", years = 1981:2005, aggr.m = "mean")
#   hist2 <- loadGridData(datasets.85_cmip5_left[x], "tas", years = 2006:2010, aggr.m = "mean")
#   bindGrid(hist1, hist2, dimension = "time")
# })
# 
# setnames <- gsub(gsub(datasets_cmip5_left, pattern = "_historical", replacement = ""),
#                  pattern = "-", replacement = ".")
# names(historical_cmip5_left) <- setnames
# save(historical_cmip5_left, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left.rda")

##########################################################################################
# earth without aggregation of rcp.85
##########################################################################################
datasets_cmip5_earth
# Only those
historical_cmip5_earth_r1r2r12 <- lapply(1:length(datasets_cmip5_earth), function(x) {
  hist1 <- loadGridData(datasets_cmip5_earth[x], "tas", years = 1981:2009, aggr.m = "mean")
  if(datasets_cmip5_extra[x] == "CMIP5_EC-EARTH_r1i1p1_historical") {
    hist2 <- loadGridData(datasets.85_cmip5_earth[x], "tas", years = 2010, aggr.m = "mean")
  } else {
    hist2 <- loadGridData(datasets_cmip5_earth[x], "tas", years = 2010, aggr.m = "mean")
  }
  bindGrid(hist1, hist2, dimension = "time")
})

names(historical_cmip5_earth_r1r2r12) <- gsub(datasets_cmip5_earth, pattern = "-", replacement = ".")
save(historical_cmip5_earth_r1r2r12, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_earth_r1r2r12.rda")


# historical_cmip5_earth_r1r2r12<- lapply(2:length(datasets_cmip5_earth), function(x) {
  # loadGridData(datasets_cmip5_earth[x], "tas", years = 1981:2010, aggr.m = "mean")
# })

################################################################################
# rcp85 2071-2100
################################################################################
#lapply(1:5, function(x)dataInventory(datasets.85_cmip5[x]))

Availablecmip5  <- UDG.datasets("CMIP5")$CMIP5_AR5
URLs.cmip5 <- (UDG.datasets("CMIP5",full.info = TRUE)$CMIP5_AR5)[,"url"]

subr1 <- grep("r1i1",Availablecmip5)
# subhist <- grep("historical",Availablecmip5)
submodels <- grep("Had|CNRM|Can|EARTH|GFDL|IPSL|MIROC|MPI|Nor",Availablecmip5)
sub85 <- grep("rcp85",Availablecmip5)
select85 <- intersect(subr1,intersect(sub85,submodels))

# Manually add r12i1p1 and r2i1p1 for earth 
subr2 <- grep("r2i1",Availablecmip5)
subr12 <- grep("r12i1",Availablecmip5)
subearth <- grep("EARTH",Availablecmip5)
earthr2rcp85 <- intersect(sub85,intersect(subr2,subearth))
earthr12rcp85 <- intersect(sub85,intersect(subr12,subearth))
select85plus <- c(select85,earthr2rcp85,earthr12rcp85)

URLs.cmip5.sub85 <- URLs.cmip5[select85plus]
datasets.85_cmip5 <- Availablecmip5[select85plus]
Availablecmip5[select85plus]

URLs.cmip5.sub85
str(URLs.cmip5.sub85)
datasets.85_cmip5

# rcp85_cmip5 <- lapply(1:length(datasets.85_cmip5), function(x) {
  # loadGridData(URLs.cmip5.sub85[x], "tas", years = 2071:2100, aggr.m = "mean")
# })
# names(rcp85_cmip5) <- gsub(datasets.85_cmip5, pattern = "-", replacement = ".")
# save(rcp85_cmip5, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5.rda")

for (i in 1:length(datasets.85_cmip5)){
  assign(gsub(datasets.85_cmip5, pattern = "-", replacement = ".")[i],loadGridData(URLs.cmip5.sub85[i], "tas", years = 2071:2100, aggr.m = "mean"))
  save(list = gsub(datasets.85_cmip5, pattern = "-", replacement = ".")[i], file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_2071_2100/",gsub(datasets.85_cmip5, pattern = "-", replacement = ".")[i],".rda"))
}

########################################
# Overige instituties rcp85 2071-2100
########################################
allmodels <- 1:length(Availablecmip5)
leftoutsubmodels <-  setdiff(allmodels,submodels)
lefoutselect85 <- intersect(subr1,intersect(sub85,leftoutsubmodels))
datasets.85_cmip5_left <- Availablecmip5[lefoutselect85]
URLs.cmip5.sub85_left <- URLs.cmip5[lefoutselect85]

# Check dates rcp85 models:
cmip5_left_85Inv <- lapply(datasets.85_cmip5_left, function(x) dataInventory(x))
names(cmip5_left_85Inv) <- datasets.85_cmip5_left
timeranges_left_85 <- sapply(cmip5_left_85Inv, function(x) x$tas$Dimensions$time$Date_range)
timeranges_left_85

cmip5_left_85Inv$`CMIP5_FGOALS-g2_r1i1p1_rcp85`
# All CMIP5_left rcp until 2100-12-31 (most common) or 2300-12-30 or 2300-12-31
# Except 
# CMIP5_bcc-csm1-1-m_r1i1p1_rcp85 until 2099-12-31
# CMIP5_FGOALS-g2_r1i1p1_rcp85 until NO AVAILABLE
# Hence we include exception for bcc: 2070-2099 instead of 2071 - 2100


# lefoutselecthist <- intersect(subr1,intersect(subhist,leftoutsubmodels))
# Availablecmip5[lefoutselecthist] --> THIS ONE IN HISTORICAL TREATHED

 # rcp85_cmip5_left <- lapply(1:length(datasets.85_cmip5_left), function(x) {
 # loadGridData(URLs.cmip5.sub85_left[x], "tas", years = 2071:2100, aggr.m = "mean")
 # })
 # names(rcp85_cmip5_left) <- gsub(datasets.85_cmip5_left, pattern = "-", replacement = ".")
 # save(rcp85_cmip5_left, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_left_2071_2100/rcp85_cmip5_left.rda")

for (i in 3:length(datasets.85_cmip5_left)){
 # for (i in 1:length(datasets.85_cmip5_left)){
   if(datasets.85_cmip5_left[i] == "CMIP5_bcc-csm1-1-m_r1i1p1_rcp85"){
     assign(gsub(datasets.85_cmip5_left, pattern = "-", replacement = ".")[i],loadGridData(URLs.cmip5.sub85_left[i], "tas", years = 2070:2099, aggr.m = "mean"))
     save(list = gsub(datasets.85_cmip5_left, pattern = "-", replacement = ".")[i], file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_left_2071_2100/",gsub(datasets.85_cmip5_left, pattern = "-", replacement = ".")[i],".rda"))
   } else{
     assign(gsub(datasets.85_cmip5_left, pattern = "-", replacement = ".")[i],loadGridData(URLs.cmip5.sub85_left[i], "tas", years = 2071:2100, aggr.m = "mean"))
     save(list = gsub(datasets.85_cmip5_left, pattern = "-", replacement = ".")[i], file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_left_2071_2100/",gsub(datasets.85_cmip5_left, pattern = "-", replacement = ".")[i],".rda"))
   }
}

datasets_cmip5_left
###############################################################################
# CMIP6 models 
###############################################################################

# download.file("https://raw.githubusercontent.com/SantanderMetGroup/ATLAS/master/AtlasHub-inventory/Hub/CMIP6Amon_Hub_20191028.csv", 
              # destfile = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP6Amon_Hub_20191028.csv", method = "curl")
cmip6inv <- read.csv("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP6Amon_Hub_20191028.csv")

# requisitos
c6subtas <- which(cmip6inv$tas == TRUE,arr.ind = TRUE)
c6subr1 <- grep("r1i1",cmip6inv$X)
c6subr1 <- grep("r1",cmip6inv$X)
c6subr2 <- grep("r2i1",cmip6inv$X)
c6subhist <- grep("historical",cmip6inv$X)
c6submodels <- grep("Had|CNRM|Can|Earth|GFDL|IPSL|MIROC|MPI|Nor",cmip6inv$X)
c6sub85 <- grep("85",cmip6inv$X)

# datasets with requistitos
c6select85 <- intersect(intersect(c6subtas,c6subr1),intersect(c6sub85,c6submodels))
c6select85 <- intersect(c6subtas,intersect(c6sub85,c6submodels))
datasets.85_cmip6 <- cmip6inv$X[c6select85]

c6selecthist <- intersect(intersect(c6subtas,c6subr1),intersect(c6subhist,c6submodels)) # CHECK IF REALMENTE R1 has all 
datasets_cmip6 <- cmip6inv$X[c6selecthist]

cmip6inv$X[intersect(intersect(c6submodels,c6subhist),c6subtas)]

# # Manually add r12i1p1 and r2i1p1 for earth CHECK 
# c6subr2 <- grep("r2i1",cmip6inv$X)
# c6subr12 <- grep("r12i1",cmip6inv$X)
# c6subearth <- grep("Earth",cmip6inv$X)
# c6earthr2h <- intersect(intersect(c6subtas,c6subhist),intersect(c6subr2,c6subearth))
# c6earthr12h <- intersect(intersect(c6subtas,c6subhist),intersect(c6subr12,c6subearth))
# c6selecthistplus <- c(c6selecthist,c6earthr2h,c6earthr12h)
# cmip6inv$X[c6selecthistplus]

#########################################################################
# jaren tot 2010 zitten in historical. 
#########################################################################
historical_cmip6 <- lapply(1:length(datasets_cmip6), function(x) {
  loadGridData(datasets_cmip6[x], "tas", years = 1981:2010, aggr.m = "mean")
})

names(historical_cmip6) <- gsub(datasets_cmip6, pattern = "-", replacement = ".")
save(historical_cmip6, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip6.rda")
##########################################################################
# Jaren met rcp85 = ssp585 zitten in new. 
##########################################################################
Availablecmip6  <- UDG.datasets("CMIP6")$CMIP6Amon
download.file("https://github.com/SantanderMetGroup/ATLAS/blob/master/ATLAS-inventory/Hub/CMIP6_mon_Hub.csv", 
destfile = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP6_mon_Hub.csv", method = "curl")
# download.file("https://raw.githubusercontent.com/SantanderMetGroup/ATLAS/master/AtlasHub-inventory/Hub/CMIP6_mon_Hub.csv", 
# destfile = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP6_mon_Hub.csv", method = "curl")
cmip6inv2 <- read.csv("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/CMIP6_mon_Hub.csv")

# requisitos
c6subtas <- which(cmip6inv2$tas == TRUE,arr.ind = TRUE)
cmip6inv2[c6subtas,]
c6subr1 <- grep("r1i1",cmip6inv2$X)
c6subhist <- grep("historical",cmip6inv2$X)
c6submodels <- grep("Had|CNRM|Can|Earth|GFDL|IPSL|MIROC|MPI|Nor",cmip6inv2$X)
c6sub85 <- grep("85",cmip6inv2$X)


cmip6inv2[intersect(c6subtas,c6sub85),]
# datasets with requistitos
c6select85 <- intersect(intersect(c6subtas,c6subr1),intersect(c6sub85,c6submodels))
# c6select85 <- intersect(intersect(c6subtas,c6subr1),c6sub85)

wanted <- cmip6inv2$X[c6select85]
todownload <- is.element(as.character(wanted),Availablecmip6)
datasets.85_cmip6 <- wanted[todownload]
# dataInventory(datasets.85_cmip6[1])

ssp585_cmip6 <- lapply(1:length(datasets.85_cmip6), function(x) {
  loadGridData(as.character(datasets.85_cmip6)[x], "tas", years = 2071:2100, aggr.m = "mean")
})

names(ssp585_cmip6) <- gsub(datasets.85_cmip6, pattern = "-", replacement = ".")
save(ssp585_cmip6, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/ssp585_cmip6.rda")
#################################################################
# Upscaling JRA55
#################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/JRA55.rda")

tas_JRA55_10d_akima_cubic <- interpGrid(JRA55,new.coordinates = getGrid(tas_interim_10dnew),
                                        method = 'bilinear',
                                        bilin.method = 'akima',
                                        linear = FALSE,
                                        extrap = TRUE)
save(tas_JRA55_10d_akima_cubic, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_JRA55_10d_akima_cubic.rda")                                    

#################################################################
# Upscaling CMIP5 extra
#################################################################

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_extra.rda")

# Lo interpolo todo a la malla del objeto tas_interim_10dnew (utilizado en el paper con interim).
# las coordenadas de los datasets de cmip5 en general no cubren todo el mundo (por ejemplo uno va de  -87, 87, y -177, 177 etc)
# solo el método cubic (method = bilinear, linear = FALSE) de akima permite extrapolación a las longitudes/latitudes externas a la
# región original del GCM. 
# Anteriormente apliqué el método lineal de 'fields' (lo cual, por cierto, no coincide con el metodo de extrapolación lineal de 'akima')

tas_historical_cmip5_extra_10d_akima_cubic <- lapply(historical_cmip5_extra, function(x){interpGrid(x,
                                                                            new.coordinates = getGrid(tas_interim_10dnew),
                                                                            method = 'bilinear',
                                                                            bilin.method = 'akima',
                                                                            linear = FALSE,
                                                                            extrap = TRUE
)
})

names(tas_historical_cmip5_extra_10d_akima_cubic) <- names(historical_cmip5_extra)
save(tas_historical_cmip5_extra_10d_akima_cubic, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_historical_cmip5_extra_10d_akima_cubic.rda")

#################################################################
# Upscaling CMIP5 earth
#################################################################

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_earth_r1r2r12.rda")

# Lo interpolo todo a la malla del objeto tas_interim_10dnew (utilizado en el paper con interim).
# las coordenadas de los datasets de cmip5 en general no cubren todo el mundo (por ejemplo uno va de  -87, 87, y -177, 177 etc)
# solo el método cubic (method = bilinear, linear = FALSE) de akima permite extrapolación a las longitudes/latitudes externas a la
# región original del GCM. 
# Anteriormente apliqué el método lineal de 'fields' (lo cual, por cierto, no coincide con el metodo de extrapolación lineal de 'akima')

tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic <- lapply(historical_cmip5_earth_r1r2r12, function(x){interpGrid(x,
                                                                                                    new.coordinates = getGrid(tas_interim_10dnew),
                                                                                                    method = 'bilinear',
                                                                                                    bilin.method = 'akima',
                                                                                                    linear = FALSE,
                                                                                                    extrap = TRUE
)
})

names(tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic) <- names(historical_cmip5_earth_r1r2r12)
save(tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda")

################################################################
# Upscaling CMIP5 left 
################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
dataset_list <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left", full.names = T)
dataset_names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left")
dataset_names <- gsub(".rda", "", dataset_names)

historical_cmip5_left <- lapply(dataset_list, function(x){get(load(x))})
names(historical_cmip5_left) <- dataset_names

tas_historical_cmip5_left_10d_akima_cubic <- lapply(historical_cmip5_left, function(x){interpGrid(x,
                                                                                                  new.coordinates = getGrid(tas_interim_10dnew),
                                                                                                  method = 'bilinear',
                                                                                                  bilin.method = 'akima',
                                                                                                  linear = FALSE,
                                                                                                  extrap = TRUE
)
})

names(tas_historical_cmip5_left_10d_akima_cubic) <- names(historical_cmip5_left)
save(tas_historical_cmip5_left_10d_akima_cubic, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_historical_cmip5_left_10d_akima_cubic.rda")
#################################################################
# Upscaling CMIP6 
#################################################################

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/historical_cmip6.rda")

# Lo interpolo todo a la malla del objeto tas_interim_10dnew (utilizado en el paper con interim).
# las coordenadas de los datasets de cmip5 en general no cubren todo el mundo (por ejemplo uno va de  -87, 87, y -177, 177 etc)
# solo el método cubic (method = bilinear, linear = FALSE) de akima permite extrapolación a las longitudes/latitudes externas a la
# región original del GCM. 
# Anteriormente apliqué el método lineal de 'fields' (lo cual, por cierto, no coincide con el metodo de extrapolación lineal de 'akima')

tas_historical_cmip6_10d_akima_cubic <- lapply(historical_cmip6, function(x){interpGrid(x,
                                                                                        new.coordinates = getGrid(tas_interim_10dnew),
                                                                                        method = 'bilinear',
                                                                                        bilin.method = 'akima',
                                                                                        linear = FALSE,
                                                                                        extrap = TRUE
)
})


save(tas_historical_cmip6_10d_akima_cubic, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_historical_cmip6_10d_akima_cubic.rda")

#################################################################
# Upscaling CMIP5 rcp85 2071-2100
#################################################################

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
rcp85_cmip5_list <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_2071_2100", full.names = T)
rcp85_cmip5_names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_2071_2100")
rcp85_cmip5_names <- gsub(".rda", "", rcp85_cmip5_names)

rcp85_cmip5 <- lapply(rcp85_cmip5_list, function(x){get(load(x))})
names(rcp85_cmip5) <- rcp85_cmip5_names


# Lo interpolo todo a la malla del objeto tas_interim_10dnew (utilizado en el paper con interim).
# las coordenadas de los datasets de cmip5 en general no cubren todo el mundo (por ejemplo uno va de  -87, 87, y -177, 177 etc)
# solo el método cubic (method = bilinear, linear = FALSE) de akima permite extrapolación a las longitudes/latitudes externas a la
# región original del GCM. 
# Anteriormente apliqué el método lineal de 'fields' (lo cual, por cierto, no coincide con el metodo de extrapolación lineal de 'akima')

tas_rcp85_cmip5_10d_akima_cubic <- lapply(rcp85_cmip5, function(x){interpGrid(x,
                                                                                new.coordinates = getGrid(tas_interim_10dnew),
                                                                                method = 'bilinear',
                                                                                bilin.method = 'akima',
                                                                                linear = FALSE,
                                                                                extrap = TRUE
)
})


save(tas_rcp85_cmip5_10d_akima_cubic, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_rcp85_cmip5_10d_akima_cubic.rda")
#################################################################
# Upscaling CMIP5 left rcp85 2071-2100
#################################################################

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
rcp85_cmip5_left_list <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_left_2071_2100", full.names = T)
rcp85_cmip5_left_names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/rcp85_cmip5_left_2071_2100")
rcp85_cmip5_left_names <- gsub(".rda", "", rcp85_cmip5_left_names)

rcp85_cmip5_left <- lapply(rcp85_cmip5_left_list, function(x){get(load(x))})
names(rcp85_cmip5_left) <- rcp85_cmip5_left_names


# Lo interpolo todo a la malla del objeto tas_interim_10dnew (utilizado en el paper con interim).
# las coordenadas de los datasets de cmip5 en general no cubren todo el mundo (por ejemplo uno va de  -87, 87, y -177, 177 etc)
# solo el método cubic (method = bilinear, linear = FALSE) de akima permite extrapolación a las longitudes/latitudes externas a la
# región original del GCM. 
# Anteriormente apliqué el método lineal de 'fields' (lo cual, por cierto, no coincide con el metodo de extrapolación lineal de 'akima')

tas_rcp85_cmip5_left_10d_akima_cubic <- lapply(rcp85_cmip5_left, function(x){interpGrid(x,
                                                                              new.coordinates = getGrid(tas_interim_10dnew),
                                                                              method = 'bilinear',
                                                                              bilin.method = 'akima',
                                                                              linear = FALSE,
                                                                              extrap = TRUE
)
})


save(tas_rcp85_cmip5_left_10d_akima_cubic, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_rcp85_cmip5_left_10d_akima_cubic.rda")



#################################################################
# Upscaling CMIP6 ssp85
#################################################################

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/ssp585_cmip6.rda")

# Lo interpolo todo a la malla del objeto tas_interim_10dnew (utilizado en el paper con interim).
# las coordenadas de los datasets de cmip5 en general no cubren todo el mundo (por ejemplo uno va de  -87, 87, y -177, 177 etc)
# solo el método cubic (method = bilinear, linear = FALSE) de akima permite extrapolación a las longitudes/latitudes externas a la
# región original del GCM. 
# Anteriormente apliqué el método lineal de 'fields' (lo cual, por cierto, no coincide con el metodo de extrapolación lineal de 'akima')

tas_ssp585_cmip6_10d_akima_cubic <- lapply(ssp585_cmip6, function(x){interpGrid(x,
                                                                                        new.coordinates = getGrid(tas_interim_10dnew),
                                                                                        method = 'bilinear',
                                                                                        bilin.method = 'akima',
                                                                                        linear = FALSE,
                                                                                        extrap = TRUE
)
})


save(tas_ssp585_cmip6_10d_akima_cubic, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/data/tas_ssp585_cmip6_10d_akima_cubic.rda")
