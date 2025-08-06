#########################
#
#########################
library(loadeR)
library(transformeR)
library(visualizeR)
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
source("../../../Trabajo/creds")

load("../../R_practice/Data/interim/tas_interim_10dnew.rda")
load("/data/Untitled/Trabajo/R_practice/Data/interim/TP_interim_5d_bil.rda")
load("/data/Untitled/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left/CMIP5_ACCESS1.0_r1i1p1.rda")
load("/data/Untitled/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left/CMIP5_ACCESS1.3_r1i1p1.rda")
load("/data/Untitled/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left/CMIP5_CMCC.CMS_r1i1p1.rda")
load("/data/Untitled/Trabajo/R_practice/exp_GCMs/data/historical_cmip5_left/CMIP5_BNU.ESM_r1i1p1.rda")
load("/data/Untitled/Trabajo/R_practice/Data/interim/tas_2T_Interim_orig.rda")
CMIP5_ACCESS1.0_r1i1p1
historical_cmip5_extra_exsubset <- list(CMIP5_BNU.ESM_r1i1p1,CMIP5_ACCESS1.0_r1i1p1,CMIP5_ACCESS1.3_r1i1p1,CMIP5_CMCC.CMS_r1i1p1,tas_2T_Interim_orig) 
ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","tas_2T_Interim")

tas_historical_cmip5_exsubset_5d_akima_cubic <- lapply(historical_cmip5_extra_exsubset, 
                                                    function(x){interpGrid(x,
                                                                           new.coordinates = getGrid(TP_interim_5d_bil),
                                                                           method = 'bilinear',
                                                                           bilin.method = 'akima',
                                                                          linear = FALSE,
                                                                          extrap = TRUE
)
})

names(tas_historical_cmip5_exsubset_5d_akima_cubic) <- ex.subset 
save(tas_historical_cmip5_exsubset_5d_akima_cubic, file = "data/tas_historical_cmip5_exsubset_5d_akima_cubic.rda")
##########################
# upscale to 20 
##########################
rm(list = ls())
library(transformeR)
library(visualizeR)
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
load("data/tas_historical_cmip5_exsubset_5d_akima_cubic.rda")

tas_historical_cmip5_exsubset_20d_akima_cubic <- lapply(tas_historical_cmip5_exsubset_5d_akima_cubic, 
                                                       function(x){upscaleGrid(x, times = 4, aggr.fun = list(FUN = mean, na.rm = TRUE))})

tas_historical_cmip5_exsubset_20d_akima_cubic <- lapply(tas_historical_cmip5_exsubset_20d_akima_cubic, redim,drop = TRUE)
save(tas_historical_cmip5_exsubset_20d_akima_cubic, file = "data/tas_historical_cmip5_exsubset_20d_akima_cubic.rda")

#########################
# CSIRO 10 
#########################
library(loadeR)
library(transformeR)
library(visualizeR)
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
source("../../../Trabajo/creds")
# load CSIRO 10 iterations
datasets_cmip5_CSIRO <- UDG.datasets("CMIP5.*CSIRO.*historical")$CMIP5_AR5
datasets.85_cmip5_CSIRO <- UDG.datasets("CMIP5.*CSIRO.*rcp85")$CMIP5_AR5

load("../../R_practice/Data/interim/tas_interim_10dnew.rda")

historical_cmip5_CSIRO <- lapply(1:length(datasets_cmip5_CSIRO), function(x) {
  pre <- loadGridData(datasets_cmip5_CSIRO[x], "tas", years = 1981:2004, aggr.m = "mean")
  pre2 <- loadGridData(datasets_cmip5_CSIRO[x], "tas", years = 2005, aggr.m = "mean")
  hist1 <- bindGrid(pre, pre2, dimension = "time")
  hist2 <- loadGridData(datasets.85_cmip5_CSIRO[x], "tas", years = 2006:2010, aggr.m = "mean")
  bindGrid(hist1, hist2, dimension = "time")
})

setnames <- gsub(gsub(datasets_cmip5_CSIRO, pattern = "_historical", replacement = ""),
                 pattern = "-", replacement = ".")
names(historical_cmip5_CSIRO) <- setnames

save(historical_cmip5_CSIRO, file = "data/historical_cmip5_CSIRO.rda")

load("data/historical_cmip5_CSIRO.rda")
tas_historical_cmip5_CSIRO_10d_akima_cubic <- lapply(historical_cmip5_CSIRO, 
                                                    function(x){interpGrid(x,
                                                                           new.coordinates = getGrid(tas_interim_10dnew),
                                                                           method = 'bilinear',
                                                                           bilin.method = 'akima',
                                                                           linear = FALSE,
                                                                           extrap = TRUE
                                                    )
                                                    })

names(tas_historical_cmip5_CSIRO_10d_akima_cubic) <- names(datasets_cmip5_CSIRO)
save(tas_historical_cmip5_CSIRO_10d_akima_cubic, file = "data/tas_historical_cmip5_CSIRO_10d_akima_cubic.rda")
