library(loadeR)
library(transformeR)
library(visualizeR)

source("/oceano/gmeteo/WORK/lisette/Trabajo/creds")

# (With old version of UDG.datasets)
# Now the below code equivalent to 
# prueba <- UDG.datasets("CMIP5.*historical")
# prueba$CMIP5_subset
UDG.datasets("CMIP5")
datasets <- UDG.datasets("CMIP5.*historical")$name[-6]
datasets.85 <- UDG.datasets("CMIP5.*rcp85")$name[-6]

# Los datasets historicos del CMIP5 llegan hasta diciembre 2005, 
# salvo "CMIP5_HadGEM2-ES_r1i1p1_historical" que llega hasta noviembre 2005 (diciembre lo cojo del rcp85).
# El periodo considerado con Interim es 1981:2010.
# Por lo tanto complemento la serie histórica de cmip5 con los años 2006:2010 del rcp85

historical <- lapply(1:length(datasets), function(x) {
  pre <- loadGridData(datasets[x], "tas", years = 1981:2004, aggr.m = "mean")
  if(datasets[x] == "CMIP5_HadGEM2-ES_r1i1p1_historical") {
    pre21 <- loadGridData(datasets[x], "tas", years = 2005, aggr.m = "mean", season = 1:11)
    pre22 <- loadGridData(datasets.85[x], "tas", years = 2005, aggr.m = "mean", season = 12)
    pre2 <- bindGrid(pre21, pre22, dimension = "time")
  } else {
    pre2 <- loadGridData(datasets[x], "tas", years = 2005, aggr.m = "mean")
  }
  hist1 <- bindGrid(pre, pre2, dimension = "time")
  hist2 <- loadGridData(datasets.85[x], "tas", years = 2006:2010, aggr.m = "mean")
  bindGrid(hist1, hist2, dimension = "time")
})

save(historical, file = "data/historical.rda")

### upscaling------------

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load("data/historical.rda")

# Lo interpolo todo a la malla del objeto tas_interim_10dnew (utilizado en el paper con interim).
# las coordenadas de los datasets de cmip5 en general no cubren todo el mundo (por ejemplo uno va de  -87, 87, y -177, 177 etc)
# solo el método cubic (method = bilinear, linear = FALSE) de akima permite extrapolación a las longitudes/latitudes externas a la
# región original del GCM. 
# Anteriormente apliqué el método lineal de 'fields' (lo cual, por cierto, no coincide con el metodo de extrapolación lineal de 'akima')

tas_historical_10d_akima_cubic <- lapply(historical, function(x){interpGrid(x,
                                          new.coordinates = getGrid(tas_interim_10dnew),
                                          method = 'bilinear',
                                          bilin.method = 'akima',
                                          linear = FALSE,
                                          extrap = TRUE
                                          )
  })


save(tas_historical_10d_akima_cubic, file = "data/tas_historical_10d_akima_cubic.rda")

######################################################################################################
# Interim data cubic akima
######################################################################################################
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_2T_Interim_orig.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")
load(file = "/data/Untitled/Trabajo/R_practice/Data/interim/tas_2T_Interim_orig.rda")
load("/data/Untitled/Trabajo/R_practice/Data/interim/tas_interim_10dnew.rda")


tas_interim_10d_akima_cubic <- interpGrid(tas_2T_Interim_orig, 
                                  new.coordinates = getGrid(tas_interim_10dnew),
                                  method = 'bilinear',
                                  bilin.method = 'akima',
                                  linear = FALSE,
                                  extrap = TRUE
)


save(tas_interim_10d_akima_cubic, file = "data/tas_interim_10d_akima_cubic.rda")

ens <- bindGrid(list(tas_interim_10dnew,tas_interim_10d_akima_cubic), dimension = "member")
spatialPlot(climatology(ens),backdrop.theme = "coastline",lonCenter = 180)

