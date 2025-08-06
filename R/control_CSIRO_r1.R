narcs(CSIRO.Mk3.6.0_r1i1p1_3_1700_1800i)
narcs(CMIP5_CSIRO.Mk3.6.0_r1i1p1_3_1700_1800i)

tas_historical_cmip5_left_10d_akima_cubic$CMIP5_CSIRO.Mk3.6.0_r1i1p1
cubicsetsCSIRO <- tas_historical_cmip5_CSIRO_10d_akima_cubic
namescubicsCSIRO <- gsub("-",".",gsub(".ncml","",gsub("https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip5/historical/historical_day_","",sapply(cubicsetsCSIRO,function(x)attr(x,"dataset")))))
names(cubicsetsCSIRO) <- namescubicsCSIRO
all.equal(tas_historical_cmip5_left_10d_akima_cubic$CMIP5_CSIRO.Mk3.6.0_r1i1p1,cubicsetsCSIRO$CSIRO.Mk3.6.0_r1i1p1)

tas_historical_cmip5_left_10d_akima_cubic$CMIP5_CSIRO.Mk3.6.0_r1i1p1$Data
cubicsetsCSIRO$CSIRO.Mk3.6.0_r1i1p1$Data
all.equal(historical_cmip5_CSIRO$CMIP5_CSIRO.Mk3.6.0_r1i1p1$Data,CMIP5_CSIRO.Mk3.6.0_r1i1p1$Data)

historical_cmip5_CSIRO$CMIP5_CSIRO.Mk3.6.0_r1i1p1

CMIP5_CSIRO.Mk3.6.0_r1i1p1
