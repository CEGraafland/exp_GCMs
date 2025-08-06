#################################################################################
# Principal Components Analisis. 
#################################################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"
cubicsets <- c(tas_historical_10d_akima_cubic_corrected,listinterim)

namescubics <- names(cubicsets)
shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
                   pattern = "_r12i1p1", replacement = "")
shortnames

data_gcms_out <- lapply(cubicsets,function(x) TimeCoordsAnom_from_Grid_rms(x, rms = TRUE))
data_gcms <- lapply(data_gcms_out, as.data.frame)

##############################################################################
#
prinCompInt <- prinComp(cubicsets[["interim_10d_akima_cubic"]],v.exp = 0.95)
plotClimatology()
plotEOF(prinCompInt,backdrop.theme = "coastline")
data_Int_anom_out <- mat2Dto3Darray(data_gcms_out[[10]], attr(data_gcms_out[[10]],"Xcoords"), attr(data_gcms_out[[10]],"Ycoords"))
grid_Int_anom <- cubicsets[["interim_10d_akima_cubic"]]
grid_Int_anom$Data <- data_Int_anom_out
prinCompInt <- prinComp(grid_Int_anom,n.eofs = 10)
prinCompInt$`2T`[[1]]$EOFs
plotEOF
transformeR::EOF2clim(prinCompInt, )


EOF2clim <- function(prinCompObj, ind.var, member, n.eofs) {
  varNames <- attributes(prinCompObj)$names
  levs <- attributes(prinCompObj)$level
  x <- attributes(prinCompObj)$xCoords
  y <- attributes(prinCompObj)$yCoords
  start <- attr(prinCompObj, "dates_start") %>% head(1)
  end <- attr(prinCompObj, "dates_end") %>% tail(1)
  season <- attr(prinCompObj, "season")
  mu <- attributes(prinCompObj[[ind.var]][[member]][["orig"]])$"scaled:center"
  sigma <- attributes(prinCompObj[[ind.var]][[member]][["orig"]])$"scaled:scale"
  eofs <- prinCompObj[[ind.var]][[member]]$EOFs[, 1:n.eofs, drop = FALSE]
  prinCompObj <- NULL
  # Rescale EOFs
  aux <- (eofs * sigma + mu) %>% t() %>% mat2Dto3Darray(x, y) %>% list()
  # Recover grids structure (EOFS are treated as members, time = 1 like a climatology
  Data <- do.call("abind", c(aux, along = -1)) %>% unname()
  attr(Data, "dimensions") <- c("time", "member", "lat", "lon")
  attr(Data, "climatology:fun") <- "transformeR::prinComp"
  xyCoords = list("x" = x, "y" = y)
  Dates <- list("start" = start, "end" = end)
  attr(Dates, "season") <- season
  start <- end <- NULL
  Variable = list("varName" = paste(varNames[ind.var], "EOFs", sep = "_"), "level" = levs[ind.var])
  out <- list("Variable" = Variable,
              "Data" = Data,
              "xyCoords" = xyCoords,
              "Dates" = Dates)
  return(out)
}