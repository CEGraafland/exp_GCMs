####################################################################################################
# Creation of precision network by conversion of HC CMIP BNs 
# one for one
####################################################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(igraph)
library(bnlearn)
library(sparsebn)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/backpermutations.rda")
load("data/tas_interim_10d_akima_cubic.rda")
############################################################################
# Function to load hciterations of models in a list 
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
###########################################################################################
# initial parameters
###########################################################################################
itfirst <- 1700
itlast <- 1800
# it <- "1700_1800|1800_1900"
it <- "1700_1800"
permk <- 3
###############################################################################
# Load hcs; choose one map
#################################################################################
string <- "data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda"
shortstring <- "CMIP5_EARTH_ONLYHIST"

string <- "tas_ncep_10d"
shortstring <- "ncep_10d"

string <- "data/tas_interim_10d_akima_cubic.rda"
shortstring <- "interim_10d_akima_cubic"

string <- "data/tas_JRA55_10d_akima_cubic.rda"
shortstring <- "JRA55_10d_akima_cubic"

string <- "data/tas_historical_cmip5_left_10d_akima_cubic.rda"
shortstring <- "CMIP5_left"

string <- "data/tas_historical_cmip5_extra_10d_akima_cubic.rda"
shortstring <- "CMIP5_extra"

string <- "data/tas_historical_cmip6_10d_akima_cubic_corrected.rda"
shortstring <- "CMIP6"

string <- "data/tas_historical_10d_akima_cubic_corrected.rda"
shortstring <- ""

string <- "data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda"
shortstring <- ""

string <- "data/tas_rcp85_cmip5_left_10d_akima_cubic.rda"
shortstring <- "FUTURE_CMIP5_left"

string <- "data/tas_ssp585_cmip6_10d_akima_cubic_corrected.rda"
shortstring <- "FUTURE_CMIP6"

string <- "data/tas_ssp585_cmip6_left_10d_akima_cubic_corrected.rda"
shortstring <- "FUTURE_CMIP6_left"

if(string == "data/tas_historical_cmip6_10d_akima_cubic_corrected.rda"){
  cubicsets <- get(load(paste0(string)))
  hc_gcms <- c(lapply(paste0(shortstring,"/",gsub(names(cubicsets), pattern = "_historical",replacement ="")), loadIterations, permused = permk, it = it))
  names(hc_gcms) <- gsub(names(cubicsets), pattern = "_historical",replacement ="")
} else if (string == "data/tas_interim_10d_akima_cubic.rda" ||string == "data/tas_JRA55_10d_akima_cubic.rda"){
  cubicsets <- list(get(load(paste0(string))))
  hc_gcms <-  c(lapply(paste0(shortstring), loadIterations, permused = permk, it = it))
  names(hc_gcms) <- shortstring
} else if(string == "tas_ncep_10d") {
  cubicsets <- list(get(load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")))
  hc_gcms <- lapply(c("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 
  names(hc_gcms) <- shortstring
  } else {cubicsets <- get(load(paste0(string)))
  if(!isTRUE(grepl("rcp85", string, fixed = TRUE)) & !isTRUE(grepl("ssp585", string, fixed = TRUE))){hc_gcms <- c(lapply(paste0(shortstring,"/",names(cubicsets)), loadIterations, permused = permk, it = it))
  names(hc_gcms) <- names(cubicsets)}
  if(isTRUE(grepl("rcp85", string, fixed = TRUE)) | isTRUE(grepl("ssp585", string, fixed = TRUE))){hc_gcms <- c(lapply(paste0(shortstring,"/FUTURE_",names(cubicsets)), loadIterations, permused = permk, it = it))
  names(hc_gcms) <- names(cubicsets)}
  }

hc_sizes <- lapply(hc_gcms,function(x) sapply(x,narcs))
hc_gcms_ordered <- mapply(function (x,y)x[order(y)], x = hc_gcms, y = hc_sizes,SIMPLIFY = FALSE) #analog hc_interims
##################################################################
# Make al datasets standardized and select hcs in loaded hcs
##################################################################
permk <- 3
# Anomalies and standarization over seasons (networks are learned with this set)
data_gcms_out <- lapply(cubicsets,function(x) TimeCoordsAnom_from_Grid_rms(x, rms = TRUE))
data_gcms_out_df <- lapply(data_gcms_out, function(x) as.data.frame(x))

# Global standarazation over the above: FIELD
data_gcms_anom_out <- lapply(data_gcms_out, function(x) mat2Dto3Darray(x,attr(x,"Xcoords"), attr(x,"Ycoords")))
grid_gcms_anom <- mapply(function(x,y) {x$Data <- y ;return(x)}, x = cubicsets, y = data_gcms_anom_out,SIMPLIFY = FALSE)
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
# assign(paste0("data_gcms_anom_scaled_",permk),lapply(data_gcms_anom_scaled, function(x) x[,permutations[[permk]]]))

# Choose between 'own optimums' or constant magnitude
# modelsize <- own_optimums
modelsize <- 18

selection_hc_gcms <- lapply(hc_gcms, function(x,y) {x[grep((y),names(x))]}, y = modelsize)
narcs_selection_hc_gcms <- lapply(selection_hc_gcms,function(x) lapply(x,narcs))

# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms <- mapply(function(x,y) lapply(x,function(x,y)bn.fit(x = x, data = y),y = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function(x,y) lapply(x,function(x,y) bn.fit(x = x, data = y),y = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
######################################################################
# Transition to precision matrix
######################################################################
f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))

NELShc3 <- lapply(selection_fits_gcms,function(x) lapply(x,as.graphNEL))
edgelistssparsehc3 <- lapply(NELShc3,function(x) lapply(x,as.edgeList))
sparsedatas <- lapply(data_gcms_out_df,function(x) sparsebnData(x[,permutations[[permk]]], type = "continuous"))

# Precisions
inv_cov_mat <- function(coefs, vars){
  # Checks: nrow = ncol
  
  pp <- nrow(coefs)
  identity_mat <- Matrix::Diagonal(pp, rep(1, pp))
  (identity_mat - coefs) %*% Matrix::solve(vars) %*% Matrix::t(identity_mat - coefs)
}

fitteddags <- mapply(function(x,y) lapply(x, function(x,y)estimate.parameters(fit = x, data = y),y = y), x = edgelistssparsehc3, y = sparsedatas,SIMPLIFY  = FALSE)
PRECShc3 <- lapply(fitteddags,function(x) lapply(x, function(x) inv_cov_mat(Matrix::Matrix(x$coefs), Matrix::Matrix(x$vars))))
precmats <- lapply(PRECShc3,function(x) lapply(x, as.matrix))

f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
PartVarsHCtoPN <- lapply(precmats, function(x)lapply(x, function(m) outer(1:nrow(m),1:ncol(m),Vecf,m)))

# edgesprecsIT <- lapply(precmats, function(x) lapply(x,function(x){ diag(x) <- 0 
# x[abs(x) != 0] <- 1
# return(x)}
# ))
# precsITgraphs <- lapply(edgesprecsIT, function(x) lapply(x, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE))
precsITgraphs <- lapply(PartVarsHCtoPN, function(x) lapply(x, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE,diag = FALSE))
nedgesprecsIT<- sapply(precsITgraphs, function(x)lapply(x,function(x) length(E(x))))
# names(precsITgraphs) <- mapply(function(x,z,y)paste0(x,"_BN",z,"_PN",y),x = names(precsITgraphs),y = nedgesprecsIT,z = narcs_selection_hc_gcms)
# precIT_graphs <- precgraphs[order(nedgesprecsIT)] 

# HIER VERDER GAAN.
init_graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_interim_10d_akima_cubic)
graphObjects <- lapply(precsITgraphs, 
       function(x) {
         lapply(x, function(x){
           graphObject <- init_graphObject;
           graphObject$graph <- x
           graphObject$adjacency <- as_adjacency_matrix(x)
           return(graphObject)
         })
       }
)
#######################################################################################
# Calculate distances HCtoPN
########################################################################################
glassoigraphs <- lapply(graphObjects,function(x) lapply(x, function(x) x$graph))
if(shortstring == ""){
  if(!isTRUE(grepl("rcp85", string, fixed = TRUE))){assign(paste0("graphs_CMIP5_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk),lapply(glassoigraphs, function(x) lapply(x, igraph.distances, perm = permutations[[permk]], data.igraph = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))))
  save(list = paste0("graphs_CMIP5_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/struct_edgeattributes/perm",permk,"/graphs_CMIP5_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk,".rda"))}
  if(isTRUE(grepl("rcp85", string, fixed = TRUE))){assign(paste0("graphs_FUTURE_CMIP5_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk),lapply(glassoigraphs, function(x) lapply(x, igraph.distances, perm = permutations[[permk]], data.igraph = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))))
    save(list = paste0("graphs_FUTURE_CMIP5_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/struct_edgeattributes/perm",permk,"/graphs_FUTURE_CMIP5_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk,".rda"))}
  } else {
  assign(paste0("graphs_",shortstring,"_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk),lapply(glassoigraphs, function(x) lapply(x, igraph.distances, perm = permutations[[permk]], data.igraph = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))))
  save(list = paste0("graphs_",shortstring,"_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/struct_edgeattributes/perm",permk,"/graphs_",shortstring,"_HC_",itfirst,"_",itlast,"_toPN_dist_perm",permk,".rda"))
}

