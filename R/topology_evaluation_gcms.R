##########################################################################
# Topology evaluation 
##########################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
rm(list = ls())
library(igraph)
library(maps)
library(geosphere)
library(broom)
library(ggplot2)
library(visualizeR)
library(bnlearn)
library(gridExtra)
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("../R/Functions/BasicNetworkFunctions.R")
load(file = "../Data/Struct_learn/permutations.rda")
load(file = "../Data/Struct_learn/backpermutations.rda")

load("data/tas_rcp85_cmip5_left_10d_akima_cubic.rda")
load("data/tas_ssp585_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_ssp585_cmip6_left_10d_akima_cubic_corrected.rda")
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
listinterim <-list(get(load("data/tas_interim_10d_akima_cubic.rda")))
names(listinterim) <- "interim_10d_akima_cubic"

listjra55 <- list(get(load("data/tas_JRA55_10d_akima_cubic.rda")))
names(listjra55) <- "JRA55_10d_akima_cubic"

listncep <- list(get(load("../Data/tas_ncep_10d.rda")))
names(listncep) <- "ncep_10d"

cubicsets5 <- get(load("data/tas_historical_10d_akima_cubic_corrected.rda"))
cubicsets5extra <- get(load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda"))
cubicsetsearth <- get(load("data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda"))
cubicsets5left <- get(load("data/tas_historical_cmip5_left_10d_akima_cubic.rda"))
cubicsets6 <- get(load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda"))

cubicsets5.rcp85 <- load("data/tas_rcp85_cmip5_10d_akima_cubic_corrected.rda")
# Leave out bcc "CMIP5_bcc.csm1.1.m_r1i1p1_rcp85" because does not reach 2100 (one year atrasado)
xout <- which(names(tas_rcp85_cmip5_left_10d_akima_cubic) == "CMIP5_bcc.csm1.1.m_r1i1p1_rcp85")
cubicsets5left.rcp85 <- tas_rcp85_cmip5_left_10d_akima_cubic[-xout]
cubicsets6.ssp585 <- tas_ssp585_cmip6_10d_akima_cubic_corrected
cubicsets6left.ssp585 <- tas_ssp585_cmip6_left_10d_akima_cubic_corrected



cubicsets <- c(tas_historical_10d_akima_cubic_corrected,listinterim)
namescubics <- names(cubicsets)
# With CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,cubicsets6,listinterim,listjra55,listncep)
names(cubicsets)
# Without CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,listinterim,listjra55,listncep)
names(cubicsets)
# Without earth
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,listinterim,listjra55,listncep)
names(cubicsets)
# Future CMIP5
cubicsets <- c(cubicsets5.rcp85,cubicsets5left.rcp85)
# Future CMIP6
cubicsets <- c(cubicsets6.ssp585,cubicsets6left.ssp585)

namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
namescubicsearth <- names(cubicsetsearth)
namescubics <- names(cubicsets)
namescubics5.rcp85 <- names(cubicsets5.rcp85)
namescubics5left.rcp85 <- names(cubicsets5left.rcp85)
namescubics6.ssp585 <- names(cubicsets6.ssp585)
namescubics6left.ssp585 <- names(cubicsets6left.ssp585)


shortnames <- gsub(gsub(names(cubicsets), pattern = "_r1i1p1", replacement = ""),
                   pattern = "_r12i1p1", replacement = "")
shortnames <- #gsub(
  gsub(names(cubicsets), pattern = "_r1i1p1", replacement = "")
#,
#pattern = "_r12i1p1", replacement = "")
shortnames2 <- names(cubicsets)
abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")
############################################################################
# Functions to load hciterations of models in a list 
############################################################################
# loadIterations <- function(pattern,permused, hctype = NULL){
#   if(is.null(hctype)){hctype <- ""}
#   hc_interim_list <- list.files(paste0("data/hciterations/",pattern,"/perm",permused,hctype), full.names = T)
#   hc_interim_names <- list.files(paste0("data/hciterations/",pattern,"/perm",permused,hctype))
#   hc_interim_names <- gsub(".rda", "", hc_interim_names)
#   
#   hc_interim_networks <- list()
#   
#   hc_interim_networks <- lapply(hc_interim_list, function(x){get(load(x))})
#   names(hc_interim_networks) <- hc_interim_names
#   interimsizes <- sapply(hc_interim_networks,narcs)
#   hc_interims <- hc_interim_networks[order(interimsizes)]
#   return(hc_interims)
# }
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
it <- "1700_1800"
permk <- 3

hc_gcms5 <- lapply(paste0(namescubics5), loadIterations, permused = permk,it = it)
hc_gcms5left<- lapply(paste0("CMIP5_left/",namescubics5left), loadIterations, permused = permk, it = it)
hc_gcms5extra <- lapply(paste0("CMIP5_extra/",namescubics5extra), loadIterations, permused = permk, it = it)
hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk, it = it)
hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_jra55 <- lapply(c("JRA55_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_ncep <- lapply(c("../Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 


hc_gcms5.rcp85 <- lapply(paste0("FUTURE_",namescubics5.rcp85), loadIterations, permused = permk, it = it)
hc_gcms5left.rcp85<- lapply(paste0("FUTURE_CMIP5_left/FUTURE_",namescubics5left.rcp85), loadIterations, permused = permk, it = it)
hc_gcms6.ssp585 <- lapply(paste0("FUTURE_CMIP6/FUTURE_",namescubics6.ssp585), loadIterations, permused = permk, it = it)
hc_gcms6left.ssp585 <- lapply(paste0("FUTURE_CMIP6_left/FUTURE_",namescubics6left.ssp585), loadIterations, permused = permk, it = it)

hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_gcms6,hc_interim)
names(hc_gcms) <- shortnames2

hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_gcms6,hc_interim)

hc_gcms <- c(hc_gcms5.rcp85,hc_gcms5left.rcp85)
names(hc_gcms) <- shortnames2

hc_gcms <- c(hc_gcms6.ssp585,hc_gcms6left.ssp585)
names(hc_gcms) <- shortnames2

hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_interim,hc_jra55,hc_ncep)
names(hc_gcms) <- shortnames2

modelsize <- 18
selection_hc_gcms <-  mapply(function(x,y) x[[grep((y),x)]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
# 
# selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
# 
# selection_hc_gcms$CMIP5_ACCESS1.0_r1i1p1
# selection_hc_gcms$CMIP5_BNU.ESM_r1i1p1
abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")
##########################################################################
# Test
##########################################################################
data.network <- TimeCoordsAnom_from_Grid_rms(cubicsets$interim_10d_akima_cubic, rms = TRUE)
x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)
# x <- (x+360)%%360
p <- expand.grid(y, x)[2:1][permutations[[3]],]

permstrenght <- bn_to_igraph.strengths(hc_gcms$CMIP5_NorESM1.M$CMIP5_NorESM1.M_3_1700_1800i, perm = permutations[[3]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))
permdist <- igraph.distances(permstrenght, perm = permutations[[3]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))
permweight <- igraph.weights(permdist, type = "bn", fromdist = 2000)

graphBN <- permweight
redtype <- "BN"

assign("graph", eval(parse(text = paste0("graph",redtype))))
# graph <- graphCN
edgelist <- as_edgelist(graph) 
nrow(edgelist)
row.names(p) <- names(V(graph))
gridpoints <- names(V(graph))

dists <- edge.attributes(graphBN)$distances
distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
distsc <- as.character(distsv)

col.1 <- adjustcolor("light blue", alpha=0.8)
col.2 <- adjustcolor("navy blue", alpha=0.8)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

arcslist <- list()
for(i in 1:nrow(edgelist))  {
  node1 <- gridpoints[gridpoints == edgelist[i,1]]
  
  node2 <- gridpoints[gridpoints == edgelist[i,2]]
  
  arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                        c(p[node2,]$Var2, p[node2,]$Var1),
                        n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
  
  edge.ind <- ceiling(edge_attr(graph,"distances",E(graph)[i])*100 / max(dists))
  col <- as.character(edge.col[edge.ind])
  negs <- which(arc[,'lon']<0)
  
  if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
  arc <- cbind(arc,edge.ind)
  # lines(arc, col=edge.col[edge.ind], lwd=edge.ind) interesting for edgewidth!
  # points(x = arc[1,1], y = arc[1,2],col = "pink", pch = 19)
  arcslist[[i]]<- arc
  
}

naampjes <- character()
for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
names(arcslist) <- naampjes

mapworld2 <- map("world2", fill = TRUE)
world_pol_df <- tidy(mapworld2, IDs = "region")

do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
  cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
})) -> all_paths



bla <- ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = edge.ind)) +
  scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
  ggtitle(paste0(redtype,": ",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

bla



plot_long_distances(hc_gcms$CMIP5_EC.EARTH$CMIP5_EC.EARTH_3_1700_1800i, data.network, minimdist = 10000, smallcol = NA, perm = permutations[[3]])
plot_long_distances(hc_gcms$interim_10d_akima_cubic$interim_10d_akima_cubic_3_1700_1800i, data.network, minimdist = 10000, smallcol = NA, perm = permutations[[3]])
plot_long_distances(hc_gcms$CMIP5_NorESM1.M$CMIP5_NorESM1.M_3_1700_1800i, data.network, minimdist = 10000, smallcol = NA, bigcol = "red", perm = permutations[[3]])

#################################################################################
# Plots of networks with ggplot
#################################################################################
data.network <- TimeCoordsAnom_from_Grid_rms(cubicsets[[1]], rms = TRUE)
x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)
# x <- (x+360)%%360
p <- expand.grid(y, x)[2:1][permutations[[3]],]

# modelsize <- 18
# selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)

permstrengths <- lapply(selection_hc_gcms, bn_to_igraph.strengths, perm = permutations[[3]],data.dag = TimeCoordsAnom_from_Grid_rms(listinterim$interim_10d_akima_cubic, rms = TRUE))
permdists <- lapply(permstrengths, igraph.distances, perm = permutations[[3]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(listinterim$interim_10d_akima_cubic, rms = TRUE))
permweights <- lapply(permdists, igraph.weights, type = "bn", fromdist = 0)


edgelists <- lapply(permweights, as_edgelist) 
row.names(p) <- names(V(permweights[[1]]))
gridpoints <- names(V(permweights[[1]]))

# row.names(p) <- names(V(graph))
# gridpoints <- names(V(graph))


col.1 <- adjustcolor("light blue", alpha=0.8)
col.2 <- adjustcolor("navy blue", alpha=0.8)
col.3 <- adjustcolor("dark green", alpha=0.8)
col.2 <- adjustcolor("orange", alpha=0.8)
col.1 <- adjustcolor("red", alpha=0.8)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

dists4plot <- function(permweight,edgelist){
arcslist <- list()
for(i in 1:nrow(edgelist))  {
  
  graphBN <- permweight
  redtype <- "BN"
  
  assign("graph", eval(parse(text = paste0("graph",redtype))))
  
  dists <- edge.attributes(graphBN)$distances
  distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
  distsc <- as.character(distsv)
  
  node1 <- gridpoints[gridpoints == edgelist[i,1]]
  
  node2 <- gridpoints[gridpoints == edgelist[i,2]]
  
  arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                        c(p[node2,]$Var2, p[node2,]$Var1),
                        n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
  
  edge.dist <- edge_attr(graph,"distances",E(graph)[i])
  edge.ind <- ceiling(edge_attr(graph,"distances",E(graph)[i])*100 / max(dists))
  col <- as.character(edge.col[edge.ind])
  negs <- which(arc[,'lon']<0)
  
  if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
  arc <- cbind(arc,edge.ind)
  # lines(arc, col=edge.col[edge.ind], lwd=edge.ind) interesting for edgewidth!
  # points(x = arc[1,1], y = arc[1,2],col = "pink", pch = 19)
  arcslist[[i]]<- arc
  
}

naampjes <- character()
for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
names(arcslist) <- naampjes

do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
  cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
})) -> all_paths

return(all_paths)
}

# all_paths_test <- dists4plot(edgelists$CMIP5_CanESM2,permweight = permweights$CMIP5_CanESM2)
all_paths_list <- mapply(dists4plot, edgelist = edgelists, permweight = permweights, SIMPLIFY = FALSE)


mapworld2 <- map("world2", fill = TRUE)
world_pol_df <- tidy(mapworld2, IDs = "region")


ggplotfunction <- function(all_paths,edgelist,shortname,col = NULL) {
  if (is.null(col)){ ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = edge.ind)) +
  scale_color_gradient2(low = col.1, mid = col.2, high = col.3,midpoint = 37.5,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
  ggtitle(paste0(redtype,": ",nrow(edgelist)," ",length(unique(all_paths$route))," ",shortname)) +
  theme(plot.title = element_text(hjust = 0.5))
  } else {  ggplot() +
      geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
      geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = "black")) +
      # scale_color_gradient2(low = col.1, mid = col.2, high = col.3,midpoint = 37.5,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
      ggtitle(paste0(redtype,": ",nrow(edgelist)," ",length(unique(all_paths$route))," ",shortname)) +
      theme(plot.title = element_text(hjust = 0.5))}
    
  }

distsc <- as.character(c(0,5000,10000,15000,20000))
redtype <- "HC"
ggplots <- mapply(ggplotfunction, all_paths_list, edgelists, abrev, SIMPLIFY = FALSE)
do.call("grid.arrange",ggplots)
names(all_paths_list)
names(edgelists)

########################################################################################
# Sub example with part of edges.  
########################################################################################
all_paths_list1 <- all_paths_list
all_paths_list2 <- all_paths_list
all_paths_list3 <- all_paths_list
all_paths_list4 <- all_paths_list
all_paths_list5 <- all_paths_list

for (x in 1:length(all_paths_list2)){
  all_paths_list1[[x]] <- all_paths_list[[x]][which(all_paths_list[[x]]$edge.ind < 25),]
  all_paths_list2[[x]] <- all_paths_list[[x]][which(all_paths_list[[x]]$edge.ind > 25 & all_paths_list[[x]]$edge.ind < 50),]
  all_paths_list3[[x]] <- all_paths_list[[x]][which(all_paths_list[[x]]$edge.ind > 50),]
}

for (x in 1:length(all_paths_list2)){
  all_paths_list4[[x]] <- all_paths_list[[x]][which(all_paths_list[[x]]$edge.ind <= 50),]
  all_paths_list5[[x]] <- all_paths_list[[x]][which(all_paths_list[[x]]$edge.ind > 50),]
}

names(ggplots)
subbies <- c("CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_BNU.ESM_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")
subbies <- c("CMIP5_NorESM1.M_r1i1p1_rcp85","CMIP5_GFDL.ESM2M_r1i1p1_rcp85","CMIP5_GFDL.ESM2G_r1i1p1_rcp85" ,"CMIP5_MPI.ESM.MR_r1i1p1_rcp85","CMIP5_MPI.ESM.LR_r1i1p1_rcp85" ,"CMIP5_CNRM.CM5_r1i1p1_rcp85","CMIP5_MIROC5_r1i1p1_rcp85","CMIP5_MIROC.ESM_r1i1p1_rcp85","CMIP5_MIROC.ESM.CHEM_r1i1p1_rcp85","CMIP5_MRI.CGCM3_r1i1p1_rcp85","CMIP5_CanESM2_r1i1p1_rcp85" ,"CMIP5_EC.EARTH_r1i1p1_rcp85","CMIP5_CESM1.BGC_r1i1p1_rcp85"   , "CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85"  ,"CMIP5_IPSL.CM5A.LR_r1i1p1_rcp85")
subbies <- c("CMIP6Amon_NorESM2.MM_ssp585_r1i1p1f1","CMIP6Amon_NorESM2.LM_ssp585_r1i1p1f1","CMIP6Amon_GFDL.ESM4_ssp585_r1i1p1f1","CMIP6Amon_MPI.ESM1.2.HR_ssp585_r1i1p1f1", "CMIP6Amon_MIROC.ES2L_ssp585_r1i1p1f2"   ,"CMIP6Amon_MIROC6_ssp585_r1i1p1f1" ,"CMIP6Amon_MRI.ESM2.0_ssp585_r1i1p1f1","CMIP6Amon_CanESM5_ssp585_r1i1p1f1","CMIP6Amon_CanESM5_ssp585_r1i1p2f1"   ,"CMIP6Amon_EC.Earth3.Veg_ssp585_r1i1p1f1","CMIP6Amon_EC.Earth3_ssp585_r1i1p1f1", "CMIP6Amon_CESM2.WACCM_ssp585_r1i1p1f1"  ,"CMIP6Amon_CESM2_ssp585_r1i1p1f1","CMIP6Amon_CNRM.ESM2.1_ssp585_r1i1p1f2","CMIP6Amon_CNRM.CM6.1_ssp585_r1i1p1f2","CMIP6Amon_IPSL.CM6A.LR_ssp585_r1i1p1f1")
#subbies <- c("CMIP5_NorESM1.M_r1i1p1_rcp85","CMIP5_GFDL.ESM2M_r1i1p1_rcp85","CMIP5_GFDL.ESM2G_r1i1p1_rcp85" ,"CMIP5_MPI.ESM.MR_r1i1p1_rcp85","CMIP5_MPI.ESM.LR_r1i1p1_rcp85" ,"CMIP5_CNRM.CM5_r1i1p1_rcp85","CMIP5_MIROC5_r1i1p1_rcp85","CMIP5_MIROC.ESM_r1i1p1_rcp85","CMIP5_MIROC.ESM.CHEM_r1i1p1_rcp85","CMIP5_MRI.CGCM3_r1i1p1_rcp85","CMIP5_CanESM2_r1i1p1_rcp85" ,"CMIP5_EC.EARTH_r1i1p1_rcp85","CMIP5_CESM1.BGC_r1i1p1_rcp85"   , "CMIP5_IPSL.CM5A.MR_r1i1p1_rcp85"  ,"CMIP5_IPSL.CM5A.LR_r1i1p1_rcp85")
subbies.ind <- sapply(subbies, function(x)which(x==names(ggplots)))

pdf(paste0('figs/future_overlap/CMIP6_overlapsubset_gglongplots_25.pdf'), width = 20, height = 20)
pdf(paste0('figs/future_overlap/CMIP5_overlapsubset_gglongplots_25.pdf'), width = 20, height = 20)
pdf(paste0('figs/exampleclustering/subset_gglongplots_25.pdf'), width = 10, height = 10)
ggplots_subbies <- mapply(ggplotfunction, all_paths_list1[subbies], edgelists[subbies], abrev[subbies.ind], SIMPLIFY = FALSE)
do.call("grid.arrange",ggplots_subbies)
dev.off()

pdf(paste0('figs/future_overlap/CMIP6_overlapsubset_gglongplots_25_50.pdf'), width = 20, height = 20)
pdf(paste0('figs/future_overlap/CMIP5_overlapsubset_gglongplots_25_50.pdf'), width = 20, height = 20)
pdf(paste0('figs/exampleclustering/subset_gglongplots_25_50.pdf'), width = 10, height = 10)
ggplots_subbies <- mapply(ggplotfunction, all_paths_list2[subbies], edgelists[subbies], abrev[subbies.ind], SIMPLIFY = FALSE)
do.call("grid.arrange",ggplots_subbies)
dev.off()

pdf(paste0('figs/future_overlap/CMIP6_overlapsubset_gglongplots_50.pdf'), width = 20, height = 20)
pdf(paste0('figs/future_overlap/CMIP5_overlapsubset_gglongplots_50.pdf'), width = 20, height = 20)
pdf(paste0('figs/exampleclustering/subset_gglongplots_50.pdf'), width = 10, height = 10)
ggplots_subbies <- mapply(ggplotfunction, all_paths_list3[subbies], edgelists[subbies], abrev[subbies.ind], SIMPLIFY = FALSE)
do.call("grid.arrange",ggplots_subbies)
dev.off()

navyblue<- adjustcolor("navy blue", alpha=0.8)
pdf(paste0('figs/exampleclustering/subset_gglongplots_alldist.pdf'), width = 10, height = 10)
ggplots_subbies <- mapply(ggplotfunction, all_paths_list[subbies], edgelists[subbies], abrev[subbies.ind], MoreArgs = list(col = navyblue) ,SIMPLIFY = FALSE)
do.call("grid.arrange",ggplots_subbies)
dev.off()

ggplots_subbies[[1]]
pdf(paste0('figs/exampleclustering/subset_gglongplots_0_50.pdf'), width = 10, height = 10)
ggplots_subbies <- mapply(ggplotfunction, all_paths_list4[subbies], edgelists[subbies], abrev[subbies.ind], SIMPLIFY = FALSE)
do.call("grid.arrange",ggplots_subbies)
dev.off()

pdf(paste0('figs/exampleclustering/subset_gglongplots_51_100.pdf'), width = 10, height = 10)
ggplots_subbies <- mapply(ggplotfunction, all_paths_list5[subbies], edgelists[subbies], abrev[subbies.ind], SIMPLIFY = FALSE)
do.call("grid.arrange",ggplots_subbies)
dev.off()
############################################################################
# Long edge plots with plot_long_distances
############################################################################
library("raster")
pdf('figs/longplots.pdf', width = 10, height = 20)
par(mfrow = c(4,3))
lapply(selection_hc_gcms, plot_long_distances, minimdist = 10000, smallcol = NA, bigcol = "red", data.dag = data.network, perm = permutations[[3]])
dev.off()
#####################################################################
# Plot networks with distances larger than X Km
#####################################################################
pdf('figs/longSplots.pdf', width = 10, height = 13)
par(mfrow = c(4,3))
mapply(plotS.spatial.igraph,network = permweights, title = abrev, MoreArgs = list(data.network = data.network,type = "bn",th = 0, smallcol = NA,th.type = "zoomweighted", from.dist = 10000,by.dist = 10000, perm = permutations[[3]],shift = TRUE))
dev.off()

#####################################################################
# Strength (loglikelihood lose) vs distance plot
#####################################################################

distplotfunction <- function(g,ab) {
  plot(E(g)$distances,E(g)$strengths,
      xlab = "distance", ylab = "strengths",
      main = paste0(length(E(g))," from  = ",min(E(g)$distances)," ",ab),
      xlim = c(0,19000) 
      #,ylim = c(-1200,0)
      )
}



pdf('figs/arcplots.pdf', width = 10, height = 13)
par(mfrow = c(4,3))
mapply(distplotfunction, g = permweights, ab = abrev)
dev.off()

############################################################################
# Complex network measures
############################################################################
############################################################################
# inverse weighted betweenness
############################################################################
weights3bnlist <- lapply(permweights, function(x)E(x)$weights)
invweights3bnlist <- lapply(weights3bnlist, function(x)1/x)

betweennessbn <- mapply(betweenness,permweights, weights = invweights3bnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
logbetweennessbn <- lapply(betweennessbn, function(x)log(1+x))
############################################################################
# unweighted betweenness
############################################################################
weights1bnlist <- lapply(permweights, function(x)E(x)$weight)
betweennessbn <- mapply(betweenness,permweights, weights = weights1bnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
logbetweennessbn <- lapply(betweennessbn, function(x)log(1+x))

###########################################################################
movmed_betw_bn <- lapply(logbetweennessbn, movingmedias)
movmed_betw_bn_mean <- lapply(movmed_betw_bn, function(x) apply(x, MARGIN = 3, FUN = mean))
climbetwbn <- lapply(movmed_betw_bn_mean, quantity2clim, what = "betweenness",ref.grid = tas_interim_10d_akima_cubic, backperm = backpermutations[[1]])


plotsbn <- list()
measbn <- climbetwbn
for (i in 1:length(climbetwbn)){
  plot <- spatialPlot(climbetwbn[[i]], main = list(paste0("BN:",abrev[i]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      lonCenter = 180,
                      col = "Reds",
                      colorkey = list(width = 0.6, lables = list(cex = 0.5))
  )
  plotsbn[[i]] <- plot
}
# First method: Not equal colorscale
do.call(grid.arrange, c(plotsbn, nrow = 2, as.table = TRUE, top = paste0("Moving Medias 1st perm BN: Log(1 + ",attr(measbn[[1]]$Data, "climatology:fun"),")")))
# Second method: Equal colorscale
ensbetwclims <- bindGrid(climbetwbn, dimension = c("member"))
ensbetwclims$Members <- abrev
spatialPlot(ensbetwclims,backdrop.theme = "coastline", as.table = TRUE, lonCenter = 180, col = "Reds", rev.colors = FALSE,main = attr(climbetwbn$interim_10d_akima_cubic$Data,"climatology:fun"))
###########################################################################
# edge betweenness
###########################################################################
weights1bnlist <- lapply(permweights, function(x)E(x)$weight)
edgebetweennessbn <- mapply(edge.betweenness,permweights, weights = weights1bnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
commu_betsbn_unw_gcms_perm3_1800$CMIP5_CanESM2$edge.betweenness

cluster_edge_betweenness()
plot(E(permweights$CMIP5_CanESM2)$strengths,commu_betsbn_unw_gcms_perm3_1800$CMIP5_CanESM2$edge.betweenness)
plot(E(permweights$CMIP5_CanESM2)$strengths,commu_betsbn_unw_gcms_perm3_1800$CMIP5_CanESM2$removed.edges)

# ,
#      xlab = "distance", ylab = "strengths",
#      main = paste0(length(E(g))," from  = ",min(E(g)$distances)," ",ab)
#      # xlim = c(0,19000) 
#      #,ylim = c(-1200,0)
# )

betwplotfunction <- function(g,ab) {
  plot(edge,E(g)$strengths,
       xlab = "distance", ylab = "strengths",
       main = paste0(length(E(g))," from  = ",min(E(g)$distances)," ",ab),
       xlim = c(0,19000) 
       #,ylim = c(-1200,0)
  )
}



pdf('figs/arcplots.pdf', width = 10, height = 13)
par(mfrow = c(4,3))
mapply(distplotfunction, g = permweights, ab = abrev)
dev.off()

#########################################################################
# unweighted closeness
#########################################################################
closenessbn <-  mapply(closeness, graph = permweights, weights = weights1bnlist, normalized = TRUE, SIMPLIFY = FALSE)
climclosenessbn <- lapply(closenessbn, quantity2clim, what = "closeness",ref.grid = tas_interim_10d_akima_cubic, backperm = backpermutations[[3]])
# First method: equal colorscale
ensclosclims <- bindGrid(climclosenessbn, dimension = c("member"))
ensclosclims$Members <- abrev
spatialPlot(ensclosclims,backdrop.theme = "coastline", as.table = TRUE, lonCenter = 180, col = "Reds", rev.colors = FALSE, set.max = 0.010,main = attr(climclosenessbn$interim_10d_akima_cubic$Data,"climatology:fun"))
# Second method: not equal colorscale
plotsbn <- list()
measbn <- climclosenessbn
for (i in 1:length(climbetwbn)){
  plot <- spatialPlot(climclosenessbn[[i]], main = list(paste0("BN:",abrev[i]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      lonCenter = 180,
                      col = "Reds",
                      colorkey = list(width = 0.6, lables = list(cex = 0.5))
  )
  plotsbn[[i]] <- plot
}
do.call(grid.arrange, c(plotsbn, nrow = 2, as.table = TRUE, top = paste0("Moving Medias 1st perm BN: Log(1 + ",attr(measbn[[1]]$Data, "climatology:fun"),")")))
#################################################################################
# area weighted connectivity
#################################################################################

awconnectivity <- function(graph,data,backperm){
  adjmat <- as.matrix(as_adjacency_matrix(graph))
  adjmat <- adjmat[backperm,backperm]
  ys <- attr(data,"VertexCoords")$y
  sumArea <- sum(cos(ys/(180)*pi))
  awconnectivity <- as.vector(cos(ys/(180)*pi)%*%adjmat) / sumArea
  return(awconnectivity)
  }

awsbn <-  mapply(awconnectivity , graph = permweights, MoreArgs = list(data = data.network, backperm = backpermutations[[3]]), SIMPLIFY = FALSE)
climawsbn <- lapply(awsbn, quantity2clim, what = "awconnectivity",ref.grid = tas_interim_10d_akima_cubic, backperm = backpermutations[[1]])
ensawsclims <- bindGrid(climawsbn, dimension = c("member"))
ensawsclims$Members <- abrev
spatialPlot(ensawsclims,backdrop.theme = "coastline", as.table = TRUE, lonCenter = 180, 
            col = "Reds", 
            rev.colors = FALSE, main = attr(climawsbn$interim_10d_akima_cubic$Data,"climatology:fun"))

#########################################################################
# Comunidades. 
#########################################################################
library(RColorBrewer)
##############################################################
# Unweighted 
##############################################################
commu_bet <- lapply(permweights, function (x) cluster_edge_betweenness(as.undirected(x), weights = E(x)$weight))

name <- paste0("commu_betsbn_unw_gcms_perm3_",modelsize*100)
assign(paste0("commu_betsbn_unw_gcms_perm3_",modelsize*100), commu_bet)
save(list = paste0("commu_betsbn_unw_gcms_perm3_",modelsize*100),file = paste0("results/communities/",name,".rda"))

load(paste0("results/communities/",name,".rda"))
commu_bet <- commu_betsbn_unw_gcms_perm3_1800

plotsbn <- list()

edgeweights <- "unw"

for(i in 1:length(commu_bet)){
  ncom <- length(levels(factor(commu_bet[[i]]$membership)))
  mem <- membership(commu_bet[[i]])
  memClim <- quantity2clim(mem, "membership", tas_interim_10d_akima_cubic, backperm = backpermutations[[3]])
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotsbn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                              set.min = NULL, set.max = NULL, 
                              lonCenter = 180, 
                              regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                              main = paste0("Comunities BN ",abrev[i]),
                              colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                              lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}
do.call(grid.arrange, c(plotsbn, ncol = 3))

commu_bet$CMIP5_CanESM2$bridges
#######################################################################
# Less merges / cut communities: Used to generate evolution plots
# Check backpermutations!
#######################################################################

plots <- list()
cuts <- 12

for(j in 1:length(commu_bet)){
  i <- cuts
  commu_bet12 <- cut_at(commu_bet[[j]],i)
  g <- permweights[[j]]
  
  
  ncom <- length(levels(factor(commu_bet12)))
  mem <- commu_bet12
  memClim <- quantity2clim(mem, "membership", tas_interim_10d_akima_cubic, backperm = backpermutations[[3]])
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plots[[j]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                            set.min = NULL, set.max = NULL, 
                            lonCenter = 180, 
                            regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                            main = paste0(length(E(g))," ",abrev[j]),
                            colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                            lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
  
}

plots
do.call(grid.arrange,c(plots))

###################################################################################
# b weights! BN
###################################################################################
# griddata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(cubicsets[[10]], rms = TRUE))


bmatgraph_from_fit <- function(fit1800){
parentslist <- lapply(coefficients(fit1800), function(x) {y <- names(x)[2:length(names(x))]
                                                          dim(y) <- c(1,length(y))
                                                          return(y)})
indexlist <- lapply(parentslist, function(x) {apply(x, MARGIN = 2, 
                                                    FUN = function(z) which(z == names(fit1800)))})
bmatrix <- mapply(function(x,y) {z <- array(data = 0, dim = c(1,648))
                                 z[x] <- y[2:length(y)]
                                 return(z)}, 
                  x = indexlist, y = coefficients(fit1800))
row.names(bmatrix) <- names(indexlist)
# diag(bmatrix) <- 1

bmatgraph <- graph_from_adjacency_matrix(bmatrix, mode = "directed", weighted = TRUE)
return(bmatgraph)
}
# set edge attributes

timecoords <- TimeCoordsAnom_from_Grid_rms(cubicsets[[10]], rms = TRUE)
set_edgeatr_bmatgraph <- function(bmatgraph, timecoords){
edgelist <- as_edgelist(bmatgraph)
names <- paste0("V",1:648)

indfrom <- c()
for(i in 1:length(edgelist[,1])){ indfrom[i] <- which(edgelist[i,1] == names)}
indto <- c()
for(i in 1:length(edgelist[,2])){ indto[i] <- which(edgelist[i,2] == names)}

# Make coordinates for every variable
longitude <- attributes(timecoords)$VertexCoords$x
lattitude <- attributes(timecoords)$VertexCoords$y

# estimate distance of all edges in igraph-edgelist.
distances <- c()
for (i in 1:nrow(edgelist)){
  x1Lat <- lattitude[indfrom[i]]
  x1Lon <- longitude[indfrom[i]]
  x2Lat <- lattitude[indto[i]]
  x2Lon <- longitude[indto[i]]
  
  disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
  distances[i] <- disti
  }

# set attributes
E(bmatgraph)$distances <- distances
E(bmatgraph)$strength <- E(bmatgraph)$weight
return(bmatgraph)}

# plot
network <- bmatgraph 
th <- 0.3
from.dist <- 5000
shift <- TRUE
curvature <- TRUE
perm <- permutations[[3]]

plot_bmatgraph <- function(network,th,from.dist,shift = TRUE,curvature = FALSE,perm = NULL,title){


E(network)[abs(strength) <= th]$color <- NA
E(network)[abs(strength) > th & strength > 0]$color <- "red"
E(network)[abs(strength) > th & strength < 0]$color <- "blue"
E(network)[distances <= from.dist]$color <- NA

# draw the graph
x <- attr(timecoords, "Xcoords", exact = FALSE)
y <- attr(timecoords, "Ycoords", exact = FALSE)

if(shift == TRUE){x <- (x+360)%%360}
points <- expand.grid(y, x)[2:1]
if (!is.null(perm)){
  points <- points[perm,]
}
if(shift == TRUE){map('world2',interior = FALSE,resolution = 0)
} else {raster::plot(lisworld)}

plot.igraph(network,
            vertex.size = 100,
            vertex.color = NA,
            vertex.label = NA,
            edge.curved = curvature,
            edge.arrow.size = 0.2,
            # edge.width = E(igraph)$width,
            lty = "dots",
            layout = as.matrix(points), add = TRUE, rescale = FALSE)

title(main = title)

}


bmatgraphs <- lapply(selection_fits_gcms, bmatgraph_from_fit)
bmatgraphse <- lapply(bmatgraphs, set_edgeatr_bmatgraph, timecoords = timecoords)
dev.off()

th <- 0.3
from.dist <- 10000
pdf(paste0("figs/bmatgraph_",th,"_",from.dist,"plots.pdf"), width = 10, height = 13)

par(mfrow = c(4,3))
mapply(plot_bmatgraph, network = bmatgraphse,  title = names(bmatgraphs), MoreArgs = list(th = th, from.dist = from.dist, curvature = TRUE,perm = permutations[[3]]))

dev.off()
############################################
# distances vs loglikes vs betas
############################################
bgr <- bmatgraphse$interim_10d_akima_cubic
bgr
gr <- bmatgraphse$CMIP5_EC.EARTH
E(gr)$strength
par(mfrow = c(4,3))
for(i in 1:length(bmatgraphse)){
# plot(E(bmatgraphse$interim_10d_akima_cubic)$distances,E(bmatgraphse$interim_10d_akima_cubic)$strength, 
  
  
  plot(E(bmatgraphse[[i]])$strength, E(bmatgraphse[[i]])$distances,
       col = "grey",
     xlab = "beta",
     ylab = "distance",
     main = abrev[[i]])
# points(E(bmatgraphse[[i]])$distances,E(bmatgraphse[[i]])$strength)
abline(h = 0.3)
abline(h = -0.3)

hist(E(bmatgraphse[[i]])$strength,
     breaks = seq(-35,35,0.2),
     freq = TRUE,
     xlab = "beta",
     xlim = c(-2,2),
     ylab = "frequency",
     ylim = c(0,400),
     main = abrev[[i]],
     add = TRUE)

}
dev.off()
E(gr)[E(gr)$strength > 15]
E(bgr)[E(bgr)$strength > 15]

par(mfrow = c(4,3))
for(i in 1:length(bmatgraphse)){
  # plot(E(bmatgraphse$interim_10d_akima_cubic)$distances,E(bmatgraphse$interim_10d_akima_cubic)$strength, 
  
  hist(E(bmatgraphse[[i]])$strength,
       breaks = seq(-35,35,0.2),
       freq = TRUE,
       xlab = "beta",
       xlim = c(-2,2),
       ylab = "frequency",
       ylim = c(0,400),
       main = abrev[[i]])
  
  par(new = T)
  plot(E(bmatgraphse[[i]])$strength, E(bmatgraphse[[i]])$distances,
       col = "grey",
       xlim = c(-2,2),
       ylim = c(0,20000),
       axes = F,
       xlab = NA,
       ylab = NA,
       main = abrev[[i]])
  axis(side = 4)
  mtext(side = 4, line = 2, 'distance',cex = 0.7)
  

  
}


i <- 4
par(mfrow = c(4,3))
for(i in 1:length(bmatgraphse)){
  big <- E(bmatgraphse[[i]])$distances>2000
  E(bmatgraphse[[i]])[big]$color <- "red"
  E(bmatgraphse[[i]])[!big]$color <- "grey"
  plot(E(bmatgraphse[[i]])$strength,-E(permstrengths[[i]])$strengths, 
       col = E(bmatgraphse[[i]])$color,
       xlim = c(-5,5),
       xlab = "beta",
       ylab = "loglik",
       main = abrev[[i]])
  # abline(h = 0.3)
  # abline(h = -0.3)
}


par(mfrow = c(4,3))
for(i in 1:length(bmatgraphse)){
  big <- E(bmatgraphse[[i]])$distances>2000
  E(bmatgraphse[[i]])[big]$color <- "red"
  E(bmatgraphse[[i]])[!big]$color <-"yellow"
  plot(E(bmatgraphse[[i]])$strength,-E(permstrengths[[i]])$strengths, 
       col = E(bmatgraphse[[i]])$color,
       xlim = c(-5,5),
       xlab = "beta",
       ylab = "loglik",
       main = abrev[[i]])
  # abline(h = 0.3)
  # abline(h = -0.3)
}

par(mfrow = c(4,3))
for(i in 1:length(bmatgraphse)){

  hist(E(bmatgraphse[[i]])$strength,
       breaks = seq(-35,35,0.2),
       freq = TRUE,
       xlab = "beta",
      xlim = c(-2,2),
       ylab = "frequency",
      ylim = c(0,400),
       main = abrev[[i]])
  # abline(h = 0.3)
  # abline(h = -0.3)
}
dev.off()

par(mfrow = c(4,3))
boxplot(E(bmatgraphse[[i]])$strength)
for(i in 1:length(bmatgraphse)){
  
  boxplot(E(bmatgraphse[[i]])$strength,
          ylim = c(-2,2),
       main = abrev[[i]])
  # abline(h = 0.3)
  # abline(h = -0.3)
}
dev.off()

#############################################################################
# How many steps are needed for CMIP5 models to go from inicial
# to end vertices for large links interim ? 
#############################################################################
medir <- ends(bmatgraphse$interim_10d_akima_cubic,E(bmatgraphse$interim_10d_akima_cubic)[distances > 10000])

distance_th <-  10000
beta_th <- 0
medirs <- lapply(bmatgraphse, function(x) ends(x,E(x)[distances > distance_th & abs(strength) > beta_th]))

#medirs$CMIP5_CanESM2
#shortest_paths(bmatgraphse$CMIP5_EC.EARTH, from = medir[,1], to = medir[,2], mode = "all", weights = NA, output = "epath")

testf <- function(x,y) {
apply(x, MARGIN = 1, function(z) shortest_paths(y, from = z[1], to = z[2], mode = "all", weights = NA, output = "epath"))}
testg <- function(x,y) {
apply(x, MARGIN = 1, function(z) distances(y, v = z[1], to = z[2], mode = "all", weights = NA))}

#testf(medirs$CMIP5_CanESM2,bmatgraphse$CMIP5_CNRM.CM5)
#testg(medirs$CMIP5_CanESM2,bmatgraphse$CMIP5_CNRM.CM5)
pathsmedirsrel <- sapply(medirs, function(z) sapply(bmatgraphse, function(w) sum(testg(z,w))/length(testg(z,w))))
pathsmedirs
pathsmedirs$CMIP5_CanESM2$epath

pathsmedirsreldata <- melt(pathsmedirsrel) 

lowcolor <- "white"
relplot <- ggplot(pathsmedirsreldata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low=lowcolor, high="darkgreen") +
  geom_text(aes(label = round(value,1))) +
  labs(x="Model 2", y="Model 1", title=paste0("Average pathlength of Model 1 to reach edges of\nlength > ",distance_th," and |beta| > ", beta_th," of Model 2")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

pathsmedirs<- sapply(medirs, function(z) sapply(bmatgraphse, function(w) sum(testg(z,w))))
pathsmedirs
pathsmedirs$CMIP5_CanESM2$epath

pathsmedirsdata <- melt(pathsmedirs) 


totalplot <- ggplot(pathsmedirsdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low=lowcolor, high="darkgreen") +
  geom_text(aes(label = round(value,1))) +
  labs(x="Model 2", y="Model 1", title=paste0("Total pathlength of Model 1 to reach edges of\nlength > ",distance_th," and |beta| > ", beta_th," of Model 2")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/strength_vs_dist/shortest_distances_",distance_th,"_",beta_th,".pdf")
plotname
pdf(plotname, width = 13, height = 5)
joep <- arrangeGrob(totalplot,relplot,nrow = 1)
grid.arrange(joep, newpage = FALSE)
dev.off()
################################################################################
# Plot tp / fp / fn
################################################################################
comparelist_tp <- sapply(selection_hc_gcms, function(y) sapply(selection_hc_gcms, function(x){ z <- compare(y,x, arcs = FALSE); return(z$tp)}))
tp_data <- melt(comparelist_tp)
comparelist_fp <- sapply(selection_hc_gcms, function(y) sapply(selection_hc_gcms, function(x){ z <- compare(y,x, arcs = FALSE); return(z$fp)}))
fp_data <- melt(comparelist_fp)
comparelist_fn <- sapply(selection_hc_gcms, function(y) sapply(selection_hc_gcms, function(x){ z <- compare(y,x, arcs = FALSE); return(z$fn)}))
fn_data <- melt(comparelist_fn)

ggplot(tp_data, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  geom_text(aes(label = round(value,0))) +
  scale_fill_gradient(limits = c(min(tp_data$value),800)) +
  # labs(x="data", y="models", title="optimum model evaluation") +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

# which of x (tweede element) also present in y (eerste element). Dus tweede element tekenen. 
comparelist_tp_arcs <- lapply(selection_hc_gcms, function(y) sapply(selection_hc_gcms, function(x){ z <- compare(y,x, arcs = TRUE); return(z$tp)}))
# fP: which in x (tweede element) not present in y (eerste element). Dus tweede element tekenen.
comparelist_fp_arcs <- lapply(selection_hc_gcms, function(y) sapply(selection_hc_gcms, function(x){ z <- compare(y,x, arcs = TRUE); return(z$fp)}))
# fn: which not in x (tweede element) but present in y (eerste element). Dus eerste element tekenen.
comparelist_fn_arcs <- lapply(selection_hc_gcms, function(y) sapply(selection_hc_gcms, function(x){ z <- compare(y,x, arcs = TRUE); return(z$fn)}))

shortnames
# To continue: Plots of tp/fn/fp edges networks with ggplot
# See paragraph below. 

# CHECKING CONDITIONS
# comparelist <- comparelist_fp_arcs
# permweights <- permweights
# xgraph <- "interim_10d_akima_cubic"
# ygraph <- "CMIP5_CanESM2"
# edgetype <- "fp"

all.equal(comparelist_fp_arcs,comparelist_fn_arcs)

cols4arcs2plot <- function(comparelist, permweights, xgraph, ygraph, edgetype = c("tp","fp","fn")) {
  
  if (edgetype == "fp"||edgetype == "tp") {
    permweight <- permweights[[ygraph]]
    edgelist <- as_edgelist(permweight)
    comparegraphlist <- comparelist[[xgraph]][[ygraph]]
  } else {permweight <- permweights[[xgraph]]
    edgelist <- as_edgelist(permweight)
    comparegraphlist <- comparelist[[ygraph]][[xgraph]]}
  
  
  indsame <- c()
  for (i in 1:nrow(comparegraphlist)){
    int <- intersect(which(edgelist[,1] == comparegraphlist[i,1]),
                     which(edgelist[,2] == comparegraphlist[i,2]))
    if(length(int)>0) {indsame[i] <- int} else {indsame[i]<- NA}
  }
  length(which(!is.na(indsame)))
  
  arcslist <- list()
  for(i in 1:nrow(edgelist))  {
    
    dists <- edge.attributes(permweight)$distances
    distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
    distsc <- as.character(distsv)
    
    node1 <- gridpoints[gridpoints == edgelist[i,1]]
    
    node2 <- gridpoints[gridpoints == edgelist[i,2]]
    
    arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                          c(p[node2,]$Var2, p[node2,]$Var1),
                          n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
    
    edge.ind <- ceiling(edge_attr(permweight,"distances",E(permweight)[i])*100 / max(dists))
    
    if (!i %in% indsame){edge.ind <- NA}

    if(node1!=node2){
      negs <- which(arc[,'lon']<0)}
    
    if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
    arc <- cbind(arc,edge.ind)
    arcslist[[i]]<- arc
  }
  
  naampjes <- character()
  for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
  names(arcslist) <- naampjes
  
  do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
    cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
  })) -> all_paths
  
  return(all_paths)
}


ggplot_edgecompare <- function(all_paths,edgelist,shortname) {ggplot() +
    geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
    geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = as.numeric(edge.ind))) +
    scale_color_gradient(low = col.1, high = col.2,limits = c(0,100),na.value = rgb(0, 0, 1, alpha = 0), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
    ggtitle(paste0(nrow(edgelist)," ",shortname)) +
    theme(plot.title = element_text(hjust = 0.5))}
##################################################################
# plot single case
##################################################################
if (edgetype == "fp"){title <- paste0("in ", ygraph," not in ",xgraph)}
if (edgetype == "fn"){title <- paste0("not in ", ygraph," but in ",xgraph)}
if (edgetype == "tp"){title <- paste0("in ", ygraph," and in ", xgraph)}
ggplot_edgecompare(all_paths,edgelist,title)

################################################################
# Plot all cases
#################################################################
shortnames
xgraph <- "CMIP5_HadGEM2.ES"
edgetype <- "fn"
assign("comparelist", eval(get(paste0("comparelist_",edgetype,"_arcs"))))
all_paths_list_compare <- lapply(X = shortnames[!grepl(xgraph,unlist(shortnames))], FUN = cols4arcs2plot, comparelist = comparelist, permweights = permweights,xgraph = xgraph,edgetype = edgetype)
if (edgetype == "fp"){titles <- lapply(shortnames[!grepl(xgraph,unlist(shortnames))],function(ygraph) paste0("in ", ygraph,"\n not in ",xgraph))}
if (edgetype == "fn"){titles <- lapply(shortnames[!grepl(xgraph,unlist(shortnames))],function(ygraph) paste0("not in ", ygraph,"\n but in ",xgraph))}
if (edgetype == "tp"){titles <- lapply(shortnames[!grepl(xgraph,unlist(shortnames))],function(ygraph) paste0("in ", ygraph,"\n and in ", xgraph))}

ggplotsedgecompare <- mapply(ggplot_edgecompare, all_paths_list_compare, edgelists[shortnames[!grepl(xgraph,unlist(shortnames))]], titles, SIMPLIFY = FALSE)

pdf(paste0("figs/tp_fp_fn/",xgraph,"_",edgetype,".pdf"), width = 15, height = 20)
do.call("grid.arrange",ggplotsedgecompare)
dev.off()

#############################################################################
# How many steps are needed for CMIP5 models to go from inicial
# to end vertices for false negatives with respect to interim ? 
# NOT YET FINISHED
#############################################################################
# mot in y graph but in interim graph 
comparelist_fn_arcs[[ygraph]][["interim_10d_akima_cubic"]]

medir <- ends(bmatgraphse$interim_10d_akima_cubic,E(bmatgraphse$interim_10d_akima_cubic)[distances > 10000])

distance_th <-  10000
beta_th <- 0
medirs <- lapply(bmatgraphse, function(ygraph) ends(comparelist_fn_arcs[[ygraph]][["interim_10d_akima_cubic"]],E(x)[distances > distance_th & abs(strength) > beta_th]))

#medirs$CMIP5_CanESM2
#shortest_paths(bmatgraphse$CMIP5_EC.EARTH, from = medir[,1], to = medir[,2], mode = "all", weights = NA, output = "epath")

testf <- function(x,y) {
  apply(x, MARGIN = 1, function(z) shortest_paths(y, from = z[1], to = z[2], mode = "all", weights = NA, output = "epath"))}
testg <- function(x,y) {
  apply(x, MARGIN = 1, function(z) distances(y, v = z[1], to = z[2], mode = "all", weights = NA))}

#testf(medirs$CMIP5_CanESM2,bmatgraphse$CMIP5_CNRM.CM5)
#testg(medirs$CMIP5_CanESM2,bmatgraphse$CMIP5_CNRM.CM5)
pathsmedirsrel <- sapply(medirs, function(z) sapply(bmatgraphse, function(w) sum(testg(z,w))/length(testg(z,w))))
pathsmedirs
pathsmedirs$CMIP5_CanESM2$epath

pathsmedirsreldata <- melt(pathsmedirsrel) 

lowcolor <- "white"
relplot <- ggplot(pathsmedirsreldata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low=lowcolor, high="darkgreen") +
  geom_text(aes(label = round(value,1))) +
  labs(x="Model 2", y="Model 1", title=paste0("Average pathlength of Model 1 to reach edges of\nlength > ",distance_th," and |beta| > ", beta_th," of Model 2")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

pathsmedirs<- sapply(medirs, function(z) sapply(bmatgraphse, function(w) sum(testg(z,w))))
pathsmedirs
pathsmedirs$CMIP5_CanESM2$epath

pathsmedirsdata <- melt(pathsmedirs) 


totalplot <- ggplot(pathsmedirsdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low=lowcolor, high="darkgreen") +
  geom_text(aes(label = round(value,1))) +
  labs(x="Model 2", y="Model 1", title=paste0("Total pathlength of Model 1 to reach edges of\nlength > ",distance_th," and |beta| > ", beta_th," of Model 2")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))

#################################################################################
#  shd
#################################################################################
shd(selection_hc_gcms$CMIP5_CanESM2,selection_hc_gcms$CMIP5_CNRM.CM5)


#################################################################################
# Plots of long strong edges networks with ggplot
#################################################################################
data.network <- TimeCoordsAnom_from_Grid_rms(cubicsets$interim_10d_akima_cubic, rms = TRUE)
x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)
# x <- (x+360)%%360
p <- expand.grid(y, x)[2:1][permutations[[3]],]

modelsize <- 18
selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
selection_hc_gcms$CMIP5_CanESM2

permstrengths <- lapply(selection_hc_gcms, bn_to_igraph.strengths, perm = permutations[[3]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))
permdists <- lapply(permstrengths, igraph.distances, perm = permutations[[3]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))
permweights <- lapply(permdists, igraph.weights, type = "bn", fromdist = 0)


edgelists <- lapply(bmatgraphse, as_edgelist) 
row.names(p) <- names(V(bmatgraphse[[10]]))
gridpoints <- names(V(bmatgraphse[[10]]))


col.1 <- adjustcolor("light blue", alpha=0.8)
col.2 <- adjustcolor("navy blue", alpha=0.8)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

col.3 <- adjustcolor("orange", alpha=0.8)
col.4 <- adjustcolor("red", alpha=0.8)
edge.pal2 <- colorRampPalette(c(col.3, col.4), alpha = TRUE)
edge.col2 <- edge.pal2(100)

edgelist <- edgelists$CMIP5_NorESM1.M
permweight <- bmatgraphse$CMIP5_NorESM1.M
cols4plot <- function(permweight,edgelist, th = 0.3,from.dist = 5000){
  arcslist <- list()

  for(i in 1:nrow(edgelist))  {
    
    graphBN <- permweight
    redtype <- "BN"
    
    assign("graph", eval(parse(text = paste0("graph",redtype))))
    
    dists <- edge.attributes(graphBN)$distances
    distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
    distsc <- as.character(distsv)
    
    node1 <- gridpoints[gridpoints == edgelist[i,1]]
    
    node2 <- gridpoints[gridpoints == edgelist[i,2]]
    
    arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                          c(p[node2,]$Var2, p[node2,]$Var1),
                          n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
    
    edge.ind <- ceiling(edge_attr(graph,"distances",E(graph)[i])*100 / max(dists))
    col <- as.character(edge.col[edge.ind])
    if(abs(edge_attr(graph,"strength",E(graph)[i])) < th) {edge.ind <- NA}
    if(abs(edge_attr(graph,"distances",E(graph)[i])) < from.dist) {edge.ind <- NA}
    if(node1!=node2){
    negs <- which(arc[,'lon']<0)}
    
    if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
    arc <- cbind(arc,edge.ind)
    # lines(arc, col=edge.col[edge.ind], lwd=edge.ind) interesting for edgewidth!
    # points(x = arc[1,1], y = arc[1,2],col = "pink", pch = 19)
    arcslist[[i]]<- arc
    
  }

  naampjes <- character()
  for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
  names(arcslist) <- naampjes
  
  do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
    cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
  })) -> all_paths
  
  return(all_paths)
}

# all_paths_test <- dists4plot(edgelists$CMIP5_CanESM2,permweight = permweights$CMIP5_CanESM2)
th <- 0.30
from.dist <- 3000
all_paths_list <- mapply(cols4plot, edgelist = edgelists, permweight = bmatgraphse, MoreArgs = list(th = th, from.dist = from.dist), SIMPLIFY = FALSE)

cols4plot(edgelist = edgelists$CMIP5_CanESM2,permweight = bmatgraphse$CMIP5_CanESM2, th = th, from.dist = from.dist)
mapworld2 <- map("world2", fill = TRUE)
world_pol_df <- tidy(mapworld2, IDs = "region")


ggplotfunction <- function(all_paths,edgelist,shortname) {ggplot() +
    geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
    geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = as.numeric(edge.ind))) +
    scale_color_gradient(low = col.1, high = col.2,limits = c(0,100),na.value = rgb(0, 0, 1, alpha = 0), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
    ggtitle(paste0(redtype,": ",nrow(edgelist)," ",shortname)) +
    theme(plot.title = element_text(hjust = 0.5))}

ggplotfunction(all_paths_list$CMIP5_CanESM2,edgelists$CMIP5_CanESM2,abrev[[1]])
ggplots <- mapply(ggplotfunction, all_paths_list, edgelists, abrev, SIMPLIFY = FALSE)
# CHECKEN WELKE ZIJN RECHT
E(bmatgraphse)


pdf(paste0("figs/ggplots_",th,"_",from.dist,".pdf"), width = 15, height = 20)
do.call("grid.arrange",ggplots)
dev.off()




graph1800 <- graph_from_cov(cubicsets$CMIP5_CanESM2_r1i1p1, bmatrix, th = 0)
graph1800$graph <- permute(graph1800$graph,permutations[[3]])
V(graph1800$graph)

dis1800 <- str_vs_dis_2(graph1800$graph,data.dag = data_gcms$CMIP5_CanESM2_r1i1p1, perm = permutations[[1]])
edge.attributes(dis1800)

dis1800 <- igraph.distances(graph1800$graph, data.igraph = data_gcms$CMIP5_CanESM2_r1i1p1,  perm = permutations[[1]])

diswei1800 <- igraph.weights(dis1800, type = "cn", perm = permutations[[1]])
diswei1800 <- cn.strengths_2(diswei1800, data.igraph = data_gcms$CMIP5_CanESM2_r1i1p1, perm = permutations[[3]])
edge.attributes(diswei1800)

g <- diswei1800
ind <- which(E(g)$weights==0)
E(g)$weights[ind] <- 0.000000001
commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                      weights = 1/E(g)$weights
)
ncom <- length(levels(factor(commu_bet$membership)))
mem <- membership(commu_bet)
memClim <- quantity2clim(mem, "membership", cubicsets$CMIP5_CanESM2_r1i1p1, backperm = backpermutations[[3]])
if(ncom <= 9){
  colRainbow <- brewer.pal(ncom,"Set1")
} else {
  colRainbow <- brewer.pal(9,"Set1")
  colRainbow<- colorRampPalette(colRainbow)(ncom)
}
if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
plotscom <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                        set.min = NULL, set.max = NULL, 
                        lonCenter = 180, 
                        regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                        # main = paste0("Comunities BN ",whichsizesBN[i]),
                        colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                        lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
)
plotscom
closeness1800 <- closeness(graph = g, weights = E(g)$weights, normalized = FALSE)

climcloseness1800 <- quantity2clim(closeness1800, what = "closeness",ref.grid = cubicsets$CMIP5_CanESM2_r1i1p1, backperm = backpermutations[[3]])

cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(80)
attr(climcloseness1800$Data, "climatology:fun")


plot <- spatialPlot(climcloseness1800, 
                    # main = list(paste0("CN:",whichsizesCN[[i]]), cex = 0.5),
                    backdrop.theme = "coastline", 
                    region = TRUE, col.regions= cb,
                    # set.max = 20000, 
                    colorkey = list(width = 0.6, lables = list(cex = 0.5)))
plot

cuts <- c(4,8,12,16,20,24)
plots<- list()
for(j in 1:length(cuts)){
  i <- cuts[[j]]
  commu_bet12 <- cut_at(commu_bet,i)
  
  
  ncom <- length(levels(factor(commu_bet12)))
  mem <- commu_bet12
  memClim <- quantity2clim(mem, "membership", cubicsets$CMIP5_CanESM2_r1i1p1, backperm = backpermutations[[3]])
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plots[[j]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                            set.min = NULL, set.max = NULL, 
                            lonCenter = 180, 
                            regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                            main = paste0("Comunities BN ",length(E(g))),
                            colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                            lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
  
}
plots
do.call(grid.arrange,c(plots))




betweenness1800 <- betweenness(g, weights = 1/E(g)$weights)
logbetweenness1800 <- log(1+betweenness1800)
climbetw1800 <- quantity2clim(logbetweenness1800, what = "betweenness", ref.grid = cubicsets$CMIP5_CanESM2_r1i1p1, backperm = backpermutations[[3]])


plot <- spatialPlot(climbetw1800, 
                    # main = list(paste0("CN:",whichsizesCN[[i]]), cex = 0.5),
                    backdrop.theme = "coastline", 
                    lonCenter = 180,
                    col = "Reds",
                    colorkey = list(width = 0.6, lables = list(cex = 0.5))
)
plot


plot_long_distances(data.dag = TimeCoordsAnom_from_Grid_rms(cubicsets$CMIP5_CanESM2_r1i1p1,rms = TRUE),dag = graph1800$graph, minimdist = 2000,perm = permutations[[3]])
plot_long_strong_distances(data.dag = TimeCoordsAnom_from_Grid_rms(gridused,rms = TRUE),dag = graph1800$graph,minimdist = 2000, perm = permutations[[3]])
plotS.spatial.igraph(network = g, data.network = TimeCoordsAnom_from_Grid_rms(cubicsets$CMIP5_CanESM2_r1i1p1, rms = TRUE), type = "cn", th = 0, th.type = "zoom", from.dist = 2000, by.dist = NULL, perm = permutations[[3]], title = "NA", remove = TRUE, shift = TRUE, curvature = FALSE)

plotS.spatial.igraph

edge.attributes(g)
