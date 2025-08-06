##########################################################################
# Topology evaluation 
##########################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
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

load("../Data/tas_ncep_10d.rda")
load("data/tas_JRA55_10d_akima_cubic.rda")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_left_10d_akima_cubic.rda")
load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda")

# library(magrittr)
# library(reshape2)
# library(gridExtra)
# library(sparsebn)
# library(RColorBrewer)
# library(gaussDiff)
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
cubicsets6 <- tas_historical_cmip6_10d_akima_cubic_corrected

# With CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,cubicsets6,listinterim,listjra55,listncep)
names(cubicsets)
# Without CMIP6
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,cubicsetsearth,listinterim,listjra55,listncep)
names(cubicsets)
# Without earth
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsets5left,listinterim,listjra55,listncep)
names(cubicsets)

namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
namescubics5left <- names(cubicsets5left)
namescubics6 <- gsub(names(cubicsets6), pattern = "_historical",replacement ="") 
namescubicsearth <- names(cubicsetsearth)
namescubics <- names(cubicsets)
shortnames <- #gsub(
  gsub(names(cubicsets), pattern = "_r1i1p1", replacement = "")
#,
#pattern = "_r12i1p1", replacement = "")
shortnames2 <- names(cubicsets)
abrev <-  gsub("_akima_cubic",gsub("CMIP5_",shortnames,replacement = ""),replacement = "")
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


#################################################################################
# Plots of networks with ggplot
#################################################################################
permk <- 3

# Anomalies and standarization over seasons (networks are learned with this set)
data_gcms_out <- lapply(cubicsets,function(x) TimeCoordsAnom_from_Grid_rms(x, rms = TRUE))
data_gcms_out_df <- lapply(data_gcms_out, function(x) as.data.frame(x))
data_gcms_anom_out <- lapply(data_gcms_out, function(x) mat2Dto3Darray(x,attr(x,"Xcoords"), attr(x,"Ycoords")))
grid_gcms_anom <- mapply(function(x,y) {x$Data <- y 
return(x)}, 
x = cubicsets, 
y = data_gcms_anom_out,SIMPLIFY = FALSE)
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

it <- "1700_1800"
hc_gcms5 <- lapply(paste0(namescubics5), loadIterations, permused = permk, it = it)
hc_gcms5extra <- lapply(paste0("CMIP5_extra/",namescubics5extra), loadIterations, permused = permk, it = it)
hc_gcms5left<- lapply(paste0("CMIP5_left/",namescubics5left), loadIterations, permused = permk, it = it)
hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk, it = it)
hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk, it = it)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_jra55 <- lapply(c("JRA55_10d_akima_cubic"), loadIterations, permused = permk, it = it) 
hc_ncep <- lapply(c("../Data/interim_struct/hciterations"), loadIterations, permused = permk,ncep = TRUE, it = it) 
# With cmip6
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_gcms6,hc_interim,hc_jra55,hc_ncep)
# Without cmip6
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_gcmsearth,hc_interim,hc_jra55,hc_ncep)
length(hc_gcms)
# Without cmip6 without earth
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcms5left,hc_interim,hc_jra55,hc_ncep)
length(hc_gcms)
# names(hc_gcms) <- shortnames
names(hc_gcms) <- shortnames2


# Choose between 'own optimums' or constant magnitude
# modelsize <- own_optimums
modelsize <- 18
selection_hc_gcms <-  mapply(function(x,y) x[[grep((y),x)]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)

#####################################################################################
# distance matrixes maken en clusteren
#####################################################################################
graphlist <- list.files(paste0("results/struct_edgeattributes/perm",permk), full.names = T)
graphlist_names <- list.files(paste0("results/struct_edgeattributes/perm",permk))
graphlist_names <- gsub(".rda", "", graphlist_names)

graphs <- lapply(graphlist, function(x){get(load(x))})
names(graphs) <- graphlist_names  
graphsizes <- list()

for(i in 1:length(graphs)){
  graph <- graphs[[i]]
  assign(graphlist_names[[i]],graphs[[i]])
  graphsizes[[i]] <- sapply(graph,function(x)lapply(x,function(x)length(E(x))))
}


#############################################################################
# How many steps are needed for CMIP5 models to go from inicial
# to end vertices for large links interim ? 
#############################################################################
# medir <- ends(bmatgraphse$interim_10d_akima_cubic,E(bmatgraphse$interim_10d_akima_cubic)[distances > 10000])

distance_sign <- "greater"
distance_th <-  10000
beta_sign <- "greater"
beta_th <- 0


cmip5list <- c(graphs_CMIP5_extra_HC_1700_1800_toPN_dist_perm3,graphs_CMIP5_HC_1700_1800_toPN_dist_perm3,graphs_CMIP5_left_HC_1700_1800_toPN_dist_perm3,graphs_interim_10d_akima_cubic_HC_1700_1800_toPN_dist_perm3,graphs_JRA55_10d_akima_cubic_HC_1700_1800_toPN_dist_perm3,graphs_ncep_10d_HC_1700_1800_toPN_dist_perm3)
cmip5graphs <- unlist(cmip5list, recursive = FALSE, use.names = FALSE)
names(cmip5graphs) <- names(cmip5list)

if(distance_sign == "greater" & beta_sign == "greater"){
  medirs <- lapply(cmip5graphs, function(x) ends(x,E(x)[distances > distance_th & abs(weight) > beta_th]))
} else if (distance_sign == "less" & beta_sign == "greater") {
  medirs <- lapply(cmip5graphs, function(x) ends(x,E(x)[distances < distance_th & abs(weight) > beta_th]))
} else if (distance_sign == "less" & beta_sign == "less") {
  medirs <- lapply(cmip5graphs, function(x) ends(x,E(x)[distances < distance_th & abs(weight) < beta_th]))
} else if (distance_sign == "greater" & beta_sign == "less") {
  medirs <- lapply(cmip5graphs, function(x) ends(x,E(x)[distances > distance_th & abs(weight) < beta_th]))
}


# Distance calculation without weights
testg <- function(x,y) {
  apply(x, MARGIN = 1, function(z) distances(y, v = z[1], to = z[2], mode = "all", weights = NA))}

distancesg <- lapply(medirs, function(z) sapply(cmip5graphs,function(w) testg(z,w)))

assign(paste0("arcs_",distance_sign,distance_th,"_",beta_sign,beta_th),distancesg)
save(list = paste0("arcs_",distance_sign,distance_th,"_",beta_sign,beta_th),file = paste0("results/distance_weights_arcs/arcs_",distance_sign,distance_th,"_",beta_sign,beta_th,".rda"))

# # Distance calculation with weights sum
# testg <- function(x,y) {
#   apply(x, MARGIN = 1, function(z) distances(y, v = z[1], to = z[2], mode = "all", weights = abs(edge_attr("weight", graph = y))))}
# 
# distancesg <- lapply(medirs, function(z) sapply(cmip5graphs,function(w) testg(z,w)))
# 
# # Distance caluclation with weights multiplication. 
# testglog <- function(x,y) {
#   apply(x, MARGIN = 1, function(z) distances(y, v = z[1], to = z[2], mode = "all", weights = -log(abs(edge_attr("weight", graph = y)))))}
# 
# distancesglog <- lapply(medirs, function(z) sapply(cmip5graphs,function(w) testglog(z,w)))
# distancesg <- lapply(distancesglog, function(x)exp(x*-1))

distancesg <-get(load(paste0("results/distance_weights_arcs/arcs_",distance_sign,distance_th,"_",beta_sign,beta_th,".rda")))
dimens <- lapply(distancesg, dim)

pathsmedirs<- sapply(distancesg, function(z) if (!is.null(dim(z))) {colSums(z)} else {rep(NA,length(z))})
# pathsmedirs<- sapply(distancesg, function(z) apply(z,MARGIN = 2,FUN = sum))
# without earth;
indearth <- grep("EC.EARTH_r1i|EC.EARTH_r2i",rownames(pathsmedirs))
pathsmedirs <- pathsmedirs[-indearth,-indearth]

pathsmedirsrel <- t(t(pathsmedirs)/diag(pathsmedirs))
pathsmedirsrel[which(pathsmedirsrel == Inf)] <- NA


###############################################################################
# correlation hellinger distance and geographical distances
###############################################################################
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
helcoef_eq <- hellinger_coefficients[rownames(pathsmedirsrel),colnames(pathsmedirsrel)]
# row: coonnexion between hellinger distance and how you reaching others
diag(apply(-log(helcoef_eq), MARGIN = c(1), function(x){apply(pathsmedirsrel, MARGIN = c(1), function(y) cor(x = x,y = y, use = "pairwise.complete.obs", method = "pearson") )}))
# column: conexion between hellinger distance and how others reaching you.
diag(cor(-log(helcoef_eq),pathsmedirsrel,use = "pairwise.complete.obs")) # column

model <- "CMIP5_HadGEM2.ES_r1i1p1"
model <- "CMIP5_bcc.csm1.1.m_r1i1p1"
model <- "CMIP5_EC.EARTH_r1i1p1"
model <- "CMIP5_BNU.ESM_r1i1p1" 
model <- "CMIP5_ACCESS1.0_r1i1p1"
model <- "interim_10d_akima_cubic"
model <- "JRA55_10d_akima_cubic"
model <- "ncep_10d"

# how others reach you; without auto correaltion
cor(-log(helcoef_eq)[model,-which(colnames(helcoef_eq) == model)],pathsmedirsrel[-which(rownames(pathsmedirsrel) == model),model],use = "pairwise.complete.obs", method = "pearson")
# how others reach you; with auto correlation
cor(-log(helcoef_eq)[model,],pathsmedirsrel[,model],use = "pairwise.complete.obs")
plot(-log(helcoef_eq)[model,],pathsmedirsrel[,model])
# how you reach others: without auto correlation
cor(-log(helcoef_eq)[model,-which(colnames(helcoef_eq) == model)],pathsmedirsrel[model,-which(colnames(pathsmedirsrel) == model)],use = "pairwise.complete.obs")
plot(-log(helcoef_eq)[model,-which(colnames(helcoef_eq) == model)],pathsmedirsrel[model,-which(colnames(pathsmedirsrel) == model)])
# how you reach others: with auto correlation
cor(-log(helcoef_eq)[model,],pathsmedirsrel[model,],use = "pairwise.complete.obs")
plot(-log(helcoef_eq)[model,],pathsmedirsrel[model,])

plot(-log(helcoef_eq[model,-which(colnames(helcoef_eq) == model)]), pathsmedirsrel[model,-which(colnames(pathsmedirsrel) == model)])
plot(-log(helcoef_eq[model,-which(colnames(helcoef_eq) == model)]), pathsmedirsrel[-which(colnames(pathsmedirsrel) == model),model])


####################################################################################
# Plot and cluster distance matrices
####################################################################################
institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor","ACCESS","bcc","BNU","CCSM4","CESM1","CMCC","CSIRO","inmcm4","MRI")

# with cmip6
institutions <- rep(institutions1,each =2)
CMIPS <- rep(c("CMIP5","CMIP6"),length(institutions1))
# without cmip6
institutions <- rep(institutions1,each =1)
CMIPS <- rep(c("CMIP5"),length(institutions1))
combinations <- cbind(CMIPS,institutions)
# without earth
combinations


namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(pathsmedirs))))
namesort <- c(namesort,30,31,32)
Order <- "Hellinger"
length(namesort)
# OLD method:
#pathsmedirs <- sapply(medirs, function(z) sapply(bmatgraphse, function(w) sum(testg(z,w))))
#pathsmedirsrel <- sapply(medirs, function(z) sapply(bmatgraphse, function(w) sum(testg(z,w))/length(testg(z,w))))
if(Order == "Hellinger"){
  pathsmedirsreldata <- melt(pathsmedirsrel[mclust$labels[mclust$order],mclust$labels[mclust$order]]) 
  pathsmedirsdata <- melt(pathsmedirs[mclust$labels[mclust$order],mclust$labels[mclust$order]]) 
} else if(Order == TRUE){
  pathsmedirsreldata <- melt(pathsmedirsrel[namesort,namesort]) 
  pathsmedirsdata <- melt(pathsmedirs[namesort,namesort]) 
} else {  pathsmedirsreldata <- melt(pathsmedirsrel) 
pathsmedirsdata <- melt(pathsmedirs) }

lowcolor <- "white"
relplot <- ggplot(pathsmedirsreldata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low=lowcolor, high="darkgreen",limits = c(2.5,4.5),midpoint = 2) +
  geom_text(aes(label = round(value,1))) +
  labs(x="Model 2", y="Model 1", title=paste0("Average pathlength of Model 1 to reach edges of\nlength > ",distance_th," and |rho| > ", beta_th," of Model 2")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")



totalplot <- ggplot(pathsmedirsdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient2(low=lowcolor, high="darkgreen") +
  geom_text(aes(label = round(value,1))) +
  labs(x="Model 2", y="Model 1", title=paste0("Total pathlength of Model 1 to reach edges of\nlength > ",distance_th," and |rho| > ", beta_th," of Model 2")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17))+ 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/strength_vs_dist/shortest_hellingercluster_distances_",distance_th,"_",beta_th,"_cmip5_PNs_left_extra.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/strength_vs_dist/shortest_distances_",distance_th,"_",beta_th,"_cmip5_PNs_left_extra.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/strength_vs_dist/shortest_weighted_distances_",distance_th,"_",beta_th,"_cmip5_PNs_left_extra.pdf")
plotname
pdf(plotname, width = 20, height = 30)
joep <- arrangeGrob(totalplot,relplot,nrow = 2)
grid.arrange(joep, newpage = FALSE)
dev.off()

##################################################################
# Partitional clustering
##################################################################
library("cluster")
library("factoextra")

clusterselection <- pathsmedirsrel
clusterselection[which(clusterselection == Inf)] <- 10

earthgrep <- unlist(apply(matrix(c("CMIP5_EC.EARTH_r1i1p1_historical","CMIP5_EC.EARTH_r2i1p1_historical","CMIP5_EC.EARTH_r1i1p1","CMIP5_EC.EARTH_r2i1p1","CMIP5_EC.EARTH_r12i1p1_historical"),nrow = 5, ncol = 1),MARGIN = 1,function(x) grep(x,rownames(clusterselection))))
pathsmedirs_without_earth <- clusterselection[-earthgrep,-earthgrep]
namesorthc <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(pathsmedirs_without_earth))))
reangrep <- unlist(apply(matrix(c("interim","ncep","JRA55"),nrow = 3, ncol = 1),MARGIN = 1,function(x) grep(x,rownames(pathsmedirs_without_earth))))

selecteddists <- pathsmedirs_without_earth[c(namesorthc,reangrep),c(namesorthc,reangrep)]
which(rownames(selecteddists) == "CMIP5_MIROC.ESM_r1i1p1.CMIP5_MIROC.ESM_3_1700_1800i")
which(colnames(selecteddists) == "CMIP5_MIROC.ESM_r1i1p1.CMIP5_MIROC.ESM_3_1700_1800i")
selecteddists <- selecteddists[-14,-14]
nrow(selecteddists)

klsort <- sort(selecteddists,index.return = TRUE)
# select optimal amount of clusters
fviz_nbclust(t(selecteddists), 
             FUNcluster = cluster::pam,
             method = "silhouette", 
             k.max =25 ,medoids = NULL, diss = t(selecteddists))
fviz_nbclust(t(selecteddists),FUNcluster = cluster::pam,method = "gap_stat",k.max = 10, diss = TRUE)
# Use PAM to perform k-medoids
pamcluster <- pam(t(selecteddists), k=3, diss = t(selecteddists),
                  # medoids = c(30,31), 
                  do.swap = FALSE)

# info
pamcluster$silinfo
names(selecteddists)
pamcluster$medoids
pamcluster$id.med
# visualize k-medoids
par(mar = c(5,15,3,4)+0.1)
plot(silhouette(pamcluster),max.strlen = 200, nmax.lab =40,cex = 0.7)
names(pamcluster)
fviz_cluster(pamcluster$clusinfo, repel = TRUE, axes = c(1,2),show.clust.cent = TRUE)
fviz_cluster(list(cluster = pamcluster$clustering, data = selecteddists),repel = TRUE)
#################################################################################
# Partitional Clustering
#################################################################################
dev.off()
library(proxy)
#selecteddists <- pathsmedirsrel
mclust <-hclust(dist(t(selecteddists),method = "minkowski", p = 2),method = "complete")
mclust <-hclust(dist(selecteddists,method = "euclidean"))
plot(mclust)
rect.hclust(mclust,k=11)
mclust$labels[mclust$order]


###############################################################################
# distances weights
###############################################################################
edge.attributes(cmip5graphs$CMIP5_GFDL.CM3_r1i1p1.CMIP5_GFDL.CM3_r1i1p1_3_1700_1800i)

########################################################
# Visualize distances example set 
########################################################
library(RColorBrewer)
plotname <- paste0("figs/exampleclustering/shortest_distances_ex_",distance_th,"_",beta_th,"_ACCESS_BNU_CMCC.CMS_interim.png")
png(plotname,units = "in", width = 9, height = 7,  res =180)
plotname <- paste0("figs/exampleclustering/shortest_distances_ex_",distance_th,"_",beta_th,"_ACCESS_BNU_CMCC.CMS_interim.pdf")
pdf(plotname, width = 9, height =7)
# n <-7

rownames(pathsmedirsrel)
ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")
ex.longdata <- melt(pathsmedirsrel[ex.subset,ex.subset])


library(RColorBrewer)

# n <- 5
# data.vec <- ex.longdata$value
 quant = seq(2,5,0.5)

quant.discreter <- function(n,data.vec,quant = NULL){
  if (is.null(quant)){br.lims <- quantile(data.vec, seq(0,1,1/n))} else {br.lims <- quant}
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

ex.longdata$d.value <- quant.discreter(n,ex.longdata$value, quant = seq(0,5,1))
                                                                        

lowcolor <- "white"
ex.relplot <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value),interpolate = FALSE) +
  scale_fill_gradientn(aesthetics = c("fill"), colours = brewer.pal(length(quant),"Greens"), guide = "legend")+
  #scale_fill_gradient2(low=lowcolor, high="darkgreen",limits = c(2.5,4.5),midpoint = 2,na.value = "white") +
  geom_text(aes(label = round(value,1)),size = 7) +
  labs(x="Model 2", y="Model 1", title=paste0("Average pathlength of Model 1 to reach edges of\nlength > ",distance_th," and |rho| > ", beta_th," of Model 2")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=17, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=17),
        plot.title=element_text(size=17)) + 
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

ex.relplot
dev.off()
# ex.longdata$d.value <- quant.discreter(n,ex.longdata$value)
# br.lims.var <- quantile(ex.longdata$value, seq(0,1,1/n))
# 
# b <- ggplot(ex.longdata, aes(x = Var2, y = Var1)) + 
#   geom_raster(aes(fill=d.value),interpolate = FALSE) + 
#   geom_text(aes(label = round(-value/10^4,0)), size = 7) +
#   scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Reds")), guide = "legend",name = "log P(data|model)\nx -10^4",n.breaks =n,labels = round(-br.lims.var[2:(n+1)]/10^4,0))+
#   labs(x="data", y="models", title=paste0(modelsize*100," model-dataset evaluation")) +
#   theme_bw() + 
#   theme(axis.text.x=element_text(size=15, angle=270, vjust=1, hjust = 1),
#         axis.text.y=element_text(size=15),
#         plot.title=element_text(size=12)) + 
#   scale_x_discrete(position = "top")
# 
