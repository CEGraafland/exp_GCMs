##########################################################################
# Topology evaluation 
##########################################################################
rm(list = ls())
library(igraph)
library(maps)
library(geosphere)
library(broom)
library(ggplot2)
library(visualizeR)
library(bnlearn)
library(gridExtra)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/backpermutations.rda")
load("data/tas_historical_10d_akima_cubic_corrected.rda")
load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda")
load("data/tas_historical_cmip6_10d_akima_cubic_corrected.rda")
load("data/tas_interim_10d_akima_cubic.rda")
load("data/tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic.rda")
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
listinterim <- list(tas_interim_10d_akima_cubic)
names(listinterim) <- "interim_10d_akima_cubic"

cubicsets5 <- tas_historical_10d_akima_cubic_corrected
cubicsets5extra <- tas_historical_cmip5_extra_10d_akima_cubic
cubicsetsearth <- tas_historical_cmip5_earth_r1r2r12_10d_akima_cubic
cubicsets6 <- tas_historical_cmip6_10d_akima_cubic_corrected
cubicsets <- c(cubicsets5,cubicsets5extra,cubicsetsearth,cubicsets6,listinterim)

namescubics5 <- names(cubicsets5)
namescubics5extra <- names(cubicsets5extra)
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
# Make al datasets standardized
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


hc_gcms5 <- lapply(paste0(namescubics5), loadIterations, permused = permk)
hc_gcms5extra <- lapply(paste0("CMIP5_extra/",namescubics5extra), loadIterations, permused = permk)
hc_gcmsearth <- lapply(paste0("CMIP5_EARTH_ONLYHIST/",namescubicsearth), loadIterations, permused = permk)
hc_gcms6 <- lapply(paste0("CMIP6/",namescubics6), loadIterations, permused = permk)
hc_interim <- lapply(c("interim_10d_akima_cubic"), loadIterations, permused = permk) 
hc_gcms <- c(hc_gcms5,hc_gcms5extra,hc_gcmsearth,hc_gcms6,hc_interim)
# names(hc_gcms) <- shortnames
names(hc_gcms) <- shortnames2

modelsize <- 18
selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)
# Make fits
# unscaled fit (uses data_gcms_out)
selection_fits_gcms <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_out_df, SIMPLIFY = FALSE)
# scaled fit (uses data_gcms_anom_scaled)
selection_fits_gcms_scaled <- mapply(function (x,y) bn.fit(x = x, data = y), x = selection_hc_gcms, y = data_gcms_anom_scaled, SIMPLIFY = FALSE)
# sel
names(selection_fits_gcms_scaled)

data.network <- TimeCoordsAnom_from_Grid_rms(cubicsets$interim_10d_akima_cubic, rms = TRUE)
x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)
# x <- (x+360)%%360
p <- expand.grid(y, x)[2:1][permutations[[3]],]

# modelsize <- 18
# selection_hc_gcms <-  mapply(function(x,y) x[[y]], x = hc_gcms, y = modelsize, SIMPLIFY = FALSE)

permstrengths <- lapply(selection_hc_gcms, bn_to_igraph.strengths, perm = permutations[[3]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))
permdists <- lapply(permstrengths, igraph.distances, perm = permutations[[3]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE))
permweights <- lapply(permdists, igraph.weights, type = "bn", fromdist = 0)
# opslaan??

edgelists <- lapply(permweights, as_edgelist) 
row.names(p) <- names(V(permweights$CMIP5_CanESM2))
gridpoints <- names(V(permweights$CMIP5_CanESM2))

# row.names(p) <- names(V(graph))
# gridpoints <- names(V(graph))


col.1 <- adjustcolor("light blue", alpha=0.8)
col.2 <- adjustcolor("navy blue", alpha=0.8)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

dists <- edge.attributes(permweights$CMIP5_CanESM2_r1i1p1)$distances
distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
distsc <- as.character(distsv)

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


ggplotfunction <- function(all_paths,edgelist,shortname) {ggplot() +
    geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
    geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = edge.ind)) +
    scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
    ggtitle(paste0(": ",nrow(edgelist)," ",shortname)) +
    theme(plot.title = element_text(hjust = 0.5))}

ggplots <- mapply(ggplotfunction, all_paths_list, edgelists, abrev, SIMPLIFY = FALSE)

do.call("grid.arrange",ggplots)

pdf('figs/longplots.pdf', width = 10, height = 20)
par(mfrow = c(4,3))
lapply(selection_hc_gcms, plot_long_distances, minimdist = 10000, smallcol = NA, bigcol = "red", data.dag = data.network, perm = permutations[[3]])
dev.off()

#####################################################################
# Plot networks with distances larger than X Km
#####################################################################
pdf('figs/longSplotsall.pdf', width = 30, height = 30)
par(mfrow = c(6,6))
mapply(plotS.spatial.igraph,network = permweights, title = abrev, MoreArgs = list(data.network = data.network,type = "bn",th = 0, smallcol = NA,th.type = "zoomweighted", from.dist = 10000,by.dist = 10000, perm = permutations[[3]],shift = TRUE))
dev.off()
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

timecoords <- TimeCoordsAnom_from_Grid_rms(cubicsets[[36]], rms = TRUE)
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
pdf(paste0("figs/bmatgraph_",th,"_",from.dist,"plotsall.pdf"), width = 10, height = 13)

par(mfrow = c(6,6))
mapply(plot_bmatgraph, network = bmatgraphse,  title = names(bmatgraphs), MoreArgs = list(th = th, from.dist = from.dist, curvature = TRUE,perm = permutations[[3]]))

dev.off()


#############################################################################
# How many steps are needed for CMIP5 models to go from inicial
# to end vertices for large links interim ? 
#############################################################################
# medir <- ends(bmatgraphse$interim_10d_akima_cubic,E(bmatgraphse$interim_10d_akima_cubic)[distances > 10000])

distance_th <-  5000
beta_th <- 0.3
medirs <- lapply(bmatgraphse, function(x) ends(x,E(x)[distances > distance_th & abs(strength) > beta_th]))

#medirs$CMIP5_CanESM2
#shortest_paths(bmatgraphse$CMIP5_EC.EARTH, from = medir[,1], to = medir[,2], mode = "all", weights = NA, output = "epath")

# testf <- function(x,y) {
  # apply(x, MARGIN = 1, function(z) shortest_paths(y, from = z[1], to = z[2], mode = "all", weights = NA, output = "epath"))}
testg <- function(x,y) {
  apply(x, MARGIN = 1, function(z) distances(y, v = z[1], to = z[2], mode = "all", weights = NA))}

#testf(medirs$CMIP5_CanESM2,bmatgraphse$CMIP5_CNRM.CM5)
#testg(medirs$CMIP5_CanESM2,bmatgraphse$CMIP5_CNRM.CM5)

distancesg <- lapply(medirs, function(z) sapply(bmatgraphse, function(w) testg(z,w)))
pathsmedirs<- sapply(distancesg, function(z) colSums(z))
pathsmedirsrel <- t(t(pathsmedirs)/diag(pathsmedirs))


institutions1 <- c("Had","CNRM","Can","EC","GFDL","IPSL","MIROC","MPI","Nor")
institutions <- rep(institutions1,each =2)
CMIPS <- rep(c("CMIP5","CMIP6"),length(institutions1))
combinations <- cbind(CMIPS,institutions)
namesort <- unlist(apply(combinations,MARGIN = 1,function(x) grep(paste0(x[1],".*",x[2],"."),rownames(pathsmedirs))))
namesort <- c(namesort,36)
Order <- TRUE

# OLD method:
#pathsmedirs <- sapply(medirs, function(z) sapply(bmatgraphse, function(w) sum(testg(z,w))))
#pathsmedirsrel <- sapply(medirs, function(z) sapply(bmatgraphse, function(w) sum(testg(z,w))/length(testg(z,w))))

if(Order == TRUE){
  pathsmedirsreldata <- melt(pathsmedirsrel[namesort,namesort]) 
  pathsmedirsdata <- melt(pathsmedirs[namesort,namesort]) 
} else {  pathsmedirsreldata <- melt(pathsmedirsrel) 
pathsmedirsdata <- melt(pathsmedirs) }

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



totalplot <- ggplot(pathsmedirsdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low=lowcolor, high="darkgreen") +
  geom_text(aes(label = round(value,1))) +
  labs(x="Model 2", y="Model 1", title=paste0("Total pathlength of Model 1 to reach edges of\nlength > ",distance_th," and |beta| > ", beta_th," of Model 2")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/strength_vs_dist/shortest_distances_",distance_th,"_",beta_th,"_cmip5_6_earth.pdf")
plotname
pdf(plotname, width = 20, height = 30)
joep <- arrangeGrob(totalplot,relplot,nrow = 2)
grid.arrange(joep, newpage = FALSE)
dev.off()


