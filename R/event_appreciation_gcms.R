#############################################################################
# Real Data Events Appreciation by GCMs 
#############################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")

library(magrittr)
library(visualizeR)
library('RColorBrewer')
library(gridExtra)

source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("../R/Functions/propagationFunctions.R")

load("../Data/Struct_learn/datapermutations.rda")
load("../Data/Struct_learn/permutations.rda")
load("../Data/tas_ncep_10d.rda")
load("data/tas_interim_10d_akima_cubic.rda")

dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_interim_10d_akima_cubic, rms = TRUE)

col.blue <- rev(brewer.pal(4,"Blues"))
col.blue
col.red <- brewer.pal(4,"Reds")

col.div <- c(col.blue, col.red)
col.div[4:5] <- "white"
col.b <- c(col.blue,col.red)
length(col.b)

##########################################################################
# select simple / combine with paperScript single data evaluation. 
##########################################################################
dataRMS2 <- as.data.frame(dataRMS)

# worst 1:.....  or best .....:360
type <- "worst"
howmany <- 16
if(type == "best"){whichdates <- (360 - howmany +1):360} else if (type == "worst"){whichdates <- 1:howmany}
# no
whichBN <- "CMIP5_EC.EARTH"
whichgcm <- orderBNs[,whichBN]

# The event indices:
whichgcm[whichdates]

plots <- list()
climBig2 <- list()
for(i in 1:length(whichdates)){
  anomaliesdate <- dataRMS2[whichgcm[whichdates[i]],]
  
  climBig2[[i]] <- quantity2clim(as.matrix(anomaliesdate), what = "Big 2",tas_interim_10d_akima_cubic)
  plots[[i]] <- spatialPlot(climBig2[[i]], lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                            col.regions = col.div,set.min = -2.5, set.max = 2.5,
                            at = seq(-3,3,1))
}
ensclimBig <- bindGrid(climBig2, dimension = c("member"))
# ensclimBig$Members <- abrev
spatialPlot(ensclimBig, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1), as.table = TRUE,
            main = paste0(type," events ",whichBN))

plotname <- paste0("figs/events/",type,"_events_",whichBN,".pdf")
pdf(plotname)
do.call(grid.arrange,c(plots, top = paste0(type," events ",whichBN)))
dev.off()

#############################################################################
# With orders (logliks_gcms_vs_int_event is from model_evaluation_gcms.R)
#############################################################################
orderBNs <- matrix(data = NA, nrow = nrow(dataRMS), ncol = length(selection_fits_gcms))
colnames(orderBNs) <- names(selection_fits_gcms)
for(i in 1:length(selection_fits_gcms)){
  orderBNs[,i] <- order(logliks_gcms_vs_int_event[,i])
}


difgcm <- numeric(length = 360)
BN1 <-2
BN2 <-8
Method <- "logdif"
if (Method == "ranking") {for(i in 1:360){difgcm[i] <- which(orderBNs[,BN1] == i) - which(orderBNs[,BN2] == i)}
} else if (Method == "logdif") {difgcm <- order(logliks_gcms_vs_int_event[,BN1] - logliks_gcms_vs_int_event[,BN2])}

# negative is: 8 beter than 9, positive is 9 beter than 8)
(logliks_gcms_vs_int_event[,BN1] - logliks_gcms_vs_int_event[,BN2])[(difgcm)]
type <- "best"
howmany <- 16
if(type == "best"){whichdates <- (360 - howmany +1):360} else if (type == "worst"){whichdates <- 1:howmany}
plots <- list()
for(i in 1:length(whichdates)){
  anomaliesdate <- dataRMS2[order(difgcm)[whichdates[i]],]
  
  climBig2 <- quantity2clim(as.matrix(anomaliesdate), what = "Big 2",tas_interim_10d_akima_cubic)
  plots[[i]] <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                            col.regions = col.div,set.min = -2.5, set.max = 2.5,
                            at = seq(-3,3,1))
}

plotname <- paste0("figs/events/",type,"events",abrev[BN1],"vs",abrev[BN2],"_",Method,".pdf")
pdf(plotname)
do.call(grid.arrange,c(plots, top = paste0(type," events ",abrev[BN1]," versus ",abrev[BN2])))
dev.off()

where <- c()
for(i in 1:length(orderBN)){
  where[i] <-  which(orderBN == orderCM[i])
}
orderBN[where]

#####################################################################
#
#####################################################################
dataRMS2 <- as.data.frame(dataRMS)

Big2ind <-  which(dataRMS2$V81 >1.5 & dataRMS2$V227 > 1.5)
Big2ind
event1 <- 1:5
event2 <- 128:130
Big2ind <-65
dataRMS2[,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_interim_10d_akima_cubic)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))


Big2ind <- which(eval(parse(text = paste0("dataRMS2$",node," >=",value))))
Big2ind <- which(eval(parse(text = paste0("dataRMS2$",node," <=",-value))))
Big2ind <-  which(dataRMS2$V81 > 1.5 & dataRMS2$V227 > 1)
Big2ind <-  which(dataRMS2$V81 > 0.5)
events <- list()
events[[1]] <- c(Big2ind[1])
m <- 1
if (!length(Big2ind)==1){
  for (k in 2:length(Big2ind)){
    if ((Big2ind[k]-1) == Big2ind[k-1]){events[[m]] <- append(events[[m]], Big2ind[k])
    } else {events[[m+1]] <- c(Big2ind[k])
    m <- m+1}
  }
} 
plots <- list()
for(i in 1:length(events)){
  extremes <- colMeans(dataRMS2[events[[i]],])
  
  climBig2 <- quantity2clim(extremes, what = "Big 2",tas_interim_10d_akima_cubic)
  plots[[i]] <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                            col.regions = col.div,set.min = -2.5, set.max = 2.5,
                            at = seq(-3,3,1),
                            main = list(paste0(length(events[[i]]))))
}

do.call(grid.arrange,c(plots))

#####################################################################
# CMIP5 subset example empirical distribution. 
#####################################################################
listinterim <- list(tas_interim_10d_akima_cubic)
cubicsets5 <- get(load("data/tas_historical_10d_akima_cubic_corrected.rda"))
cubicsets5extra <- get(load("data/tas_historical_cmip5_extra_10d_akima_cubic.rda"))
cubicsets5left <- get(load("data/tas_historical_cmip5_left_10d_akima_cubic.rda"))

subbies <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1")
listsubs <- cubicsets5left[subbies]
listexample <- c(listsubs,listinterim)
names(listexample)<- c(subbies,c("tas_interim_10d_akima_cubic"))
listdataRMS <- lapply(listexample,TimeCoordsAnom_from_Grid_rms,rms = TRUE)
listdataRMS2 <- lapply(listdataRMS,as.data.frame)

listBig2ind <- lapply(listdataRMS2, function(x)x$V81>2)

listextremes <- mapply(function(x,y) colMeans(x[y,]), 
                       x = listdataRMS2,y = listBig2ind, SIMPLIFY = FALSE)
listclimBig2 <- mapply(function(x,y) quantity2clim(x,what = "Big 2",y),x = listextremes,y = listexample, SIMPLIFY = FALSE)

ensclimBig2 <- bindGrid(listclimBig2, dimension = c("member"))
ensclimBig2$Members <- names(listexample)
ensclimBig2$Members

plotname <- "figs/events/example_empirical.pdf"
pdf(plotname,width = 10,height = 3 )
spatialPlot(ensclimBig2, names.attr = names(listexample), layout = c(5,1), as.table = TRUE,lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,0.75), main = "V81 >= 2")
dev.off()

