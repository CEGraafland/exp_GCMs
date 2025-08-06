#########################################################################
# Analysis Bhattacharya distances full CMIP5 and CMIP6 Historical
#########################################################################
# Partitional clustering
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
library("cluster")
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
selecteddists <- hellinger_coefficients
###############################################################
# Hierarchical clustering batachari
###############################################################
library(proxy)
which(rownames(selecteddists) == "CMIP5_EC.EARTH_r1i1p1_historical")
selecteddists <-selecteddists[-32,-32]
d <- dist(-log(selecteddists),method = "minkowski", p = 2)
mclust <-hclust(dist(-log(selecteddists),method = "minkowski", p = 2),method = "complete")
plot(mclust, hang = -1)
rect.hclust(mclust,k=11)
dev.off()

orderby <- "ncep_10d"
orderby <- "interim_10d_akima_cubic"
ordersort <- sort(-log(selecteddists)[,orderby],index.return = TRUE,decreasing = TRUE)

ordersets <- TRUE
if (ordersets == TRUE){ longdata <- melt(d[ordersort$ix,ordersort$ix])
} else {longdata <- melt(-log(selecteddists))}

# order by clustering.

get_upper_tri <- function(data2){
  data2[upper.tri(data2)]<- NA
  return(data2)
}
orderedmat<- -log(selecteddists)[mclust$order,mclust$order]
upper_tri <- get_upper_tri(orderedmat)
longdata <- melt(upper_tri, na.rm = TRUE)

##############################################################################
# Visualize  Hellinger Distance global set full matrix
##############################################################################
code <- FALSE
if(code == TRUE) {r.i. <- sapply(rownames(orderedmat),function(x)which(x == codes$rn.orderedmat))
toChange <- as.logical(sapply(r.i., length))
rownames(orderedmat)[toChange]<- codes$abrev.y.atm[unlist(r.i.[toChange])]
colnames(orderedmat)[toChange]<- codes$abrev.y.atm[unlist(r.i.[toChange])]}

rownames(orderedmat) <- gsub("CMIP6Amon_","",gsub("historical","",rownames(orderedmat)))
colnames(orderedmat) <- gsub("CMIP6Amon_","",gsub("historical","",colnames(orderedmat)))
longdata <- melt(orderedmat)
n <- 7
br.lims <- quantile(longdata$value, seq(0,1,1/n))
data.vec <- longdata$value

quant.discreter <- function(n,data.vec){
  br.lims <- quantile(data.vec, seq(0,1,1/n))
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}

longdata$d.value <- quant.discreter(n,longdata$value)
modelsize <- 18

b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[2:(n+1)],0))+
  geom_text(aes(label = round(value,0)),size =8) +
  labs(title=paste0("Bhattacharyi distance",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=25, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=25),
        plot.title=element_text(size=25)) +
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

b

plotname <- "figs/Bhattacharya_CMIP5_CMIP6/Bhattacharya_CMIP5_CMIP6.pdf"
pdf(plotname,width = 15, height = 15)
b             
dev.off()

plotname <- paste0("figs/Bhattacharya_CMIP5_CMIP6/Bhattacharya_CMIP5_CMIP6.png")
png(plotname,units = "in", width = 30, height = 30, res =180, pointsize =12)
b             
dev.off()


# example subset GFDL
subr <-grep("GFDL",rownames(orderedmat),fixed = FALSE)
rownames(orderedmat)[subr]

longdata <- melt(orderedmat[subr,])
n <- 7
br.lims <- quantile(longdata$value, seq(0,1,1/n))
data.vec <- longdata$value

quant.discreter <- function(n,data.vec){
  br.lims <- quantile(data.vec, seq(0,1,1/n))
  d.vec <- numeric(length = length(data.vec))
  # j <- 1
  for(j in 1:length(data.vec)){
    y <- data.vec[j]
    i <- 1
    while(i <= (length(br.lims)-1)){
      if(y >= br.lims[i] & y <= br.lims[i+1]){
        d.vec[j] <- i
        break
      } else {i <- i +1}
    } 
  }
  return(d.vec)
}

longdata$d.value <- quant.discreter(n,longdata$value)
modelsize <- 18

b <- ggplot(longdata, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=d.value)) + 
  scale_fill_gradientn(aesthetics = c("fill"), colours = rev(brewer.pal(n,"Blues")),guide = "legend",name = "value",labels = round(br.lims[2:(n+1)],0))+
  geom_text(aes(label = round(value,0)),size =8) +
  labs(x="", y="",title=paste0("Bhattacharyi distance",modelsize*100," CMIP models")) +
  theme_bw() + 
  theme(axis.text.x=element_text(size=25, angle=270, vjust=1, hjust = 1),
        axis.text.y=element_text(size=25),
        plot.title=element_text(size=25)) +
  scale_y_discrete(position = "left") + 
  scale_x_discrete(position = "top")

b

plotname <- "figs/Bhattacharya_CMIP5_CMIP6/Bhattacharya_CMIP5_CMIP6_GFDL.pdf"
pdf(plotname,width = 15, height =10)
b             
dev.off()

plotname <- paste0("figs/Bhattacharya_CMIP5_CMIP6/Bhattacharya_CMIP5_CMIP6_GFDL.png")
png(plotname,units = "in", width = 30, height = 12, res =180, pointsize =12)
b             
dev.off()


grep()