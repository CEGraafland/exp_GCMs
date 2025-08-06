###########################################################################
# give weights by hellinger distance
###########################################################################
load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
ex.subset <- c("CMIP5_BNU.ESM_r1i1p1","CMIP5_ACCESS1.0_r1i1p1","CMIP5_ACCESS1.3_r1i1p1","CMIP5_CMCC.CMS_r1i1p1","interim_10d_akima_cubic")
bdists <- hellinger_coefficients

bexp <- -log(hellinger_coefficients[ex.subset,ex.subset])[-2,-2]
ref.perf <- "interim_10d_akima_cubic"
refID <- which(rownames(bdists) == ref.perf)
D <- bdists[ref.perf,-refID]
S <- bdists[-refID,-refID]
namei <- "CMIP5_BNU.ESM_r1i1p1"
sS <- 30
sD <- 25


IPweight <- function(bdists,ref.perf,sD,sS){
  refID <- which(rownames(bdists) == ref.perf)
  D <- bdists[ref.perf,-refID]
  S <- bdists[-refID,-refID]
  W <- numeric(length=length(D))
  for(i in 1:length(names(D))){
    name.i <-names(D)[i]
    nameID.S <- which(rownames(S) == name.i)
    nameID.D <- which(names(D) == name.i)
    Di <- D[nameID.D]
    Sij <- S[nameID.S,][-nameID.S]
    sumS <- sum(exp(-Sij^2/sS^2))
    W[i] <- exp(-Di^2/sD^2)/(1+sumS)
  }
  W <- W/sum(W)
  names(W) <- names(D)
  return(W)
}


##########################################
# Model weights cmip5 set
##########################################
matri <--log(hellinger_coefficients)[1:(nrow(hellinger_coefficients)-2),1:(ncol(hellinger_coefficients)-2)]
modelweights <- IPweight(matri,"interim_10d_akima_cubic",40,25)

par(mar = c(8, 4, 0, 2))
plot(modelweights,xaxt =  "n",xlab = "")
xlabels=gsub("r1i1p1","",gsub("historical","",names(modelweights)))
axis(1, at=1:length(modelweights),labels = xlabels, col.axis="red", las=2)
abline(h = 1/length(modelweights), lty = 2, col = "red")
dev.off()
###########################################
# Model weights subset example 
###########################################
modelweights.ex <- IPweight(-log(hellinger_coefficients[ex.subset,ex.subset]),"interim_10d_akima_cubic",sD = 40, sS = 25)
bdists <- hellinger_coefficients
par(mar = c(5,4,4,2)+0.1)
plot(modelweights.ex,xaxt =  "n", ylim = c(0,0.6), ylab = expression("modelweights w"[i]),xlab = "", pch = 16, col = "blue")

abline(h = 0.25, lty = 2, col = "red")
xlabels=gsub("r1i1p1","",gsub("historical","",names(modelweights.ex)))
axis(1, at=1:length(modelweights.ex),labels = xlabels, col.axis="red", las=2)
