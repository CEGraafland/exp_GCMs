#####################################################################################
# Make correlation figures models. 
#####################################################################################
distance_sign <- "less"
distance_th <-  7500
beta_sign <- "greater"
beta_th <- 0.15

distancesgreater <- get(load(paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/distance_weights_arcs/arcs_","greater",distance_th,"_",beta_sign,beta_th,".rda")))
distancesless <- get(load(paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/results/distance_weights_arcs/arcs_","less",distance_th,"_",beta_sign,beta_th,".rda")))

dimens <- lapply(distancesg, dim)
pathsmedirsgreater<- sapply(distancesgreater, function(z) if (!is.null(dim(z))) {colSums(z)} else {rep(NA,length(z))})
pathsmedirsless <- sapply(distancesless, function(z) if (!is.null(dim(z))) {colSums(z)} else {rep(NA,length(z))})

pathsmedirsrelgreater <- t(t(pathsmedirsgreater)/diag(pathsmedirsgreater))
pathsmedirsrelgreater[which(pathsmedirsrelgreater == Inf)] <- NA

pathsmedirsrelless <- t(t(pathsmedirsless)/diag(pathsmedirsless))
pathsmedirsrelless[which(pathsmedirsrelless == Inf)] <- NA

load(file ="results/hellinger_coefficient/hellinger_coefficients.rda")
helcoef_eq <- hellinger_coefficients[rownames(pathsmedirsrelgreater),colnames(pathsmedirsrelgreater)]



 model <- "CMIP5_HadGEM2.ES_r1i1p1"
# model <- "CMIP5_bcc.csm1.1.m_r1i1p1"
# model <- "CMIP5_EC.EARTH_r1i1p1"
# model <- "CMIP5_BNU.ESM_r1i1p1" 
# model <- "CMIP5_ACCESS1.0_r1i1p1"
# model <- "interim_10d_akima_cubic"
# model <- "JRA55_10d_akima_cubic"

for (model in rownames(pathsmedirsrelgreater)[-13]){
# greater
# how others reach you; without auto correaltion
  cor_greater_tomodel <- cor(-log(helcoef_eq)[model,-which(colnames(helcoef_eq) == model)],pathsmedirsrelgreater[-which(rownames(pathsmedirsrelgreater) == model),model],use = "pairwise.complete.obs")
# how you reach others: without auto correlation
cor_greater_frommodel <- cor(-log(helcoef_eq)[model,-which(colnames(helcoef_eq) == model)],pathsmedirsrelgreater[model,-which(colnames(pathsmedirsrelgreater) == model)],use = "pairwise.complete.obs")

# less
# how others reach you; without auto correaltion
cor_less_tomodel <- cor(-log(helcoef_eq)[model,-which(colnames(helcoef_eq) == model)],pathsmedirsrelless[-which(rownames(pathsmedirsrelless) == model),model],use = "pairwise.complete.obs")
# how you reach others: without auto correlation
cor_less_frommodel <- cor(-log(helcoef_eq)[model,-which(colnames(helcoef_eq) == model)],pathsmedirsrelless[model,-which(colnames(pathsmedirsrelless) == model)],use = "pairwise.complete.obs")

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/figs/hellinger_wrt_distance/",model,"_distance",distance_th,"_weight_",beta_th,".pdf")
pdf(plotname)
par(mfrow = c(2,2))
plot(-log(helcoef_eq[model,-which(colnames(helcoef_eq) == model)]), pathsmedirsrelgreater[model,-which(colnames(pathsmedirsrelless) == model)], xlab = "Hellinger distance", ylab = paste0("mean shortest path to edges > ",distance_th, " other models"), main = paste0(model,"\n cor = ",cor_greater_frommodel))
if(all(is.na(pathsmedirsrelgreater[-which(colnames(pathsmedirsrelless) == model),model])) == TRUE){plot(0,type='n',axes=FALSE,ann=FALSE)} else {plot(-log(helcoef_eq[model,-which(colnames(helcoef_eq) == model)]), pathsmedirsrelgreater[-which(colnames(pathsmedirsrelless) == model),model], xlab = "Hellinger distance", ylab = paste0("mean shortest path other models to edges > ",distance_th,"\n ",model), main = paste0(model,"\n cor = ",cor_greater_tomodel))}
plot(-log(helcoef_eq[model,-which(colnames(helcoef_eq) == model)]), pathsmedirsrelless[model,-which(colnames(pathsmedirsrelless) == model)], xlab = "Hellinger distance", ylab = paste0("mean shortest path to edges < ",distance_th, " other models"), main = paste0(model,"\n cor = ",cor_less_frommodel))
plot(-log(helcoef_eq[model,-which(colnames(helcoef_eq) == model)]), pathsmedirsrelless[-which(colnames(pathsmedirsrelless) == model),model], xlab = "Hellinger distance", ylab = paste0("mean shortest path other models to edges < ",distance_th,"\n ",model), main = paste0(model,"\n cor = ",cor_less_tomodel))
dev.off()
}
model


# row: coonnexion between hellinger distance and how you reaching others
diag(apply(-log(helcoef_eq), MARGIN = c(1), function(x){apply(pathsmedirsrel, MARGIN = c(1), function(y) cor(x = x,y = y, use = "pairwise.complete.obs", method = "pearson") )}))
# column: conexion between hellinger distance and how others reaching you.
diag(cor(-log(helcoef_eq),pathsmedirsrel,use = "pairwise.complete.obs")) # column
