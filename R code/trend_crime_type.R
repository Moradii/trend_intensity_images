library(spatstat)
library(remote)
library(raster)
library(gimms)
library(trend)
library(maptools)

## Import Data
load("crime_type_lpp.RData")

## Import Network
load("london_data.RData")

#########################################################
#########################################################
######################################################### Antisocial
#########################################################
#########################################################

antisoc.lpp <- antisocial.lpp
plot(antisoc.lpp[[1]],col=1,pch=20,cex=2,main="")

n_antisoc <- length(antisoc.lpp)


#########################################################
######################################################### common bw
#########################################################
sigma <- as.vector(unlist(lapply(antisoc.lpp, bw.scott.iso)))
sigma1 <- exp(mean(log(sigma)))

#########################################################
######################################################### Intensity estimate
#########################################################
d_antisoc <- lapply(X=1:n_antisoc, function(i){
  # s <- bw.scott.iso(antisoc.lpp[[i]])
  # out <- density.lpp(antisoc.lpp[[i]],sigma=s,distance = "euclidean",positive=TRUE)
  out <- density.lpp(antisoc.lpp[[i]],sigma=sigma1,distance = "euclidean",positive=TRUE)
  raster(as.im(out))
})

d_antisoc_brick <- brick(d_antisoc)
month <- rep(array(month.name),5)
year <- rep(c(2013:2017),each=12)
names(d_antisoc_brick) <- paste(month, year, sep="_")

##########################################################
########################################################## trend detection
##########################################################

d_antisoc_brick_deseason <- deseason(d_antisoc_brick)
d_antisoc_brick_deseason_agg <- aggregate(d_antisoc_brick_deseason,fac=2)
spplot(d_antisoc_brick_deseason_agg[[1:12]])
matplot(t(as.data.frame(d_antisoc_brick_deseason_agg)),type="l",ylab = "")

pacfs_antisoc <- lapply(X=1:ncell(d_antisoc_brick_deseason_agg), function(i){
  print(i)
  x <- d_antisoc_brick_deseason_agg[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_antisoc)
pacf_coef_antisoc <- numeric()
for (i in 1:length(pacfs_antisoc)) {
  if(class(pacfs_antisoc[[i]])=="numeric"){pacf_coef_antisoc[i] <- NA}
  else{pacf_coef_antisoc[i] <- pacfs_antisoc[[i]]$acf[1]}
}

mean(pacf_coef_antisoc[!is.na(pacf_coef_antisoc)])
summary(pacf_coef_antisoc[!is.na(pacf_coef_antisoc)])
par(mar=rep(4,4))
boxplot(pacf_coef_antisoc[!is.na(pacf_coef_antisoc)])

pacf.antisoc.image <- d_antisoc_brick_deseason_agg[[1]]
pacf.antisoc.image@data@values <- pacf_coef_antisoc
spplot(pacf.antisoc.image,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

d_antisoc_brick_deseason_agg_mult <- t(as.matrix(d_antisoc_brick_deseason_agg))
d_antisoc_brick_deseason_agg_mult <- d_antisoc_brick_deseason_agg_mult [,!apply(d_antisoc_brick_deseason_agg_mult,2, anyNA)]
mult.mk.test(ts(d_antisoc_brick_deseason_agg_mult))
# Multivariate Mann-Kendall Trend Test
# 
# data:  ts(d_antisoc_brick_deseason_agg_mult)
# z = 1.9229, p-value = 0.0545
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S        varS 
# 297400 23921544144

##########################################################
########################################################## Uni MK
##########################################################

antisoc_trend <- significantTau(d_antisoc_brick_deseason_agg,prewhitening=FALSE,p=0.05)
plot(antisoc_trend)


#########################################################
#########################################################
######################################################### Shoplifting
#########################################################
#########################################################


plot(shoplifting.lpp[[1]],col=1,pch=20,cex=2,main="")

n_shoplifting <- length(shoplifting.lpp)

#########################################################
######################################################### common bw
#########################################################
sigmashop <- as.vector(unlist(lapply(shoplifting.lpp, bw.scott.iso)))
sigmashop1 <- exp(mean(log(sigmashop)))

#########################################################
######################################################### Intensity estimate
#########################################################
d_shoplifting <- lapply(X=1:n_shoplifting, function(i){
  # s <- bw.scott.iso(shoplifting.lpp[[i]])
  # out <- density.lpp(shoplifting.lpp[[i]],sigma=s,distance = "euclidean",positive=TRUE)
  out <- density.lpp(shoplifting.lpp[[i]],sigma=sigmashop1,distance = "euclidean",positive=TRUE)
  raster(as.im(out))
})

d_shoplifting_brick <- brick(d_shoplifting)
month <- rep(array(month.name),5)
year <- rep(c(2013:2017),each=12)
names(d_shoplifting_brick) <- paste(month, year, sep="_")

##########################################################
########################################################## trend detection
##########################################################

d_shoplifting_brick_deseason <- deseason(d_shoplifting_brick)
d_shoplifting_brick_deseason_agg <- aggregate(d_shoplifting_brick_deseason,fac=2)
spplot(d_shoplifting_brick_deseason_agg[[1:12]])
matplot(t(as.data.frame(d_shoplifting_brick_deseason_agg)),type="l",ylab = "")

pacfs_shoplifting <- lapply(X=1:ncell(d_shoplifting_brick_deseason_agg), function(i){
  print(i)
  x <- d_shoplifting_brick_deseason_agg[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_shoplifting)
pacf_coef_shoplifting <- numeric()
for (i in 1:length(pacfs_shoplifting)) {
  if(class(pacfs_shoplifting[[i]])=="numeric"){pacf_coef_shoplifting[i] <- NA}
  else{pacf_coef_shoplifting[i] <- pacfs_shoplifting[[i]]$acf[1]}
}

mean(pacf_coef_shoplifting[!is.na(pacf_coef_shoplifting)])
summary(pacf_coef_shoplifting[!is.na(pacf_coef_shoplifting)])
par(mar=rep(4,4))
boxplot(pacf_coef_shoplifting[!is.na(pacf_coef_shoplifting)])

pacf.shoplifting.image <- d_shoplifting_brick_deseason_agg[[1]]
pacf.shoplifting.image@data@values <- pacf_coef_shoplifting
spplot(pacf.shoplifting.image,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

d_shoplifting_brick_deseason_agg_mult <- t(as.matrix(d_shoplifting_brick_deseason_agg))
d_shoplifting_brick_deseason_agg_mult <- d_shoplifting_brick_deseason_agg_mult [,!apply(d_shoplifting_brick_deseason_agg_mult,2, anyNA)]
mult.mk.test(ts(d_shoplifting_brick_deseason_agg_mult))
# Multivariate Mann-Kendall Trend Test
# 
# data:  ts(d_shoplifting_brick_deseason_agg_mult)
# z = 2.5993, p-value = 0.00934
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S       varS 
# 216636 6945994774


##########################################################
########################################################## Uni MK
##########################################################

shoplifting_trend <- significantTau(d_shoplifting_brick_deseason_agg,prewhitening=FALSE,p=0.05)
plot(shoplifting_trend)

#########################################################
#########################################################
######################################################### Drug Dealings
#########################################################
#########################################################


plot(drug.lpp[[1]],col=1,pch=20,cex=2,main="")
n_drug <- length(drug.lpp)


#########################################################
######################################################### common bw
#########################################################
sigmadrug <- as.vector(unlist(lapply(drug.lpp, bw.scott.iso)))
sigmadrug1 <- exp(mean(log(sigmadrug)))

#########################################################
######################################################### Intensity estimate
#########################################################
d_drug <- lapply(X=1:n_drug, function(i){
  # s <- bw.scott.iso(drug.lpp[[i]])
  # out <- density.lpp(drug.lpp[[i]],sigma=s,distance = "euclidean",positive=TRUE)
  out <- density.lpp(drug.lpp[[i]],sigma=sigmadrug1,distance = "euclidean",positive=TRUE)
  raster(as.im(out))
})

d_drug_brick <- brick(d_drug)
month <- rep(array(month.name),5)
year <- rep(c(2013:2017),each=12)
names(d_drug_brick) <- paste(month, year, sep="_")

##########################################################
########################################################## trend detection
##########################################################

d_drug_brick_deseason <- deseason(d_drug_brick)
d_drug_brick_deseason_agg <- aggregate(d_drug_brick_deseason,fac=2)
spplot(d_drug_brick_deseason_agg[[1:12]])
matplot(t(as.data.frame(d_drug_brick_deseason_agg)),type="l",ylab = "")

pacfs_drug <- lapply(X=1:ncell(d_drug_brick_deseason_agg), function(i){
  print(i)
  x <- d_drug_brick_deseason_agg[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_drug)
pacf_coef_drug <- numeric()
for (i in 1:length(pacfs_drug)) {
  if(class(pacfs_drug[[i]])=="numeric"){pacf_coef_drug[i] <- NA}
  else{pacf_coef_drug[i] <- pacfs_drug[[i]]$acf[1]}
}

mean(pacf_coef_drug[!is.na(pacf_coef_drug)])
summary(pacf_coef_drug[!is.na(pacf_coef_drug)])
par(mar=rep(4,4))
boxplot(pacf_coef_drug[!is.na(pacf_coef_drug)])

pacf.drug.image <- d_drug_brick_deseason_agg[[1]]
pacf.drug.image@data@values <- pacf_coef_drug
spplot(pacf.drug.image,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

d_drug_brick_deseason_agg_mult <- t(as.matrix(d_drug_brick_deseason_agg))
d_drug_brick_deseason_agg_mult <- d_drug_brick_deseason_agg_mult [,!apply(d_drug_brick_deseason_agg_mult,2, anyNA)]
mult.mk.test(ts(d_drug_brick_deseason_agg_mult))
# Multivariate Mann-Kendall Trend Test
# 
# data:  ts(d_drug_brick_deseason_agg_mult)
# z = -3.9133, p-value = 9.104e-05
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S       varS 
# -315112 6483962237 

##########################################################
########################################################## Uni MK
##########################################################

drug_trend <- significantTau(d_drug_brick_deseason_agg,prewhitening=FALSE,p=0.05)
plot(drug_trend)


#########################################################
#########################################################
######################################################### Vehicle Crime
#########################################################
#########################################################

plot(vehicle.lpp[[1]],col=1,pch=20,cex=2,main="")
n_vehicle <- length(vehicle.lpp)

#########################################################
######################################################### common bw
#########################################################
sigmavehicle <- as.vector(unlist(lapply(vehicle.lpp, bw.scott.iso)))
sigmavehicle1 <- exp(mean(log(sigmavehicle)))

#########################################################
######################################################### Intensity estimate
#########################################################
d_vehicle <- lapply(X=1:n_vehicle, function(i){
  # s <- bw.scott.iso(vehicle.lpp[[i]])
  # out <- density.lpp(vehicle.lpp[[i]],sigma=s,distance = "euclidean",positive=TRUE)
  out <- density.lpp(vehicle.lpp[[i]],sigma=sigmavehicle1,distance = "euclidean",positive=TRUE)
  raster(as.im(out))
})

d_vehicle_brick <- brick(d_vehicle)
month <- rep(array(month.name),5)
year <- rep(c(2013:2017),each=12)
names(d_vehicle_brick) <- paste(month, year, sep="_")

##########################################################
########################################################## trend detection
##########################################################

d_vehicle_brick_deseason <- deseason(d_vehicle_brick)
d_vehicle_brick_deseason_agg <- aggregate(d_vehicle_brick_deseason,fac=2)
spplot(d_vehicle_brick_deseason_agg[[1:12]])
matplot(t(as.data.frame(d_vehicle_brick_deseason_agg)),type="l",ylab = "")

pacfs_vehicle <- lapply(X=1:ncell(d_vehicle_brick_deseason_agg), function(i){
  print(i)
  x <- d_vehicle_brick_deseason_agg[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_vehicle)
pacf_coef_vehicle <- numeric()
for (i in 1:length(pacfs_vehicle)) {
  if(class(pacfs_vehicle[[i]])=="numeric"){pacf_coef_vehicle[i] <- NA}
  else{pacf_coef_vehicle[i] <- pacfs_vehicle[[i]]$acf[1]}
}

mean(pacf_coef_vehicle[!is.na(pacf_coef_vehicle)])
summary(pacf_coef_vehicle[!is.na(pacf_coef_vehicle)])
par(mar=rep(4,4))
boxplot(pacf_coef_vehicle[!is.na(pacf_coef_vehicle)])

pacf.vehicle.image <- d_vehicle_brick_deseason_agg[[1]]
pacf.vehicle.image@data@values <- pacf_coef_vehicle
spplot(pacf.vehicle.image,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

d_vehicle_brick_deseason_agg_mult <- t(as.matrix(d_vehicle_brick_deseason_agg))
d_vehicle_brick_deseason_agg_mult <- d_vehicle_brick_deseason_agg_mult [,!apply(d_vehicle_brick_deseason_agg_mult,2, anyNA)]
mult.mk.test(ts(d_vehicle_brick_deseason_agg_mult))
# Multivariate Mann-Kendall Trend Test
# 
# data:  ts(d_vehicle_brick_deseason_agg_mult)
# z = 0.59152, p-value = 0.5542
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S        varS 
# 76368 16667763064 

##########################################################
########################################################## Uni MK
##########################################################

vehicle_trend <- significantTau(d_vehicle_brick_deseason_agg,prewhitening=FALSE,p=0.05)
plot(vehicle_trend)


##########################################################
##########################################################
##########################################################    All plots
##########################################################
##########################################################

North <- list("SpatialPolygonsRescale", layout.north.arrow(type=1),
              offset = c(530500, 181500), scale = 800, which = 1)
################################################

## Combined PACF Plot
trend_detected <- stack(pacf.antisoc.image,pacf.shoplifting.image,pacf.drug.image, pacf.vehicle.image)
names(trend_detected) <- c("Antisocial_Behaviour","Shoplifting","Drugs", "Vehicle_Crime")

png("pacfcrimetype.png",width = 940, height = 340)
spplot(trend_detected,sp.layout=list(list(roads),North), colorkey=list(labels=list(cex=1.5)),
       col.regions = rev(topo.colors(20, alpha = 1)),
       scales=list(draw=T,cex=1.5),
       par.strip.text=list(cex=1.5),
       cex = 1.5)
dev.off()

## Combined Trend Plot
trend_detected <- stack(antisoc_trend, shoplifting_trend, drug_trend, vehicle_trend)
names(trend_detected) <- c("Antisocial_Behaviour","Shoplifting","Drugs", "Vehicle_Crime")

png("trendcrimetype.png",width = 940, height = 340)
spplot(trend_detected,sp.layout=list(list(roads),North), colorkey=list(labels=list(cex=1.5)),
       col.regions = rev(topo.colors(20, alpha = 1)),
       scales=list(draw=T,cex=1.5),
       par.strip.text=list(cex=1.5),
       cex = 1.5)
dev.off()

save.image("results(crimetype).RData")
