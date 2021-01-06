library(dplyr)
library(ggplot2)
library(gimms)
library(maptools)
library(raster)
library(remote)
library(sparr)
library(spatstat)
library(trend)

## Import Crime Data
load("Data/crime_type_lpp.RData")

## Import Network
load("Data/london_data.RData")



#########################################################
#########################################################
######################################################### relative risk (antisocial - drug dealings)
#########################################################
#########################################################

#########################################################
######################################################### common bw
#########################################################
n <- length(antisocial.lpp)

sigma_ratio2 <- lapply(X=1:n, function(i){
  x2 <- as.ppp(antisocial.lpp[[i]])
  y2 <- as.ppp(drug.lpp[[i]])
  s2 <- LSCV.risk(x2,y2,method = "davies",edge = F,hlim = c(0,2000))
  return(s2)
})
sigma_ratio2 <- unlist(sigma_ratio2)
sigma_ratio2 <- exp(mean(log(sigma_ratio2)))

#########################################################
######################################################### risk estimate
#########################################################

r2 <- lapply(X=1:n, function(i){
  # x <- as.ppp(acc.ppp[[i]])
  # y <- as.ppp(crime.ppp[[i]])
  # s <- LSCV.risk(x,y,method = "davies",edge = F,hlim = c(0,2000))
  # print(s)
  out1_2 <- density.lpp(antisocial.lpp[[i]],sigma=sigma_ratio2,distance = "euclidean",positive=TRUE)
  out2_2 <- density.lpp(drug.lpp[[i]],sigma=sigma_ratio2,distance = "euclidean",positive=TRUE)
  return(raster(as.im(out1_2/out2_2)))
})

r_brick2 <- brick(r2)
spplot(r_brick2[[1:12]])

##########################################################
########################################################## trend detection
##########################################################


r_deseason2 <- deseason(r_brick2)
r_deseason_agg2 <- aggregate(r_deseason2, fac=2, na.rm=TRUE, na.action=na.pass) 



spplot(r_deseason_agg2[[1:12]])
matplot(t(as.data.frame(r_deseason_agg2)),type="l")

pacfs_r2 <- lapply(X=1:ncell(r_deseason_agg2), function(i){
  print(i)
  x <- r_deseason_agg2[i]
  #x <- s(x)  # to omit NaN values
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_r2)
pacf_coef_r2 <- numeric()
for (i in 1:length(pacfs_r2)) {
  if(class(pacfs_r2[[i]])=="numeric"){pacf_coef_r2[i] <- NA}
  else{pacf_coef_r2[i] <- pacfs_r2[[i]]$acf[1]}
}

mean(pacf_coef_r2[!is.na(pacf_coef_r2)])
summary(pacf_coef_r2[!is.na(pacf_coef_r2)])
par(mar=rep(4,4))
boxplot(pacf_coef_r2[!is.na(pacf_coef_r2)])

pacf.image_r2 <- r_deseason_agg2[[1]]
pacf.image_r2@data@values <- pacf_coef_r2
spplot(pacf.image_r2,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

r_deseason_agg_mult2 <- t(as.matrix(r_deseason_agg2))
r_deseason_agg_mult2 <- r_deseason_agg_mult2 [,!apply(r_deseason_agg_mult2  ,2, anyNA)]

mv_mk_2 <- mult.mk.test(ts(r_deseason_agg_mult2)) 


#####################################################
#####################################################
#####################################################
# Multivariate Mann-Kendall Trend Test             ##
#                                                  ##
# data:  ts(r_deseason_agg_mult2)                  ##
# z = 2.1472, p-value = 0.03177                    ##
# alternative hypothesis: true S is not equal to 0 ##
# sample estimates:                                ##
#  S        varS                                   ##
# 258067 14444422845                               ##
#####################################################
#####################################################
#####################################################


##########################################################
########################################################## Uni MK
##########################################################

r_trend2 <- significantTau(r_deseason_agg2,prewhitening=FALSE,p=0.05)
plot(r_trend2)

r_trend_pixels2 <- Which(r_trend2,cells=T)
#part_r <- d_crime_brick_deseason_agg[r_trend_pixels]

part_r2 <- r_deseason_agg2[r_trend_pixels2]
matplot(colMeans(part_r2),type = "l")

# Relative PACF Summary
summary(pacf_coef_r2[r_trend_pixels2])

###########################################################
#     Min.  1st Qu.   Median    Mean    3rd Qu.   Max.   ##
# -0.17121 -0.02066  0.01056  0.03554  0.08427  0.33387  ## 
###########################################################


#########################################################
#########################################################
######################################################### relative risk (drug dealings - shoplifting)
#########################################################
#########################################################

#########################################################
######################################################### common bw
#########################################################


n <- length(shoplifting.lpp)

sigma_ratio4 <- lapply(X=1:n, function(i){
  x4 <- as.ppp(shoplifting.lpp[[i]])
  y4 <- as.ppp(drug.lpp[[i]])
  s4 <- LSCV.risk(x4,y4,method = "davies",edge = F,hlim = c(0,2000))
  return(s4)
})
sigma_ratio4 <- unlist(sigma_ratio4)
sigma_ratio4 <- exp(mean(log(sigma_ratio4)))

#########################################################
######################################################### risk estimate
#########################################################

r4 <- lapply(X=1:n, function(i){
  # x <- as.ppp(acc.ppp[[i]])
  # y <- as.ppp(crime.ppp[[i]])
  # s <- LSCV.risk(x,y,method = "davies",edge = F,hlim = c(0,2000))
  # print(s)
  out1_4 <- density.lpp(shoplifting.lpp[[i]],sigma=sigma_ratio4,distance = "euclidean",positive=TRUE)
  out2_4 <- density.lpp(drug.lpp[[i]],sigma=sigma_ratio4,distance = "euclidean",positive=TRUE)
  return(raster(as.im(out2_4/out1_4)))
})

r_brick4 <- brick(r4)
spplot(r_brick4[[1:12]])

##########################################################
########################################################## trend detection
##########################################################

r_deseason4 <- deseason(r_brick4)
r_deseason_agg4 <- aggregate(r_deseason4,fac=2)
spplot(r_deseason_agg4[[1:12]])
matplot(t(as.data.frame(r_deseason_agg4)),type="l")

pacfs_r4 <- lapply(X=1:ncell(r_deseason_agg4), function(i){
  print(i)
  x <- r_deseason_agg4[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_r4)
pacf_coef_r4 <- numeric()
for (i in 1:length(pacfs_r4)) {
  if(class(pacfs_r4[[i]])=="numeric"){pacf_coef_r4[i] <- NA}
  else{pacf_coef_r4[i] <- pacfs_r4[[i]]$acf[1]}
}

mean(pacf_coef_r4[!is.na(pacf_coef_r4)])
summary(pacf_coef_r4[!is.na(pacf_coef_r4)])
par(mar=rep(4,4))
boxplot(pacf_coef_r4[!is.na(pacf_coef_r4)])

pacf.image_r4 <- r_deseason_agg4[[1]]
pacf.image_r4@data@values <- pacf_coef_r4
spplot(pacf.image_r4,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

r_deseason_agg_mult4 <- t(as.matrix(r_deseason_agg4))
r_deseason_agg_mult4 <- r_deseason_agg_mult4 [,!apply(r_deseason_agg_mult4  ,2, anyNA)]

mv_mk_4 <- mult.mk.test(ts(r_deseason_agg_mult4))


#####################################################
#####################################################
#####################################################
# Multivariate Mann-Kendall Trend Test             ##
#                                                  ##
# data:  ts(r_deseason_agg_mult4)                  ##
# z = -4.1479, p-value = 3.355e-05                 ##
# alternative hypothesis: true S is not equal to 0 ##
# sample estimates:                                ##
#  S       varS                                    ##
# -403666 9470690756                               ##
#####################################################
#####################################################
#####################################################


##########################################################
########################################################## Uni MK
##########################################################

r_trend4 <- significantTau(r_deseason_agg4,prewhitening=FALSE,p=0.05)
plot(r_trend4)

r_trend_pixels4 <- Which(r_trend4,cells=T)
#part_r <- d_crime_brick_deseason_agg[r_trend_pixels]

part_r4 <- r_deseason_agg4[r_trend_pixels4]
matplot(colMeans(part_r4),type = "l")

# Relative PACF Summary
summary(pacf_coef_r4[r_trend_pixels4])

###########################################################
#     Min.  1st Qu.   Median    Mean    3rd Qu.   Max.   ##
# -0.09728  0.07246  0.21162  0.20764  0.32144  0.53077  ## 
###########################################################


#########################################################
#########################################################
######################################################### relative risk (vehicle crime - shoplifting)
#########################################################
#########################################################

#########################################################
######################################################### common bw
#########################################################


n <- length(shoplifting.lpp)

sigma_ratio5 <- lapply(X=1:n, function(i){
  x5 <- as.ppp(shoplifting.lpp[[i]])
  y5 <- as.ppp(vehicle.lpp[[i]])
  s5 <- LSCV.risk(x5,y5,method = "davies",edge = F,hlim = c(0,2000))
  return(s5)
})
sigma_ratio5 <- unlist(sigma_ratio5)
sigma_ratio5 <- exp(mean(log(sigma_ratio5)))

#########################################################
######################################################### risk estimate
#########################################################

r5 <- lapply(X=1:n, function(i){
  # x <- as.ppp(acc.ppp[[i]])
  # y <- as.ppp(crime.ppp[[i]])
  # s <- LSCV.risk(x,y,method = "davies",edge = F,hlim = c(0,2000))
  # print(s)
  out1_5 <- density.lpp(shoplifting.lpp[[i]],sigma=sigma_ratio5,distance = "euclidean",positive=TRUE)
  out2_5 <- density.lpp(vehicle.lpp[[i]],sigma=sigma_ratio5,distance = "euclidean",positive=TRUE)
  return(raster(as.im(out2_5/out1_5)))
})

r_brick5 <- brick(r5)
spplot(r_brick5[[1:12]])

##########################################################
########################################################## trend detection
##########################################################

r_deseason5 <- deseason(r_brick5)
r_deseason_agg5 <- aggregate(r_deseason5,fac=2)
spplot(r_deseason_agg5[[1:12]])
matplot(t(as.data.frame(r_deseason_agg5)),type="l")

pacfs_r5 <- lapply(X=1:ncell(r_deseason_agg5), function(i){
  print(i)
  x <- r_deseason_agg5[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_r5)
pacf_coef_r5 <- numeric()
for (i in 1:length(pacfs_r5)) {
  if(class(pacfs_r5[[i]])=="numeric"){pacf_coef_r5[i] <- NA}
  else{pacf_coef_r5[i] <- pacfs_r5[[i]]$acf[1]}
}

mean(pacf_coef_r5[!is.na(pacf_coef_r5)])
summary(pacf_coef_r5[!is.na(pacf_coef_r5)])
par(mar=rep(4,4))
boxplot(pacf_coef_r5[!is.na(pacf_coef_r5)])

pacf.image_r5 <- r_deseason_agg5[[1]]
pacf.image_r5@data@values <- pacf_coef_r5
spplot(pacf.image_r5,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

r_deseason_agg_mult5 <- t(as.matrix(r_deseason_agg5))
r_deseason_agg_mult5 <- r_deseason_agg_mult5 [,!apply(r_deseason_agg_mult5  ,2, anyNA)]

mv_mk_5 <- mult.mk.test(ts(r_deseason_agg_mult5))


#####################################################
#####################################################
#####################################################
# Multivariate Mann-Kendall Trend Test             ##
#                                                  ##
# data:  ts(r_deseason_agg_mult5)                  ##
# z = -0.77913, p-value = 0.4359                   ##
# alternative hypothesis: true S is not equal to 0 ##
# sample estimates:                                ##
#  S        varS                                   ##
# -82000 11076632771                               ##
#####################################################
#####################################################
#####################################################


##########################################################
########################################################## Uni MK
##########################################################

r_trend5 <- significantTau(r_deseason_agg5,prewhitening=FALSE,p=0.05)
plot(r_trend5)

r_trend_pixels5 <- Which(r_trend5,cells=T)
#part_r <- d_crime_brick_deseason_agg[r_trend_pixels]

part_r5 <- r_deseason_agg5[r_trend_pixels5]
matplot(colMeans(part_r5),type = "l")

# Relative PACF Summary
summary(pacf_coef_r5[r_trend_pixels5])

###########################################################
#     Min.  1st Qu.   Median    Mean    3rd Qu.   Max.   ##
#   0.09081 0.15372  0.28746  0.28527   0.40579 0.47209  ## 
###########################################################

#########################################################
#########################################################
######################################################### relative risk (vehicle crime - drug dealings)
#########################################################
#########################################################

#########################################################
######################################################### common bw
#########################################################

n <- length(drug.lpp)

sigma_ratio6 <- lapply(X=1:n, function(i){
  x6 <- as.ppp(drug.lpp[[i]])
  y6 <- as.ppp(vehicle.lpp[[i]])
  s6 <- LSCV.risk(x6,y6,method = "davies",edge = F,hlim = c(0,2000))
  return(s6)
})
sigma_ratio6 <- unlist(sigma_ratio6)
sigma_ratio6 <- exp(mean(log(sigma_ratio6)))

#########################################################
######################################################### risk estimate
#########################################################

r6 <- lapply(X=1:n, function(i){
  # x <- as.ppp(acc.ppp[[i]])
  # y <- as.ppp(crime.ppp[[i]])
  # s <- LSCV.risk(x,y,method = "davies",edge = F,hlim = c(0,2000))
  # print(s)
  out1_6 <- density.lpp(drug.lpp[[i]],sigma=sigma_ratio6,distance = "euclidean",positive=TRUE)
  out2_6 <- density.lpp(vehicle.lpp[[i]],sigma=sigma_ratio6,distance = "euclidean",positive=TRUE)
  return(raster(as.im(out2_6/out1_6)))
})

r_brick6 <- brick(r6)
spplot(r_brick6[[1:12]])

##########################################################
########################################################## trend detection
##########################################################

r_deseason6 <- deseason(r_brick6)
r_deseason_agg6 <- aggregate(r_deseason6,fac=2)
spplot(r_deseason_agg6[[1:12]])
matplot(t(as.data.frame(r_deseason_agg6)),type="l")

pacfs_r6 <- lapply(X=1:ncell(r_deseason_agg6), function(i){
  print(i)
  x <- r_deseason_agg6[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_r6)
pacf_coef_r6 <- numeric()
for (i in 1:length(pacfs_r6)) {
  if(class(pacfs_r6[[i]])=="numeric"){pacf_coef_r6[i] <- NA}
  else{pacf_coef_r6[i] <- pacfs_r6[[i]]$acf[1]}
}

mean(pacf_coef_r6[!is.na(pacf_coef_r6)])
summary(pacf_coef_r6[!is.na(pacf_coef_r6)])
par(mar=rep(4,4))
boxplot(pacf_coef_r6[!is.na(pacf_coef_r6)])

pacf.image_r6 <- r_deseason_agg6[[1]]
pacf.image_r6@data@values <- pacf_coef_r6
spplot(pacf.image_r6,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

r_deseason_agg_mult6 <- t(as.matrix(r_deseason_agg6))
r_deseason_agg_mult6 <- r_deseason_agg_mult6 [,!apply(r_deseason_agg_mult6  ,2, anyNA)]

mv_mk_6 <- mult.mk.test(ts(r_deseason_agg_mult6))


#####################################################
#####################################################
#####################################################
# Multivariate Mann-Kendall Trend Test             ##
#                                                  ##
# data:  ts(r_deseason_agg_mult6)                  ##
# z = 1.5158, p-value = 0.1296                     ##
# alternative hypothesis: true S is not equal to 0 ##
# sample estimates:                                ##
#  S        varS                                   ##
# 151874 10039121143                               ##
#####################################################
#####################################################
#####################################################


##########################################################
########################################################## Uni MK
##########################################################

r_trend6 <- significantTau(r_deseason_agg6,prewhitening=FALSE,p=0.05)
plot(r_trend6)

r_trend_pixels6 <- Which(r_trend6,cells=T)
#part_r <- d_crime_brick_deseason_agg[r_trend_pixels]

part_r6 <- r_deseason_agg6[r_trend_pixels6]
matplot(colMeans(part_r6),type = "l")

# Relative PACF Summary
summary(pacf_coef_r6[r_trend_pixels6])

###########################################################
#     Min.  1st Qu.   Median    Mean    3rd Qu.   Max.   ##
# -0.02262  0.03290  0.08194  0.08688  0.14051  0.22728  ## 
###########################################################


##########################################################
##########################################################
##########################################################    All plots
##########################################################
##########################################################


############################################
############################################
## Relative Trends                       ###
                                         ###
## Antisocial vs. Drug  (a)              ###
## Drug vs. Shoplifting (b)              ###
## Vehicle_Crime vs. Shoplifting (c)     ###
## Vehicle_Crime vs. Drug (d)            ###
############################################
############################################


North <- list("SpatialPolygonsRescale", layout.north.arrow(type=1), 
              offset = c(530500, 181500), scale = 800, which = 1)


## Combined Trend Plot
trend_detected <- stack(r_trend2, r_trend4, r_trend5, r_trend6)
#names(trend_detected) <- c("Antisocial vs. Drug","Drug vs. Shoplifting", "Vehicle_Crime vs. Shoplifting", "Vehicle_Crime vs. Drug")
names(trend_detected) <- c("a","b", "c", "d")

png("trend_pattern_crime_type.png",width = 940, height = 340)
spplot(trend_detected,sp.layout=list(list(roads),North), colorkey=list(labels=list(cex=1.5)),
       col.regions = rev(topo.colors(20, alpha = 1)),
       scales=list(draw=T,cex=1.5),
       par.strip.text=list(cex=1.5),
       cex = 1.5)
dev.off()


save.image(file = "relative_trend(crime_type).RData")