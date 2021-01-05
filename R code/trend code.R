library(spatstat)
library(remote)
library(raster)
library(gimms)
library(sparr)
library(trend)
#########################################################
#########################################################
######################################################### accident
#########################################################
#########################################################
plot(acc.ppp[[1]],col=1,pch=20,cex=2,main="")
n.acc <- unlist(lapply(acc.ppp, npoints))
ts.plot(.Last.value[1:12])

n <- length(acc.ppp)

sp_acc <- lapply(X=1:n, function(i){
  SpatialPoints(coords(acc.ppp[[i]]))
})

#########################################################
######################################################### common bw
#########################################################
sigma <- as.vector(unlist(lapply(acc.ppp, bw.scott.iso)))
sigma1 <- exp(mean(log(sigma)))

#########################################################
######################################################### Intensity estimate
#########################################################
d <- lapply(X=1:n, function(i){
  # s <- bw.scott.iso(acc.ppp[[i]])
  # out <- density.lpp(acc.ppp[[i]],sigma=s,distance = "euclidean",positive=TRUE)
  out <- density.lpp(acc.ppp[[i]],sigma=sigma1,distance = "euclidean",positive=TRUE)
  raster(as.im(out))
})

d_brick <- brick(d)
month <- rep(array(month.name),5)
year <- rep(c(2013:2017),each=12)
names(d_brick) <- paste(month, year, sep="_")

##########################################################
########################################################## trend detection
##########################################################

d_brick_deseason <- deseason(d_brick)
d_brick_deseason_agg <- aggregate(d_brick_deseason,fac=2)
spplot(d_brick_deseason_agg[[1:12]])
matplot(t(as.data.frame(d_brick_deseason_agg)),type="l",ylab = "")

pacfs <- lapply(X=1:ncell(d_brick_deseason_agg), function(i){
  print(i)
  x <- d_brick_deseason_agg[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs)
pacf_coef <- numeric()
for (i in 1:length(pacfs)) {
  if(class(pacfs[[i]])=="numeric"){pacf_coef[i] <- NA}
  else{pacf_coef[i] <- pacfs[[i]]$acf[1]}
}

mean(pacf_coef[!is.na(pacf_coef)])
summary(pacf_coef[!is.na(pacf_coef)])
par(mar=rep(4,4))
boxplot(pacf_coef[!is.na(pacf_coef)])

pacf.acc.image <- d_brick_deseason_agg[[1]]
pacf.acc.image@data@values <- pacf_coef
spplot(pacf.acc.image,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

d_brick_deseason_agg_mult <- t(as.matrix(d_brick_deseason_agg))
d_brick_deseason_agg_mult <- d_brick_deseason_agg_mult [,!apply(d_brick_deseason_agg_mult  ,2, anyNA)]
mult.mk.test(ts(d_brick_deseason_agg_mult ))
# Multivariate Mann-Kendall Trend Test
# 
# data:  ts(d_brick_deseason_agg_mult)
# z = -0.064862, p-value = 0.9483
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S       varS 
# -6316 9481939763
#####################################################
#####################################################
#####################################################
# Multivariate Mann-Kendall Trend Test             ##
#                                                  ##
# data:  ts(d_brick_deseason_agg_mult)             ##
# z = -0.089269, p-value = 0.9289                  ##
# alternative hypothesis: true S is not equal to 0 ## using different bandwidth
# sample estimates:                                ##
#   S       varS                                   ##
# -8190 8417264108                                 ##
#####################################################
#####################################################
#####################################################


##########################################################
########################################################## Uni MK
##########################################################

acc_trend <- significantTau(d_brick_deseason_agg,prewhitening=FALSE,p=0.05)
plot(acc_trend)

acc_trend_pixels <- Which(acc_trend,cells=T)
part_acc <- d_brick_deseason_agg[acc_trend_pixels]
matplot(colMeans(part_acc),type = "l")


########################################################## average behavior of not-deseasoned data in the detected pixels
d_brick_agg <- aggregate(d_brick,fac=2)
acc_trend_pixels <- Which(acc_trend,cells=T)
part_1 <- d_brick_agg [acc_trend_pixels]
matplot(colMeans(part_1),type = "l")

#########################################################
#########################################################
######################################################### Crime
#########################################################
#########################################################

plot(crime.ppp[[1]],col=1,pch=20,cex=2,main="")
n.crime <- unlist(lapply(crime.ppp, npoints))
ts.plot(.Last.value[1:12])

#########################################################
######################################################### common bw
#########################################################

sigma_crime <- as.vector(unlist(lapply(crime.ppp, bw.scott.iso)))
sigma_crime <- exp(mean(log(sigma_crime)))

#########################################################
######################################################### Intensity estimate
#########################################################

n_crime <- length(crime.ppp)
d_crime <- lapply(X=1:n_crime, function(i){
  # s <- bw.scott.iso(crime.ppp[[i]])
  # out <- density.lpp(crime.ppp[[i]],sigma=s,distance = "euclidean",positive=TRUE)
  out <- density.lpp(crime.ppp[[i]],sigma=sigma_crime,distance = "euclidean",positive=TRUE)
  raster(as.im(out))
})

d_crime_brick <- brick(d_crime)
names(d_crime_brick) <- names(d_brick)

##########################################################
########################################################## trend detection
##########################################################

d_crime_brick_deseason <- deseason(d_crime_brick)
d_crime_brick_deseason_agg <- aggregate(d_crime_brick_deseason,fac=2)
spplot(d_crime_brick_deseason_agg[[1:12]])
matplot(t(as.data.frame(d_crime_brick_deseason_agg)),type="l",ylab = "")

pacfs_crime <- lapply(X=1:ncell(d_crime_brick_deseason_agg), function(i){
  print(i)
  x <- d_crime_brick_deseason_agg[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_crime)
pacf_coef_crime <- numeric()
for (i in 1:length(pacfs_crime)) {
  if(class(pacfs_crime[[i]])=="numeric"){pacf_coef_crime[i] <- NA}
  else{pacf_coef_crime[i] <- pacfs_crime[[i]]$acf[1]}
}

mean(pacf_coef_crime[!is.na(pacf_coef_crime)])
summary(pacf_coef_crime[!is.na(pacf_coef_crime)])
par(mar=rep(4,4))
boxplot(pacf_coef_crime[!is.na(pacf_coef_crime)])

pacf.image_crime <- d_crime_brick_deseason_agg[[1]]
pacf.image_crime@data@values <- pacf_coef_crime
spplot(pacf.image_crime,scales=list(draw=T,cex=1.4))


##########################################################
########################################################## Multi MK
##########################################################

d_crime_brick_deseason_agg_mult <- t(as.matrix(d_crime_brick_deseason_agg))
d_crime_brick_deseason_agg_mult <- d_crime_brick_deseason_agg_mult [,!apply(d_crime_brick_deseason_agg_mult  ,2, anyNA)]
mult.mk.test(ts(d_crime_brick_deseason_agg_mult))
# Multivariate Mann-Kendall Trend Test
# 
# data:  ts(d_crime_brick_deseason_agg_mult)
# z = 4.6263, p-value = 3.723e-06
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S        varS 
# 555152 14399945040 
#####################################################
#####################################################
#####################################################
# Multivariate Mann-Kendall Trend Test             ##
#                                                  ##
# data:  ts(d_crime_brick_deseason_agg_mult)       ##
# z = 4.2605, p-value = 2.04e-05                   ##
# alternative hypothesis: true S is not equal to 0 ## using different bandwidth
# sample estimates:                                ##
#   S       varS                                   ##
# 467229 12026670384                               ##
#####################################################
#####################################################
#####################################################


##########################################################
########################################################## Uni MK
##########################################################

crime_trend <- significantTau(d_crime_brick_deseason_agg,prewhitening=FALSE,p=0.05)
plot(crime_trend)

# crime_trend1 <- significantTau(d_crime_brick_deseason_agg,prewhitening=FALSE,p=0.01)
# plot(crime_trend1)

crime_trend_pixels <- Which(crime_trend,cells=T)
part_crime <- d_crime_brick_deseason_agg[crime_trend_pixels]
matplot(colMeans(part_crime),type = "l")


########################################################## average behavior of not-deseasoned data in the detected pixels
d_crime_brick_agg <- aggregate(d_crime_brick,fac=2)
crime_trend_pixels <- Which(crime_trend,cells=T)
part_2 <- d_crime_brick_agg[crime_trend_pixels]
matplot(colMeans(part_2),type = "l")



#########################################################
#########################################################
######################################################### relative risk
#########################################################
#########################################################

#########################################################
######################################################### common bw
#########################################################

sigma_ratio <- lapply(X=1:n, function(i){
  x <- as.ppp(acc.ppp[[i]])
  y <- as.ppp(crime.ppp[[i]])
  s <- LSCV.risk(x,y,method = "davies",edge = F,hlim = c(0,2000))
  return(s)
})
sigma_ratio <- unlist(sigma_ratio)
sigma_ratio <- exp(mean(log(sigma_ratio)))

#########################################################
######################################################### risk estimate
#########################################################

r <- lapply(X=1:n, function(i){
  # x <- as.ppp(acc.ppp[[i]])
  # y <- as.ppp(crime.ppp[[i]])
  # s <- LSCV.risk(x,y,method = "davies",edge = F,hlim = c(0,2000))
  # print(s)
  out1 <- density.lpp(acc.ppp[[i]],sigma=sigma_ratio,distance = "euclidean",positive=TRUE)
  out2 <- density.lpp(crime.ppp[[i]],sigma=sigma_ratio,distance = "euclidean",positive=TRUE)
  return(raster(as.im(out2/out1)))
})

r_brick <- brick(r)
spplot(r_brick[[1:12]])

##########################################################
########################################################## trend detection
##########################################################

r_deseason <- deseason(r_brick)
r_deseason_agg <- aggregate(r_deseason,fac=2)
spplot(r_deseason_agg[[1:12]])
matplot(t(as.data.frame(r_deseason_agg)),type="l")

pacfs_r <- lapply(X=1:ncell(r_deseason_agg), function(i){
  print(i)
  x <- r_deseason_agg[i]
  if(!anyNA(x)){return(pacf(x))}
  else{return(0)}
})

length(pacfs_r)
pacf_coef_r <- numeric()
for (i in 1:length(pacfs_r)) {
  if(class(pacfs_r[[i]])=="numeric"){pacf_coef_r[i] <- NA}
  else{pacf_coef_r[i] <- pacfs_r[[i]]$acf[1]}
}

mean(pacf_coef_r[!is.na(pacf_coef_r)])
summary(pacf_coef_r[!is.na(pacf_coef_r)])
par(mar=rep(4,4))
boxplot(pacf_coef_r[!is.na(pacf_coef_r)])

pacf.image_r <- r_deseason_agg[[1]]
pacf.image_r@data@values <- pacf_coef_r
spplot(pacf.image_r,scales=list(draw=T,cex=1.4))

##########################################################
########################################################## Multi MK
##########################################################

r_deseason_agg_mult <- t(as.matrix(r_deseason_agg))
r_deseason_agg_mult <- r_deseason_agg_mult [,!apply(r_deseason_agg_mult  ,2, anyNA)]
mult.mk.test(ts(r_deseason_agg_mult))
# Multivariate Mann-Kendall Trend Test
# 
# data:  ts(r_deseason_agg_mult)
# z = 2.6572, p-value = 0.00788
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S       varS 
# 232288 7642210923 
#####################################################
#####################################################
#####################################################
# Multivariate Mann-Kendall Trend Test             ##
#                                                  ##
# data:  data:  ts(r_deseason_agg_mult)            ##
# z = 2.7506, p-value = 0.005948                   ##
# alternative hypothesis: true S is not equal to 0 ## using different bandwidth
# sample estimates:                                ##
#   S       varS                                   ##
# 237064 7427819051                                ##
#####################################################
#####################################################
#####################################################

##########################################################
########################################################## Uni MK
##########################################################

r_trend <- significantTau(r_deseason_agg,prewhitening=FALSE,p=0.05)
plot(r_trend)

r_trend_pixels <- Which(r_trend,cells=T)
part_r <- d_crime_brick_deseason_agg[r_trend_pixels]
matplot(colMeans(part_r),type = "l")


########################################################## average behavior of not-deseasoned data in the detected pixels
r_brick_agg <- aggregate(r_brick,fac=2)
r_trend_pixels <- Which(r_trend,cells=T)
part_3 <- r_brick_agg[r_trend_pixels]
matplot(colMeans(part_3),type = "l")

##########################################################
##########################################################
##########################################################    All plots
##########################################################
##########################################################

North <- list("SpatialPolygonsRescale", layout.north.arrow(type=1), 
              offset = c(530500, 181500), scale = 800, which = 1)
# scale1 <- list("SpatialPolygonsRescale", layout.scale.bar(), 
#                offset = c(530500, 180000), scale = 500, fill=c("transparent","black"), which = 1)
# s1_text0 <- list("sp.text", c(530500, 180000 + 150), "0", cex = .5, which = 1)
# s1_text1 <- list("sp.text", c(530500 + 500, 180000 + 150), "500 m", cex = .5, which = 1)

png("acc_int.png",width = 1480, height = 880)
spplot(d_brick[[1:12]],sp.layout=list(North),colorkey=list(labels=list(cex=1.5)),
       scales=list(draw=T,cex=1.5),
       col.regions = rev(topo.colors(20, alpha = 1)),
       par.strip.text=list(cex=1.5),
       cex = 1.5
)
dev.off()

png("crime_int.png",width = 1480, height = 880)
spplot(d_crime_brick[[1:12]],sp.layout=list(North),colorkey=list(labels=list(cex=1.5)),
       scales=list(draw=T,cex=1.5),
       col.regions = rev(topo.colors(20, alpha = 1)),
       par.strip.text=list(cex=1.5),
       cex = 1.5
)
dev.off()



pacf_all <- stack(pacf.acc.image,pacf.image_crime,pacf.image_r)
names(pacf_all) <- c("Traffic_Accident","Street_Crime","Relative_Risk")
# spplot(pacf_all,sp.layout=list(list(roads),North,scale1,s1_text0,s1_text1),colorkey=list(labels=list(cex=1.5)),
#        col.regions = rev(terrain.colors(20)),scales=list(draw=T,cex=1.4))


trend_detected <- stack(acc_trend,crime_trend,r_trend)
names(trend_detected) <- c("Traffic_Accident","Street_Crime","Relative_Risk")
# spplot(trend_detected,sp.layout=list(list(roads),North),
#        col.regions = rev(terrain.colors(20)),scales=list(draw=T,cex=1.4))

png("pacf.png",width = 940, height = 340)
spplot(pacf_all,sp.layout=list(list(roads),North),colorkey=list(labels=list(cex=1.5)),
       col.regions = rev(topo.colors(20, alpha = 1)),
       scales=list(draw=T,cex=1.5),
       par.strip.text=list(cex=1.5),
       cex = 1.5)
dev.off()

png("trend.png",width = 940, height = 340)
spplot(trend_detected,sp.layout=list(list(roads),North),colorkey=list(labels=list(cex=1.5)),
       col.regions = rev(topo.colors(20, alpha = 1)),
         scales=list(draw=T,cex=1.5),
       par.strip.text=list(cex=1.5),
       cex = 1.5)
dev.off()


df1 <- data.frame(x=c(1:60),y=as.numeric(colMeans(part_acc)))
df1 <- cbind(df1,ID=rep("Traffic Accident",60))
df2 <- data.frame(x=c(1:60),y=as.numeric(colMeans(part_crime)))
df2 <- cbind(df2,ID=rep("Street Crime",60))
df3 <- data.frame(x=c(1:60),y=as.numeric(colMeans(part_r)))
df3 <- cbind(df3,ID=rep("Relative risk",60))
df <- rbind(df1,df2,df3)

library(ggplot2)
png("trendlines.png",width = 740, height = 540)
ggplot(df,aes(x,y))+geom_line()+facet_wrap(~ID,nrow=1)+
  geom_smooth(method = "loess", se = FALSE,col=2,lwd=1.2)+
  labs(x = "time-period",y="")+
  theme(legend.position="bottom",aspect.ratio=2,axis.text=element_text(size=20),
        axis.title=element_text(size=20),strip.text = element_text(size = 20)
        ,legend.text=element_text(size=rel(3)),legend.title=element_text(size=20))
dev.off()

save.image("results(general).RData")
