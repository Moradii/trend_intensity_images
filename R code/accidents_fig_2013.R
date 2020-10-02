library(dplyr)
library(sp)
library(maptools)
library(spatstat)
library(zoo)



## Data Import
load("london_data.RData")  # "acc.uk" "crime.street" "roads"     
attach(acc.uk)
glimpse(acc.uk)
acc.uk$day_of_week_name <- ordered(acc.uk$day_of_week_name, levels=c("Sunday","Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))
acc.uk$month <- ordered(acc.uk$month, levels=c(1,2,3,4,5,6,7,8,9,10,11,12))



acc.uk$week_end_night<-ifelse((acc.uk$day_of_week==6 & acc.uk$time_slot %in% c(22:23)) | 
                                acc.uk$day_of_week==7 & acc.uk$time_slot %in% c(22:23,0:6) |
                                acc.uk$day_of_week==1 & acc.uk$time_slot %in% c(0:6),"Yes","No") 
attach(acc.uk)



### Get x and y coords in two vectors
### Get coordinates 
coords <- as.matrix(acc.uk[,3:4])


lon <- acc.uk$location_e
lat <- acc.uk$location_n



### Create two vectors xrange and yrange with dimensions of triangle that contain all points

xrange <- range(lon, na.rm=T)
yrange <- range(lat, na.rm=T)


### Get OSM Network of the study region
proj4string(roads) <- CRS("+init=epsg:27700")

# Create linnet for study region network
nw <- as.linnet.SpatialLines(roads,fuse=TRUE,trace=TRUE)

## Density Plot (Region)
acc.uk.ppp<-ppp(x=acc.uk$location_e, y=acc.uk$location_n, xrange, yrange, checkdup=FALSE)

## Density Plot (Network)
X=acc.uk.ppp
acc.nw.lpp=lpp(X,nw)

## Point Patterns 

png("accident_2013.png",width = 1150, height = 670)
par(mfrow = c(3, 4))
par(mar = c(0, 0, 1, 0))
a <- NULL
ln_ppp <- NULL
k <- 1
aa <- dplyr::filter(acc.uk, year == 2013)
mths <- sort(unique(aa$month))
n.ppp <- list(length(mths))
nw.ppp <- list(length(mths))
j<- NULL
for (i in seq_along(mths)) {
  n.ppp[[i]] <- ppp(
    x = aa$location_e[aa$month == mths[i]],
    y = aa$location_n[aa$month == mths[i]],
    xrange, yrange, checkdup = FALSE
  )
  # plot(n.ppp[[i]])
  nw.ppp[[i]] <- lpp(n.ppp[[i]], nw)
}



#a<- st_as_sf(nw.ppp[[1]])
#plot(a)

arrow1 <-  layout.north.arrow(type = 1) 
Narrow1 <- maptools::elide(arrow1)
#Narrow1 <- maptools::elide(arrow1, shift = c(extent(530320.3, 534141.8, 179782.4, 182375.05)))


names(nw.ppp) <- as.yearmon(seq(as.Date("2013/1/1"), as.Date("2013/12/31"), "month"))
for (i in seq_along(mths)) {
  
  plot(nw.ppp[[i]], col = "black", lwd = 0.15, cols = "brown2", cex = 2, pch = 19, main = "")
  title(main = names(nw.ppp[i]), cex.main = 1.45)
  #plot(Narrow1, add = TRUE,col = "black")
  box(which = "figure", lty = "solid")
  
  if(i==1) {
    arrow.x = c(-0.3, 0.3, 0.3, 0.8, 0, -0.8, -0.3, -0.4)
    arrow.y = c( -4.5,    -4.5,  2,    2,  4,   2,   2,      0)
    polygon(530600 + arrow.x * 75, 182100 + arrow.y * 75, col = "black")
    text(x = 530600, y = 182100, labels = "N", cex = 2.4)
    #addnortharrow(scale = 0.3, padin = c(0.40, 0.30), pos = "bottomleft")
  }
}
dev.off()
