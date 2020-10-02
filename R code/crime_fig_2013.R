library(dplyr)
library(sp)
library(maptools)
library(spatstat)
library(zoo)



## Data Import
load("london_data.RData")  # "acc.uk" "crime.street" "roads"  
attach(crime.street)
glimpse(crime.street)

## Process to convert to linear point pattern

# Get x and y coords in two vectors
lon <- crime.street$lon
lat <- crime.street$lat

# Create two vectors xrange and yrange with dimensions of triangle that contain all points
xrange <- range(lon, na.rm=T)
yrange <- range(lat, na.rm=T)

# Create ppp
regionN <- ppp(lon, lat, xrange, yrange, checkdup=FALSE)


## Get coordinates of crime sites

min_lon = min(crime.street[,10])
max_lon = max(crime.street[,10])
min_lat = min(crime.street[,11])
max_lat = max(crime.street[,11])

### Get coordinates 
coords <- as.matrix(crime.street[,10:11])


### Get OSM Network of the study region
proj4string(roads) <- CRS("+init=epsg:27700")

# Create linnet for study region network
nw <- as.linnet.SpatialLines(roads,fuse=TRUE,trace=TRUE)
#plot(nw)
nw_sp <- roads

# Transform points to same projection
# EPSG 27700
#proj4string(nw_sp) <- CRS("+init=epsg:27700")



## Density Plot (Region)

crime.street.ppp<-ppp(x=crime.street$lon, y=crime.street$lat, xrange, yrange, checkdup=FALSE)




## Density Plot (Network)


X=crime.street.ppp
crime.nw.lpp=lpp(X,nw)

## Point Patterns 

png("crime_2013.png",width = 1150, height = 670)

par(mfrow = c(3, 4))
par(mar = c(0, 0, 1, 0))

aa <- dplyr::filter(crime.street, year == 2013)
mths <- sort(unique(aa$month))
n.ppp <- list(length(mths))
nw.ppp <- list(length(mths))



n.ppp  <- list(length(mths))
nw.ppp <- list(length(mths))

for (i in seq_along(mths)){
  n.ppp[[i]] <- ppp(x=aa$lon[aa$month==mths[i]],
                    y=aa$lat[aa$month==mths[i]],
                    xrange,yrange, checkdup=FALSE)
  #plot(n.ppp[[i]])
  nw.ppp[[i]] <- lpp(n.ppp[[i]],nw)
  
}
names(nw.ppp) <- as.yearmon(seq(as.Date("2013/1/1"), as.Date("2013/12/31"), "month"))
for (i in seq_along(mths)) {
  plot(nw.ppp[[i]], col = "black", lwd = 0.15, cols = "brown2", cex = 1.7, pch = 19, main = "")
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
