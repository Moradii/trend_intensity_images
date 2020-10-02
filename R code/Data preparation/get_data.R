## Load libraries
library(dplyr)
library(osmdata)
library(maptools)
library(spatstat)


setwd("~/Point Pattern Analysis/Accident Crime")

############################################
## Import Accident Data (2013-2017)
acc.uk <-read.csv("Accident/Acc Data/acc_uk.csv")
#acc.uk$day_of_week_name <- ordered(acc.uk$day_of_week_name, levels=c("Sunday","Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))
#acc.uk$month <- ordered(acc.uk$month, levels=c(1,2,3,4,5,6,7,8,9,10,11,12))
glimpse(acc.uk)
attach(acc.uk)


### Get coordinates 
#coords.a <- as.matrix(acc.uk[,5:6])
coords.a <- as.matrix(acc.uk[,3:4])
# Get x and y coords in two vectors
lon.a <- acc.uk$location_e
lat.a <- acc.uk$location_n

# Create two vectors xrange and yrange with dimensions of triangle that contain all points
xrange.a <- range(lon.a, na.rm=T)
yrange.a <- range(lat.a, na.rm=T)


#############################################
## Import Crime Data (2013-2017)
crime.street <- read.csv("Crime/Crime Data/crime_uk_new.csv")
glimpse(crime.street)
attach(crime.street)

### Get coordinates 
coords.c <- as.matrix(crime.street[,4:5])

# Get x and y coords in two vectors
lon.c <- crime.street$lon
lat.c <- crime.street$lat

# Create two vectors xrange and yrange with dimensions of triangle that contain all points
xrange.c <- range(lon.c, na.rm=T)
yrange.c <- range(lat.c, na.rm=T)



################################################
## Import OSM Network

# Extent of Bounding Box for OSM road network
min_lon = min(acc.uk[,5], crime.street[,4])
max_lon = max(acc.uk[,5], crime.street[,4])
min_lat = min(acc.uk[,6], crime.street[,5])
max_lat = max(acc.uk[,6], crime.street[,5])

h <- opq(bbox=c(min_lon,min_lat,max_lon,max_lat))%>%
  add_osm_feature(key = 'highway', 
                  value = c('motorway', 'trunk','primary','secondary','tertiary', 
                            'residential', 'motorway_link',  'trunk_link',
                            'primary_link', 'secondary_link','tertiary_link', 
                            'pedestrian', 'track', 'cycleway', 'road','living_street' 
                            )
                  ) 

h_sp <- osmdata_sp(h)
road <- h_sp$osm_lines
# Export OSM network as shapefile
# writeOGR(road, "osm", "osm_london", driver="ESRI Shapefile")


##########################################################
## Create linnet for study region network

# Access OSM network
roads<-readShapeSpatial("Network/osm/osm_nw_proj")
plot(roads)

#Create linnet
nw_ln <- as.linnet.SpatialLines(roads,fuse=TRUE,trace=TRUE)


##########################################################
## Accident Point pattern on linear network (Accident LPP)

yrs <-  sort(unique(acc.uk$year))
acc.yrs <- list(length(yrs))
acc.lpp <- NULL
ln_ppp <- NULL
k=1
for(j in yrs)
{
  acc.yrs <- dplyr::filter(acc.uk, year==j)
  mths <- sort(unique(acc.yrs$month))
  n.ppp  <- list(length(mths))
  nw.ppp <- list(length(mths))
  nw.ker <- list(length(mths))
  nw.im <- list(length(mths))
  nw.ras <- list(length(mths))
  for (i in seq_along(mths)){
    n.ppp[[i]] <- ppp(x=acc.yrs$location_e[acc.yrs$month==mths[i]],
                      y=acc.yrs$location_n[acc.yrs$month==mths[i]],
                      xrange.a,yrange.a, checkdup=FALSE)

    nw.ppp[[i]] <- lpp(n.ppp[[i]],nw_ln)
    acc.lpp[k] <- nw.ppp[i]
    k=k+1
  }
  
}



##########################################################
## Crime Point pattern on linear network (Crime LPP)

yrs <-  sort(unique(crime.street$year))
crime.yrs <- list(length(yrs))
crime.lpp <- NULL
k=1
for(j in yrs)
{
  crime.yrs <- dplyr::filter(crime.street, year==j)
  mths <- sort(unique(crime.yrs$month))
  n.ppp  <- list(length(mths))
  nw.ppp <- list(length(mths))
  nw.ker <- list(length(mths))
  nw.im <- list(length(mths))
  nw.ras <- list(length(mths))
  for (i in seq_along(mths)){
    n.ppp[[i]] <- ppp(x=crime.yrs$lon[crime.yrs$month==mths[i]],
                      y=crime.yrs$lat[crime.yrs$month==mths[i]],
                      xrange.c,yrange.c, checkdup=FALSE)
    
    
    nw.ppp[[i]] <- lpp(n.ppp[[i]],nw_ln)
    crime.lpp[k] <- nw.ppp[i]
    k=k+1
  }
  
}

##############################################
# Save LPP and Network as Rdata file

save(list = c("acc.lpp", "crime.lpp", "roads"), file = "london_df.Rdata")
