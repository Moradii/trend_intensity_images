library(dplyr)
library(spatstat)
library(maptools)



## Data Import
load("london_data.RData")  # "acc.uk" "crime.street" "roads"  
attach(crime.street)
glimpse(crime.street)

unique(crime.street$crime_type)
######################################

# "Anti-social behaviour"   "Violent crime"     "Public disorder and weapons"  
# "Shoplifting"   "Drugs"   "Vehicle crime"     "Robbery"   "Theft from the person"       
# "Violence and sexual offences" "Bicycle theft"  "Public order"    
#######################################


########################################
########################################  Anti Social
########################################
########################################

crime.antisoc <- dplyr::filter(crime.street, crime_type == "Anti-social behaviour")


## Process to convert to linear point pattern

# Get x and y coords in two vectors
lon <- crime.antisoc$lon
lat <- crime.antisoc$lat

# Create two vectors xrange and yrange with dimensions of triangle that contain all points
xrange <- range(lon, na.rm=T)
yrange <- range(lat, na.rm=T)

# Create ppp
regionN <- ppp(lon, lat, xrange, yrange, checkdup=FALSE)

## Get coordinates of crime sites
min_lon = min(crime.antisoc[,10])
max_lon = max(crime.antisoc[,10])
min_lat = min(crime.antisoc[,11])
max_lat = max(crime.antisoc[,11])

### Get coordinates 
coords <- as.matrix(crime.antisoc[,10:11])
#xy <- as.matrix(crime.street[,8:9])

## Get OSM Network of the study region
glimpse(roads)

# Create linnet for study region network
nw <- as.linnet.SpatialLines(roads,fuse=TRUE,trace=TRUE)
#plot(nw)
nw_sp <- roads

## LPP objects for each crime type

par(mfrow = c(3, 4))
par(mar = c(0, 0, 1, 0))

yrs <-  sort(unique(crime.antisoc$year))

crime.yrs <- list(length(yrs))
res <- NULL
ln_ppp <- NULL
k=1
for(j in yrs)
{
  crime.yrs <- dplyr::filter(crime.antisoc, year==j)
  mths <- sort(unique(crime.yrs$month))
  n.ppp  <- list(length(mths))
  nw.ppp <- list(length(mths))
  nw.ker <- list(length(mths))
  
  for (i in seq_along(mths)){
    n.ppp[[i]] <- ppp(x=crime.yrs$lon[crime.yrs$month==mths[i]],
                      y=crime.yrs$lat[crime.yrs$month==mths[i]],
                      xrange,yrange, checkdup=FALSE)
    
    nw.ppp[[i]] <- lpp(n.ppp[[i]],nw)
   
    # Create raster stack
    res[k] <- nw.ppp[i]
    k=k+1
  }
  
}
antisoc.lpp <- res


########################################
########################################
########################################  Shoplifting
########################################
########################################

crime.shoplifting <- dplyr::filter(crime.street, crime_type == "Shoplifting")


## Process to convert to linear point pattern

# Get x and y coords in two vectors
lon <- crime.shoplifting$lon
lat <- crime.shoplifting$lat

# Create two vectors xrange and yrange with dimensions of triangle that contain all points
xrange <- range(lon, na.rm=T)
yrange <- range(lat, na.rm=T)

# Create ppp
regionN <- ppp(lon, lat, xrange, yrange, checkdup=FALSE)

## Get coordinates of crime sites
min_lon = min(crime.shoplifting[,10])
max_lon = max(crime.shoplifting[,10])
min_lat = min(crime.shoplifting[,11])
max_lat = max(crime.shoplifting[,11])

### Get coordinates 
coords <- as.matrix(crime.shoplifting[,10:11])
#xy <- as.matrix(crime.street[,8:9])

## Get OSM Network of the study region
glimpse(roads)

# Create linnet for study region network
nw <- as.linnet.SpatialLines(roads,fuse=TRUE,trace=TRUE)
#plot(nw)
nw_sp <- roads

## LPP objects for each crime type

par(mfrow = c(3, 4))
par(mar = c(0, 0, 1, 0))

yrs <-  sort(unique(crime.shoplifting$year))

crime.yrs <- list(length(yrs))
res <- NULL
ln_ppp <- NULL
k=1
for(j in yrs)
{
  crime.yrs <- dplyr::filter(crime.shoplifting, year==j)
  mths <- sort(unique(crime.yrs$month))
  n.ppp  <- list(length(mths))
  nw.ppp <- list(length(mths))
  nw.ker <- list(length(mths))
  
  for (i in seq_along(mths)){
    n.ppp[[i]] <- ppp(x=crime.yrs$lon[crime.yrs$month==mths[i]],
                      y=crime.yrs$lat[crime.yrs$month==mths[i]],
                      xrange,yrange, checkdup=FALSE)
    
    nw.ppp[[i]] <- lpp(n.ppp[[i]],nw)
    
    # Create raster stack
    res[k] <- nw.ppp[i]
    k=k+1
  }
  
}
shoplifting.lpp <- res



########################################
########################################
########################################  Drugs
########################################
########################################

crime.drug <- dplyr::filter(crime.street, crime_type == "Drugs")


## Process to convert to linear point pattern

# Get x and y coords in two vectors
lon <- crime.drug$lon
lat <- crime.drug$lat

# Create two vectors xrange and yrange with dimensions of triangle that contain all points
xrange <- range(lon, na.rm=T)
yrange <- range(lat, na.rm=T)

# Create ppp
regionN <- ppp(lon, lat, xrange, yrange, checkdup=FALSE)

## Get coordinates of crime sites
min_lon = min(crime.drug[,10])
max_lon = max(crime.drug[,10])
min_lat = min(crime.drug[,11])
max_lat = max(crime.drug[,11])

### Get coordinates 
coords <- as.matrix(crime.drug[,10:11])
#xy <- as.matrix(crime.street[,8:9])

## Get OSM Network of the study region
glimpse(roads)

# Create linnet for study region network
nw <- as.linnet.SpatialLines(roads,fuse=TRUE,trace=TRUE)
#plot(nw)
nw_sp <- roads

## LPP objects for each crime type

par(mfrow = c(3, 4))
par(mar = c(0, 0, 1, 0))

yrs <-  sort(unique(crime.drug$year))

crime.yrs <- list(length(yrs))
res <- NULL
ln_ppp <- NULL
k=1
for(j in yrs)
{
  crime.yrs <- dplyr::filter(crime.drug, year==j)
  mths <- sort(unique(crime.yrs$month))
  n.ppp  <- list(length(mths))
  nw.ppp <- list(length(mths))
  nw.ker <- list(length(mths))
  
  for (i in seq_along(mths)){
    n.ppp[[i]] <- ppp(x=crime.yrs$lon[crime.yrs$month==mths[i]],
                      y=crime.yrs$lat[crime.yrs$month==mths[i]],
                      xrange,yrange, checkdup=FALSE)
    
    nw.ppp[[i]] <- lpp(n.ppp[[i]],nw)
    
    # Create raster stack
    res[k] <- nw.ppp[i]
    k=k+1
  }
  
}
drug.lpp <- res





########################################
########################################
########################################  Vehicle Crime
########################################
########################################

crime.vehicle <- dplyr::filter(crime.street, crime_type == "Vehicle crime")


## Process to convert to linear point pattern

# Get x and y coords in two vectors
lon <- crime.vehicle$lon
lat <- crime.vehicle$lat

# Create two vectors xrange and yrange with dimensions of triangle that contain all points
xrange <- range(lon, na.rm=T)
yrange <- range(lat, na.rm=T)

# Create ppp
regionN <- ppp(lon, lat, xrange, yrange, checkdup=FALSE)

## Get coordinates of crime sites
min_lon = min(crime.vehicle[,10])
max_lon = max(crime.vehicle[,10])
min_lat = min(crime.vehicle[,11])
max_lat = max(crime.vehicle[,11])

### Get coordinates 
coords <- as.matrix(crime.vehicle[,10:11])
#xy <- as.matrix(crime.street[,8:9])

## Get OSM Network of the study region
glimpse(roads)

# Create linnet for study region network
nw <- as.linnet.SpatialLines(roads,fuse=TRUE,trace=TRUE)
#plot(nw)
nw_sp <- roads

## LPP objects for each crime type

par(mfrow = c(3, 4))
par(mar = c(0, 0, 1, 0))

yrs <-  sort(unique(crime.vehicle$year))

crime.yrs <- list(length(yrs))
res <- NULL
ln_ppp <- NULL
k=1
for(j in yrs)
{
  crime.yrs <- dplyr::filter(crime.vehicle, year==j)
  mths <- sort(unique(crime.yrs$month))
  n.ppp  <- list(length(mths))
  nw.ppp <- list(length(mths))
  nw.ker <- list(length(mths))
  
  for (i in seq_along(mths)){
    n.ppp[[i]] <- ppp(x=crime.yrs$lon[crime.yrs$month==mths[i]],
                      y=crime.yrs$lat[crime.yrs$month==mths[i]],
                      xrange,yrange, checkdup=FALSE)
    
    nw.ppp[[i]] <- lpp(n.ppp[[i]],nw)
    
    # Create raster stack
    res[k] <- nw.ppp[i]
    k=k+1
  }
  
}
vehicle.lpp <- res


## For LPP Creation
#antisoc.lpp 
#shoplifting.lpp 
#drug.lpp 
#vehicle.lpp 

save(antisoc.lpp, shoplifting.lpp, drug.lpp, vehicle.lpp, file = "crime_type_lpp.RData")




