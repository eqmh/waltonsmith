library(fields)
library(raster)
library(rerddap)

### these two methods create the same bathymetry that can be used in underway and bottom plotting

################## geographic scope
### west Florida Shelf
lonbox_e <- -79 ### Florida Bay
lonbox_w <- -88 ### mouth of Mississippi River
latbox_n <- 31 ### northern coast
latbox_s <- 24 ### southern edge of Key West

### ERDDAP method
etopo <- info('etopo180')
# define region of interest
latitude = c(latbox_s, latbox_n)
longitude = c(lonbox_w, lonbox_e)
# ERDDAP extract
etopo_grab <- griddap(etopo, latitude=latitude, longitude=longitude, stride=c(1,1) ,fields='altitude')
# reshape data as matrix
topo_lon <- sort(unique(etopo_grab$data$longitude))
topo_lat <- sort(unique(etopo_grab$data$latitude))
topo <- matrix(etopo_grab$data$altitude,
                 length(topo_lon),
                 length(topo_lat))
# plot
imagePlot(topo_lon,
          topo_lat,
          topo,
          asp=1)


### .grd Method download
### .grd file downloaded from: https://www.ngdc.noaa.gov/mgg/global/
setwd('~/Downloads')
etopo_r <- raster('ETOPO1_Bed_c_gdal.grd')
# reshape data as matrix; this is the whole world so it is alot of data
topo2 <- matrix(getValues(etopo_r),
               etopo_r@ncols,
               etopo_r@nrows)
# define lon/lat
lons <- seq(etopo_r@extent@xmin,
            etopo_r@extent@xmax,
            length.out=etopo_r@ncols)
lats <- seq(etopo_r@extent@ymax,
            etopo_r@extent@ymin,
            length.out=etopo_r@nrows)
# define region of interest to extract
# lats
row <- which.min(abs(latbox_s-lats))
row_end <- which.min(abs(latbox_n-lats))
nrows <- length(row:row_end)
#lons
col <- which.min(abs(lonbox_w-lons))
col_end <- which.min(abs(lonbox_e-lons))
ncols <- length(col:col_end)
### extract region of interest
topo2_lon <- lons[col:col_end]
topo2_lat <- lats[row:row_end]
topo2 <- topo2[col:col_end,row:row_end]
# plot
imagePlot(topo2_lon,
          topo2_lat,
          topo2,
          asp=1)
