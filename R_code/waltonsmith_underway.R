rm(list=ls())

library(fields)
library(lubridate)
library(ncdf4)
library(NISTunits)
library(raster)
library(rgdal)

### load bathymetric data for plotting and masking shallow waters
### Can be downloaded from: https://www.ngdc.noaa.gov/mgg/global/
### there's also an ERDDAP version if your more comfortable: https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.graph
setwd("~/Desktop/professional/biblioteca/data")
bathy <- nc_open('etopo1.nc')
topo <- ncvar_get(bathy, 'Band1')
topo_lat <- ncvar_get(bathy, 'lat')
topo_lon <- ncvar_get(bathy, 'lon')
nc_close(bathy)


################## geographic scope
lonbox_e <- -79 ### Florida Bay
lonbox_w <- -84 ### mouth of Mississippi River
latbox_n <- 38 ### northern coast
latbox_s <- 24 ### remove the Keys

ind_lat <- which(topo_lat<latbox_n & topo_lat>latbox_s)
ind_lon <- which(topo_lon<lonbox_e & topo_lon>lonbox_w)

topo_lat <- topo_lat[ind_lat]
topo_lon <- topo_lon[ind_lon]
topo <- topo[ind_lon,ind_lat]

### load map
### coastline map can be downloaded from: https://www.ngdc.noaa.gov/mgg/shorelines/
setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
world <- readOGR('GSHHS_h_L1.shp')
### crop to make shapefile smaller and use less RAM
world <- crop(world, extent(-84, -79, 24.5, 28))


### load underway data
setwd('~/Downloads/')
underway <- read.csv('WS22072_out.csv')
### cruise name for file naming
cruise <- 'WS22072'
### make date times
underway$time <- dmy_hms(underway$time)
### rename out of convenience
underway$lon <- underway$lon.dd
underway$lat <- underway$lat.dd
### copy to plot tracklines later
orig <- underway

### remove NA values; krigging hates NAs (but who doesn't)
underway <- underway[-which(is.na(underway$chl.c3p.raw)),]
### remove low salinity and when ship not moving very fast; salinity threshold is arbitrary can plot it out to see if this value makes sense
underway$sal.tsg[which(underway$sal.tsg<20)] <- NA
underway <- underway[which(underway$sog>2),]

### the convolution filters smooths the data; excessively noisy data can either slow down the algorithm and make hard to interpret krigged surfaces; smoothing span is arbitrary and may be superfluous depending upon individual cruises
### chlorophyll
underway$chl <- filter(underway$chl.c3p.raw,rep(1/5,5),'convolution')
# plot(underway$chl.c3p.raw,typ='l')
# points(underway$chl,col=2,typ='l')
# plot(underway$chl.c3p.raw,underway$chl)

### scaling with 1 as the max; not used anymore
# underway$chl <- underway$chl/max(underway$chl,na.rm=T)

### salinity
underway$sal <- filter(underway$sal.tsg,rep(1/5,5),'convolution')
# plot(underway$sal.tsg,typ='l')
# points(underway$sal,col=2,typ='l')
# plot(underway$sal.tsg,underway$sal)

### temperature
underway$tempC <- filter(underway$temp.c3p,rep(1/5,5),'convolution')
# plot(underway$temp.c3p,typ='l')
# points(underway$tempC,col=2,typ='l')
# plot(underway$temp.c3p,underway$tempC)

### convert from Celsius to Fahrenheit
underway$temp <- NISTdegCtOdegF(underway$tempC)


### subsample to make krigging go faster; also subsample span arbitrary
underway <- underway[seq(1,nrow(underway),8),]


### colorpalettes
### colors are modified versions of recommended oceanography color schemes: https://cran.r-project.org/web/packages/cmocean/vignettes/cmocean.html
temp_col <- colorRampPalette(c('gray20','purple','darkorange','gold'))
sal_col <- colorRampPalette(c('purple4','dodgerblue4','seagreen3','khaki1'))
chl_col <- colorRampPalette(c('honeydew2','darkseagreen3','forestgreen','darkslategrey'))


#### -------------- kriging --------------
### krigging resolution
resolution <- .01
### krigging locations
underway.loc <- cbind(lon=underway$lon,lat=underway$lat)
loc.grid <- list(lon=seq(min(underway$lon,na.rm=T)-adj, max(underway$lon,na.rm=T)+adj,resolution),
                 lat=seq(min(underway$lat,na.rm=T)-adj, max(underway$lat,na.rm=T)+adj,resolution))

### ----------------- Chlorophyll krig -----------------
my.krig <- spatialProcess(underway.loc, underway$chl)
chl_kriged <- predictSurface(my.krig, loc.grid, extrap=TRUE)
chl_SE <- predictSurfaceSE(my.krig, loc.grid, extrap=TRUE)

### corrects the krigged min and max to observed min and max
if(max(chl_kriged$z,na.rm=T)>max(underway$chl,na.rm=T)){
  chl_kriged$z[which(chl_kriged$z>max(underway$chl,na.rm=T))] <- max(underway$chl,na.rm=T)
}
if(min(chl_kriged$z,na.rm=T)<min(underway$chl,na.rm=T)){
  chl_kriged$z[which(chl_kriged$z<min(underway$chl,na.rm=T))] <- min(underway$chl,na.rm=T)
}

### color and contour breaks
chl_breaks <- pretty(underway$chl,n=15)
chl_cols <- chl_col(length(chl_breaks)-1)


### ----------------- Salinity krig -----------------
my.krig <- spatialProcess(underway.loc, underway$sal)
sal_kriged <- predictSurface(my.krig, loc.grid, extrap=TRUE)
# sal_SE <- predictSurfaceSE(my.krig, loc.grid, extrap=TRUE)

if(max(sal_kriged$z,na.rm=T)>max(underway$sal,na.rm=T)){
  sal_kriged$z[which(sal_kriged$z>max(underway$sal,na.rm=T))] <- max(underway$sal,na.rm=T)
}
if(min(sal_kriged$z,na.rm=T)<min(underway$sal,na.rm=T)){
  sal_kriged$z[which(sal_kriged$z<min(underway$sal,na.rm=T))] <- min(underway$sal,na.rm=T)
}

sal_breaks <- pretty(underway$sal,n=15)
sal_cols <- sal_col(length(sal_breaks)-1)


### ----------------- Temperature krig -----------------
my.krig <- spatialProcess(underway.loc, underway$temp)
temp_kriged <- predictSurface(my.krig, loc.grid, extrap=TRUE)
# temp_SE <- predictSurfaceSE(my.krig, loc.grid, extrap=TRUE)

if(max(temp_kriged$z,na.rm=T)>max(underway$temp,na.rm=T)){
  temp_kriged$z[which(temp_kriged$z>max(underway$temp,na.rm=T))] <- max(underway$temp,na.rm=T)
}
if(min(temp_kriged$z,na.rm=T)<min(underway$temp,na.rm=T)){
  temp_kriged$z[which(temp_kriged$z<min(underway$temp,na.rm=T))] <- min(underway$temp,na.rm=T)
}

temp_breaks <- pretty(underway$temp,n=15)
temp_cols <- temp_col(length(temp_breaks)-1)


### ----------------- plots -----------------
### limits for plotting
adj <- .1
xlims <- c(min(orig$lon,na.rm=T)-adj, max(orig$lon,na.rm=T)+adj)
ylims <- c(min(orig$lat,na.rm=T)-adj, max(orig$lat,na.rm=T)+adj)

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png(paste0(cruise,'_underway.png'), height = 11, width = 4, units = 'in', res=300)
par(mfrow=c(3,1),mar=c(4.5,4,2,1),oma=c(4,1,4,1))
imagePlot(temp_kriged$x,
          temp_kriged$y,
          temp_kriged$z,
          col=temp_cols,breaks=temp_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(cols),legend.width=.7,legend.mar=3)
contour(temp_kriged$x,
        temp_kriged$y,
        temp_kriged$z,
        levels=temp_breaks,add=T)
# image(temp_SE,add=T,breaks=quantile(temp_SE$z,c(.6,1),na.rm=T),col='white')
image(chl_SE,add=T,breaks=quantile(chl_SE$z,c(.3,1),na.rm=T),col='white') # chlorophyll standard error used to mask surface with SE greater than arbitrary threshold; can be adjusted for individual cruises
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T) # mask surface values that are in less than 1 meter depth
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
lines(orig$lon,orig$lat,col='green',lwd=1) ### plot track line
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext(expression(paste('Surface Temperature (',degree,'F)')),adj=1,cex=.75)


imagePlot(sal_kriged$x,
          sal_kriged$y,
          sal_kriged$z,
          col=sal_cols,breaks=sal_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(cols),legend.width=.7,legend.mar=3)
contour(sal_kriged$x,
        sal_kriged$y,
        sal_kriged$z,
        levels=sal_breaks,add=T)
# image(sal_SE,add=T,breaks=quantile(sal_SE$z,c(.4,1),na.rm=T),col='white')
image(chl_SE,add=T,breaks=quantile(chl_SE$z,c(.3,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
lines(orig$lon,orig$lat,col='magenta',lwd=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('Surface Salinity (PSU)',adj=1,cex=.75)


imagePlot(chl_kriged$x,
          chl_kriged$y,
          chl_kriged$z,
          col=chl_cols,breaks=chl_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(cols),legend.width=.7,legend.mar=3)
contour(chl_kriged$x,
        chl_kriged$y,
        chl_kriged$z,
        levels=chl_breaks,add=T)
image(chl_SE,add=T,breaks=quantile(chl_SE$z,c(.3,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
# plot(FL,col='gray70',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
lines(orig$lon,orig$lat,col='purple',lwd=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('Surface Chlorophyll Fluorescence (RFU)',adj=1,cex=.75)

mtext('Walton Smith Bulletin',
      outer=T,line=1,side=3,font=2,at=.05,adj=0,cex=1.25)
mtext(paste('Collected:',
            paste(
              paste(month.abb[month(underway$time[1])],
                    day(underway$time[1])),
              paste(month.abb[month(underway$time[nrow(underway)])],
                    day(underway$time[nrow(underway)])),
              sep='-')),
      outer=T,line=-.1,side=3,at=.05,adj=0)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov with concerns \nProcessed: ',as.Date(Sys.time())),
      outer=T,line=2,side=1,col='red',font=2,at=.01,adj=0,cex=.75)
dev.off()


### Code written using:
# sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.6.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] rgdal_1.5-29      ncdf4_1.19        raster_3.5-15     sp_1.4-6          NISTunits_1.0.1   lubridate_1.8.0   fields_13.3      
# [8] viridis_0.6.2     viridisLite_0.4.0 spam_2.8-0       
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.8.3     pillar_1.7.0     compiler_4.1.3   tools_4.1.3      dotCall64_1.0-1  lifecycle_1.0.1  tibble_3.1.6    
# [8] gtable_0.3.0     lattice_0.20-45  pkgconfig_2.0.3  rlang_1.0.2      cli_3.2.0        gridExtra_2.3    terra_1.5-21    
# [15] dplyr_1.0.8      rgeos_0.5-9      generics_0.1.2   vctrs_0.3.8      maps_3.4.0       grid_4.1.3       tidyselect_1.1.2
# [22] glue_1.6.2       R6_2.5.1         fansi_1.0.2      ggplot2_3.3.5    purrr_0.3.4      magrittr_2.1.2   scales_1.1.1    
# [29] codetools_0.2-18 ellipsis_0.3.2   colorspace_2.0-3 utf8_1.2.2       munsell_0.5.0    crayon_1.5.0 