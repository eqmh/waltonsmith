rm(list=ls())

library(fields)
library(lubridate)
library(NISTunits)
library(raster)
library(ncdf4)
library(rgdal)


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
# setwd("C:/Users/brendan.turley/Desktop/FL_habs/ne_10m_admin_0_countries")
# setwd("~/Desktop/professional/biblioteca/data/shapefiles/ne_10m_admin_0_countries")
# world <- readOGR('ne_10m_admin_0_countries.shp')
setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
world <- readOGR('GSHHS_h_L1.shp')
world <- crop(world, extent(-84, -79, 24.5, 28))
# setwd("~/Desktop/professional/biblioteca/data/shapefiles/Florida_Shoreline__1_to_40%2C000_Scale_-shp")
# FL <- readOGR('Florida_Shoreline__1_to_40%2C000_Scale_.shp')

### colorpalettes
### breaks and colors
temp_col <- colorRampPalette(c(1,'purple','darkorange','gold'))
sal_col <- colorRampPalette(c('purple4','dodgerblue4','seagreen3','khaki1'))
chl_col <- colorRampPalette(c('honeydew2','darkseagreen3','forestgreen','darkslategrey'))
ox.col1 <- colorRampPalette(c(1,'firebrick4','red'))
ox.col2 <- colorRampPalette(c('darkgoldenrod4','goldenrod2','gold'))
ox.col3 <- colorRampPalette(c('dodgerblue4','deepskyblue2','cadetblue1'))



setwd('~/Downloads')
data <- read.csv('WS22022_Sample_log.csv')
ind <- which(data$Depth!=0)
data2 <- data[ind,]
st_rm <- c('2','MR','9','9.5','10','12','18','21/LK','WS','KW1','KW2')
data3 <- data2[!is.element(data2$Station,st_rm),]
data3$Date.GMT <- mdy(data3$Date.GMT)


adj <- .1
resolution <- .01
ctd.loc <- cbind(lon=data3$Longitude.Decimal,lat=data3$Latitude.Decimal)
loc.grid <- list(lon=seq(min(data3$Longitude.Decimal,na.rm=T)-adj, max(data3$Longitude.Decimal,na.rm=T)+adj,resolution),
                 lat=seq(min(data3$Latitude.Decimal,na.rm=T)-adj, max(data3$Latitude.Decimal,na.rm=T)+adj,resolution))
xlims <- range(loc.grid$lon)
ylims <- range(loc.grid$lat)

### ----------------- Temperature krig -----------------
data3$tempF <- NISTdegCtOdegF(data3$Temperature..C..CTD.data)
my.krig <- spatialProcess(ctd.loc,data3$tempF)
temp_kriged <- predictSurface(my.krig, loc.grid, extrap=F)
temp_se <- predictSurfaceSE(my.krig, loc.grid, extrap=F)

if(max(temp_kriged$z,na.rm=T)>max(data3$tempF,na.rm=T)){
  temp_kriged$z[which(temp_kriged$z>max(data3$tempF,na.rm=T))] <- max(data3$tempF,na.rm=T)
}
if(min(temp_kriged$z,na.rm=T)<min(data3$tempF,na.rm=T)){
  temp_kriged$z[which(temp_kriged$z<min(data3$tempF,na.rm=T))] <- min(data3$tempF,na.rm=T)
}

temp_breaks <- pretty(data3$tempF,n=10)
temp_cols <- temp_col(length(temp_breaks)-1)


### ----------------- Salinity krig -----------------
my.krig <- spatialProcess(ctd.loc,data3$Salinity..PSU..CTD.data)
sal_kriged <- predictSurface(my.krig, loc.grid, extrap=F)
sal_se <- predictSurfaceSE(my.krig, loc.grid, extrap=F)

if(max(sal_kriged$z,na.rm=T)>max(data3$Salinity..PSU..CTD.data,na.rm=T)){
  sal_kriged$z[which(sal_kriged$z>max(data3$Salinity..PSU..CTD.data,na.rm=T))] <- max(data3$Salinity..PSU..CTD.data,na.rm=T)
}
if(min(sal_kriged$z,na.rm=T)<min(data3$Salinity..PSU..CTD.data,na.rm=T)){
  sal_kriged$z[which(sal_kriged$z<min(data3$Salinity..PSU..CTD.data,na.rm=T))] <- min(data3$Salinity..PSU..CTD.data,na.rm=T)
}

# sal_breaks <- pretty(sal_kriged$z,n=10)
sal_breaks <- pretty(data3$Salinity..PSU..CTD.data,n=10)
sal_cols <- sal_col(length(sal_breaks)-1)


### ----------------- DO krig -----------------
my.krig <- spatialProcess(ctd.loc,data3$Oxygen.mg.l..CTD.data)
do_kriged <- predictSurface(my.krig, loc.grid, extrap=F)
do_se <- predictSurfaceSE(my.krig, loc.grid, extrap=F)

if(max(do_kriged$z,na.rm=T)>max(data3$Oxygen.mg.l..CTD.data,na.rm=T)){
  do_kriged$z[which(do_kriged$z>max(data3$Oxygen.mg.l..CTD.data,na.rm=T))] <- max(data3$Oxygen.mg.l..CTD.data,na.rm=T)
}
if(min(do_kriged$z,na.rm=T)<min(data3$Oxygen.mg.l..CTD.data,na.rm=T)){
  do_kriged$z[which(do_kriged$z<min(data3$Oxygen.mg.l..CTD.data,na.rm=T))] <- min(data3$Oxygen.mg.l..CTD.data,na.rm=T)
}

breaks <- pretty(data3$Oxygen.mg.l..CTD.data,n=10)
cols <- c(ox.col1(length(breaks[breaks<2])),
          ox.col2(length(breaks[breaks>=2 & breaks<3.5])),
          ox.col3(length(breaks[breaks>=3.5])-1))


### ----------------- plots -----------------
setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('WS22022_bottom.png', height = 11, width = 4, units = 'in', res=300)
par(mfrow=c(3,1),mar=c(4.5,4,2,1),oma=c(4,1,4,1))
imagePlot(temp_kriged$x,
          temp_kriged$y,
          temp_kriged$z,
          col=temp_cols,breaks=temp_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(temp_cols),legend.width=.7,legend.mar=3)
contour(temp_kriged$x,
        temp_kriged$y,
        temp_kriged$z,
        levels=temp_breaks,add=T)
image(temp_se,add=T,breaks=quantile(temp_se$z,c(.9,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=20,col='green')
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext(expression(paste('Bottom Temperature (',degree,'F)')),adj=1,cex=.75)


imagePlot(sal_kriged$x,
          sal_kriged$y,
          sal_kriged$z,
          col=sal_cols,breaks=sal_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(sal_cols),legend.width=.7,legend.mar=3)
contour(sal_kriged$x,
        sal_kriged$y,
        sal_kriged$z,
        levels=sal_breaks,add=T)
image(sal_se,add=T,breaks=quantile(sal_se$z,c(.9,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=20,col='magenta')
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('Salinity (PSU)',adj=1,cex=.75)

imagePlot(do_kriged$x,
          do_kriged$y,
          do_kriged$z,
          col=cols,breaks=breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(cols),legend.width=.7,legend.mar=3)
contour(do_kriged$x,
        do_kriged$y,
        do_kriged$z,
        levels=breaks,add=T)
image(do_se,add=T,breaks=quantile(do_se$z,c(.9,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-2,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=20,col='orange')
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext(expression(paste('Bottom DO (mg l'^-1,')')),adj=1,cex=.75)

mtext('Walton Smith Bulletin',
      outer=T,line=1,side=3,font=2,at=.05,adj=0,cex=1.25)
mtext(paste('Collected:',
            paste(
              paste(month.abb[month(data3$Date.GMT[1])],
                    day(data3$Date.GMT[1])),
              paste(month.abb[month(data3$Date.GMT[nrow(data3)])],
                    day(data3$Date.GMT[nrow(data3)])),
              sep='-')),
      outer=T,line=-.1,side=3,at=.05,adj=0)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov with concerns \nProcessed: ',as.Date(Sys.time())),
      outer=T,line=2,side=1,col='red',font=2,at=.01,adj=0,cex=.75)
dev.off()


