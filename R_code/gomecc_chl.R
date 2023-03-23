rm(list=ls())

library(fields)
library(lubridate)
library(ncdf4)
library(NISTunits)
library(raster)
library(rgdal)
library(magrittr)

### load bathymetric data for plotting and masking shallow waters
### Can be downloaded from: https://www.ngdc.noaa.gov/mgg/global/
### NOAA NCEI updated this website and the new version is a .grd which can be parsed using the raster package
### there's also an ERDDAP version if your more comfortable: https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.graph
### There is an additional code to extract the ERDDAP version
setwd("G:/My Drive/ongoing_research/pigment_SF")
bathy <- nc_open('etopo1.nc')
topo <- ncvar_get(bathy, 'Band1')
topo_lat <- ncvar_get(bathy, 'lat')
topo_lon <- ncvar_get(bathy, 'lon')
nc_close(bathy)


################## geographic scope
lonbox_e <- -79 ### Florida Bay
lonbox_w <- -86 ### mouth of Mississippi River
latbox_n <- 31 ### northern coast
latbox_s <- 24 ### include the Keys

ind_lat <- which(topo_lat<latbox_n & topo_lat>latbox_s)
ind_lon <- which(topo_lon<lonbox_e & topo_lon>lonbox_w)

topo_lat <- topo_lat[ind_lat]
topo_lon <- topo_lon[ind_lon]
topo <- topo[ind_lon,ind_lat]

### load map
### coastline map can be downloaded from: https://www.ngdc.noaa.gov/mgg/shorelines/
setwd("G:/My Drive/ongoing_research/pigment_SF/GSHHS_shp/h")
world <- readOGR('GSHHS_h_L1.shp')
### crop to make shapefile smaller and use less RAM
world <- crop(world, extent(-84, -79, 24, 28))

### load HPLC data
setwd("G:/My Drive/ongoing_research/pigment_SF")
hplc <- read.csv('HPLC_reports_consolidated.csv')

### transform date format
hplc <- hplc %>% dplyr::mutate(Date2 =ymd(Date), .before=Date)
# cruise <- 'WS16074'

### select pigment vectors, remove bad pigment values, and subset by date range
min_date <- "2016-01-01"
max_date <-  "2022-01-31"

season <- c("Apr", "May", "Jun")

sel_data <- hplc[, c("Original.PI.Sample.Label", "Station", "Gregorian.Month", "Longitude", "Latitude", "Sampling.Depth..meters.", "Date2", "X.Tot_Chl_a.", "X.Fuco.", "X.Zea.", "X.DV_Chl_a.")]
colnames(sel_data) <- c("sample.label", "station","month",  "longitude", "latitude", "depth", "date", "chla", "fuco", "zea", "dvchla")

### FUCO
fuco_data <- sel_data[sel_data$chla > 0, ] ### remove bad pigment values (smaller than zero, e.g. -9999 or -8888)
fuco_data <- fuco_data[fuco_data$fuco > 0, ]
fuco_data_filt <- fuco_data[, 1:9] %>% 
  dplyr::mutate(fuco_ratio = fuco_data$fuco/fuco_data$chla, .after=fuco) %>%
  dplyr::filter(date >= min_date & date <= max_date) %>% dplyr::filter(month %in% season)

### ZEA
zea_data <- sel_data[sel_data$chla > 0, ] ### remove bad pigment values (smaller than zero, e.g. -9999 or -8888)
zea_data <- zea_data[zea_data$fuco > 0, ]
zea_data_filt <- zea_data[, c(1:8, 10)] %>% 
  dplyr::mutate(zea_ratio = zea_data$zea/zea_data$chla, .after=zea) %>%
  dplyr::filter(date >= min_date & date <= max_date) %>% dplyr::filter(month %in% season)

### DVCHLA
dv_data <- sel_data[sel_data$chla > 0, ] ### remove bad pigment values (smaller than zero, e.g. -9999 or -8888)
dv_data <- dv_data[dv_data$dvchla > 0, ]
dv_data_filt <- dv_data[, c(1:8, 11)] %>% 
  dplyr::mutate(dvchla_ratio = dv_data$dvchla/dv_data$chla, .after=dvchla) %>%
  dplyr::filter(date >= min_date & date <= max_date) %>% dplyr::filter(month %in% season)

### generate a single pigment ratio value per station (the interpolation function below doesn't like more than one value per location)
### FUCO
sta_list_fuco <- unique(fuco_data_filt$station)
df_fuco <- matrix(ncol = 6, nrow = length(sta_list_fuco)) # creates an empty data frame
for (i in 1:length(sta_list_fuco)) {
  sta_val <- fuco_data_filt[fuco_data_filt$station == sta_list_fuco[i], ] 
  chla_avg <- mean(sta_val$chla)
  fuco_avg <- mean(sta_val$fuco)
  fuco_ratio_avg <- mean(sta_val$fuco_ratio)
  
  df_fuco[i,] <- cbind(sta_list_fuco[i], sta_val$longitude[1], sta_val$latitude[1], chla_avg, fuco_avg, fuco_ratio_avg) # concatenates rows
}

df_fuco <- as.data.frame(df_fuco)
colnames(df_fuco) <- c("station", "longitude", "latitude", "chla", "fuco", "fuco_ratio") # change column names
cols.num_fuco <- c("longitude", "latitude", "chla", "fuco", "fuco_ratio") # make these columns numeric vectors (double-check with "sapply(df, class)")
df_fuco[cols.num_fuco] <- sapply(df_fuco[cols.num_fuco],as.numeric)

### ZEA
sta_list_zea <- unique(zea_data_filt$station)
df_zea <- matrix(ncol = 6, nrow = length(sta_list_zea)) # creates an empty data frame
for (i in 1:length(sta_list_zea)) {
  sta_val <- zea_data_filt[zea_data_filt$station == sta_list_zea[i], ] 
  chla_avg <- mean(sta_val$chla)
  zea_avg <- mean(sta_val$zea)
  zea_ratio_avg <- mean(sta_val$zea_ratio)
  
  df_zea[i,] <- cbind(sta_list_zea[i], sta_val$longitude[1], sta_val$latitude[1], chla_avg, zea_avg, zea_ratio_avg) # concatenates rows
}

df_zea <- as.data.frame(df_zea)
colnames(df_zea) <- c("station", "longitude", "latitude", "chla", "zea", "zea_ratio") # change column names
cols.num_zea <- c("longitude", "latitude", "chla", "zea", "zea_ratio") # make these columns numeric vectors (double-check with "sapply(df, class)")
df_zea[cols.num_zea] <- sapply(df_zea[cols.num_zea],as.numeric)

### DVChla
sta_list_dv <- unique(dv_data_filt$station)
df_dv <- matrix(ncol = 6, nrow = length(sta_list_dv)) # creates an empty data frame
for (i in 1:length(sta_list_dv)) {
  sta_val <- dv_data_filt[dv_data_filt$station == sta_list_dv[i], ] 
  chla_avg <- mean(sta_val$chla)
  dv_avg <- mean(sta_val$dvchla)
  dv_ratio_avg <- mean(sta_val$dvchla_ratio)
  
  df_dv[i,] <- cbind(sta_list_dv[i], sta_val$longitude[1], sta_val$latitude[1], chla_avg, dv_avg, dv_ratio_avg) # concatenates rows
}

df_dv <- as.data.frame(df_dv)
colnames(df_dv) <- c("station", "longitude", "latitude", "chla", "dvchla", "dvchla_ratio") # change column names
cols.num_dv <- c("longitude", "latitude", "chla", "dvchla", "dvchla_ratio") # make these columns numeric vectors (double-check with "sapply(df, class)")
df_dv[cols.num_dv] <- sapply(df_dv[cols.num_dv],as.numeric)  

### colorpalettes
### colors are modified versions of recommended oceanography color schemes: https://cran.r-project.org/web/packages/cmocean/vignettes/cmocean.html
chl_col <- colorRampPalette(c('honeydew2','darkseagreen3','forestgreen','darkslategrey'))


#### -------------- kriging FUCO --------------
### krigging resolution
resolution <- .01
adj <- .1
### krigging locations
df_fuco.loc <- cbind(lon=df_fuco$longitude,lat=df_fuco$latitude)
loc.grid_fuco <- list(lon=seq(min(df_fuco$longitude, na.rm=T) - adj, max(df_fuco$longitude, na.rm=T) + adj, resolution),
                 lat=seq(min(df_fuco$latitude, na.rm=T) - adj, max(df_fuco$latitude, na.rm=T) + adj, resolution))
### limits for plotting
xlims_fuco <- range(loc.grid_fuco$lon)
ylims_fuco <- range(loc.grid_fuco$lat)


### ----------------- Pigment FUCO ratio krig -----------------
my.krig_fuco <- spatialProcess(df_fuco.loc, df_fuco$fuco_ratio)
ratio_kriged_fuco <- predictSurface(my.krig_fuco, loc.grid_fuco, extrap=F)
ratio_SE_fuco <- predictSurfaceSE(my.krig_fuco, loc.grid_fuco, extrap=F)

### corrects the krigged min and max to observed min and max
if(max(ratio_kriged_fuco$z, na.rm=T) > max(df_fuco$fuco_ratio, na.rm=T)){
  ratio_kriged_fuco$z[which(ratio_kriged_fuco$z > max(df_fuco$fuco_ratio, na.rm=T))] <- max(df_fuco$fuco_ratio, na.rm=T)
}
if(min(ratio_kriged_fuco$z, na.rm=T) < min(df_fuco$fuco_ratio, na.rm=T)){
  ratio_kriged_fuco$z[which(ratio_kriged_fuco$z < min(df_fuco$fuco_ratio, na.rm=T))] <- min(df_fuco$fuco_ratio, na.rm=T)
}

### color and contour breaks FUCO
breaks_fuco <- pretty(df_fuco$fuco_ratio, n=10)
breaks_fuco <- seq(0, 0.6, by=.05) ### use for fixed breaks.
cols_fuco <- chl_col(length(breaks_fuco)-1)

#### -------------- kriging ZEA --------------
### krigging resolution
resolution <- .01
adj <- .1
### krigging locations
df_zea.loc <- cbind(lon=df_zea$longitude,lat=df_zea$latitude)
loc.grid_zea <- list(lon=seq(min(df_zea$longitude, na.rm=T) - adj, max(df_zea$longitude, na.rm=T) + adj, resolution),
                      lat=seq(min(df_zea$latitude, na.rm=T) - adj, max(df_zea$latitude, na.rm=T) + adj, resolution))
### limits for plotting
xlims_zea <- range(loc.grid_zea$lon)
ylims_zea <- range(loc.grid_zea$lat)


### ----------------- Pigment ZEA ratio krig -----------------
my.krig_zea <- spatialProcess(df_zea.loc, df_zea$zea_ratio)
ratio_kriged_zea <- predictSurface(my.krig_zea, loc.grid_zea, extrap=F)
ratio_SE_zea <- predictSurfaceSE(my.krig_zea, loc.grid_zea, extrap=F)

### corrects the krigged min and max to observed min and max
if(max(ratio_kriged_zea$z, na.rm=T) > max(df_zea$zea_ratio, na.rm=T)){
  ratio_kriged_zea$z[which(ratio_kriged_zea$z > max(df_zea$zea_ratio, na.rm=T))] <- max(df_zea$zea_ratio, na.rm=T)
}
if(min(ratio_kriged_zea$z, na.rm=T) < min(df_zea$zea_ratio, na.rm=T)){
  ratio_kriged_zea$z[which(ratio_kriged_zea$z < min(df_zea$zea_ratio, na.rm=T))] <- min(df_zea$zea_ratio, na.rm=T)
}

### color and contour breaks ZEA
breaks_zea <- pretty(df_zea$zea_ratio, n=10)
breaks_zea <- seq(0, 0.6, by=.05) ### use for fixed breaks.
cols_zea <- chl_col(length(breaks_zea)-1)


#### -------------- kriging DVCHLA --------------
### krigging resolution
resolution <- .01
adj <- .1
### krigging locations
df_dv.loc <- cbind(lon=df_dv$longitude,lat=df_dv$latitude)
loc.grid_dv <- list(lon=seq(min(df_dv$longitude, na.rm=T) - adj, max(df_dv$longitude, na.rm=T) + adj, resolution),
                     lat=seq(min(df_dv$latitude, na.rm=T) - adj, max(df_dv$latitude, na.rm=T) + adj, resolution))
### limits for plotting
xlims_dv <- range(loc.grid_dv$lon)
ylims_dv <- range(loc.grid_dv$lat)


### ----------------- Pigment ZEA ratio krig -----------------
my.krig_dv <- spatialProcess(df_dv.loc, df_dv$dvchla_ratio)
ratio_kriged_dv <- predictSurface(my.krig_dv, loc.grid_dv, extrap=F)
ratio_SE_dv <- predictSurfaceSE(my.krig_dv, loc.grid_dv, extrap=F)

### corrects the krigged min and max to observed min and max
if(max(ratio_kriged_dv$z, na.rm=T) > max(df_dv$dvchla_ratio, na.rm=T)){
  ratio_kriged_dv$z[which(ratio_kriged_dv$z > max(df_dv$dvchla_ratio, na.rm=T))] <- max(df_dv$dvchla_ratio, na.rm=T)
}
if(min(ratio_kriged_dv$z, na.rm=T) < min(df_dv$dvchla_ratio, na.rm=T)){
  ratio_kriged_dv$z[which(ratio_kriged_dv$z < min(df_dv$dvchla_ratio, na.rm=T))] <- min(df_dv$dvchla_ratio, na.rm=T)
}

### color and contour breaks ZEA
breaks_dv <- pretty(df_dv$dvchla_ratio, n=10)
breaks_dv <- seq(0, 0.4, by=.05) ### use for fixed breaks.
cols_dv <- chl_col(length(breaks_dv)-1)

### ----------------- plot ALL -----------------
setwd('~/waltonsmith/figures')
png('MBON_pigment_ratios_spring.png', height = 11, width = 4, units = 'in', res=300)

### FUCO ###
par(mfrow=c(3,1),mar=c(4.5,4,2,1),oma=c(4,1,4,1))
imagePlot(ratio_kriged_fuco$x,
          ratio_kriged_fuco$y,
          ratio_kriged_fuco$z,
          col=cols_fuco,breaks=breaks_fuco,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims_fuco,ylim=ylims_fuco,
          nlevel=length(cols_fuco),legend.width=1,legend.mar=3)
contour(ratio_kriged_fuco$x,
        ratio_kriged_fuco$y,
        ratio_kriged_fuco$z,
        levels=breaks_fuco,add=T)
image(ratio_SE_fuco,add=T,breaks=quantile(ratio_SE_fuco$z,c(1,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
# plot(FL,col='gray70',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(df_fuco$longitude,df_fuco$latitude,pch=20,col='yellow')
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('Fucoxanthin:TChla: spring',adj=1,cex=.75)

### ZEA ###
imagePlot(ratio_kriged_zea$x,
          ratio_kriged_zea$y,
          ratio_kriged_zea$z,
          col=cols_zea,breaks=breaks_zea,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims_zea,ylim=ylims_zea,
          nlevel=length(cols_zea),legend.width=1,legend.mar=3)
contour(ratio_kriged_zea$x,
        ratio_kriged_zea$y,
        ratio_kriged_zea$z,
        levels=breaks_zea,add=T)
image(ratio_SE_zea,add=T,breaks=quantile(ratio_SE_zea$z,c(1,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
# plot(FL,col='gray70',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(df_zea$longitude,df_zea$latitude,pch=20,col='yellow')
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('Zeaxanthin:TChla: spring',adj=1,cex=.75)

### DVCHLA ###
imagePlot(ratio_kriged_dv$x,
          ratio_kriged_dv$y,
          ratio_kriged_dv$z,
          col=cols_dv,breaks=breaks_dv,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims_dv,ylim=ylims_dv,
          nlevel=length(cols_dv),legend.width=1,legend.mar=3)
contour(ratio_kriged_dv$x,
        ratio_kriged_dv$y,
        ratio_kriged_dv$z,
        levels=breaks_dv,add=T)
image(ratio_SE_dv,add=T,breaks=quantile(ratio_SE_dv$z,c(1,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
# plot(FL,col='gray70',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(df_dv$longitude,df_dv$latitude,pch=20,col='yellow')
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('DV_Chla:TChla: spring',adj=1,cex=.75)

mtext('MBON Pigment Ratios Report',
      outer=T,line=1,side=3,font=2,at=.05,adj=0,cex=1.25)
mtext(paste('Sampling period:',
            paste(
              paste(month.abb[month(hplc$Date2[1])],
                    day(hplc$Date2[1]), year(hplc$Date2[1])),
              paste(month.abb[month(hplc$Date2[nrow(hplc)])],
                    day(hplc$Date2[nrow(hplc)]), year(hplc$Date2[nrow(hplc)])),
              sep='-')),
      outer=T,line=-.1,side=3,at=.05,adj=0)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact enrique.montes@noaa.gov with concerns \nProcessed: ',as.Date(Sys.time())),
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