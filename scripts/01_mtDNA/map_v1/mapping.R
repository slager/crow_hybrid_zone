library(OpenStreetMap)
library(raster)
library(mapplots)
library(rgdal)
library(plotrix)

## CUSTOM FUNCTION DEPENDENCIES:
source("add.pie.ask.asp.R")
source("draw.circle.ask.asp.R")

## FILE DEPENDENCIES:
# pop_haps.csv
# crow_rangemap_good_1feature.shp and associated files

my.projection <- "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=45 +lon_0=-105 +x_0=6200000 +y_0=3000000 +datum=NAD83 +units=m +no_defs  <>"
pie.size <- .03
pop_haps <- read.csv("pop_haps.csv",stringsAsFactors=F)

## Get AMCRNOCR range map
amcrnocr <- readOGR("rangemaps",layer="crow_rangemap_good_1feature")
spTransform(amcrnocr, CRS(my.projection)) -> amcrnocr   ## Reproject

## PNW ONLY MAP (Run first because need some values to draw Inset mini box)

## Latlong bounds for Lambert projection
minlat <- 45
maxlat <- 54
minlong <- -135
maxlong <- -119

## Get basemap - Longish
pnw <- openmap(c(lat= maxlat,lon= minlong),c(lat= minlat,lon= maxlong),minNumTiles=36,type="nps")
pnw <- openproj(pnw,my.projection)
pnwprojbackup <-pnw

## User specified values
pnw.final.w.in <- 1.375
#plot(pnw)
#drawExtent()  # Manually get extent of what to use as Inset
pnw.xmin <- 4560000
pnw.xmax <- 5040000
pnw.ymin <- 3240000
pnw.ymax <- 3930000

## Autocalc'd and some console outputs
pnw$bbox$p2[1]-pnw$bbox$p1[1] -> pnw.bbox.width
pnw$bbox$p1[2]-pnw$bbox$p2[2] -> pnw.bbox.height
pnw.raw.w.in <- pnw.final.w.in/(pnw.xmax-pnw.xmin)*pnw.bbox.width
pnw.raw.h.in <- pnw.raw.w.in*(pnw.bbox.height/pnw.bbox.width)
#left crop inches
pnw.raw.w.in*(pnw.xmin-pnw$bbox$p1[1])/pnw.bbox.width
#right crop inches
pnw.raw.w.in*(pnw$bbox$p2[1]-pnw.xmax)/pnw.bbox.width
#top crop inches
pnw.raw.h.in*(pnw$bbox$p1[2]-pnw.ymax)/pnw.bbox.height
#bottom crop inches
pnw.raw.h.in*(pnw.ymin-pnw$bbox$p2[2])/pnw.bbox.height
#Final width PNW
pnw.final.w.in
#Final height PNW
pnw.final.w.in/(pnw.xmax-pnw.xmin)*(pnw.ymax-pnw.ymin)

pdf("pnw.pdf",pnw.raw.w.in,pnw.raw.h.in)
plot(pnw)
plot(amcrnocr,border=NA,col=rgb(0,0,0,alpha=0.3),add=T)
#abline(v=pnw.xmin)
#abline(v=pnw.xmax)
#abline(h=pnw.ymin)
#abline(h=pnw.ymax)

##Inset pies

subset(pop_haps,pop_haps$plot=='inset') -> a
SpatialPoints(a[,c('x','y')],CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) -> coords
spTransform(coords,my.projection) -> coords
for (i in 1:nrow(a)){
  if (min(c(a[i,'A'],a[i,'N'])) > 0){
    add.pie.ask.asp(asp=1,z=c(a[i,'A'],a[i,'N']),x=coords@coords[i,'x'],y=coords@coords[i,'y'],radius=pie.size/pnw.final.w.in*(pnw.xmax-pnw.xmin)*sqrt(a[i,'A']+a[i,'N']),col=c('red','blue'),init.angle=a[i,'ia'],clockwise=a[i,'cw'],labels="")
  }
  else {
    draw.circle.ask.asp(asp=1,x=coords@coords[i,'x'],y=coords@coords[i,'y'],radius=pie.size/pnw.final.w.in*(pnw.xmax-pnw.xmin)*sqrt(a[i,'A']+a[i,'N']),border='black',col=ifelse(a[i,'A']==0,'blue','red'))
  }
}

dev.off()

## ALL NA MAP

##Latlong Bounds for ALL NA Lambert Projection
minlat <- 17
maxlat <- 70
minlong <- -170
maxlong <- -35

## Download ALL NA basemap (long)
allna <- openmap(c(lat= maxlat,lon= minlong),c(lat= minlat,lon= maxlong),minNumTiles=36,type="nps")
allnabackup <- allna
allnabackup -> allna
allna <- openproj(allna,my.projection)
allna -> allnaprojbackup

## User specified values
allna.final.w.in <- 7
#plot(allna)
#drawExtent()  #To manually get plot lim values right below
# This is the desired extent of the final, 7-inch wide PDF
allna.xmin <-3027801
allna.xmax <-9913361
allna.ymin <-1019102
allna.ymax <-5658132

## Autocalc'd and some console outputs for ALL NA
allna$bbox$p2[1]-allna$bbox$p1[1] -> allna.bbox.width
allna$bbox$p1[2]-allna$bbox$p2[2] -> allna.bbox.height
allna.raw.w.in <- allna.final.w.in/(allna.xmax-allna.xmin)*allna.bbox.width
allna.raw.h.in <- allna.raw.w.in*(allna.bbox.height/allna.bbox.width)
#left crop inches
allna.raw.w.in*(allna.xmin-allna$bbox$p1[1])/allna.bbox.width
#right crop inches
allna.raw.w.in*(allna$bbox$p2[1]-allna.xmax)/allna.bbox.width
#top crop inches
allna.raw.h.in*(allna$bbox$p1[2]-allna.ymax)/allna.bbox.height
#bottom crop inches
allna.raw.h.in*(allna.ymin-allna$bbox$p2[2])/allna.bbox.height
#Final width all.na
allna.final.w.in
#Final height all.na
allna.final.w.in/(allna.xmax-allna.xmin)*(allna.ymax-allna.ymin)

## Actually plot ALL NA map

pdf("allna.pdf",allna.raw.w.in,allna.raw.h.in)
plot(allna)
plot(amcrnocr,border=NA,col=rgb(0,0,0,alpha=0.3),add=T)
#abline(v=plotxlim1)
#abline(v=plotxlim2)
#abline(h=plotylim1)
#abline(h=plotylim2)

# Draw small InsetBox on ALL NA map
rect(pnw.xmin, pnw.ymin, pnw.xmax, pnw.ymax)
## Main
subset(pop_haps,pop_haps$plot=='main') -> a
SpatialPoints(a[,c('x','y')],CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) -> coords
spTransform(coords,my.projection) -> coords
for (i in 1:nrow(a)){
  if (min(c(a[i,'A'],a[i,'N'])) > 0){
    add.pie.ask.asp(asp=1,z=c(a[i,'A'],a[i,'N']),x=coords@coords[i,'x'],y=coords@coords[i,'y'],radius=pie.size/allna.final.w.in*(allna.xmax-allna.xmin)*sqrt(a[i,'A']+a[i,'N']),col=c('red','blue'),init.angle=a[i,'ia'],clockwise=a[i,'cw'],labels="")
  }
  else {
    draw.circle.ask.asp(asp=1,x=coords@coords[i,'x'],y=coords@coords[i,'y'],radius=pie.size/allna.final.w.in*(allna.xmax-allna.xmin)*sqrt(a[i,'A']+a[i,'N']),border='black',col=ifelse(a[i,'A']==0,'blue','red'))
  }
}
##Louisiana (hashed out b/c got ND2s for those now)
#subset(pop_haps,pop_haps$plot=='la') -> a
#SpatialPoints(a[,c('x','y')],CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) -> coords
#spTransform(coords,my.projection) -> coords
#draw.circle.ask.asp(asp=1,x=coords@coords[1,'x'],y=coords@coords[1,'y'],radius=pie.size/allna.final.w.in*(allna.xmax-allna.xmin)*sqrt(2),border='black',col=NA)

dev.off()


## OBSOLETE CODE

#usa <- getData("GADM",country="USA",level=1)
#us <- us[us$NAME_1 %in% c("Washington","Oregon","Alaska"),]
#can <- getData("GADM",country="CAN",level=1)
#mex <- getData("GADM",country="MEX",level=1)
#bhs <- getData("GADM",country="BHS",level=1)


#crop(usa,extent(minlong,maxlong,minlat,maxlat)) -> usa_c
#crop(can,extent(minlong,maxlong,minlat,maxlat)) -> can_c
#crop(mex,extent(minlong,maxlong,minlat,maxlat)) -> mex_c
#crop(bhs,extent(minlong,maxlong,minlat,maxlat)) -> bhs_c

#spTransform(usa_c, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")) -> usa_c   ## OpenStreetMap mercator
#spTransform(can_c, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")) -> can_c   ## OpenStreetMap mercator
#spTransform(mex_c, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")) -> mex_c   ## OpenStreetMap mercator
#spTransform(bhs_c, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")) -> bhs_c   ## OpenStreetMap mercator


#crop(us,extent(par()$usr)) -> us
#plot(usa_c,add=T,lwd=0.5)
#plot(can_c,add=T,lwd=0.5)
#plot(mex_c,add=T,lwd=0.5)
#plot(bhs_c,add=T,lwd=0.5)



#ca.provinces <- canada[canada$NAME_1 %in% provinces,]
#mx.states <- mexico[mexico$NAME_1 %in% 

#us.bbox <- bbox(us)
#ca.bbox <- bbox(ca.provinces)
#xlim <- c(min(us.bbox[1,1],ca.bbox[1,1]),max(us.bbox[1,2],ca.bbox[1,2]))
#ylim <- c(min(us.bbox[2,1],ca.bbox[2,1]),max(us.bbox[2,2],ca.bbox[2,2]))
#plot(ca.provinces, xlim=xlim, ylim=ylim, add=T)

#crop(us,extent(-135,-120,45,56)) -> uscrop
#spTransform(us, CRS("+init=epsg:3857")) -> us   ## OpenStreetMap mercator

#crop(usa,extent(minlong,maxlong,minlat,maxlat)) -> usa_c
#crop(can,extent(minlong,maxlong,minlat,maxlat)) -> can_c
#spTransform(usa_c, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")) -> usa_c   ## OpenStreetMap mercator
#spTransform(can_c, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")) -> can_c   ## OpenStreetMap mercator
#crop(us,extent(par()$usr)) -> us
#plot(usa_c,add=T)
#plot(can_c,add=T)