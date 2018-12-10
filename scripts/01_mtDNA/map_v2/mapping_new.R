#library(OpenStreetMap)
library(raster)
#library(mapplots)
library(rgdal)
#library(plotrix)
library(maps)
library(sp)
library(raster) 
library(RColorBrewer)
library(rgeos)



## CUSTOM FUNCTION DEPENDENCIES for Pies:
source("add.pie.ask.asp.R")
source("draw.circle.ask.asp.R")

## Transparent color function
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  invisible(t.col)
}

t_percent <- 0 # transparency
pie.size <- .03
pop_haps <- read.csv("pop_haps.csv",stringsAsFactors=F)



## Get Shapefiles
amcr <- readOGR("rangemaps/amcr_only",layer="amcr_only")
nocr <- readOGR("rangemaps/nocr_only_final",layer="nocr_only_final")
overlap_zone <- readOGR("rangemaps/overlap_zone_final",layer="overlap_zone_final")

#can <- readOGR("gis_new/lpr_000b16a_e",layer="lpr_000b16a_e")
can <- readOGR("gis_new/province",layer="province")
#can <- readRDS("gadm36_CAN_1_sp.rds")
usa <- readOGR("gis_new/gz_2010_us_040_00_500k",layer="gz_2010_us_040_00_500k")
#usa <- readRDS("gadm36_USA_1_sp.rds")

#usa <- map("state", fill = TRUE)
#IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
#usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

## Simplify shapefiles


#can <-gSimplify(can,tol=0.1, topologyPreserve=T)
#> England3 = SpatialPolygonsDataFrame(England2, data=England@data)

## Define UTM projection
my.projection <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

## Manual encode UTM-zone 10 boundaries for 2/3 aspect ratio Western US map
x0 <- -446329
x1 <- 1394510
y0 <- 4010916
y1 <- 6772176

## Reproject from lat-long to UTM
amcr <- spTransform(amcr,my.projection)
nocr <- spTransform(nocr,my.projection)
overlap_zone <- spTransform(overlap_zone,my.projection)
can <- spTransform(can,my.projection)
usa <- spTransform(usa,my.projection)

## Crop to desired extent
amcr <- crop(amcr, extent(x0, x1, y0, y1))
nocr <- crop(nocr, extent(x0, x1, y0, y1))
overlap_zone <- crop(overlap_zone, extent(x0, x1, y0, y1))
can <- crop(can, extent(x0, x1, y0, y1))
usa <- crop(usa, extent(x0, x1, y0, y1))

'#FDC086' -> orange
'#BEAED4' -> purple
'#7FC97F' -> green

'#CC3311' -> red
'#0077BB' -> blue

#png("trial.png",w=13/3,h=6.5,units="in",res=1000)
pdf("crow_map_final.pdf",13/3,6.5)
par(mar=rep(0,4))
plot(1, type="n", xlab="", ylab="",xaxs = "i", yaxs = "i",xlim=c(x0, x1), ylim=c(y0, y1))
plot(nocr,add=T,lty=0,col=t_col(green,t_percent))
plot(amcr,add=T,lty=0,col=t_col(orange,t_percent))
plot(overlap_zone,add=T,lty=0,col=t_col(purple,t_percent))
plot(can,add=T)
plot(usa,add=T)

#dev.off()

## User specified values
pnw.final.w.in <- 13/3  ## Hard-coded for now
pnw.xmin <- x0  ## Lazy redefinition of variable names for old code
pnw.xmax <- x1
pnw.ymin <- y0
pnw.ymax <- y1

## Autocalc'd and some console outputs
x1-x0 -> pnw.bbox.width
y1-y0 -> pnw.bbox.height
pnw.raw.w.in <- pnw.final.w.in/(pnw.xmax-pnw.xmin)*pnw.bbox.width
pnw.raw.h.in <- pnw.raw.w.in*(pnw.bbox.height/pnw.bbox.width)
#left crop inches
#pnw.raw.w.in*(pnw.xmin-pnw$bbox$p1[1])/pnw.bbox.width
#right crop inches
#pnw.raw.w.in*(pnw$bbox$p2[1]-pnw.xmax)/pnw.bbox.width
#top crop inches
#pnw.raw.h.in*(pnw$bbox$p1[2]-pnw.ymax)/pnw.bbox.height
#bottom crop inches
#pnw.raw.h.in*(pnw.ymin-pnw$bbox$p2[2])/pnw.bbox.height
#Final width PNW
#pnw.final.w.in
#Final height PNW
#pnw.final.w.in/(pnw.xmax-pnw.xmin)*(pnw.ymax-pnw.ymin)

#pdf("pnw.pdf",pnw.raw.w.in,pnw.raw.h.in)
#plot(pnw)
#plot(amcrnocr,border=NA,col=rgb(0,0,0,alpha=0.3),add=T)
#abline(v=pnw.xmin)
#abline(v=pnw.xmax)
#abline(h=pnw.ymin)
#abline(h=pnw.ymax)

##Inset pies

#subset(pop_haps,pop_haps$plot=='inset') -> a
pop_haps -> a
SpatialPoints(a[,c('x','y')],CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) -> coords
spTransform(coords,my.projection) -> coords
for (i in 1:nrow(a)){
  if (min(c(a[i,'A'],a[i,'N'])) > 0){
    add.pie.ask.asp(asp=1,z=c(a[i,'A'],a[i,'N']),x=coords@coords[i,'x'],y=coords@coords[i,'y'],radius=pie.size/pnw.final.w.in*(pnw.xmax-pnw.xmin)*sqrt(a[i,'A']+a[i,'N']),col=c(red,blue),init.angle=a[i,'ia'],clockwise=a[i,'cw'],labels="")
  }
  else {
    draw.circle.ask.asp(asp=1,x=coords@coords[i,'x'],y=coords@coords[i,'y'],radius=pie.size/pnw.final.w.in*(pnw.xmax-pnw.xmin)*sqrt(a[i,'A']+a[i,'N']),border='black',col=ifelse(a[i,'A']==0,blue,red))
  }
}

dev.off()
