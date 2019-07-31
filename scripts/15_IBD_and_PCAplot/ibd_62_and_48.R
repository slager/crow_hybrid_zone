library(adegenet)
library(rgdal)
library(sp)
library(MASS)
library(RColorBrewer)
library(magrittr)

### See adegenet basic tutorial.  ISOLATION BY DISTANCE
### ALL SAMPLES


read.structure("Cb.s62.r0.05.maf0.02.het0.5.tsv.unlinked.nomd.str",
               n.ind=62,n.loc=738,onerowperind=F,
               col.lab=1,col.pop=0,col.others=0,
               row.marknames=1,NA.char="0") -> str

#str(str)
#class(str@tab)
#str@tab[1:4,1:4]
#rownames(str@tab)

pop <- substr(rownames(str@tab),1,nchar(rownames(str@tab))-2)
str$pop <- as.factor(pop)

read.csv("crow_geo.csv",header=T,stringsAsFactors=F) -> crowgeo
crowgeo -> coords
coordinates(coords) <- c("long","lat")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
utm <- spTransform(coords, CRS("+proj=utm +zone=9 ellps=WGS84"))
coordinates(utm) -> crowgeo[,c("x","y")] ## x and y are UTMs
crowgeo[,c("x","y")] -> xy
crowgeo[,"loc"] -> rownames(xy)
as.matrix(xy) -> str$other$xy

#toto <- genind2genpop(str,missing="NA")
toto <- genind2genpop(str)
Dgen <- dist.genpop(toto,method=2)
Dgeo <- dist(str$other$xy)/1000 # km
ibd <- mantel.randtest(Dgen,Dgeo,nrepet=1000)  ## Throws error if NAs
ibd

pdf("IBDsim.all62.pdf",6,6)
plot(ibd)
dev.off()

pdf("IBDscat.all62.pdf",6,6)
plot(Dgeo,Dgen,xlab="Geographic distance (km)",ylab="Edwards' genetic distance")
abline(lm(Dgen~Dgeo), col="red",lty=2)
dev.off()
# 
# # KDE version of plot
# dens <- kde2d(Dgeo,Dgen, n=1000)
# myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
# pdf("IBDkernel.all62.pdf",8,5)
# plot(Dgeo, Dgen, pch=20,cex=.5)
# image(dens, col=transp(myPal(300),.7), add=TRUE)
# abline(lm(Dgen~Dgeo))
# title("Isolation by distance plot")
# dev.off()

### ISOLATION BY DISTANCE -- COASTAL PACIFIC ONLY

read.structure("Cb.s62.r0.05.maf0.02.het0.5.str.tsv.coastal.unlinked.nomd.str",
               n.ind=48,n.loc=905,onerowperind=F,
               col.lab=1,col.pop=0,col.others=0,
               row.marknames=1,NA.char="0") -> str

pop <- substr(rownames(str@tab),1,nchar(rownames(str@tab))-2)
#names(pop) <- NULL
str$pop <- as.factor(pop)


read.csv("crow_geo_COASTAL.csv",header=T,stringsAsFactors=F) -> crowgeo
crowgeo -> coords
coordinates(coords) <- c("long","lat")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
utm <- spTransform(coords, CRS("+proj=utm +zone=9 ellps=WGS84"))
coordinates(utm) -> crowgeo[,c("x","y")] ## x and y are UTMs
crowgeo[,c("x","y")] -> xy
crowgeo[,"loc"] -> rownames(xy)
as.matrix(xy) -> str$other$xy

#toto <- genind2genpop(str,missing="NA")
toto <- genind2genpop(str)
Dgen <- dist.genpop(toto,method=2)
Dgeo <- dist(str$other$xy)/1000 #km
ibd <- mantel.randtest(Dgen,Dgeo,nrepet=100000)  ## Throws NA error
ibd

pdf("IBDsim.coastal48.pdf",6,6)
plot(ibd)
dev.off()

pdf("IBDscat.coastal48.pdf",6,6)
plot(Dgeo,Dgen,xlab="Geographic distance (km)",ylab="Edwards' genetic distance")
abline(lm(Dgen~Dgeo), col="red",lty=2)
dev.off()
# 
# # KDE version of plot
# dens <- kde2d(Dgeo,Dgen, n=1000)
# myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
# pdf("IBDkernel.coastal48.pdf",8,5)
# plot(Dgeo, Dgen, pch=20,cex=.5)
# image(dens, col=transp(myPal(300),.7), add=TRUE)
# abline(lm(Dgen~Dgeo))
# title("Isolation by distance plot")
# dev.off()