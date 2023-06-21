#Author: Penelope A. Hancock
library(INLA)
library(raster)
library(rgdal)
source("INLA marker frequency.R")
source("invelogit2.r")
#Set up a regular grid over Uganda
Uganda <- readOGR(dsn = ".", layer = "Uganda.shp")
stepsize <- (1/111)
nxy <- c(535,394)
rast_grid<-raster(xmn=29.70466,xmx=34.52284,ymn=-1.213042,ymx=2.339036,nrow=nxy[2],ncol=nxy[1])
ex<-extent(rast_grid)
pred_locs<-xyFromCell(rast_grid,1:ncell(rast_grid))

#Set up a window for plotting the maps
dev.new(noRStudioGD = T,width=12,height=8)
par(mfrow=c(3,5))
par(mar=c(3,3,3,3))
title_vec=list("round 1","round 2","round 3","round 4","round 5")

#Run the INLA model and plot the maps for each marker and round
for (marker in 1:3){
  i=1
  if (marker==1){
    data1<-read.csv("kdr_frequenciesF_gamb2.csv")
    res<-inla_marker_frequency(data1)
    bin.res<-res$bin.res
    mesh<-res$mesh
    mesh1d<-res$mesh1d
    breaks1=seq(0,0.5,length=21)
    }
  if (marker==2){
    data1<-read.csv("Cyp6p4_236M_frequencies_gamb2.csv")
    res<-inla_marker_frequency(data1)
    bin.res<-res$bin.res
    mesh<-res$mesh
    mesh1d<-res$mesh1d
    breaks1=seq(0.85,1,length=21)
  }
  if (marker==3){
    data1<-read.csv("Cyp4j5_43F_frequencies_gamb2.csv")
    res<-inla_marker_frequency(data1)
    bin.res<-res$bin.res
    mesh<-res$mesh
    mesh1d<-res$mesh1d
    breaks1=seq(0.3,1,length=21)
  }
  for (round in seq(1,5,1)){
    A.pred<-inla.spde.make.A(mesh,loc=cbind(pred_locs[,'x'],pred_locs[,'y']),group=rep(round,nrow(pred_locs)),group.mesh=mesh1d)
    lpmean<-bin.res$summary.fixed$mean+ drop(A.pred%*%bin.res$summary.random$one.field$mean)
    P<-rast_grid
    P[]<-invelogit(lpmean,100)
    P_mask<-mask(P,Uganda)
    plot(P_mask,col=topo.colors(21),breaks=breaks1,legend=FALSE,asp=1)
    plot(P,legend.only=TRUE,col=topo.colors(21),breaks=breaks1,axis.args=list(at=breaks1))
    title(title_vec[[i]]);i=i+1;
  }
}

