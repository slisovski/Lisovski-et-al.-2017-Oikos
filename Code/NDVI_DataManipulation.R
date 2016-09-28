## NDVI - seasonal length and productivity
library(splancs)
library(rgeos)
library(rgdal)
library(raster)
library(bbmle)
library(maptools)
  data(wrld_simpl)
source("splitYear.R")
source('~/Dropbox/Apps/getProj.R', chdir = TRUE)
load("NDVI_raster.RData")


tmp1 <- NDVI[[1]]
tmp2 <- rasterize(spTransform(wrld_simpl, CRS(proj4string(NDVI[[1]]))), tmp1)
area <- matrix(c(-180, rep(-45, 100), rep(-180, 100),
                 0, seq(0, 80, length=100), seq(80, 0, length=100)), ncol=2)
areaSP <- SpatialPolygons(list(Polygons(list(Polygon(area)),ID=1)), proj4string=CRS(proj4string(wrld_simpl)))
bbox   <- spTransform(areaSP, CRS(proj4string(NDVI[[1]])))

tmp2[] <- 1:length(tmp2[])
ind <- unlist(extract(tmp2, bbox))

## NDVI Date index
index <- matrix(as.numeric(unlist(strsplit(names(NDVI), "_"))), 
                nrow = length(NDVI), ncol = 2, byrow = T)
## NDVI coordinates
crds  <- coordinates(NDVI[[1]]) 


out <- matrix(nrow = nrow(crds), ncol = 5)

for(i in ind) {
  cat(paste(i, "\r    "))
  tmp <- data.frame(year  = rep(2007:2012, each = 52),
                    week = rep(1:52, length(2007:2012)), ndvi = NA)
  for(n in which(index[,1]==2007)[1]:nrow(index)) {  
    tmp$ndvi[which(paste(tmp[,1], tmp[,2], sep = "_")==names(NDVI)[n])] <- values(NDVI[[n]])[i]
  }  
  tmp[,3] <- ifelse(tmp[,3]<(-100), NA, tmp[,3])
  if(sum(!is.na(tmp$ndvi))>150){
    tmp2 <- splitYear(tmp$ndvi, tmp$week, tmp$year)
    tmp3 <- data.frame(week = unlist(lapply(tmp2, function(x) x[,2])),
                       ndvi = unlist(lapply(tmp2, function(x) x[,4])))
    # plot(tmp3[,2]~tmp3[,1])
    ls <- loess(tmp3[,2]~tmp3[,1], span = 0.3)
    pr <- predict(ls, newdata = 1:52)
    prtmp <- predict(loess(tmp3[,2]~tmp3[,1], span = 0.7), newdata = 1:52)
    if(length(which(prtmp>c(prtmp[-1], NA) & prtmp>c(NA, prtmp[-length(prtmp)])))>0){    
    # lines(1:52, pr, lwd=2, col="orange")
    
    tmp4 <- cbind(c(1:52), pr)
      mx <- max(which(tmp4[,2]>c(tmp4[-1,2], NA) & tmp4[,2]>c(NA, tmp4[-nrow(tmp4),2])))
      ct01  <- tmp4[min(which(!is.na(tmp4[,2]))):mx,]
      ind01 <- ct01[which(ct01[,2]<c(NA, ct01[-nrow(ct01),2]) & ct01[,2]<c(ct01[-1,2], NA)),1]
      if(length(ind01)>0) ct1 <- tmp4[ind01[which.min(tmp4[ind01,2])],1] else ct1 <- min(which(!is.na(pr)))
    # abline(v=ct1, lty = 2, lwd = 2, col = "grey20")
    
    
      ct02  <- tmp4[mx:max(which(!is.na(tmp4[,2]))),]
      ind02 <- ct02[which(ct02[,2]<c(NA, ct02[-nrow(ct02),2]) & ct02[,2]<c(ct02[-1,2], NA)),1]
      if(length(ind02)>0) ct2 <- tmp4[ind02[which.min(tmp4[ind02,2])],1] else ct2 <- max(which(!is.na(pr)))
    # abline(v=ct2, lty = 2, lwd = 2, col = "grey20")
    
    
    pr2 <- predict(ls, newdata = ct1:ct2)
    # lines(ct1:ct2, pr2, lwd=4)
    
    ## Data
    ind3 <- which.max(pr2)
    
    ## amp
    out[i,3] <- diff(range(pr2))
    
    ## sd
    out[i,4] <- sd(pr2)
    
    ## Onset
    onc <- approxfun(x = pr2[1:ind3], y = ct1:(ct1 + (ind3-1)))
    out[i,1] <- onc(pr2[1] + ((25*(pr2[ind3]-pr2[1]))/100)) 
    # points(out[i,1], pr2[1] + ((25*(pr2[ind3]-pr2[1]))/100), pch = 16, cex = 3, col = "firebrick")
    
    ## end
    end <- approxfun(x = pr2[ind3:length(pr2)], y = c(ct1:ct2)[ind3:length(pr2)])
    out[i,2] <- end(pr2[length(pr2)] + ((25*(pr2[ind3]-pr2[length(pr2)]))/100))
    # points(out[i,2], pr2[length(pr2)] + ((25*(pr2[ind3]-pr2[length(pr2)]))/100), pch = 16, cex = 3, col = "firebrick")
    
    # area
    if(sum(pr2>=50)>1){
      s <- ifelse((pr2[1] + ((25*(pr2[ind3]-pr2[1]))/100))>50, out[i,1], onc(50))
      e <- ifelse((pr2[length(pr2)] + ((25*(pr2[ind3]-pr2[length(pr2)]))/100))>50, out[i,2], end(50))
      pr3 <- c(50, predict(ls, newdata = seq(s, e, by = 1)), 50)
      x3  <- c(s, seq(s, e, by = 1), e)
      # polygon(x3, pr3, col="red")
      out[i, 5] <- areapl(cbind(x3, pr3))
    }
    }
  } # end extraction 
}

ndviUSphen <- out 
  colnames(ndviUSphen) <- c("start", "end", "max", "sd", "area")
save(ndviUSphen, file="ndviUSphen.RData")
load("ndviUSphen.RData")


  rm(NDVI)
US <- gIntersection(wrld_simpl, areaSP , byid=T)
 load("NDVI_raster.RData")
tt <- spTransform(US, CRS("+proj=moll +lon_0=-105 +lat_0=0"))
  tt1 <- crop(NDVI[[1]], bbox)
  tt1 <- projectRaster(tt1,crs="+proj=moll +lon_0=-105 +lat_0=0")
  tt2 <- rasterize(tt, tt1)


dur <- NDVI[[1]]
  dur[] <- ndviUSphen[,2] - ndviUSphen[,1]
  dur <- crop(dur, bbox)
  dur <- projectRaster(dur, crs="+proj=moll +lon_0=-105 +lat_0=0")
  dur[] <- ifelse(tt2[]>0, dur[], NA)
plot(dur)

amp <- NDVI[[1]]
  amp[] <- ndviUSphen[,3]
  amp <- crop(amp, bbox)
  amp <- projectRaster(amp, crs="+proj=moll +lon_0=-105 +lat_0=0")
  amp[] <- ifelse(tt2[]>0, amp[], NA)
plot(amp)

area <- NDVI[[1]]
  area[] <- ndviUSphen[,3]
  area <- crop(area, bbox)
area <- projectRaster(area, crs="+proj=moll +lon_0=-105 +lat_0=0")
area[] <- ifelse(tt2[]>0, area[], NA)
plot(area)



save(amp, file = "amp_moll.RData")
save(dur, file = "dur_moll.RData")
save(area, file = "area_moll.RData")
