library(maps)
library(rgdal)
library(sp)
library(maptools)
  data(wrld_simpl)
library(cluster)
library(flexclust)
library(mgcv)
library(rgeos)
library(Bolstad2)

source('~/Dropbox/Science/Code/getProj.R', chdir = TRUE)
source("~/Dropbox/Science/Code/projMap.R", chdir = TRUE)

  

sDat <- read.csv("~/Dropbox/Science/Projects/SeasonalityAndAIV/Results_27092016.csv")


xlim <- c(-170, -50)
ylim <- c(5, 80)


## extract coordinates
lon <- apply(cbind(as.character(sDat$Longitude)), 1, FUN=function(x) as.numeric(substring(as.character(x), 1, 6)))
lat <- apply(cbind(as.character(sDat$Latitude)), 1, FUN=function(x) as.numeric(substring(as.character(x), 1, 6)))

map(xlim = xlim, ylim = ylim)
points(lon, lat, pch = 16, col = "red")



## get dates
tmp1 <- strsplit(as.character(sDat$Collection.Date), split = "/")
date <- as.POSIXct(strptime(sapply(tmp1, FUN = function(x) paste(as.numeric(x[3]), as.numeric(x[1]), as.numeric(x[2]), sep = "-")),
                           format = "%y-%m-%d"), tz = "GMT")


## get test results
test <- ifelse(sDat$Flu.Test.Status=="Negative", 0, 1)


## create new table
survTab <- data.frame(institute=sDat$Collector.Institution, date=date, country=sDat$Country,
                      species=sDat$Host.Species, name=sDat$Host.Common.Name, test=test,
                      lon=lon, lat=lat, capture=sDat$Sample.Capture.Status, 
                      type=sDat$Bird.Behavior, age=sDat$Age)

#	table(survTab$type)
#	Captive Wild     Domestic      Unknown         Wild 
# 257              137           666             134015  
survTab <- survTab[survTab$type=="Wild",]


# only data with location
survTab <- survTab[!is.na(survTab$lon) & !is.na(survTab$lat),]


# exclude institutions only providing positive results
tmp1 <- table(survTab[,which(names(survTab)%in%c("institute", "test"))])		  
tmp2 <- row.names(tmp1)[which((tmp1[,2]/tmp1[,1])<=0.95)]
if(length(tmp2)>0) survTab <- subset(survTab, institute%in%tmp2)


## only dubbling ducks
ind <- as.data.frame(table(survTab$species))
# ind <- ind[ind[,2]>50,] 
ind <- ind[order(ind[,2], decreasing = T),]
# unique(ind)


### Dabbling ducks
ducks <- c("Anas platyrhynchos", "Anas discors", "Anas acuta", "Anas clypeata",
           "Anas carolinensis", "Anas americana", "Anas strepera", "Anas crecca",
           "Anas crecca carolinensis","Anas rubripes", "Anas cyanoptera", "Anas sp.",
           "Anas fulvigula", "Anas sp", "Anas [platyrhynchos x sp]", "Anas [platyrhynchos x rubripes]",
           "Anas penelope", "Anas platyrhynchos / Anas rubripes", "Anas platyrhynchos domesticus",
           "Anas [rupripes x sp]", "Anas rubripes/Anas platyrhynchos")


## only dabbling ducks
tab <- subset(survTab, species%in%ducks)


points(tab$lon, tab$lat, pch = 16, col = "green")


# nrow(tab)
# [1] 69233 # 27-September-2016
# hist(as.numeric(format(tab$date, "%Y")), breaks = 1981:2016)
tab <- tab[as.numeric(format(tab$date, "%Y"))%in%c(2005:2013),]
# nrow(tab)



proj <- getProj(range(tab$lon), range(tab$lat))

crds <- project(cbind(tab$lon, tab$lat), proj)
crdsU <- unique(cbind(crds[,1], crds[,2]))


opar <- par(mfrow = c(1,1))
plot(crdsU, pch=20, cex=0.5, col="magenta")
wrldp <- spTransform(wrld_simpl, CRS(proj))
plot(wrldp, add=T)
par(opar)

radius <- 5000

# spatial cluster
cl1 <- qtclust(crds, radius = radius)	


medoids <- aggregate(data.frame(lon=crds[,1], lat=crds[,2]), by=list(clusters(cl1)), FUN = "mean")
names(medoids) <- c("cluster", "lon.med", "lat.med")

med.table <- merge(data.frame(lon=crds[,1], lat=crds[,2], cluster=clusters(cl1)), medoids, by="cluster", all.x=T)

t <- as.data.frame(table(data.frame(sum=rep(1, nrow(med.table)), cluster=med.table$cluster)))

barplot(t(t[,3]))

tab$cluster <- clusters(cl1)

c.lon <- medoids[match(tab$cluster, medoids$cluster), 2]
c.lat <- medoids[match(tab$cluster, medoids$cluster), 3]

center <- project(cbind(c.lon, c.lat), proj, inv = T)

tab$lon.center <- center[,1]
tab$lat.center <- center[,2]



# load('~/Dropbox/Science/Projects/NDVIProjects/CodeChunks/NDVI_raster.RData')
load('~/Dropbox/Science/Projects/SeasonalityAndAIV/Analysis/amp_moll.RData')
load('~/Dropbox/Science/Projects/SeasonalityAndAIV/Analysis/dur_moll.RData')


area <- matrix(c(-180, rep(-45, 100), rep(-180, 100),
                 0, seq(0, 80, length=100), seq(80, 0, length=100)), ncol=2)
areaSP <- SpatialPolygons(list(Polygons(list(Polygon(area)),ID=1)), proj4string=CRS(proj4string(wrld_simpl)))
bbox   <- spTransform(areaSP, CRS(proj4string(amp)))
US <- gIntersection(wrld_simpl, areaSP , byid=T)
US <- spTransform(US, CRS("+proj=moll +lon_0=-105 +lat_0=0"))

tmp <- as.data.frame(table(tab$cluster))
tab <- tab[tab$country!="Guatemala", ]
locs <- project(cbind(tab$lon.center, tab$lat.center), proj4string(amp))

tmp2 <- aggregate(locs, by=list(tab$cluster), mean)

loc.amp <- extract(amp, tmp2[,2:3], buffer=350000, fun=mean, na.rm=T)
loc.dur <- extract(dur, tmp2[,2:3], buffer=350000, fun=mean, na.rm=T) 

## clustering
y.pam <- pam(matrix(c(loc.dur,loc.amp),ncol=2),3, 
             stand=TRUE)
pam.cluster <- y.pam$clustering


## No seasonality = but in closest cluster (3)
pam.cluster[loc.amp<100] <- 2


pam <- data.frame(cluster=tmp2[,1], lon.center=tmp2[,2], lat.center=tmp2[,3],
                  pam = pam.cluster)
tab$pam <- pam$pam[match(tab$cluster, pam$cluster)]
fig <- list(tab, locs, pam.cluster, loc.amp, loc.dur, pam)



cols <- c("#003A7D", "#D2960C", "#7DAF00") 
pch.cl <- approxfun(y = seq(2,5, length = 10), x = seq(50, max(table(fig[[1]]$cluster)), length = 10))

pdf("~/Dropbox/Science/Projects/SeasonalityAndAIV/Figures/Fig02_raw1.pdf", useDingbats=FALSE, width = 12, height = 8)
layout(matrix(c(1,2,3,1,2,3,1,4,5,8,4,5,8,6,7,8,6,7), ncol = 3, byrow = T))

mp <- projMap(c(-168, -45), c(12, 80), proj = "+proj=moll +lon_0=-105 +lat_0=0", map.col="grey95", 
              grid=TRUE, grid.lty=3, grid.lwd=0.8, grid.col="grey60", lab.offset=0.8, asp = 1, font = 3, lab.cex = 1.3,
              mar = c(5,5,4,3))

clL  <- aggregate(fig[[2]], by=list(fig[[1]]$cluster), mean)
n <- aggregate(rep(1, length=nrow(fig[[2]])), by=list(fig[[1]]$cluster), sum)$x


points(clL[rev(order(n)),2:3], pch=21, cex=pch.cl(n)[rev(order(n))],
       bg=cols[fig[[3]][rev(order(n))]], col = "black")


# par(mar=c(4, 10, 0.5, 2), las = 1)
# plot(fig[[5]][rev(order(n))], fig[[4]][rev(order(n))], pch=16, cex=1.4, col="grey90",
#      xlab="", ylab="", xaxt="n", yaxt = "n", type = "n", ylim = c(100, 650), xlim = c(5, 45))
# points(ifelse(amp[]>100, dur[], NA), ifelse(amp[]>100, amp[], NA), pch = 16, cex = 1.5, col = "grey90")     
# mtext("Seasonal amplitude \n(NDVI)", 2, line = 4.6, las = 3, cex = 1.4)
# mtext("Seasonal duration (weeks)", 1, line = 4.5, las = 1, cex = 1.4)
# points(fig[[5]][rev(order(n))], fig[[4]][rev(order(n))], pch=21, cex=pch.cl(n)[rev(order(n))], 
#        bg=cols[fig[[3]][rev(order(n))]], col = "black")
# 
# 
# axis(1, cex.lab = 2, lwd = 1.4, labels = FALSE)
# axis(1, cex.lab = 2, lwd = 0, labels = TRUE, line = 0.7)
# axis(2, cex.lab = 2, lwd = 1.4)
# box(lwd = 1.7)



### NDVI time series
# 
# clusterTab <- data.frame(cluster = tab$pam, lon = tab$lon.center, lat = tab$lat.center)
# ind <- data.frame(year = as.numeric(substring(names(NDVI), 1, 4)), week = as.numeric(substring(names(NDVI), 6, 7)))
# 
# ts <- array(dim = c(nrow(clusterTab), 52, 9))
# 
# for(i in 2005:2012) {
#  for(j in 1:52) {
#    ind01 <- which(ind[,1]==i & ind[,2]==j)[1]
#    if(!is.na(ind01)){
#    tmp01 <- NDVI[[which(ind[,1]==i & ind[,2]==j)[1]]]
#    ts[,j,which(2005:2012==i)] <- extract(tmp01, as.matrix(project(cbind(clusterTab$lon, clusterTab$lat), proj = "+proj=moll +over +ellps=WGS84")), radius = 35000, na.rm = T, fun = mean)
#    }
#  } 
# }

out <- array(dim=c(3,3,3))

par(mar = c(2,2,1,2))

for(i in c(1,3,2)) {
col.t <- col2rgb(cols[i])/256

plot(1,1)
### NDVi

# ndvi1 <- ts[which(clusterTab$cluster==i),,]
# ndvi1[ndvi1<0] <- NA 

# a1 <- apply(ndvi1, c(1,2), FUN = median, na.rm = T)

# l1 <- apply(a1, 2, median, na.rm = T)/1000
# l2 <- apply(a1, 2, function(x) quantile(x, probs = 0.975, na.rm = T))/1000
# l3 <- apply(a1, 2, function(x) quantile(x, probs = 0.025, na.rm = T))/1000

# plot(NA, xlim = c(1,52), ylim = c(0,0.7), xlab = "", ylab = "", yaxt = "n", xaxt = "n")

# polygon(c(c(1:52)[!is.na(l1)], rev(c(1:52)[!is.na(l1)])), c(l2[!is.na(l1)], rev(l3[!is.na(l1)])), 
#        col = "grey", border = NA)

# lines(1:52, l1, type = "l", lwd = 2)
# axis(2, at = seq(0, 0.6, length = 4))

tm1 <- seq(as.POSIXct("2012-01-01"), as.POSIXct("2012-12-01"), by = "month")
wk1 <- as.numeric(format(tm1, "%U"))
lab <- format(tm1, "%b")

# axis(1, at = wk1[c(TRUE, FALSE)], labels = lab[c(TRUE, FALSE)])  


# out[which(c(1,3,2)==i), 1, 3] <- quantile(apply(a1, 1, function(x) abs(diff(range(x, na.rm = T))))/1000, probs = 0.5)
# out[which(c(1,3,2)==i), 2, 3] <- quantile(apply(a1, 1, function(x) abs(diff(range(x, na.rm = T))))/1000, probs = 0.975)
# out[which(c(1,3,2)==i), 3, 3] <- quantile(apply(a1, 1, function(x) abs(diff(range(x, na.rm = T))))/1000, probs = 0.025)



## Infection curves

tab01 <- subset(tab, as.numeric(format(date, "%U"))%in%24:47)

with(tab01[tab01$pam==i,], plot(jitter(as.numeric(format(date, "%U")), 0.8), 
                            jitter(test, 0.2),
                            pch = 16, col = rgb(col.t[1], col.t[2], col.t[3], alpha = 0.1),
                            xlim = c(24, 48), yaxt = "n", xaxt = "n"))

axis(4, at = c(0,1))

tmp <- aggregate(tab01$test[tab01$pam==i], by = list(as.numeric(format(tab01$date[tab01$pam==i], "%U"))),
                 FUN = function(x) length(x))
names(tmp) <- c("Week", "N")
tmp$p <- aggregate(tab01$test[tab01$pam==i], by = list(as.numeric(format(tab01$date[tab01$pam==i], "%U"))),
                   FUN = function(x) sum(x))[,2]


tmp2 <- merge(data.frame(Week = 1:52), tmp, all.x = T)
tmp2$N <- ifelse(tmp2$N<5, NA, tmp2$N)
tmp2$upper <- apply(tmp2, 1, function(x) ifelse(!is.na(x[2]), binom.test(x[3], x[2])$conf.int[1], NA))
tmp2$lower <- apply(tmp2, 1, function(x) ifelse(!is.na(x[2]), binom.test(x[3], x[2])$conf.int[2], NA))


par(new = T)

plot(NA, xlim = c(24, 48), ylim = c(-0.1, 0.6), xaxt = "n", yaxt = "n", bty = "n")

arrows(tmp2$Week, tmp2$upper,
       tmp2$Week, tmp2$lower, length = 0)

points(tmp2$Week,
       tmp2$p/tmp2$N, pch = 22, bg = "white", type = "p")


kk <- tab01[tab01$pam==i,]
kk$Week <- as.numeric(format(kk$date, "%U"))
model <- gam(test~s(Week), data = kk, family = binomial)

x <- seq(24, 47, by = 0.25)

p <- predict(model, newdata=list(Week=x), type = "link", se.fit = TRUE)
upr <- p$fit + (2 * p$se.fit)
lwr <- p$fit - (2 * p$se.fit)

upr <- model$family$linkinv(upr)
lwr <- model$family$linkinv(lwr)

fit <- predict(model,newdata=list(Week=x),type="response")

lines(x, fit,lwd=4, col = rgb(col.t[1], col.t[2], col.t[3], alpha = 0.6))
lines(x, upr, lty = 3)
lines(x, lwr, lty = 3)
axis(2, at = c(0, 0.5))
axis(2, at = c(0.1, 0.2, 0.3, 0.4), labels = FALSE)
  
axis(1, at = wk1[c(TRUE, FALSE)], labels = lab[c(TRUE, FALSE)])

# using sintegral in Bolstad2

out[which(c(1,3,2)==i), 1, 1] <- sintegral(x,fit)$int
out[which(c(1,3,2)==i), 2, 1] <- sintegral(x,upr)$int
out[which(c(1,3,2)==i), 3, 1] <- sintegral(x,lwr)$int


# out[which(c(1,3,2)==i), 1, 2] <- sintegral(1:52,l1)$int
# out[which(c(1,3,2)==i), 2, 2] <- sintegral(1:52,l2)$int
# out[which(c(1,3,2)==i), 3, 2] <- sintegral(1:52,l3)$int


}

out[,,3] <- out[,,3]*10


par(mfrow = c(1,3), mar = c(4,3,4,3), las = 2)

bp <- barplot(out[c(2,3,1),1,c(1)], ylim = c(0,4), col = cols[c(3,1,2)])
arrows(bp, out[c(2,3,1),2,1], bp, out[c(2,3,1),3,1], length = 0, lwd = 1.2)

temp <- data.frame(cl = y.pam$clustering, amp = loc.amp, dur = loc.dur)
nd.out1 <- aggregate(temp[,2:3], by = list(temp[,1]), 
                     FUN = function(x) quantile(x, probs = c(0.5)))
nd.out2 <- aggregate(temp[,2:3], by = list(temp[,1]), 
                     FUN = function(x) quantile(x, probs = c(0.975)))
nd.out3 <- aggregate(temp[,2:3], by = list(temp[,1]), 
                     FUN = function(x) quantile(x, probs = c(0.0275)))


bp <- barplot(nd.out1$amp[c(3,2,1)]/100, ylim = c(0,6))
arrows(bp, nd.out2$amp[c(3,2,1)]/100, bp, nd.out3$amp[c(3,2,1)]/100, length = 0, lwd = 1.2)

bp <- barplot(nd.out1$dur[c(3,2,1)], ylim = c(0, 40))
arrows(bp, nd.out2$dur[c(3,2,1)], bp, nd.out3$dur[c(3,2,1)], length = 0, lwd = 1.2)



axis(4, at = seq(0, 6), labels = seq(0,6)/10, las = 2)

# dev.off()







## Species list
species <- as.data.frame(table(tab$species))
species <- species[species[,2]>0,]
species <- species[rev(order(species[,2])),]

cl1 <- subset(tab, pam == 1)
cl1s <- as.data.frame(table(cl1$species))

species$cl1 <- merge(species, cl1s, by = "Var1", all.x = T, sort = F)[,3]


cl2 <- subset(tab, pam == 2)
cl2s <- as.data.frame(table(cl2$species))

species$cl2 <- merge(species, cl2s, by = "Var1", all.x = T, sort = F)[,4]


cl3 <- subset(tab, pam == 3)
cl3s <- as.data.frame(table(cl3$species))

species$cl3 <- merge(species, cl3s, by = "Var1", all.x = T, sort = F)[,5]

sum <- apply(species[,-1], 2, sum)


figS2 <- aggregate(rep(1, nrow(tab)), by = list(format(tab$date, "%Y-%m"), tab$pam), sum)
  figS2 <- data.frame(year = as.numeric(substring(figS2$Group.1, 1, 4)),
                      month = as.numeric(substring(figS2$Group.1, 6,7)),
                      cluster = figS2$Group.2,
                      n = figS2$x)
  
opar <- par(mfrow = c(3,1), oma = c(1,4,1,1), mar = c(3,3,0,0))

for(i in c(3,1,2)) {
  
  tmp01 <- data.frame(year = rep(2006:2013, each = 12), month = rep(1:12, length(2006:2013)))
    tmp01$n <- merge(tmp01, figS2[figS2$cluster==i,], all.x = T)
  
  bp <- barplot(height = tmp01$n$n, col = cols[i])
  if(i!=2) {
    axis(1, at = bp, labels = NA)
    abline(v =  bp[tmp01$n$month==01], lty = 2, col = "grey90")
    } else {
    axis(1, at = bp, labels = NA)
    axis(1, at = bp[tmp01$n$month==01], labels = tmp01$n$year[tmp01$n$month==01], tick = 2)
  }
}
mtext("Sample Size", 2, outer = T, line = 1.5)  
  


seq(as.POSIXct("2006-01-01"), as.POSIXct("2014-01-01"), by = "month")

