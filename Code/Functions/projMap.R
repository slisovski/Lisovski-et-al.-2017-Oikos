projMap <- function(xlim, ylim, proj="+proj=laea", map.col="grey90", grid=TRUE, grid.lty=2, grid.lwd=1, grid.col="grey80", lab.offset=NULL, asp = 1, font = 1, lab.cex = 1, mar = mar) {
	  require(maptools)
    require(rgdal)
    require(raster)
    require(rgeos)
    if(missing(xlim)) xlim <- c(-175, 175)
    if(missing(ylim)) ylim <- c(-89, 89)
    xcentre <- round(xlim[1] + diff(xlim)/2)
    ycentre <- round(ylim[1] + diff(ylim)/2)
    if (length(grep("+lon_0", proj)) == 0) {
      proj <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", proj, xcentre, ycentre)
    }

	data(wrld_simpl)
	wrld_simpl <- gIntersection(wrld_simpl, as(extent(xlim[1], xlim[2], ylim[1], ylim[2]), "SpatialPolygons"), byid=T)
  mp <- spTransform(wrld_simpl, CRS(proj))
	
    if(grid){	
    gl <- gridlines(wrld_simpl)
gl <- spTransform(gl, CRS(proj))
ga <- gridat(wrld_simpl, side="WS", offset = 0)
ga <- spTransform(ga, CRS(proj))
## bug in gridat should use stringsAsFactors = FALSE
ga$labels <- as.character(ga$labels)
    ga$offset <- ifelse(is.null(lab.offset), ga$offset, rep(lab.offset,(length(ga$offset))))
    }
    
    par(mar=c(3,3,3,3))
    plot(mp, col = "transparent", asp = asp,  xlab = "", ylab = "", axes = FALSE)
    plot(mp, add = TRUE, col = map.col)
    
    plot(spTransform(gl, CRS(proj)), col = grid.col, lty = grid.lty, lwd= grid.lwd, add = TRUE)
text(coordinates(ga), labels=parse(text=ga$labels), pos=ga$pos, offset=ga$offset, xpd=T, font = font, cex = lab.cex)

return(mp)
}