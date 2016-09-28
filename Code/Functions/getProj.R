getProj <- function(xlim, ylim, proj="+proj=laea"){
	
	if(missing(xlim)) xlim <- c(-175, 175)
  if(missing(ylim)) ylim <- c(-89, 89)
    xcentre <- round(xlim[1] + diff(xlim)/2)
    ycentre <- round(ylim[1] + diff(ylim)/2)
  if (length(grep("+lon_0", proj)) == 0) {
    proj <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", proj, xcentre, ycentre)
  }

return(proj)
	
}