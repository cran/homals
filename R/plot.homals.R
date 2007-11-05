`plot.homals` <-
function(x, plot.dim = c(1,2), plot.type, var.subset, main, type, xlab, ylab, 
         xlim, ylim, leg.pos = "topright", ...)
{
#S3 plot method for objects of class "homals"
#Produces various 2D-plots
#plot.dim ... vector of length 2 with dimensions to be plotted against
#plot.type ... type of plot to be drawn: "catplot","graphplot","hullplot","labplot",
#              "lossplot","objplot","prjplot","spanplot","starplot", "trfplot", 
#              "vecplot","vorplot", "jointplot","loadplot".
#var.subset ... numeric vector with subset of variables




#plot.type <- plot.type[1]          #use first plot-type only

options(locatorBell = FALSE)
if (x$ndim == 1) stop("No plots can be drawn for ndim = 1 !")
if (length(plot.dim) !=  2) stop("plot.dim must be of length 2!")
if (plot.type != "trfplot") {      #plot.dim are ignored for trfplot
  pd1 <- plot.dim[1]
  pd2 <- plot.dim[2]
  if (pd2 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")
  if (missing(xlab)) xlab <- paste("Dimension",pd1)
  if (missing(ylab)) ylab <- paste("Dimension",pd2)
}

nvar <- dim(x$dframe)[2]
if (missing(var.subset)) var.subset <- 1:nvar
   

#----------------------------------loadplot-------------------------------------
if (plot.type == "loadplot") {
  xycoor <- t(sapply(x$cat.loadings, function(xy) xy[1,c(pd1,pd2)]))
  if (missing(main)) main1 <- "Loadings plot" else main1 <- main
  xlim.min <- min(xycoor[,1],0)
  xlim.max <- max(xycoor[,1],0)
  ylim.min <- min(xycoor[,2],0)
  ylim.max <- max(xycoor[,2],0)
  if (missing(xlim)) xlim <- c(xlim.min,xlim.max)*1.2
  if (missing(ylim)) ylim <- c(ylim.min,ylim.max)*1.2
  plot(xycoor,type = "p", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main1,...)

  for (i in 1:nvar) lines(rbind(xycoor[i,],c(0,0)),...)
  identify(xycoor, labels = rownames(xycoor))
}
#-------------------------------- end loadplot ---------------------------------




#----------------------------------catplot--------------------------------------
#plots the rank-restricted category quantifications for each variable

if (plot.type == "catplot") {

  if (missing(type)) type <- "b"
  if (missing(xlim)) xlim <- range(sapply(x$rank.cat, function(zz) range(zz[,pd1])))
  if (missing(ylim)) ylim <- range(sapply(x$rank.cat, function(zz) range(zz[,pd2])))

  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Category plot for",colnames(x$dframe[i]))  else main1 <- main

    par("ask"=TRUE)
    plot(x$rank.cat[[i]][,c(pd1,pd2)], type = type, xlim = xlim, ylim = ylim,
    main = main1, xlab = xlab, ylab = ylab, ...)
    text(x$rank.cat[[i]][,c(pd1,pd2)], levels(x$dframe[,i]), pos = 3, ...)
    abline(h = 0, v = 0, col = "gray", lty = 2, ...)
  }
}
#----------------------------------end catplot----------------------------------

#------------------------------------jointplot----------------------------------
if (plot.type == "jointplot") {
  xylist <- lapply(x$rank.cat, apply, 2, range)         
  xytab <- sapply(xylist, function(yy) yy[,c(pd1,pd2)])
  xmin <- min(xytab[1,], x$scores[,pd1])
  xmax <- max(xytab[2,], x$scores[,pd1])
  ymin <- min(xytab[3,], x$scores[,pd2])
  ymax <- max(xytab[4,], x$scores[,pd2])
    
  if (missing(xlim)) xlim <- c(xmin, xmax)
  if (missing(ylim)) ylim <- c(ymin, ymax)
  
  if (missing(main)) main <- "Joint Plot"
  plot(x$scores[,c(pd1,pd2)], type = "n", main = main, xlab = xlab, ylab = ylab, 
  xlim = xlim, ylim = ylim, ...)           #draw scores
  text(x$scores[,c(pd1,pd2)], labels = rownames(x$dframe), col = 1)
  catcol <- rainbow(ncol(x$dframe))
  
  catleg <- NULL
  for (j in var.subset)
  {
    text(x$rank.cat[[j]][,c(pd1,pd2)], labels = rownames(x$rank.cat[[j]]), col = catcol[j])
    catleg <- c(catleg, catcol[j])
  }
  legend(leg.pos,colnames(x$dframe)[var.subset], col = catleg, pch = 22)
}

#----------------------------------end jointplot--------------------------------

#------------------------------------graphplot----------------------------------
if (plot.type == "graphplot") {

  if (missing(main)) main <- "Graphplot"
  if (missing(xlim)) xlim <- range(x$scores[,pd1])*1.2
  if (missing(ylim)) ylim <- range(x$scores[,pd2])*1.2
  plot(x$scores[,c(pd1,pd2)], col = "GREEN", pch = 8, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,...)           #draw scores

  dmat <- NULL
  for (j in 1:ncol(x$dframe))
  {
    y <- computeY(x$dframe[,j], x$scores[,c(pd1,pd2)])
    dmat <- rbind(dmat, y)
    points(y, col = "RED", pch = 16)                             #insert points
    for (i in 1:nrow(x$dframe))
      lines(rbind(x$scores[i,c(pd1,pd2)], y[x$dframe[i,j],]))                #insert lines
  }
  repvec <- sapply(x$rank.cat, function(yy) dim(yy)[1])
  varnames <- rep(colnames(x$dframe),repvec)
  rownames(dmat) <- paste(varnames, rownames(dmat))
  xycoor <- rbind(dmat, x$scores[,c(pd1,pd2)])
  identify(xycoor, labels = rownames(xycoor))
}

#----------------------------------end graphplot--------------------------------

#------------------------------------hullplot-----------------------------------
#plots the convex hulls
if (plot.type == "hullplot") {

  for (i in var.subset) {
    
    if (missing(main)) main1 <- paste("Hullplot for",colnames(x$dframe[i]))  else main1 <- main
    
    par("ask" = TRUE) 
    plot(x$scores[,c(pd1,pd2)], col = "GREEN", pch = 8, main = main1, xlab = xlab, ylab = ylab, ...)
       
    for (j in levels(x$dframe[,i])) 
    {
      ind <- which(j==x$dframe[,i])                      #object index for convex hulls
  	  lst <- ind[chull(x$scores[ind,c(pd1,pd2)])]                  #convex hull over ind
  	  lines(x$scores[c(lst,lst[1]),c(pd1,pd2)], ...)
  	  text(x$scores[lst,c(pd1,pd2)], j, ...)
    }
  }
}
#-----------------------------------end hullplot--------------------------------

#--------------------------------------labplot----------------------------------
#plot labeled object scores (for each variable separately)

if (plot.type == "labplot") {

  if (missing(xlim)) xlim <- range(x$scores[,pd1])
  if (missing(ylim)) ylim <- range(x$scores[,pd2])
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Labplot for",colnames(x$dframe[i]))  else main1 <- main
    par("ask" = TRUE)
    plot(x$scores[,c(pd1,pd2)], type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main1, ...)
    text(x$scores[,c(pd1,pd2)], as.vector(x$dframe[,i]), ...)
  }
} 
#-----------------------------------end labplot---------------------------------

#------------------------------------lossplot-----------------------------------
if (plot.type == "lossplot") {
  
  for (i in var.subset) { 
    
    if (missing(main)) main1 <- paste("Lossplot for",colnames(x$dframe[i])) else main1 <- main
    
    z <- computeY(x$dframe[,i], x$scores[,c(pd1,pd2)])
    k <- dim(z)[1]
    
    if (missing(xlim)) xlim1 <- range(c(z[,2],x$rank.cat[[i]][,pd1])) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(c(z[,2],x$rank.cat[[i]][,pd2])) else ylim1 <- ylim
    
    par("ask" = TRUE)
    plot(x$rank.cat[[i]][,c(pd1,pd2)], type = "p", main = main1, xlab = xlab, ylab = ylab,
         xlim = xlim1, ylim = ylim1, col = "RED", ...)
    text(x$rank.cat[[i]][,c(pd1,pd2)], levels(x$dframe[,i]), pos = 3, col = "RED", ...)
    lines(x$rank.cat[[i]][,c(pd1,pd2)], col = "RED")
    
    points(z,type="p", col = "blue")
    text(z, levels(x$dframe[,i]), col="blue", pos = 3, ...)
    lines(z, col="blue")
    for (j in 1:k) lines(rbind(x$rank.cat[[i]][j,c(pd1,pd2)],z[j,]),col="lightgray", lty=3)
    
    abline(h = 0, v = 0, col = "gray", lty = 2, ...)
  }  
}
#----------------------------------end lossplot---------------------------------

#------------------------------------objplot------------------------------------
#draws labeled object score plot

if (plot.type == "objplot") {
    
  if (missing(xlim)) xlim <- range(x$scores[,pd1])
  if (missing(ylim)) ylim <- range(x$scores[,pd2])
  if (missing(main)) main1 <- "Object score plot" else main1 <- main

  plot(x$scores[,c(pd1,pd2)], type = "n", main = main1, xlab = xlab, ylab = ylab, 
       xlim = xlim, ylim = ylim, ...)
  text(x$scores[,c(pd1,pd2)], rownames(x$scores), ...)
} 

#---------------------------------end objplot-----------------------------------

#----------------------------------prjplots-------------------------------------
#draws projection plot

if (plot.type == "prjplot") {

  if (missing(xlim)) xlim <- range(x$scores[,pd1])
  if (missing(ylim)) ylim <- range(x$scores[,pd2])
  xylim <- c(min(xlim[1],ylim[1]),max(xlim[2],ylim[2]))
  
  for (i in var.subset) {
  
    if (missing(main)) main1 <- paste("Projection plot for", colnames(x$dframe[i])) else main1 <- main
     
    a <- x$cat.loadings[[i]][c(pd1,pd2)]
    
    if (x$rank.vec[i] == 1) {
      par("ask" = TRUE)
      plot(x$scores[,c(pd1,pd2)], type = "p", main = main1, xlab = xlab, ylab = ylab,
           xlim = xylim, ylim = xylim, ...)
      text(x$scores[,c(pd1,pd2)], as.vector(x$dframe[,i]), pos = 3) 
      text(x$rank.cat[[i]], levels(x$dframe[,i]), col = "RED", pos = 3)
      slope = a[2]/a[1]
      abline(coef = c(0,slope))
      slope = -a[1]/a[2]
      
      for (j in 1:(dim(x$rank.cat[[i]])[1]))
  	    abline(coef = c(x$rank.cat[[i]][j,pd2] - slope*x$rank.cat[[i]][j,pd1],slope), col = "RED")
  	  for (k in 1:(dim(x$scores)[1])) {
  	    j <- x$dframe[,i][k]
        icpt <- x$rank.cat[[i]][j,pd2] - slope*x$rank.cat[[i]][j,pd1] 
  	    u <-(x$scores[k,pd1] + slope*(x$scores[k,pd2]-icpt))/(1+slope^2)  
        lines(rbind(x$scores[k,c(pd1,pd2)],c(u,icpt+slope*u)), col = "BLUE")
  	  }
      abline(h = 0, v = 0, col = "gray", lty = 2, ...)   
    } else {
      warning("No projection plot for ", colnames(x$dframe[i]),"(rank != 1).", call. = FALSE) 
    }
  }
}
#---------------------------------end prjplot-----------------------------------

#-------------------------------------spanplot----------------------------------
if (plot.type == "spanplot") {
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Span plot for", colnames(x$dframe[i])) else main1 <- main
    
    par("ask" = TRUE)
    plot(x$scores[,c(pd1,pd2)], col = 1, pch = 21, main = main1, xlab = xlab, ylab = ylab)
    lev <- levels(x$dframe[,i])
    rb<-rainbow(length(lev))
    for (k in lev) {
  	  ind <- which(k==x$dframe[,i])
  	  n <- length(ind)
  	  mm <- mst(dist(x$scores[ind,c(pd1,pd2)]))
  	  for (j in 1:n) {
  		  jnd <- which(1 == as.vector(mm[j,]))
  		  sapply(jnd, function(r) lines(rbind(x$scores[ind[j],c(pd1,pd2)], x$scores[ind[r],c(pd1,pd2)]),
        col = rb[which(lev==k)]))
  		}
  	 legend(leg.pos,paste("Category",lev), col = rb, lty = 1, ...)
    } 
  }
}
#----------------------------------end spanplot---------------------------------


#------------------------------------starplot-----------------------------------
if (plot.type == "starplot") {
  
  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Star plot for", colnames(x$dframe[i])) else main1 <- main
    
    par("ask" = TRUE)
    plot(x$scores[,c(pd1,pd2)], col = "BLUE", main = main1, xlab = xlab, ylab = ylab, ...)
    
    z <- computeY(x$dframe[,i], x$scores[,c(pd1,pd2)])
    points(z, type="o", pch = 24, col = "RED")
    text(z, levels(x$dframe[,i]), col = "RED", pos = 3)
    for (j in 1:length(x$dframe[,i])) 
      lines(rbind(x$scores[j,c(pd1,pd2)],z[x$dframe[,i][j],]), col = "BLUE")
    identify(x$scores[,c(pd1,pd2)], labels = rownames(x$dframe), col = "BLUE") 
  }
}
#----------------------------------end starplot---------------------------------

#------------------------------------vecplot------------------------------------
#draws vector plots

if (plot.type == "vecplot") {
  
  if (missing(xlim)) xlim <- range(x$scores[,pd1])
  if (missing(ylim)) ylim <- range(x$scores[,pd2])
  xylim <- c(min(xlim[1],ylim[1]),max(xlim[2],ylim[2]))

  for (i in var.subset) {
    if (missing(main)) main1 <- paste("Vector plot for", colnames(x$dframe[i])) else main1 <- main
     
    if (x$rank.vec[i] == 1) {
      a <- x$cat.loadings[[i]][,c(pd1,pd2)]
      par("ask" = TRUE)
      plot(x$scores[,c(pd1,pd2)], type = "p", main = main1, xlab = xlab, ylab = ylab,
           xlim = xylim, ylim = xylim, ...)
      text(x$scores[,c(pd1,pd2)], as.vector(x$dframe[,i]), pos = 3, ...) 
      text(x$rank.cat[[i]][,c(pd1,pd2)], levels(x$dframe[,i]), col = "RED", ...)
      slope = a[2]/a[1]
      abline(coef = c(0,slope))
      
      for (j in 1:length(x$dframe[,i])) {
        xs <- xe <- x$scores[j,c(pd1,pd2)]  
        xe[1] <- (xs[1]+(xs[2]*slope))/(1+(slope^2))
        xe[2] <- slope*xe[1]     
  	    lines(rbind(xs,xe), col = "BLUE", ...)
  	  }
      abline(h = 0, v = 0, col = "gray", lty = 2, ...) 
    } else {
      warning("No vector plot for ", colnames(x$dframe[i]),"(rank != 1).", call. = FALSE) 
    }
  }
}
#----------------------------------end vecplot----------------------------------


#------------------------------------trfplot------------------------------------
#draws transformation plots

if (plot.type == "trfplot") {

   if (missing(type)) type <- "b"
   if (missing(xlab)) xlab <- "original scale"
   if (missing(ylab)) ylab <- "transformed scale"

   for (i in var.subset) {
     
     if (missing(main)) main1 <- paste("Transformation plot for", colnames(x$dframe[i])) else main1 <- main
     if (missing(ylim)) ylim1 <- range(x$low.rank[[i]]) else ylim1 <- ylim
     
     p <- dim(x$low.rank[[i]])[2]          #number of dimensions
     vlev <- levels(x$dframe[,i])
     
     par("ask" = TRUE)                     #first dimensions
     matplot(x$low.rank[[i]], type = type, main = main1, ylim = ylim1, xlab = xlab, 
            ylab = ylab, xaxt = "n", pch = 20, col = 1:p, lty = 1:p,...)
     if (p != 1) legend(leg.pos,paste("Solution",1:p),col = 1:p, lty = 1:p,...)
     axis(1, at = 1:length(vlev), labels = vlev)
     
   }
}
#----------------------------------end trfplot----------------------------------

#------------------------------------vorplot------------------------------------
#draws voronoi regions

if (plot.type == "vorplot") {

  for (i in var.subset) {
     
     if (missing(main)) main1 <- paste("Voronoi plot for", colnames(x$dframe[i])) else main1 <- main
     z <- rbind(x$scores[,c(pd1,pd2)],x$rank.cat[[i]][,c(pd1,pd2)])
     if (missing(xlim)) xlim1 <- range(z[,1]) else xlim1 <- xlim
     if (missing(ylim)) ylim1 <- range(z[,2]) else ylim1 <- ylim
   
     par("ask" = TRUE)
     plot(x$scores[,c(pd1,pd2)], type = "n", main = main1, xlab = xlab, ylab = ylab, 
          xlim = xlim1, ylim = ylim1, ...)
     drawEdges(x$rank.cat[[i]][,c(pd1,pd2)], far = 1000)
     text(x$scores[,c(pd1,pd2)], as.vector(x$dframe[,i]), ...)
  }
}
#----------------------------------end vorplot----------------------------------
}

