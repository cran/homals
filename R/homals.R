`homals` <-
function(dframe,   # data (in data-frame)
		sets=0, 				# list of vectors of indices
		ndim=2,              	# dimensionality (default 2)
		active=TRUE,            # which variables are active (single TRUE means all)
		rank=ndim,           	# which quantification ranks (default all ndim)
		level="nominal",		# which measurement level (default all nominal)
		eps=1e-6,           	# iteration precision eigenvalues (default 1e-6)
		itermax=100,         	# maximum number of iterations (default 100)
		verbose=0				# debugging output level
)

{

#-----------------------------set some constants--------------------------------

name <- deparse(substitute(dframe))		# frame name
nobj <- nrow(dframe)					# number of objects
nvar <- ncol(dframe)					# number of variables
vname <- names(dframe)					# variable names
rname <- rownames(dframe)				# object names

#-----------------------------convert to factors--------------------------------

for (j in 1:nvar) {
	dframe[, j] <- as.factor(dframe[, j])
    levels(dframe[, j])<-sort(levels(dframe[, j]))
	}

#-----------------------------parameter consistency-----------------------------

active<-checkPars(active,nvar)
rank<-checkPars(rank,nvar)
level<-checkPars(level,nvar)

if (length(sets) == 1) sets <- lapply(1:nvar,"c")
if (!all(sort(unlist(sets)) == (1:nvar))) {
	print(cat("sets union",sort(unlist(sets)),"\n"))
	stop("inappropriate set structure !")
	}
nset <- length(sets)

mis<-rep(0,nobj)
for (l in 1:nset) {
	lset<-sets[[l]]
	if (all(!active[lset])) next()
	jset<-lset[which(active[lset])]
	for (i in 1:nobj) {
		if (any(is.na(dframe[i,jset]))) 
			dframe[i,jset] <- NA
		else mis[i] <- mis[i] + 1
		}
	}
	
for (j in 1:nvar) {
	k<-length(levels(dframe[,j]))
	if (rank[j] > min(ndim,k-1)) rank[j]<-min(ndim,k-1)
	}
	
#----------------initialize scores and counters-----------------------------

x <- cbind(orthogonalPolynomials(mis,1:nobj,ndim))
x <- normX(centerX(x,mis),mis)$q
y <- lapply(1:nvar, function(j) computeY(dframe[,j],x))
y <- updateY(dframe,x,y,active,rank,level,sets)
sold <- totalLoss(dframe,x,y,active,rank,level,sets)
iter <- pops <- 0

#-----------------------------main computation--------------------------------

repeat {
	iter <- iter + 1
	y<-updateY(dframe,x,y,active,rank,level,sets,verbose=verbose)
	smid <- totalLoss(dframe,x,y,active,rank,level,sets)
	ssum <- totalSum(dframe,x,y,active,rank,level,sets)
	qv <- normX(centerX((1/mis)*ssum,mis),mis)
	z <- qv$q
	snew<-totalLoss(dframe,z,y,active,rank,level,sets)
	if (verbose > 0) cat("Itel:",formatC(iter,digits=3,width=3),"Loss Total: ", formatC(c(sold,smid,snew),digits=6,width=9,format="f"),"\n")
	r <- qv$r
	if (iter == itermax) {
		stop("maximum number of iterations reached")
		}
	if (snew > sold) {
		stop(cat("loss function increases in iteration ",iter,"\n"))
		}
	if ((sold - snew) < eps) break()
		else {x <- z; sold <- snew}
	}

#-----------------------------store final version--------------------------------

ylist<-alist<-clist<-ulist<-NULL
for (j in 1:nvar) {
  gg<-dframe[,j]; c<-computeY(gg,z); d<-as.vector(table(gg))
  lst<-restrictY(d,c,rank[j],level[j])
  y<-lst$y; a<-lst$a; u<-lst$z
  ylist<-c(ylist,list(y)); alist<-c(alist,list(a)); clist<-c(clist,list(c)); ulist<-c(ulist,list(u))
}

#--------------------------preparing/labeling output----------------------------

dimlab <- paste("D", 1:ndim, sep = "")
for (i in 1:nvar) {
  if (ndim == 1) {
    ylist[[i]] <- cbind(ylist[[i]])
    ulist[[i]] <- cbind(ulist[[i]])
    clist[[i]] <- cbind(clist[[i]])
    #alist[[i]] <- cbind(alist[[i]])
  }
  rownames(ylist[[i]]) <- rownames(ulist[[i]]) <- rownames(clist[[i]])
  rownames(alist[[i]]) <- paste(1:dim(alist[[i]])[1])
  colnames(clist[[i]]) <- colnames(ylist[[i]]) <- colnames(alist[[i]]) <- dimlab
  colnames(ulist[[i]]) <- paste(1:dim(as.matrix(ulist[[i]]))[2])
}
names(ylist) <- names(ulist) <- names(clist) <- names(alist) <- colnames(dframe)
rownames(z) <- rownames(dframe)
colnames(z) <- dimlab
#alist.t <- lapply(alist,t)

#--------------------------end preparing/labeling output------------------------

result <- list(datname = name, dframe = dframe, ndim = ndim, niter = iter, level = level, 
               eigenvalues = r, loss = snew, rank.vec = rank,
               scores = z, rank.cat = ylist, cat.centroids = clist,
               cat.loadings = alist, low.rank = ulist, active = active)
class(result) <- "homals"
result
}

