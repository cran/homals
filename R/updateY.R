`updateY` <-
function(dframe,x,y,active,rank,level,sets,verbose=0){
  nobj<-dim(x)[1]; ndim<-dim(x)[2]; nset<-length(sets)
  for (l in 1:nset) {
  	indi<-sets[[l]]; jndi<-indi[which(active[indi])]
    if (length(jndi) == 0) next()
  	ii<-which(!is.na(dframe[,jndi[1]]))
	if (length(ii) == 0) next()
  	ss<-sumSet(dframe,nobj,ndim,y,jndi)
    for (j in jndi) {
		gg<-dframe[ii,j]; yy<-y[[j]]; d<-as.vector(table(gg))           
		s1<-sum((x[ii,]-ss[ii,])^2)
		ss[ii,]<-ss[ii,]-yy[gg,]
		yc<-computeY(gg,x[ii,]-ss[ii,])
		yy<-restrictY(d,yc,rank[j],level[j],verbose=verbose)$y
		ss[ii,]<-ss[ii,]+yy[gg,]
		s2<-sum((x[ii,]-ss[ii,])^2)
		y[[j]]<-yy
		if (verbose > 1) cat("Set: ",formatC(l,digits=3,width=3),
			" After Variable: ",formatC(j,digits=3,width=3),
			" Loss: ", formatC(c(s1,s2),digits=6,width=9, format="f"),"\n")		
   		}
  	}
return(y)
}

