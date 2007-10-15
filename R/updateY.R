`updateY` <-
function(dframe,x,y,active,rank,level,sets,verbose=FALSE){
  nobj<-dim(x)[1]; ndim<-dim(x)[2]; nset<-length(sets)
  for (l in 1:nset) {
  	indi<-sets[[l]]; jndi<-indi[which(active[indi])]
    if (length(jndi) == 0) next()
  	ss<-sumSet(dframe,nobj,ndim,y,jndi)
    ii<-which(!is.na(dframe[,jndi[1]]))
  	for (j in jndi) {
		gg<-dframe[,j]; yy<-y[[j]]; d<-as.vector(table(gg))
		s1<-sum((x[ii,]-ss[ii,])^2)
		ss[ii,]<-ss[ii,]-yy[gg[ii],]
		yc<-computeY(gg[ii],x[ii,]-ss[ii,])
		yy<-restrictY(d,yc,rank[j],level[j])$y
		ss[ii,]<-ss[ii,]+yy[gg[ii],]
		s2<-sum((x[ii,]-ss[ii,])^2)
		y[[j]]<-yy
		if (verbose) cat("Set: ",formatC(l,digits=3,width=3),
			" After Variable: ",formatC(j,digits=3,width=3),
			" Loss: ", formatC(c(s1,s2),digits=6,width=9, format="f"),"\n")		
   		}
  	}
return(y)
}

