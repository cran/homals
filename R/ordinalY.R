`ordinalY` <-
function(d,y,r,itermax=100,eps=1e-6,verbose=FALSE) {
r<-min(r,dim(y))
if (r == 1)	return(singord(d,y,r,itermax=itermax,eps=eps,verbose=verbose))
	else return(multord(d,y,r,itermax=itermax,eps=eps,verbose=verbose))
}

