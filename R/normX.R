`normX` <-
function(x,w) {
  qq<-La.svd((1/sqrt(w))*x); list(q=(1/sqrt(w))*(qq$u)%*%(qq$vt),r=qq$d)
}

