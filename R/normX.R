`normX` <-
function(x,w) {
  qq<-La.svd((1/sqrt(w))*x)
  list(q=(1/sqrt(w))*(qq$u)%*%(qq$vt),r=qq$d)    
  #list(q=sqrt(dim(x)[1]*sqrt(w))*(qq$u)%*%(qq$vt),r=qq$d)               #variance = 1, X'MX = nmI
}

