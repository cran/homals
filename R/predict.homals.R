`predict.homals` <-
function(object, ...)
{ 
  #computes classification table and misclassification rate
 
  cl.table <- NULL
  nvec <- colnames(object$dframe)
  for (i in 1:length(object$cat.loadings)) 
  {
    ax <- object$cat.loadings[[i]]
    ay <- object$catscores[[i]]
    ag <- object$dframe[,i]
    bx <- lsfit(t(ax),t(object$objscores),intercept=FALSE)$coef
    ux <- crossprod(bx,ax)
    d <- outer(rowSums(ux^2),rowSums(ay^2),"+")-2*tcrossprod(ux,ay)
    h <- levels(ag)[apply(d,1,which.min)]
    cl.table[[i]] <- table(ag,h,dnn=list("obs","pre"))
    names(cl.table)[[i]] <- nvec[i]
  }
  cr.vec <- sapply(cl.table, function(x) (sum(diag(x)))/sum(x))
  result <- list(cl.table = cl.table, cr.vec = cr.vec)
  class(result) <- "predict.homals"
  result
}

