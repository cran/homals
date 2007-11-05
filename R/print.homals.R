`print.homals` <-
function(x, ...)
{
  nvar <- dim(x$dframe)[2]
  cat("\nLoss:",x$loss,"\n")
  cat("\nEigenvalues:\n")
  eigen.val <- round(x$eigenvalues, 4)
  names(eigen.val) <- paste("D",1:x$ndim,sep="")
  print(eigen.val)
  cat("\n")
  
  loadmat <- t(sapply(x$cat.loadings, function(xx) (xx[1,])))
  if (x$ndim == 1) loadmat <- t(loadmat)
  colnames(loadmat) <- paste("D", 1:x$ndim, sep = "")
  cat("Loadings (first solution only):\n")
  print(loadmat)
  cat("\n")
}

