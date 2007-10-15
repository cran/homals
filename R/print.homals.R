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
 
  
 # cat("\n\nRank-restricted Category Quantifications:")
 # for (i in 1:nvar)
 # {
 #   rcq <- round(x$rank.cat[[i]], 4)
 #   cat("\n",colnames(x$dframe[i]),":\n",sep="")
 #   print(rcq)
 # }  
 # cat("\n")
}

