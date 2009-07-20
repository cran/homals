normX <- function (x,w)
	{
	   qq <- qr((1/sqrt(w)) * x)
	   list(q = (1/sqrt(w)) * qr.Q(qq), r=abs(diag(qr.R(qq))))
	}


