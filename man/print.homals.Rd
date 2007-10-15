\name{print.homals}
\alias{print.homals}
\alias{summary.homals}
\title{S3 methods for homals objects}
\description{Print and summary method for objects of class \code{"homals"}
}
\usage{
\method{print}{homals}(x, ...)
\method{summary}{homals}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"homals"}.}
  \item{object}{Object of class \code{"homals"}.}
  \item{\dots}{Further arguments are ignored}
}
\details{The \code{print} methods provides the eigenvalues and the rank-restricted
category quantifications for each variable. The \code{summary} method gives the 
caterogy loadings, category centroids, rank-restricted category quantifications, and
lower rank quantifications for each variable separately.
}

\examples{

data(hartigan)
res <- homals(hartigan)
res
summary(res)

}

\keyword{ models }
