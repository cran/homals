\name{homals-package}
\alias{homals-package}
\docType{package}
\title{Homogeneity Analysis
}
\description{This package performs a homogeneity analysis and various extensions. Rank restrictions on the category quantifications can be imposed (nonlinear PCA). The categories are transformed by means of optimal scaling with options for nominal, ordinal, and numerical scale levels (for rank-1 restrictions). Variables can be grouped into sets, in order to emulate regression analysis and canonical correlation analysis. 
}
\details{
\tabular{ll}{
Package: \tab homals\cr
Type: \tab Package\cr
Version: \tab 0.9-10\cr
Date: \tab 2007-02-22\cr
License: \tab GPL\cr
URL: \tab \url{http://r-forge.r-project.org/projects/psychor/}\cr
}
The basic function of the homals package is \code{homals(data)} which performs 
homogeneity analysis and returns an object of class \code{"homals"}. It allows for 
grouping variables into sets, constraining ranks on category quantifications and levels.
Various plot functions are provided in \code{plot.homals} and \code{plot3d.homals}.
}
\author{
Jan de Leeuw, Patrick Mair

Maintainer: Jan de Leeuw <deleeuw@ucla.edu>
}
\references{
de Leeuw, J., \& Mair, P. (2007). Homogeneity analysis in R: The package homals. Preprint available at \url{http://gifi.stat.ucla.edu/homalsR.pdf}. 

Gifi, A. (1990). Nonlinear Multivariate Analysis. New York: Wiley. 
}

\keyword{ package }
\seealso{\code{\link{homals}} 
}
