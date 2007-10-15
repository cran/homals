\name{homals-package}
\alias{homals-package}
\docType{package}
\title{Homogeneity Analysis
}
\description{This package performs a homogeneity analysis, aka a multiple correspondence
  analysis, but with many additional options. Variables can be grouped into
  sets, in order to emulate regression analysis and canonical analysis. For
  each variable there are, in addition, rank constraints on the category
  quantifications (or transformations) and level constraints (which allows
  one to treat a variable as nominal, ordinal, or numerical). 
}
\details{
\tabular{ll}{
Package: \tab homals\cr
Type: \tab Package\cr
Version: \tab 0.9-0\cr
Date: \tab 2007-10-15\cr
License: \tab GPL\cr
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
Gifi, A. (1990). Nonlinear Multivariate Analysis. New York: Wiley. 
}

\keyword{ package }
\seealso{\code{\link{homals}} 
}
