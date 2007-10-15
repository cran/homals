\name{hartigan}
\alias{hartigan}

\docType{data}
\title{Hartigan's Hardware}
\description{
  A number of bolts, nails, screws, and tacks are classified according to a 
  number of criteria.
}
\usage{hartigan}
\format{
  A data frame with the following variables and categories:
  
  Thread: Y = Yes, N = No
  
  Head: F = Flat, C = Cup, O = Cone, R = Round, Y = Cylinder
  
  Head indentation: N = None, T = Star, L = Slit
  
  Bottom: S = Sharp, F = Flat
  
  Length: (in half inches)
  
  Brass: Y = Yes, N = No
}

\references{
Hartigan, J.A. (1975). Clustering Algorithms. New York: Wiley.
}


\examples{
   data(hartigan)
   hartigan
}
\keyword{datasets}
