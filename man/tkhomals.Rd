\name{tkhomals}
\alias{tkhomals}
\title{Homogeneity Analysis GUI}
\description{
  This function provides a tcl/tk gui for the homals function.
}
\usage{
tkhomals(data)
}
\arguments{
  \item{data}{data (in data-frame)}
}
\details{
  This functions requires data to be stored in a data-frame. It can produce a variety
  of graphs which are stored in a pdf file. Results are stored in an ascii file.
}
\note{
  Needs to be executed from a writeable directory. File names are determined by the
  name of the data argument passed to the function. 
}
\author{Arno Ouwehand}
\examples{
## present tcl/tk gui for data set senate
\dontrun{data(senate)}
\dontrun{tkhomals(senate)}
}
\keyword{multivariate}
