\name{predict.homals}
\alias{predict.homals}
\alias{print.predict.homals}
\alias{summary.predict.homals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Classification tables for homals solutions}
\description{Given a multiple homals solution we can reconstruct the indicator
matrix by assigning each object to the closest category point of
the variable. We can then find out how well we have reconstructed the
original data in terms of a classification table.
}
\usage{
\method{predict}{homals}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{"homals"}}
  \item{\dots}{Further arguments ignored}
}
\details{For variables with rank restrictions we first project the objects on the 
hyperplane spanned by the category quantifications, and then compute distances in that plane.
If there are sets of variables, the percentage correctly classified
for a variable is usually not the best thing to look at; except in some
special cases, in which some sets consist of a single variable (such as
regression and discriminant analysis). In any case a square table with observed vs.
predicted values is computed.
}
\value{
  \item{cl.table}{List of classification tables for each variable.}
  \item{cr.vec}{Vector with classification rates.}
}

\seealso{\code{\link{homals}}}
\examples{
  data(galo)
  res <- homals(galo, active=c(TRUE, TRUE,TRUE, TRUE, FALSE))
  pr.res <- predict(res)
  pr.res
  summary(pr.res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
