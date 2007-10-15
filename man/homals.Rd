\name{homals}
\alias{homals}
\title{Homogeneity Analysis}
\description{
  This function performs a homogeneity analysis, aka a multiple correspondence
  analysis, but with many additional options. Variables can be grouped into
  sets, in order to emulate regression analysis and canonical analysis. For
  each variable there are, in addition, rank constraints on the category
  quantifications (or transformations) and level constraints (which allows
  one to treat a variable as nominal, ordinal, or numerical). 
}
\usage{
homals(dframe, sets = 0, ndim = 2, active = TRUE, rank = ndim, level = "nominal", 
       eps1 = 1e-06, eps2 = 1e-06, itermax = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dframe}{Data in data-frame.}
  \item{sets}{List of vectors of indices.}
  \item{ndim}{Number of dimensions to be extracted.}
  \item{active}{Which variables are active for computation (\code{TRUE} means all)}
  \item{rank}{Which quantification ranks. Default is number of dimensions \code{ndim}}
  \item{level}{Which quantification levels. Possible values are \code{"nominal"},
  \code{"ordinal"}, \code{"numerical"}, and \code{"polynomial"} which can be defined as single 
  character (if all variable are of the same level) or as vector which length corresponds to the
  number of variables.} 
  \item{eps1}{Iteration precision eigenvalues.}
  \item{eps2}{Iteration precision eigenvectors.}
  \item{itermax}{Maximum number of iterations.}
}
\details{
  The general idea of homogeneity analysis is to make a joint plot in p-space of all
  objects (or individuals) and the categories of all variables. The technique, in
  its most simple form, makes the graph plot in such a way that the sum of squares
  of the length of the nm lines is a small as possible, subject to a normalization
  of the object scores (their n x p coordinate matrix must be orthonormal).
  
  Rank constraints require the category quantifications of the categories of
  a variable to lie in a subspace of p-space. Requiring rank equal to one for all
  variables reduces homogeneity analysis to principal component analysis
  (with optimal scaling of the variables). 
  
  Sets of variables are incorporated by using additivity restrictions on the
  category quantifications (i.e. we code the variables within a set interactively,
  but then use quantifications based on main effects only). 
  
  By combining the various types of restrictions we obtain far-reaching
  generalizations of principal component analysis, canonical analysis,
  discriminant analysis and regression analysis.
}

\value{
  Returns an object of class \code{"homals"}.
  
  \item{datname}{Name of the data file.}
  \item{dframe}{Data file.}
  \item{ndim}{Number of dimensions.}
  \item{niter}{Number of iterations.}
  \item{level}{Vector with variable levels.}
  \item{eigenvalues}{Final eigenvalues.}
  \item{loss}{Value of the loss function.}
  \item{rank.vec}{Vector with ranks for each variable.}
  \item{scores}{Matrix with object scores on single dimension.}
  \item{rank.cat}{List of rank restriced category quantifications.}
  \item{cat.centroids}{List of category centroids.}
  \item{cat.loadings}{List of category loadings.}
  \item{low.rank}{List of lower rank quantifications.}
}
\references{
Gifi, A. (1990). Nonlinear Multivariate Analysis. New York: Wiley. 
}
\author{Jan de Leeuw, Patrick Mair}
\seealso{\code{\link{plot.homals}}}
\examples{
  data(galo)
  
  ##simple homals 2D-homals solution (School not active)
  res <- homals(galo, active = c(rep(TRUE, 4), FALSE))
  
  ##predict IQ 
  res <- homals(galo, active = c(rep(TRUE, 4), FALSE), sets = list(c(1,3,4),2,5))
  
  ##Single ranks for each variable (non-linear PCA)
  res <- homals(galo, active = c(rep(TRUE, 4), FALSE), sets = list(c(1,3,4),2,5))
  
  ##Nonlinear canonical correlation (2 sets of variables)
  res <- homals(galo, active = c(rep(TRUE, 4), FALSE), sets = list(c(1,2),c(3,4),5))
  
  ##regression (numerical)
  data(neumann)
  res <- homals(neumann, sets = list(3, 1:2), level = "numerical", rank = 1)
  
  ##regression (monotone)
  res <- homals(neumann, sets = list(3, 1:2), level = "ordinal", rank = 1)
}
\keyword{models}
