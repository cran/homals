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
homals(data, ndim = 2, rank = ndim, level = "nominal", sets = 0, active = TRUE, 
       eps = 1e-06, itermax = 100, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data in data frame.}
  \item{ndim}{Number of dimensions to be extracted.}
  \item{rank}{Which quantification ranks. Default is number of dimensions \code{ndim}}
  \item{level}{Which quantification levels. Possible values are \code{"nominal"},  
  \code{"ordinal"}, \code{"numerical"}, and \code{"polynomial"} which can be defined as single 
  character (if all variable are of the same level) or as vector which length corresponds to the
  number of variables.} 
  \item{sets}{List of vectors of indices.}
  \item{active}{Which variables are active for computation (\code{TRUE} means all)}
  \item{eps}{Iteration precision eigenvalues.}
  \item{itermax}{Maximum number of iterations.}
  \item{verbose}{0 for no iteration output, 1 prints outer iterations, 2 prints outer and inner iterations, 3 prints outer and inner and innermost iterations.}
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
  \item{catscores}{List of category scores for each variable.}
  \item{scoremat}{Reproduced data matrix with category scores.}
  \item{objscores}{Matrix with object scores on each dimension.}
  \item{cat.centroids}{List of category centroids.}
  \item{ind.mat}{Dummy coded indicator matrix.}
  \item{cat.loadings}{List of category loadings.}
  \item{low.rank}{List of lower rank quantifications.}
  \item{niter}{Number of iterations.}
  \item{eigenvalues}{Final eigenvalues.}
  \item{loss}{Value of the loss function.}
  \item{rank.vec}{Vector with ranks for each variable.}
  \item{active}{Vector with active/inactive variables.}
}
\references{
de Leeuw, J., \& Mair, P. (2007). Homogeneity analysis in R: The package homals. Preprint available at \url{http://gifi.stat.ucla.edu/homalsR.pdf}. 

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
  
  ##3-dimensional senate solution
  data(senate)
  res <- homals(senate, active = c(FALSE, rep(TRUE, 20)), ndim = 3)
  
  ##Predicting Species in iris dataset (polynomial level constraints for predictors)
  data(iris)
  res <- homals(iris, sets = list(1:4,5), level = c(rep("polynomial",4),"nominal"), rank = 2, itermax = 1000)

}
\keyword{models}
