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
  
  The general idea of homogeneity analysis is to make a joint plot in p-space of all
  objects (or individuals) and the categories of all variables. In this plot
  we connect objects with the categories they are in, thus producing
  a graph plot. If there are m variables, then m lines depart from each
  object, and with n object the graph plot has nm lines. The technique, in
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
\usage{
homals(data, sets=0, ndim=2, active=T, rank=ndim, level="NO", starplots = FALSE,
catplots =  FALSE, trfplots =  FALSE, lossplots =  FALSE, graphplot =  FALSE,
objplot =  FALSE, objscores =  FALSE, objlabel =  FALSE, offset =  1.20,
eps1 =  -Inf, eps2 =  1e-6, itermax =  100, voronoi =  FALSE, save.me =  FALSE,
demo =  FALSE , timer =  FALSE , tk = FALSE, img1, img2, img3, name)
}
\arguments{
  \item{data}{data (in data-frame)}
  \item{sets}{list of vectors of indices}
  \item{ndim}{dimensionality (default 2)}
  \item{active}{which variables are active (single T means all)}
  \item{rank}{which quantification ranks (default all p)}
  \item{level}{which quantification levels NO (nominal), OR (orthogonal) (default all nominal)}
  \item{starplots}{which starplots, default is none (FALSE)}
  \item{catplots}{which category plots, default is none (FALSE)}
  \item{trfplots}{which transformation plots, default is none (FALSE)}
  \item{lossplots}{which loss plots, default is none (FALSE)}
  \item{graphplot}{graphplot, default is no (FALSE)}
  \item{objplot}{object score plot, default is no (FALSE)}
  \item{objscores}{object scores written to file, default is no (FALSE)}
  \item{objlabel}{object score plot labeled, default is no (FALSE)}
  \item{offset}{offset for labeled plots, default is 1.20}
  \item{eps1}{iteration precision eigenvalues, default is -Inf}
  \item{eps2}{iteration precision eigenvectors, default is 1e-6}
  \item{itermax}{maximum number of iterations, default is 100}
  \item{voronoi}{voronoi diagram, default is no (FALSE)}
  \item{save.me}{do we return the results, default is no (FALSE)}
  \item{demo}{animated iteration demo, default is no (FALSE)}
  \item{timer}{time the steps of program, default is no (FALSE)}
  \item{tk}{create tk output, default is no (FALSE), used by \code{tkhomals}}
  \item{img1}{tkrplot image placeholder, used by \code{tkhomals}} 
  \item{img2}{tkrplot image placeholder, used by \code{tkhomals}} 
  \item{img3}{tkrplot image placeholder, used by \code{tkhomals}} 
  \item{name}{dataframe name, used by \code{tkhomals}}
}
\details{
  This functions requires data to be stored in a data-frame. It can produce a variety
  of graphs which are stored in a pdf file. Results are stored in an ascii file.
}
\note{
  Needs to be executed from a writeable directory. File names are determined by the
  name of the data argument passed to the function. 
}
\author{Jan de Leeuw}
\examples{
## produce graphplot and objplot for data set mammals,
## results are save into an ascii file "mammals.out"
## all graphs are is saved into a pdf file named "mammals.pdf"
\dontrun{data(mammals)}
\dontrun{homals(mammals, graphplot = TRUE, objplot = TRUE)}
}
\keyword{multivariate}
