\name{mammals}
\alias{mammals}

\docType{data}
\title{Mammals dataset}
\description{
  Dentition of mammals.
}
\usage{mammals}
\format{
  A data frame with the tooth pattern, the number of top incisors, bottom incisors, top canines, bottom canines,
  top premolars, bottom premolars, top molars, and bottom molars.
  
  - top incisors: (1) zero incisors; (2) one incisor;  (3) two incisors; (4) three or more incisors.

  - bottom incisors: (1) zero incisors; (2) one incisor; (3) two incisors; (4) three incisors; (5) four
incisors. 

  - top canine: (1) zero canines; (2) one canine.
  
  - bottom canine: (1) zero canines; (2) one canine.
  
  - top premolar (1) zero premolars; (2) one premolar; (3) two premolars; (3) two premolars; (4)
three premolars; (5) four premolars.

  - bottom premolar: (1) zero premolars; (2) one premolar; (3) two premolars; (3) two premolars;
(4) three premolars; (5) four premolars.

  - top molar: (1) zero, one or two molars; (2) more than two molars.
  
  - bottom molar: (1) zero, one or two molars; (2) more than two molars.
}


\references{
Hartigan, J.A. (1975). Clustering Algorithms. New York: Wiley.
}
\examples{
   data(mammals)
   mammals
}
\keyword{datasets}