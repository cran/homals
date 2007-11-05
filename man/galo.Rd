\name{galo}
\alias{galo}

\docType{data}
\title{GALO dataset}
\description{
The objects (individuals) are 1290 school children in the sixth grade of elementary
school in the city of Groningen (Netherlands) in 1959.
}
\usage{galo}
\format{
Data frame with the five variables Gender, IQ, Advice, SES and School. IQ (original range 60 to 144) 
has been categorized into 9 ordered categories and the schools are enumerated from 1 to 37. 

SES:
 
LoWC = Lower white collar; MidWC = Middle white collar; Prof = Professional, Managers;
Shop = Shopkeepers; Skil = Schooled labor; Unsk = Unskilled labor. 
 
Advice: 

Agr = Agricultural; Ext = Extended primary education; Gen = General; 
Grls = Secondary school for girls; Man = Manual, including housekeeping; 
None = No further education; Uni = Pre-University. 
}

\references{
Peschar, J.L. (1975). School, Milieu, Beroep. Groningen: Tjeek Willink.
}


\examples{
   data(galo)
   galo
}
\keyword{datasets}
