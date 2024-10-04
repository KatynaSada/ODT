\name{mutations_w12}
\alias{mutations_w12}
\docType{data}
\title{mutations_w12 Data Set}
\description{A binary matrix containing mutation values, wether the mutation is present or not in the patient from patients retrieved from waves 1 and 2. Used as a toy example in \code{trainTree} and \code{predictTree} and \code{niceTree}
}
\usage{data("mutations_w12")}
\format{
  A binary matrix with,where there are 247 patients from the waves and    70 different mutations   are considered.
  
  The format is:
  num [1:247, 1:70] 0 0 0 0 0 1...
}
\examples{
data(mutations_w12)
}
\keyword{datasets}
