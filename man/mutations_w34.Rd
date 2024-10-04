\name{mutations_w34}
\alias{mutations_w34}
\docType{data}
\title{mutations_w34 Data Set}
\description{A binary matrix containing mutation values, wether the mutation is present or not in the patient from patients retrieved from waves 3 and 4. Used as a toy example in \code{trainTree} and \code{predictTree} and \code{niceTree}
}
\usage{data("mutations_w34")}
\format{
  A binary matrix with,where there are 142 patients from the waves and    70 different mutations   are considered.
  
  The format is:
  num [1:142, 1:70] 0 0 0 1 0 0...
}
\examples{
data(mutations_w34)
}
\keyword{datasets}