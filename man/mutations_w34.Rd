\name{mutations_w34}
\alias{mutations_w34}
\docType{data}
\title{mutations_w34 Data Set}
\description{
A binary matrix representing mutation status for patients from Waves 3 and 4 of the BeatAML cohort, indicating whether specific mutations are present (1) or absent (0) in each patient. This dataset serves as a toy example for demonstrating the functionality of the \code{trainTree}, \code{predictTree}, and \code{niceTree} functions.
}
\usage{data("mutations_w34")}
\format{
  A binary matrix consisting of 142 rows and 70 columns, where each row represents a different patient and each column corresponds to a specific mutation. 
  
  The format is as follows:
  num [1:247, 1:70] 0 0 0 0 0 1 ...
}

\examples{
data(mutations_w34)
}
\keyword{datasets}