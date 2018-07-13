\name{dynRB-package}
\alias{dynRB-package}
\alias{dynRB}
\docType{package}
\title{
Dynamic Range Boxes
}
\description{
The package DynRB improves the concept of multivariate range boxes, which is highly susceptible for outlines and does not consider the distribution of the data. The package uses dynamic range boxes to overcome these problems.
}
\details{
\tabular{ll}{
Package: \tab dynRB\cr
Type: \tab Package\cr
Version: \tab 0.13\cr
Date: \tab 2018-07-13\cr
}
}
\author{
Manuela Schreyer \email{manuelalarissa.schreyer@sbg.ac.at},\cr 
Wolfgang Trutschnig \email{Wolfgang.Trutschnig@sbg.ac.at},\cr 
Robert R. Junker \email{Robert.Junker@sbg.ac.at} (corresponding author),\cr 
Jonas Kuppler \email{Jonas.Kuppler@uni-duesseldorf.de},\cr 
Arne Bathke \email{Arne.Bathke@sbg.ac.at}
}
\references{
Junker RR, Kuppler J, Bathke AC, Schreyer ML, Trutschnig W (2016) Dynamic range boxes - A robust non-parametric approach to quantify size and overlap of n-dimensional hypervolumes. Methods in Ecology and Evolution doi: 10.1111/2041-210X.12611
}
\note{
A tutorial that introduces the functions implemented in the package, explains the outputs, gives examples on how to use the output for follow-up analyses, and finally shows how to visualize the results can be found here: \href{https://tinyurl.com/y7xgdeq7}{Tutorial dynRB}
}
\keyword{ package }
\examples{
# example function dynRB_VPa
# for reliable results use steps = 201
data(finch2)
r<-dynRB_VPa(finch2, steps = 101) 
r$result
}

