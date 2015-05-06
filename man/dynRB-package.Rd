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
Version: \tab 0.1\cr
Date: \tab 2015-05-06\cr
}
}
\author{
Manuela Schreyer \email{manuelalarissa.schreyer@sbg.ac.at} and
Wolfgang Trutschnig \email{Wolfgang.Trutschnig@sbg.ac.at}
}
\references{
Will be added as soon as the corresponding paper describing the approach is accepted for publication.
}
\keyword{ package }
\examples{
# example function dynRB_SQac
# for reliable results use steps = 201
data(finch2)
r<-dynRB_SQac(finch2, steps = 101)  
}
