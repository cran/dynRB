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
Version: \tab 0.16\cr
Date: \tab 2021-05-11\cr
}
}
\author{
Manuela Schreyer \email{manuelalarissa.schreyer@sbg.ac.at},\cr 
Wolfgang Trutschnig \email{Wolfgang.Trutschnig@sbg.ac.at},\cr 
Robert R. Junker \email{Robert.Junker@sbg.ac.at} (corresponding author),\cr 
Jonas Kuppler \email{Jonas.Kuppler@uni-duesseldorf.de},\cr 
Arne Bathke \email{Arne.Bathke@sbg.ac.at},\cr 
Judith H. Parkinson \email{judith.parkinson@sbg.ac.at},\cr 
Raoul Kutil \email{Raoul_Hugo.Kutil@stud.sbg.ac.at}
}
\references{
Junker RR, Kuppler J, Bathke AC, Schreyer ML, Trutschnig W (2016) Dynamic range boxes - A robust non-parametric approach to quantify size and overlap of n-dimensional hypervolumes. Methods in Ecology and Evolution doi: 10.1111/2041-210X.12611\cr\cr
Judith H. Parkinson, Raoul Kutil, Jonas Kuppler, Robert R. Junker, Wolfgang Trutschnig, Arne C. Bathke: A Fast and Robust Way to Estimate Overlap of Niches and Draw Inference, International Journal of Biostatistics (2018)
}
\keyword{ package }
\examples{
# example function dynRB_VPa
# for reliable results use steps = 201
data(finch2)
r<-dynRB_VPa(finch2, steps = 101) 
r$result
}

