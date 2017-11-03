\name{dynRB_Pn}
\alias{dynRB_Pn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise overlaps for each dimension
}
\description{
Function returns pairwise overlaps for each dimension \eqn{n}. Number of dynamic range boxes (\code{steps}) can be adjusted. Default: \code{steps = 201}}
\usage{
dynRB_Pn(A = A, steps = 201, correlogram = FALSE, row_col = c(2, 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Data frame, where the first column is a character vector containing the objects (e.g. species) and the other columns are numeric vectors
(containing measurements).
}
  \item{steps}{
Number of range boxes. Default: \code{steps = 201}
}
  \item{correlogram}{
If \code{TRUE}, the correlogram for each species is shown. If \code{FALSE}, no correlogram is shown. Default: \code{correlogram = FALSE}
}
\item{row_col}{
  Number of rows and columns of the figures (correlogram for each species). Default: \code{row_col = c(2, 2)}
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Data frame containing the summarized overlaps for each pair of objects and dimension.
}
\references{
Junker RR, Kuppler J, Bathke AC, Schreyer ML, Trutschnig W (2016) Dynamic range boxes - A robust non-parametric approach to quantify size and overlap of n-dimensional hypervolumes. Methods in Ecology and Evolution doi: 10.1111/2041-210X.12611
}
\author{
Manuela Schreyer \email{manuelalarissa.schreyer@sbg.ac.at},\cr 
Wolfgang Trutschnig \email{Wolfgang.Trutschnig@sbg.ac.at},\cr 
Robert R. Junker \email{Robert.Junker@sbg.ac.at} (corresponding author),\cr 
Jonas Kuppler \email{Jonas.Kuppler@uni-duesseldorf.de},\cr 
Arne Bathke \email{Arne.Bathke@sbg.ac.at}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# example function dynRB_Pn
# for reliable results use steps = 201
data(finch2)
r<-dynRB_Pn(finch2, steps = 101)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
