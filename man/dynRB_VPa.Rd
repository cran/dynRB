\name{dynRB_VPa}
\alias{dynRB_VPa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Size and pairwise overlap
}
\description{
Function returns size and pairwise overlaps of niches or trait-spaces. Size or overlaps of dimensions can be aggregated by using either "product", "mean" or "geometric mean" as aggregation method. The results obtained by using the product are automatically printed. Number of dynamic range boxes (\code{steps}) can be adjusted. Default: \code{steps = 201}
}
\usage{
dynRB_VPa(A = A, steps = 201, correlogram = FALSE, row_col = c(2, 2), 
          pca.corr = FALSE, var.thres = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Data frame, where the first column is a character vector and the other columns are numeric vectors.
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
\item{pca.corr}{
  If \code{TRUE}, a principal components analysis is performed. 
}
\item{var.thres}{
  Variance predicted by the PCA-axes, if \code{pca.corr = TRUE}. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Data frame containing the summarized niche overlap (and volume) for each pair of objects aggregated by all three possible
choices (i.e. product, mean, geometric mean).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# example function dynRB_VPa
# for reliable results use steps = 201
data(finch2)
r<-dynRB_VPa(finch2, steps = 101, correlogram = TRUE, row_col = c(1,1))
r$result
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
