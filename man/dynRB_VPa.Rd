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
dynRB_VPa(A = A, steps = 201, graphic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Data frame, where the first column is a character vector and the other columns are numeric vectors.
}
  \item{steps}{
Number of range boxes. Default: \code{steps = 201} 
}
  \item{graphic}{
If \code{TRUE}, graphics are shown. If \code{FALSE}, no graphics are shown. Default: \code{graphic = FALSE}
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
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Manuela Schreyer \email{manuelalarissa.schreyer@sbg.ac.at} and
Wolfgang Trutschnig \email{Wolfgang.Trutschnig@sbg.ac.at}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# example function dynRB_VPa
# for reliable results use steps = 201
data(finch2)
r<-dynRB_VPa(finch2, steps = 101)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
