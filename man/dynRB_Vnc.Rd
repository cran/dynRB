\name{dynRB_Vnc}
\alias{dynRB_Vnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Corrected Relative Dynamic Range Box size per dimension and object
}
\description{
Function returns adjusted Dynamic Range Box size of each dimension \eqn{n}. Number of dynamic range boxes (\code{steps}) can be adjusted. Default: \code{steps = 201}
}
\usage{
dynRB_Vnc(A = A, steps = 201, graphic = FALSE)
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
Data frame containing the summarized (adjusted) niche length for each object and dimension.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
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
# example function dynRB_Vnc
# for reliable results use steps = 201
data(finch2)
r<-dynRB_Vnc(finch2, steps = 101)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
