\name{ranks_OV}
\alias{ranks_OV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Overlaps for each dimension using ranks
}
\description{
Function returns the asymmetric overlaps for each dimension, calculated by the method published by Parkinson et al. (2018) using ranks. Further two confidence intervals are returned for each estimate. The confidence level, as well as the repetitions for bootstrap can be adjusted.
}
\usage{
ranks_OV(A = A, alpha = 0.05, reps4boot = 1000, digit = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Data frame, where the first column contains two objects (e.g. species) and the other columns are numeric vectors (containing measurments).
}
  \item{alpha}{
The confidence level. Default: \code{alpha = 0.05} 
}
\item{reps4boot}{
Number of repetitions for the bootstrap. . Default: \code{reps4boot = 1000}
}
\item{digit}{
  Number of digits after which the results are cut off. Default: \code{digit = 3}
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Data Frame containing the two asymmetric overlaps for each dimension together with their confidence intervals. The last row contains the d-dimensional asymmetric overlaps.
}
\references{
Judith H. Parkinson, Raoul Kutil, Jonas Kuppler, Robert R. Junker, Wolfgang Trutschnig, Arne C. Bathke: A Fast and Robust Way to Estimate Overlap of Niches and Draw Inference, International Journal of Biostatistics (2018)
}
\author{
Judith H. Parkinson \email{judith.parkinson@sbg.ac.at},\cr 
Raoul Kutil \email{Raoul_Hugo.Kutil@stud.sbg.ac.at},\cr 
Jonas Kuppler \email{Jonas.Kuppler@uni-duesseldorf.de},\cr 
Robert R. Junker \email{Robert.Junker@sbg.ac.at} (corresponding author),\cr 
Wolfgang Trutschnig \email{Wolfgang.Trutschnig@sbg.ac.at},\cr 
Arne Bathke \email{Arne.Bathke@sbg.ac.at}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# example function ranks_OV
data(finch2)
head(finch2)
ranks_OV(finch2[1:4], alpha = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
