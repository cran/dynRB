\name{finch2}
\alias{finch2}
\docType{data}
\title{
Subset of data set finch
}
\description{
To demonstrate the application of the functions for real world data, we used existing data sets on niches and trait-spaces and quantified their sizes and overlaps. The data set finch2 is a data set on morphological measurements of three Darwin finches. The data set comprises quantitative measurements of nine traits characterizing two species of finches, each trait was measured at least in 10 individuals per species.
}
\usage{data("finch2")}
\format{
  A data frame with 103 observations on the following 10 variables.
  \describe{
    \item{\code{Species}}{a character vector of the Species Geospiza fuliginosa parvula and Geospiza fortis fortis}
    \item{\code{BodyL}}{a numeric vector}
    \item{\code{WingL}}{a numeric vector}
    \item{\code{TailL}}{a numeric vector}
    \item{\code{BeakW}}{a numeric vector}
    \item{\code{BeakH}}{a numeric vector}
    \item{\code{LBeakL}}{a numeric vector}
    \item{\code{UBeakL}}{a numeric vector}
    \item{\code{N.UBkL}}{a numeric vector}
    \item{\code{TarsusL}}{a numeric vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Snodgrass R and Heller E (1904) Papers from the Hopkins-Stanford Galapagos Expedition, 1898-99. XVI. Birds. Proceedings of the Washington Academy of Sciences 5: 231-372.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(finch2)
## quick overview
head(finch2)
}
\keyword{datasets}
