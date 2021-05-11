\name{finch}
\alias{finch}
\docType{data}
\title{
Data set finch
}
\description{
To demonstrate the application of the functions for real world data, we used existing data sets on niches and trait-spaces and quantified their sizes and overlaps. The data set finch is a data set on morphological measurements of Darwin finches. The data set comprises quantitative measurements of nine traits characterizing five species of finches, each trait was measured at least in 10 individuals per species.
}
\usage{data("finch")}
\format{
  A data frame with 146 observations on the following 10 variables.
  \describe{
    \item{\code{Species}}{a character vector of the Species Geospiza heliobates, Geospiza prosthemelas prosthemelas, Geospiza fuliginosa parvula, Geospiza fortis fortis and Geospiza fortis platyrhyncha}
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
data(finch)
## quick overview
head(finch)
}
\keyword{datasets}
