\name{plot.survfitr}
\alias{plot.survfitr}
\alias{lines.survfitr}

\title{Plots estimated survival function from an object of class `survrec'.}

\description{Additional plots can be added to the same of axes using `lines.survrec'.}

\usage{
plot.survfitr(x,prob=FALSE,...)
lines.survfitr(x,prob=FALSE, ...)
}

\arguments{
 \item{x}{
    Object of class \code{survrec} (output from calling \code{survrec}
    function).}
 \item{prob}{Print of the probability or survival function if its value is 
              TRUE or FALSE respectively.} 
 \item{...} {additional arguments passed to the plot function.}
}

\value{
  Prints a plot or a line of class \code{survrec}
}

\seealso{\code{\link{psh.fit}} \code{\link{mlefrailty.fit}}  \code{\link{wc.fit}}}

\keyword{survival}
