\name{survfitr}
\alias{survfitr}
\title{ Compute a Survival Curve for Recurrent Event Data}
\description{
 Computes an estimate of a survival curve for recurrent event data using either 
the Peña-Strawderman-Hollander, Wang-Chang or MLE Frailty estimators. The resulting object of class "survfitr" 
   is plotted by `plot.survfitr', before it is returned.
}
\usage{
survfitr(formula, data, type="MLEfrailty",...) 
}

\arguments{
 \item{formula}{A formula object. If a formula object is supplied it must have
a Survr object as the response on the left of the ~ operator and a term on the
right. For a single survival curve the "~1" part of the formula is required.  }
  \item{data}{a data frame in wich to interpret the variables named in the formula.}
  \item{type}{a character string specifying the type of survival curve. Possible
value are "pena-strawderman-hollander", "wang-chang" or "MLEfrailty". The default
is "MLEfrailty". Only the first words are required, e.g "pe","wa","ML"}
\item{...} {additional arguments passed to the type of estimator.
             (see \code{psh.fit}, \code{wc.fit}, or \code{mlefrailty.fit}}  }

}
\details{
  See the help details of psh.fit, wc.fit or mlefrailty depending on the type chosen
}
\value{

a survfitr object. Methods defined for survfitr objects are provided for print,plot,lines 
and summary.
}

\references{

1. Peña, E.A., Strawderman, R. and Hollander, M. (2001). Nonparametric Estimation with Recurrent Event Data.  \emph{J. Amer. Statist. Assoc} \bold{96}, 1299-1315. \cr
2. Wang, M.-C. and Chang, S.-H. (1999). Nonparametric Estimation of a Recurrent Survival Function. \emph{J. Amer. Statist. Assoc} \bold{94}, 146-153.

}

\seealso{\code{\link{print.survfitr}},\code{\link{plot.survfitr}},
   \code{\link{lines.survfitr}}, \code{\link{summary.survfitr}},
  \code{\link{Survr}},\code{\link{psh.fit}},\code{\link{wc.fit}},
  \code{\link{mlefrailty.fit}}  }

\examples{
  data(MMC)
  # fit a pena-strawderman-hollander and plot it
  fit<-survfitr(Survr(id,time,event)~group,data=MMC,type="pena")
  plot(fit)
  # print the survival estimators
  fit
  summary(fit)

  # fit a MLE Frailty and plot it
  fit<-survfitr(Survr(id,time,event)~group,data=MMC,type="MLE")
  plot(fit)
  # print the survival estimators
  fit
  summary(fit)

}


\keyword{survival }
