% from pmg() as a template
\name{hetsar}
\alias{hetsar}
\alias{summary.hetsar}
\alias{print.summary.hetsar}
\alias{residuals.hetsar}
\title{Heterogeneous estimation of spatial lag (HSAR) panel models}
\usage{
hetsar(
  formula,
  data,
  w,
  na.action,
  index = NULL,
  ...
)

\method{summary}{hetsar}(object, ...)

\method{print}{summary.hetsar}(
  x,
  digits = max(3, getOption("digits") - 2),
  width = getOption("width"),
  ...
)

\method{residuals}{hetsar}(object, ...)
}
\arguments{
\item{formula}{a symbolic description of the model to be estimated,}

\item{data}{a \code{data.frame},}

\item{w}{a \code{matrix} or a \code{listw} object,}

\item{na.action}{see \code{\link[=lm]{lm()}},}

\item{index}{the indexes, see \code{\link[=pdata.frame]{pdata.frame()}},}

\item{\dots}{further arguments.}

\item{object, x}{an object of class \code{hetsar},}

\item{digits}{digits,}

\item{width}{the maximum length of the lines in the print output,}
}
\value{
An object of class \code{c("hetsar", "panelmodel")} containing:
\item{coefficients}{the vector of coefficients,}
\item{residuals}{the vector of residuals,}
\item{fitted.values}{the vector of fitted values,}
\item{vcov}{the covariance matrix of the coefficients,}
\item{df.residual}{degrees of freedom of the residuals,}
\item{model}{a data.frame containing the variables used for the
estimation,} \item{call}{the call,} \item{sigma}{always \code{NULL},
\code{sigma} is here only for compatibility reasons (to allow using
the same \code{summary} and \code{print} methods as \code{pggls}),}
\item{indcoef}{the matrix of individual coefficients from
separate time series regressions.}
}
\description{
Heterogeneous estimation of spatial lag (HSAR) panel models
}
\details{
\code{hetsar} is a function for the estimation of linear spatial
autoregressive panel models with
heterogeneous coefficients (HSAR) by the maximum likelihood estimator.
}
\examples{
data(HousePricesUS, package="pder")
data(usaw49, package="pder")
## HSAR estimator
hsarmod <- hetsar(log(price) ~ log(income),
                  HousePricesUS, usaw49)

summary(hsarmod)
}
\references{
ABP
}
\author{
Giovanni Millo
}
\keyword{spatial panel heterogeneous autoregression}