\name{coxsplit}
\alias{coxsplit}\alias{coxsplity}\alias{coxsplitw}
\title{
Split Data for Cross-validation
}
\description{
Specify which fold each observation is in to obtain \code{foldid}.

\code{coxsplit} randomly split data. \code{coxsplity} splits data by the `status' of \code{y}, and within each status, data is randomly split. \code{coxsplitw} splits data in the order of \code{w}.
}
\usage{
coxsplit(y, nfolds)

coxsplity(y, nfolds)

coxsplitw(w, nfolds)
}
\arguments{
  \item{y}{response variable. \code{y} should be a two-column matrix with columns named `time' and `status'. The latter is a binary variable, with `1' indicating event, and `0' indicating right censored.
}
  \item{w}{input vector, same length as \code{y}. The coefficients vary with \code{w}.
}
  \item{nfolds}{number of folds.
}
}
\details{
By default, \code{coxsplit} is incorporated in \code{\link{Coxnet}} and \code{coxsplitw} in \code{\link{loCoxnet}}. To use other data splitting method, generate \code{foldid} at the outset and supply it as an augment in \code{\link{Coxnet}} and \code{\link{loCoxnet}}.
}
\value{
a vector of values between 1 and \code{nfolds} specifying which fold each observation is in.
}
\author{
Xiang Li, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xl2473@columbia.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Coxnet}}, \code{\link{loCoxnet}}, \code{\link{print.Coxnet}}
}
\examples{
set.seed(1213)
N=100
xb=rnorm(N)
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1)  # censoring indicator
y=cbind(time=ty,status=1-tcens)
foldid=coxsplit(y,10)  # 10-fold
}
\keyword{Cox model}
\keyword{Regularization}
\keyword{Split}


