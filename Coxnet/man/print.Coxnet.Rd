\name{print.Coxnet}
\alias{print.Coxnet}
\title{
Print a Coxnet Object
}
\description{
Print a summary of results along the path of \code{lambda}.
}
\usage{
\method{print}{Coxnet}(x, digits = 4, ...)
}
\arguments{
  \item{x}{fitted Coxnet object
}
  \item{digits}{significant digits in printout
}
  \item{\dots}{additional print arguments
}
}
\details{
The performed model is printed, followed by \code{fit} and \code{fit0} (if any) from a fitted \code{Coxnet} object, or \code{fit} and \code{cvh} (if any) from a fitted \code{loCoxnet} object
}
\value{
The matrix above is silently returned
}
\author{
Xiang Li, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xl2473@columbia.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Coxnet}}, \code{\link{loCoxnet}}, \code{\link{coxsplit}}
}
\examples{
set.seed(1213)
N=1000;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]%*%beta
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1)  # censoring indicator
y=cbind(time=ty,status=1-tcens)

fiti=Coxnet(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
print(fiti)
}
\keyword{Cox model}
\keyword{Regularization}
\keyword{Print}


