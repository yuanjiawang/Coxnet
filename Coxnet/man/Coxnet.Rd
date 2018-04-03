\name{Coxnet, loCoxnet}
\alias{Coxnet}\alias{loCoxnet}
\title{
Fit a Cox Model with Various Regularization Forms
}
\description{
\code{Coxnet} fits a Cox model regularized with net, elastic-net or lasso penalty, and their adaptive forms, such as adaptive lasso and net adjusting for signs of linked coefficients.
Moreover, it treats the number of non-zero coefficients as another tuning parameter and simultaneously selects with the regularization parameter \code{lambda}.

\code{loCoxnet} fits a varying coefficient Cox model by kernel smoothing, incorporated with the aforementioned penalties.

The package uses one-step coordinate descent algorithm and runs extremely fast by taking into account the sparsity structure of coefficients.
}
\usage{

Coxnet(x, y, Omega = NULL, penalty = c("Lasso", "Enet", "Net"), 
  alpha = 1, lambda = NULL, nlambda = 50, rlambda = NULL, nfolds = 1, foldid = NULL,
  inzero = TRUE, adaptive = c(FALSE, TRUE), aini = NULL, isd = FALSE,
  ifast = TRUE, keep.beta = FALSE, thresh = 1e-6, maxit = 1e+5)

loCoxnet(x, y, w, w0 = NULL, h = NULL, hnext = NULL, Omega = NULL,
  penalty = c("Lasso", "Enet", "Net"), alpha = 1, lambda = NULL,
  nlambda = 50, rlambda = NULL, nfolds = 1, foldid = NULL,
  adaptive = c(FALSE, TRUE), aini = NULL, isd = FALSE, keep.beta = FALSE,
  thresh = 1e-6, thresh2 = 1e-10, maxit = 1e+5)
}
\arguments{
  \item{x}{input matrix. Each row is an observation vector.
}
  \item{y}{response variable. \code{y} should be a two-column matrix with columns named `time' and `status'. The latter is a binary variable, with `1' indicating event, and `0' indicating right censored.
}
  \item{w}{input vector, same length as \code{y}. The coefficients vary with \code{w}.
}
  \item{w0}{evaluation local points. The output of estimates are evaludated at these local value \code{w0}. If \code{w0 = NULL}, \code{w0} is generated as 10 equally spaced points in the range of \code{w}.
}
  \item{h}{bandwidth.
}
  \item{hnext}{an increase in bandwidth \code{h}. Default is 0.01.
}
  \item{Omega}{correlation/adjancy matrix with zero diagonal, used for \code{penalty = "Net"} to calculate Laplacian matrix.
}
  \item{penalty}{penalty type. Can choose \code{"Net"}, \code{"Enet"} and \code{"Lasso"}. For \code{"Net"}, need to specify \code{Omega}; otherwises, \code{"Enet"} is performed.
}
  \item{alpha}{ratio between L1 and Laplacian for \code{"Net"}, or between L1 and L2 for \code{"Enet"}. Can be zero and one. For \code{penalty = "Net"}, the penalty is defined as \deqn{\lambda*{\alpha*||\beta||_1+(1-\alpha)/2*(\beta^{T}L\beta)},}
where \eqn{L} is a Laplacian matrix calculated from \code{Omega}. For \code{adaptive = c(FALSE,TRUE)}, its calculation is also based on an initial estimate of \eqn{\beta} from regularized Cox model with \code{penalty = "Enet",alpha = 1} to adjust for signs of coefficents. For \code{penalty = "Enet"}, the penalty is defined as \deqn{\lambda*{\alpha*||\beta||_1+(1-\alpha)/2*||\beta||_2}}.
}
  \item{lambda}{a user supplied decreasing sequence. If \code{lambda = NULL}, a sequency of \code{lambda} is generated based on \code{nlambda} and \code{rlambda}. Supplying a value of \code{lambda} overrides this.
}  
  \item{nlambda}{number of \code{lambda} values. Default is 50.
}
  \item{rlambda}{fraction of \code{lambda.max} to determine the smallest value for \code{lambda}. The default is \code{rlambda = 0.0001} when the number of observations is larger than or equal to the number of variables; otherwise, \code{rlambda = 0.01}.
}
  \item{nfolds}{number of folds. Default is \code{nfolds = 1} and \code{foldid = NULL} and cross-validation is not performed. For cross-validation, smallest value allowable is \code{nfolds = 3}. Specifying \code{foldid} overrisdes this.
}
  \item{foldid}{an optional vector of values between 1 and \code{nfolds} specifying which fold each observation is in.
}
  \item{inzero}{logical flag for simultaneously tuning the number of non-zero coefficients with \code{lambda}. Default is \code{inzero = TRUE}.
}
  \item{adaptive}{logical flags for adaptive version. Default is \code{adaptive = c(FALSE, TRUE)}. The first element is for adaptive on \eqn{\beta} in L1 and the second for adjusting for signs of linked coefficients in Laplacian matrix.
}
  \item{aini}{a user supplied initial estimate of \eqn{\beta}. It is a list including \code{wbeta} for adaptive L1 and \code{sgn} for adjusting Laplacian matrix. \code{wbeta} is the absolute value of inverse initial estimates. If \code{aini = NULL} but \code{adaptive} is required, \code{aini} is generated from regularized Cox model with \code{penalty = "Enet", alpha = 1}.
}
  \item{isd}{logical flag for outputing standardized coefficients. \code{x} is always standardized prior to fitting the model. Default is \code{isd = FALSE}, returning \eqn{\beta} on the original scale.
}
  \item{ifast}{logical flag for efficient calculation of risk set updates. Default is \code{ifast = TRUE}.
}
  \item{keep.beta}{logical flag for returning estimates for all \code{lambda} values. For \code{keep.beta = FALSE}, only return the estimate with the largest cross-validation partial likelihood.
}
  \item{thresh}{convergence threshold for coordinate descent. Default value is \code{1E-6}.
}
  \item{thresh2}{threshold for removing very small \code{lambda} value for local methods. The algorithm computes along a sequence of \code{lambda} value until any absolute value of the second derivative is smaller than \code{thresh2}. The estimates are reported based on this smaller set of \code{lambda} value. Default value is \code{1E-10}.
}
  \item{maxit}{Maximum number of iterations for coordinate descent. Default is \code{10^5}.
}
}
\details{
  One-step coordinate descent algorithm is applied for each \code{lambda}. \code{ifast = TRUE} adopts an efficient way to update risk set and sometimes the algorithm ends before all \code{nlambda} values of \code{lambda} have been evaluated. To evaluate small values of \code{lambda}, use \code{ifast = FALSE}. The two methods only affect the efficiency of algorithm, not the estimates.
  
  Cross-validation partial likelihood is used for tuning parameters. For \code{inzero =  TRUE}, we further select the number of non-zero coefficients obtained from regularized Cox model at each \code{lambda}. This is motivated by formulating L0 variable selection in ADMM form.  
  
  For vayring coefficients methods, the bandwidth is selected by cross-validation.  We recommend to check whether a small increase of \code{h}, say \code{h+hnext}, will improve the current \code{cvm}.
}
\value{
\code{Coxnet} outputs an object with S3 class \code{"Coxnet"}.
  \item{Beta}{estimated coefficients.}
  \item{Beta0}{coefficients after tuning the number of non-zeros, for \code{inzero = TRUE}.}
  \item{fit}{a data.frame containing \code{lambda} and the number of non-zero coefficients \code{nzero}. For cross-validation, additional results are reported, such as average cross-validation partial likelihood \code{cvm} and its standard error \code{cvse}, and \code{index} with \code{max} indicating the largest \code{cvm}.}
  \item{fit0}{a data.frame containing \code{lambda}, \code{cvm} and \code{nzero} based on \code{inzero = TRUE}.}
  
  \item{lambda.max}{value of \code{lambda} that gives maximum \code{cvm}.}
  \item{lambda.opt}{value of \code{lambda} based on \code{inzero = TRUE}.}
  
  \item{cv.nzero}{\code{cvm} with length of number of non-zero components of \code{Beta0}. The kth value of \code{cv.nzero} corresponds to retaining the k largest non-zero coefficients (absolute values) in \code{Beta0}. The optimal number of non-zero is selected by the maximum value of \code{cv.nzero} at \code{lambda = lambda.opt}.}
  \item{penalty}{penalty type.}
  \item{adaptive}{logical flags for adaptive version (see above).}
  \item{flag}{convergence flag (for internal debugging). \code{flag = 0} means converged.}
  
\code{loCoxnet} outputs an object with S3 class \code{"Coxnet"} and \code{"loCoxnet"}.
  \item{Beta}{a list of estimated coefficients with length of \code{lambda}. If there are more than one \code{w0} value, each element of the list is a matrix with \code{p} rows and the number of columns is the length of \code{w0}. If there is one \code{w0}, \code{Beta} is a matrix rather than a list, with \code{p} rows and \code{nlambda} columns.}
  \item{fit}{a data.frame containing \code{lambda} and the number of non-zero coefficients \code{nzero}. For cross-validation, additional results are reported, such as average cross-validation partial likelihood \code{cvm} and its standard error \code{cvse}, and \code{index} with \code{max} indicating the largest \code{cvm}.}
  
  \item{lambda.max}{value of \code{lambda} that gives maximum \code{cvm}.}
  
  \item{cvh}{a data.frame containing bandwidth, \code{cvm} and \code{cvse}.}
    
  \item{penalty}{penalty type.}
  \item{adaptive}{logical flags for adaptive version (see above).}
  \item{flag}{convergence flag (for internal debugging). \code{flag = 0} means converged.}
}
\author{
Xiang Li, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xl2473@columbia.edu>
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent, Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
  Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011)
  \emph{Regularization Paths for Cox's Proportional Hazards Model via
    Coordinate Descent, Journal of Statistical Software, Vol. 39(5)
    1-13}\cr
  \url{http://www.jstatsoft.org/v39/i05/}\cr
  Sun, H., Lin, W., Feng, R., and Li, H. (2014)
  \emph{Network-regularized high-dimensional cox regression for analysis of genomic data, Statistica Sinica.}\cr
  \url{http://www3.stat.sinica.edu.tw/statistica/j24n3/j24n319/j24n319.html}\cr
  van Houwelingen, H. C., Bruinsma, T., Hart, A. A., van't Veer, L. J., & Wessels, L. F. (2006)
  \emph{Cross-validated Cox regression on microarray gene expression data. Statistics in medicine, 25(18), 3201-3216.}\cr
  \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.2353/full}\cr
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.Coxnet}}, \code{\link{coxsplit}}
}
\examples{
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]%*%beta
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1)  # censoring indicator
y=cbind(time=ty,status=1-tcens)

fiti=Coxnet(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)
}
\keyword{Cox model}
\keyword{Regularization}
