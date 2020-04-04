\name{CobbDouglas_boot}
\alias{CobbDouglas_boot}
\title{Bootstrap confidence intervals for a Cobb-Douglas frontier}
\description{Boostrap resampling to approximate confidence intervals for parameters and fitted values of a Cobb-Douglas production frontier.}
\usage{CobbDouglas_boot(x, nboot=500, conf=0.95)}
\arguments{
  \item{x}{An object of class \code{CobbDouglas}.}
  \item{nboot}{The number of bootstrap replications. It must be at least 50.}
  \item{conf}{The confidence level. Default is 0.95.}
}
\value{An object of class \code{CobbDouglas_boot}, that is a list with the following components:
  \item{\code{parameters}}{Bootstrap confidence intervals at level \code{conf} for the parameters.}
  \item{\code{fitted}}{Bootstrap confidence intervals at level \code{conf} for the fitted values.}
  %\item{\code{boot.parameters}}{The boostrap simulation for the parameters.}
  %\item{\code{boot.fitted}}{The bootstrap simulation for the fitted values.}
}
%\note{}
%\references{}
\seealso{\link{CobbDouglas}.}
\examples{
data(proc)

m2 <- CobbDouglas("output", c("labour","capital"), data=proc)
set.seed(123)
CobbDouglas_boot(m2, nboot=150)

m2c <- CobbDouglas("output", c("labour","capital"), data=proc, beta.sum=1)
set.seed(123)
CobbDouglas_boot(m2c, nboot=150)

m2d <- CobbDouglas("output", c("labour","capital"), data=proc, beta.sum=0.7)
set.seed(123)
CobbDouglas_boot(m2d, nboot=150)
}
