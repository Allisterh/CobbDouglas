\name{CobbDouglas_boot}
\alias{CobbDouglas_boot}
\title{Bootstrap confidence intervals for a Cobb-Douglas frontier}
\description{Boostrap resampling to approximate confidence intervals for parameters and technical efficiencies of a Cobb-Douglas production frontier.}
\usage{CobbDouglas_boot(x, nboot=500, conf=0.95)}
\arguments{
  \item{x}{An object of class \code{CobbDouglas}.}
  \item{nboot}{The number of bootstrap replications. It must be at least 50.}
  \item{conf}{The confidence level. Default is 0.95.}
}
\value{An object of class \code{CobbDouglas_boot}, that is a list with three components: \code{parameters}, \code{y.side} and \code{x.side}. Each component is a data.frame containing the estimated value and the bootstrap confidence interval at level \code{conf} for each parameter (component \code{parameters}), output-side technical efficiencies (component \code{y.side}) and input-side technical efficiencies (component \code{x.side}).}
\details{Non-parameteric bootstrap resampling is performed (Efron, 1979) and bias-corrected boostrap confidence intervals are computed (Efron, 1987).}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
B. Efron (1987). Better Bootstrap Confidence Intervals. \emph{Journal of the American Statistical Association}, 82(397): 171-185.

B. Efron (1979). Bootstrap Methods: Another Look at the Jackknife. \emph{The Annals of Statistics}, 7(1): 1-26.
}
\seealso{\link{CobbDouglas}.}
\examples{
data(production)

m2 <- CobbDouglas("output", c("labour","capital"), data=production)
set.seed(123)
CobbDouglas_boot(m2, nboot=150)

m2c <- CobbDouglas("output", c("labour","capital"), data=production, beta.sum=1)
set.seed(123)
CobbDouglas_boot(m2c, nboot=150)

m2d <- CobbDouglas("output", c("labour","capital"), data=production, beta.sum=0.7)
set.seed(123)
CobbDouglas_boot(m2d, nboot=150)
}
