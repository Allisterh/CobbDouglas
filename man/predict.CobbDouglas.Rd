\name{predict.CobbDouglas}
\alias{predict.CobbDouglas}
\title{Prediction using a Cobb-Douglas production frontier}
\description{Predict the maximum producible output or technical efficiency using a Cobb-Douglas production frontier.}
\usage{
\method{predict}{CobbDouglas}(object, newdata=NULL, type="output", ...)
}
\arguments{
  \item{object}{An object of class \code{CobbDouglas}.}
  \item{newdata}{A data.frame in which to look for variables with which to predict the maximum producible output (if \code{type="output"}) or technical efficiency (if \code{type="efficiency"}). If \code{NULL} (the default), fitted values or technical efficiencies of the sample units are returned.}
  \item{type}{The type of prediction: \code{"output"} (maximum producible output)) or \code{"efficiency"} (technical efficiency). It can be abbreviated.}
  \item{...}{Further arguments passed to the generic \code{predict} method.}
}
\value{An object of class \code{data.frame}.}
%\note{}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
C. W. Cobb and P. H. Douglas (1928). A Theory of Production. \emph{American Economic Review}, 18: 139-165.
}
\seealso{\link{CobbDouglas}.}
\examples{
data(production)
m2 <- CobbDouglas("output", c("labour","capital"), data=production)

# prediction of the maximum producible output
predict(m2, newdata=data.frame(labour=20,capital=5))

# prediction of technical efficiency
predict(m2, newdata=data.frame(output=15,labour=20,capital=5), type="eff")
}
