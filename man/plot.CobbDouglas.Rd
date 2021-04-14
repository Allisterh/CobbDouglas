\name{plot.CobbDouglas}
\alias{plot.CobbDouglas}
\title{Graphic for an estimated Cobb-Douglas frontier}
\description{Display the margin of an estimated Cobb-Douglas frontier with respect to a selected input variable.}
\usage{
\method{plot}{CobbDouglas}(x, x.name=NULL, x.fixed=NULL, add.grid=TRUE, add.points=TRUE, add.legend=TRUE,
  cex.legend=0.9, digits=3, xlab=NULL, ylab=NULL, col="red", ...)
}
\arguments{
  \item{x}{Object of class \code{CobbDouglas}.}
  \item{x.name}{Character string indicating the name of the input variable to consider. If \code{NULL} (the default), the first input variable is selected.}
  \item{x.fixed}{Named numerical vector including the values at which input variables besides the considered one should be held constant. If \code{NULL} (the default) or invalid, the empirical means are used. Ignored if there is a single input variable.}
  \item{add.grid}{Logical value indicating whether a grid should be added to the graphic. Default is \code{TRUE}.}
  \item{add.points}{Logical value indicating whether the observed points should be added to the graphic. Default is \code{TRUE}. Ignored if there are more than one input variable.}
  \item{add.legend}{Logical value indicating whether a legend showing the values of the input variables held constant should be added to the graphic. Default is \code{TRUE}.}
  \item{cex.legend}{Size of the legend. Default is 0.9.}
  \item{digits}{Integer non-negative value indicating the number of decimal places to be used in the legend. Default is 3. Ignored if \code{add.legend}=\code{FALSE}.}
  \item{ylab}{Text for y-axis (optional).}
  \item{xlab}{Text for x-axis (optional).}
  \item{col}{The color for the frontier. Default is \code{"red"}.}
  \item{...}{Further graphical parameters.}
}
%\value{}
%\details{}
%\note{}
%\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
%\references{}
\seealso{\link{CobbDouglas}.}
\examples{
data(production)


### one input variable: 'labour'

m1 <- CobbDouglas("output", "labour", data=production)
summary(m1)

# plot the estimated frontier
plot(m1, cex.axis=1.1, cex.lab=1.2)


### two input variables: 'labour' and 'capital'

m2 <- CobbDouglas("output", c("labour","capital"), data=production)
summary(m2)

## plot the estimated frontier from the side of 'labour'
# 'capital' fixed at its empirical mean
plot(m2, x.name="labour", cex.axis=1.1, cex.lab=1.2)
# 'capital' fixed at value 1
plot(m2, x.name="labour", x.fixed=c(capital=1), cex.axis=1.1, cex.lab=1.2)

## plot the estimated frontier from the side of 'capital'
# 'labour' fixed at its empirical mean
plot(m2, x.name="capital", cex.axis=1.1, cex.lab=1.2)
# 'labour' fixed at value 1
plot(m2, x.name="capital", x.fixed=c(labour=1), cex.axis=1.1, cex.lab=1.2)
}
