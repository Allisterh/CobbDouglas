\name{CobbDouglas-package}
\alias{CobbDouglas-package}
\docType{package}
\title{Estimation and efficiency analysis for the Cobb-Douglas production frontier}
\description{Functionalities to estimate a Cobb-Douglas production frontier from sample data through constrained least squares
(possibly setting the desired returns to scale),
to predict the maximum producible output or technical efficiency,
and to obtain boostrap confidence intervals for parameters and fitted values.}
\details{
  \tabular{ll}{
  Package: \tab CobbDouglas\cr
  Type: \tab Package\cr
  Version: \tab 1.0\cr
  Date: \tab 2020-04-04\cr
  License: \tab GPL-2
  }
The main functions of the package are:
\itemize{
  \item{\link{CobbDouglas}, }{to estimate the frontier from sample data;}
  \item{\link{predict.CobbDouglas}, }{to predict the maximum producible output or technical efficiency;}
  \item{\link{CobbDouglas_boot}, }{to obtain bootstrap confidence intervals for parameters and fitted values.}
  }
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
C. W. Cobb and P. H. Douglas (1928). A Theory of Production. \emph{American Economic Review}, 18: 139-165.
}
