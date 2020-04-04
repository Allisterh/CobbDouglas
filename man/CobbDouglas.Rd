\name{CobbDouglas}
\alias{CobbDouglas}
\title{Estimation of a Cobb-Douglas production frontier}
\description{Estimation of a Cobb-Douglas production frontier from sample data.}
\usage{CobbDouglas(y.name, x.names=NULL, data, beta.sum=NULL)}
\arguments{
  \item{y.name}{The name of the output variable.}
  \item{x.names}{The names of the input variables. If \code{NULL} (the default), it is set equal to the name of the variables in \code{data} besides \code{y.name}.}
  \item{data}{A data.frame containing the output and the input variables.}
  \item{beta.sum}{Constraint on the sum of beta parameters. If \code{NULL} (the default), beta parameters are freely estimated (subjected only to positive constraint).}
}
\value{An object of class \code{CobbDouglas}, that is a list with the following components:
  \item{\code{parameters}}{Parameter estimates.}
  \item{\code{efficiency}}{Technical efficiencies of the sample units.}
  \item{\code{fitted}}{Fitted values on both logarithmic and original scale.}
  \item{\code{residuals}}{Residuals on both logarithmic and original scale.}
  \item{\code{beta.sum}}{Value passed to argument \code{beta.sum}.}
  \item{\code{y.name}}{Value passed to argument \code{y.name}.}
  \item{\code{x.names}}{Value passed to argument \code{x.names}.}
  \item{\code{data}}{Data used in the estimation.}
}
\details{
Consider a sample of \eqn{n} production units, for which the quantity of the output \eqn{Y} and of \eqn{H} input variables \eqn{X_1,\ldots,X_H} is measured.
Let \eqn{y_i} be the quantity of the output for unit \eqn{i},
and \eqn{x_{hi}} be the quantity of the \eqn{h}-th input for unit \eqn{i}.
A Cobb-Douglas production frontier is defined as:
\deqn{y^*_i = \tau\prod_{h=1}^H x_{hi}^{\beta_h}
}
where \eqn{y^*_i} is the maximum producible output for unit \eqn{i},
\eqn{\tau} is a parameter representing the total factor productivity for a technically efficient unit,
and \eqn{\beta_h} (\eqn{h=1,\ldots,H}) is a parameter representing the elasticity of the output with respect to the \eqn{h}-th input.

Constant returns to scale holds if \eqn{\sum_{h=1}^H \beta_h=1} (obtained setting the argument \code{beta.sum} to value 1).
Instead, \eqn{\sum_{h=1}^H \beta_h<1} implies decreasing returns to scale, while
\eqn{\sum_{h=1}^H \beta_h>1} implies increasing returns to scale.
For \eqn{i=1,\ldots,n}, quantities \eqn{y^*_i} are denoted as \emph{fitted values},
while quantities \eqn{e_i=y_i-y^*_i} are denoted as \emph{residuals}.

Estimation of the Cobb-Douglas production frontier is performed through constrained least squares on the logarithmic scale:
\deqn{(\hat{\tau},\hat{\beta}_1,\ldots,\hat{\beta}_H)=\code{argmin}_{\tau,\beta_1,\ldots,\beta_H}\sum_{i=1}^n(\log y_i-\log y^*_i)^2}
where:
\deqn{\log y^*_i = \log\tau+\sum_{h=1}^H \beta_h \log x_{hi}}
subjected to constraints:
\deqn{\log y^*_i \geq \log y_i \hspace{.7cm} i=1,\ldots,n}
\deqn{\beta_h \geq 0 \hspace{.7cm} h=1,\ldots,H}

S3 methods available for class \code{CobbDouglas} are:
\itemize{
  \item{\code{print}: }{ to get essential information.}
  \item{\code{summary}: }{ to get summaries of estimation.}
  \item{\code{plot}: }{ to display the scatterplot with the estimated frontier (only for frontiers with a single input).}
  \item{\code{predict}: }{ to predict the maximum producible output or technical efficiency. See \link{predict.CobbDouglas}.}
  }
Also, the method \link{CobbDouglas_boot} is available to approximate confidence intervals of parameters and fitted values.
}
%note{}
\references{
C. W. Cobb and P. H. Douglas (1928). A theory of production. \emph{American Economic Review}, 18: 139-165.
}
\seealso{\link{predict.CobbDouglas}; \link{CobbDouglas_boot}.}
\examples{
data(proc)


### 1 input variable: 'labour'

m1 <- CobbDouglas("output", "labour", data=proc)
summary(m1)

# plot the estimated frontier
plot(m1, cex.axis=1.1, cex.lab=1.2)

# technical efficiencies
m1_eff <- m1$efficiency
## NOT RUN:
# m1_eff

# efficient units
m1_eff[which(m1_eff$y.side==1),]


### 1 input variable: 'labour' + constraint on beta

# beta=1 (constant returns to scale)
m1c <- CobbDouglas("output", "labour", data=proc, beta.sum=1)
m1c$parameters
m1c$efficiency[which(m1c$efficiency$y.side==1),]
plot(m1c, cex.axis=1.1, cex.lab=1.2, main="beta = 1", cex.main=1.6)

# beta=1.25 (increasing returns to scale)
m1i <- CobbDouglas("output", "labour", data=proc, beta.sum=1.25)
m1i$parameters
m1i$efficiency[which(m1i$efficiency$y.side==1),]
plot(m1i, cex.axis=1.1, cex.lab=1.2, main="beta = 1.25", cex.main=1.6)

# beta=0.3 (decreasing returns to scale)
m1d <- CobbDouglas("output", "labour", data=proc, beta.sum=0.3)
m1d$parameters
m1d$efficiency[which(m1d$efficiency$y.side==1),]
plot(m1d, cex.axis=1.1, cex.lab=1.2, main="beta = 0.3", cex.main=1.6)


### 2 input variables: 'labour' and 'capital'

# no constraints on the sum of beta parameters
m2 <- CobbDouglas("output", c("labour","capital"), data=proc)
summary(m2)
m2$efficiency[which(m2$efficiency$y.side==1),]

# beta.sum=1 (constant returns to scale)
m2c <- CobbDouglas("output", c("labour","capital"), data=proc, beta.sum=1)
summary(m2c)
m2c$efficiency[which(m2c$efficiency$y.side==1),]

# beta.sum=0.7 (decreasing returns to scale)
m2d <- CobbDouglas("output", c("labour","capital"), data=proc, beta.sum=0.7)
summary(m2d)
m2d$efficiency[which(m2d$efficiency$y.side==1),]
}
