# CobbDouglas
__Efficiency analysis using the Cobb-Douglas production frontier__

`CobbDouglas` is an R package implementing functionalities to estimate a Cobb-Douglas production frontier from sample data through constrained least squares.
It is possible to set the desired returns to scale, to predict the maximum producible
output or technical efficiency, and to obtain boostrap confidence intervals for parameters.

R (The R Project for Statistical Computing) needs to be installed on your system in order
to use the `CobbDouglas` package. R can be downloaded from https://www.r-project.org/.

To install the `CobbDouglas` package, open the console of R and type:
```
install.packages("devtools")  ## do not run if package 'devtools' is already installed
library(devtools)
install_github("alessandromagrini/CobbDouglas")
```

For any request or feedback, please write to <alessandro.magrini@unifi.it> (Alessandro Magrini)

Below, you find some examples of use of the package.
_________________________________________________________________

Load simulated data
```
data(production)
```
Frontier with 1 input variable: 'labour'
```
m1 <- CobbDouglas("output", "labour", data=production)
summary(m1)

# plot the estimated frontier
plot(m1, cex.axis=1.1, cex.lab=1.2)

# technical efficiencies
m1_eff <- m1$efficiency
## NOT RUN:
# m1_eff

# efficient units
m1_eff[which(m1_eff$y.side==1),]

# setting beta=1 (constant returns to scale) seems not so good
m1c <- CobbDouglas("output", "labour", data=production, beta.sum=1)
m1c$parameters
m1c$efficiency[which(m1c$efficiency$y.side==1),]
plot(m1c, cex.axis=1.1, cex.lab=1.2, main="beta = 1", cex.main=1.6)
```
Frontier with 2 input variables: 'labour' and 'capital'
```
# no constraints on the sum of beta parameters
m2 <- CobbDouglas("output", c("labour","capital"), data=production)
summary(m2)
m2$efficiency[which(m2$efficiency$y.side==1),]

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

# beta.sum=1 (constant returns to scale)
m2c <- CobbDouglas("output", c("labour","capital"), data=production, beta.sum=1)
summary(m2c)
m2c$efficiency[which(m2c$efficiency$y.side==1),]

# beta.sum=0.7 (decreasing returns to scale)
m2d <- CobbDouglas("output", c("labour","capital"), data=production, beta.sum=0.7)
summary(m2d)
m2d$efficiency[which(m2d$efficiency$y.side==1),]
```
Prediction of the maximum producible output
```
predict(m2, newdata=data.frame(labour=20,capital=5))
```
Prediction of technical efficiency
```
predict(m2, newdata=data.frame(output=15,labour=20,capital=5), type="eff")
```
Bootstrap confidence intervals
```
set.seed(123)
CobbDouglas_boot(m2, nboot=500)
```
