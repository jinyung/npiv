\name{calc_power}
\alias{calc_power}
\title{Calculate power of attenuation}
\description{
Given distance and speed, calculate the the attenuation power with respect to 
distance.
}
\usage{
calc_power(r, u)
}
\arguments{
  \item{r}{(num) distance, such as radius of area of influence}
  \item{u}{(num) speed}
}
\value{
return the attenuation power (slope) from log-log regression.
}
\details{
For fluid speed attenuates exponentially with distance following the equation:

\deqn{u = a \cdot r^b}

the attenuation power \code{b} can be found by utilizing log-log regression. 

By taking log both sides: \eqn{log(u) = log(a \cdot r^b)} and rewrite the RHS 
into \eqn{log(a) + log(r^b)} and into \eqn{log(a) + b \cdot log(r)}, the 
expression is in the form of a simple linear regression of \eqn{Y = a+bx}:

\deqn{log(u) = log(a) + b \cdot log(r)} 

where \code{b} can be obtained from finding the slope of the expression.
}
\seealso{
\code{\link{calc_r}}, \code{\link{plot_ur}}
}