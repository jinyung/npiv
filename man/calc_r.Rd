\name{calc_r}
\alias{calc_r}
\alias{calc_s}
\title{Calculate area of influence}
\description{
Calculate area of influence (and radius of circle with equivalent area) from 
generated flow.
}
\usage{
calc_r(dir, idlist, thre = c(0.5, 0.95), thre.len = 15, val.area)

calc_s(dir, idlist, ustar = 0.0005, val.area)
}
\arguments{
  \item{dir}{(char) directory path containing the vector field output files}
  \item{idlist}{(int) index to frames for calculations}
  \item{thre}{(num) threshold values to adjust U*, two values are required}
  \item{thre.len}{(int) number of U* values to use}
  \item{val.area}{(num) area of each velocity vector value covers, by default  
    is \eqn{16px \times 16px \times (\textnormal{scale in mm/px})^2}}
  \item{ustar}{(num) theshold value for area calculation used in \code{cal_s}}
}
\value{
  \item{r}{(num) radius of circle with equivalent area}
  \item{S}{(num) area of influence}
  \item{u}{(num) theshold values used by \code{calc_r}}
}
\details{
\code{calc_s} is intended for calculation of area of influence from flow 
exceeding a single threshold value.

\code{calc_r} is intended for calculation of area of influence over a series of 
threshold values, facilitating the calculation of flow attenuation over distance 
by using another function: \code{\link{calc_power}}. Values output from 
\code{calc_r} (\code{r} and \code{u}) are in matrix, with the second dimension 
being the values of theshold used in calculations.
}

\seealso{
\code{\link{calc_power}}, \code{\link{plot_ur}}
}