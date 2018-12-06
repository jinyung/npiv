\name{plot_ur}
\alias{plot_ur}
\title{Plot flow speed against distance}
\description{
Plot log-log regression fit of flow speed against distance (radius of circle 
with equivalent area.)
}
\usage{
plot_ur(result, col, ...)
}
\arguments{
  \item{result}{(list) returned values from \code{\link{calc_r}}}
  \item{col}{color for each time frame, default is rainbow color series if not 
    specified}
  \item{...}{arguments passed into underlying \code{\link{matplot}}function}
}

\value{
Produce a plot
}

\details{
Plot came with regression fit lines of nearest rounded values (up and down) of 
attenuation power.
}
\seealso{
\code{\link{calc_r}}, \code{\link{calc_power}}
}