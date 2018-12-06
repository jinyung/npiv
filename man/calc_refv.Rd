\name{calc_refv}
\alias{calc_refv}
\title{Calculate moving velocity}
\description{
Calculate moving velocity of nauplius
}
\usage{
calc_refv(body_land_array, type = c('top', 'side'), scal, t_stamp) 
}
\arguments{
  \item{body_land_array}{(3d array) landmarks on nauplius carapace, with 
    3rd dimension is the frames}
  \item{type}{(char) data from either top or side view}
  \item{scal}{(num) scale in length/pixel}
  \item{t_stamp}{(num) time of each frame}
}

\value{
(num) moving velocity
}

\note{
This function will be combined with \code{\link{disp}} in next release.
}