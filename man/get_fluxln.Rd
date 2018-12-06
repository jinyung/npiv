\name{get_fluxln}
\alias{get_fluxln}
\title{Create a 'flux line'}
\description{
Compute line segment for 2D flux calculation
}
\usage{
get_fluxln(body_land, vect_scal = 1, ln_width, ln_width_scal = 1, 
  type = c('top', 'side'))
}
\arguments{
  \item{body_land}{(matrix) landmarks on nauplius carapace}
  \item{vect_scal}{(num) scaling factor for centroid-tail vector, see details}
  \item{ln_width}{(num) for top view, carapace width; for side view, carapace 
    height}
  \item{ln_width_scal}{(num) scaling factor for \code{ln_width}}
  \item{type}{(char) choices of top or side view}
}
\value{
xy coordinates of two points define the line segment in \code{\link{xy.coords}} 
format.
}
\details{
For \code{type = 'top'}, the \emph{center point} of the line segment is by default 
the end point of the centroid-to-tail vector. This position of the center point 
can be adjusted by changing \code{vect_scal}, which scale the centroid-to-tail 
vector (extend when >1 and shrink when <1). The length of the line segment is by
default the width of the carapace (specified by \code{ln_width}), but can be 
scaled by \code{ln_width_scal}.

For \code{type = 'side'}, the \emph{starting point} of the line segment is by 
default the end point of the centroid-to-tail vector. Similarly, the position of 
the starting point can be adjusted by changing \code{vect_scal}. The length of 
the line segment is by default the height of the carapace (specified by 
\code{ln_width}), but can be scaled by \code{ln_width_scal}.
}
\seealso{
\code{\link{calc_flux}}, \code{\link{calc_refv}}
}