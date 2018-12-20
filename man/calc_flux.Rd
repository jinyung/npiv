\name{calc_flux}
\alias{calc_flux}
\alias{tcalc_flux}
\alias{proj_v}
\title{Calculate flux}
\description{
Calculate 2D 'flux' over a line segment.
}
\usage{
calc_flux(FUN = get_fluxln, vf_path, side = c('right', 'left'),
  refv = 0, width, scal, ...) 
  
tcalc_flux(body_land_array, vf_dir, id, vect_scal, 
  type = c('top', 'side'), relative = TRUE, ln_width, 
  ln_width_scal, side = c('right', 'left'), t_stamp, scal) 
  
proj_v(xycoord, vf_path, nseg = 8, side = c('right', 'left'))
}

\arguments{
  \item{FUN}{(function) function to give xy coordinates of line segment for flux 
    calculation, and should in R's \code{\link{xy.coords}} format}
  \item{vf_path}{(char) path to file containing the vector field}
  \item{side}{(char) there are two possible normal directions to line segment, 
    side specify which to output, see details}
  \item{refv}{(num) velocity of frame of reference, see details. Single value 
    for \code{calc_flux} but multiple values for \code{tcalc_flux}}
  \item{width}{(num) length of the flux line segment}
  \item{...}{arguments to be passed into \code{FUN}}
  \item{body_land_array}{(3d array) array of landmarks of nauplius body}
  \item{vf_dir}{(char) path to directory containing files of vector field}
  \item{id}{(int) vector of integers specifying the frames id of vector fields}
  \item{vect_scal}{(num) passed to \code{link{get_fluxln}}, single value}
  \item{type}{(char) either for top view or side view videos}
  \item{relative}{(logical) if \code{TRUE}, \code{refv} will be calculated from 
    \code{body_land_array}}
  \item{ln_width}{passed to \code{\link{get_fluxln}}, single value}
  \item{ln_width_scal}{passed to \code{\link{get_fluxln}}, single value}
  \item{t_stamp}{passed to \code{\link{calc_refv}}}
  \item{scal}{scale in length/ pixel to scale \code{width}}
  \item{xycoord}{xy coordinates of two points defining the line segment to 
    calculate flux, should in \code{R \link{xy.coords}} format}
  \item{nseg}{(int) how many segments the flux line consists of, arbitrarily set}
}
\value{
\code{calc_flux} returns single value of 2D flux

\code{tcalc_flux} returns multiple values of 2D flux over frames

\code{proj_v} returns \code{nseg} number of projected vectors 
}
\details{
\code{tcalc_flux} is a wrapper of \code{calc_flux}, and \code{calc_flux} is a 
wrapper of \code{proj_vf}.

Given xy coordinates of two points in the vector field, \code{proj_v} 
interpolate vectors on points on line segment formed by these two points, 
project them on normal direction, and returns their magnitude with direction 
indicated by sign. There are two possible unit normal vectors to the line 
segment, and determined as follows: 1) the two points forming the line segment 
is viewed as a vector with first point as the starting point and second as the 
ending point; 2) Rotate this vector 90 degrees clockwise and anticlockwise gives 
two normal vectors. Argument \code{side} of \code{'right'} outputs the projected 
values of clockwise normal direction and \code{'left'} outputs the anticlockwise 
counterpart. Changing the side only flip the sign of the values.

\code{calc_flux} takes the mean of the magnitudes of the projected vectors and 
times it with the length of the line segment, giving 2D 'flux'. To change the 
frame of reference, \code{refv} can be specified to give the velocity of frame 
of reference (e.g. moving body of nauplius). \code{refv} should be velocity that 
has only direction in one axis, i.e. in the axis normal to the flux line. Thus 
it is only a scalar value with plus minus sign, and the sign should be 
designated in the same direction as that of normal vector specified by 
\code{side}.

\code{tcalc_flux} is the wrapper function used in the manuscript to batch 
process 2D flux over many frames of videos, specifically use the flux line 
definition described in the manuscript (see \code{\link{get_fluxln}}). For 
generalized use, \code{calc_flux} is the recommended function, which allows 
user to insert custom definition of flux line.
}
\seealso{
\code{\link{get_fluxln}}, \code{\link{calc_refv}}
}