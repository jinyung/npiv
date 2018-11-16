\name{flux}
\alias{flux}
\alias{side_flux}
\title{Calculate flux}
\description{
Calculate 2D 'flux' over a line segment.
}
\usage{
flux(body_land, vf_dir, id, vect_scale, body_width, scal, 
     plot = FALSE, relative = TRUE, t_stamp)
                 
side_flux(body_land, vf_dir, id, vect_scale, body_height, scal, 
          relative = TRUE, t_stamp, plot = FALSE)
}
\arguments{
  \item{body_land}{Defined differently for \code{side_flux}, see details}
  \item{vf_dir}{(char) path to directory containing the vector field output from
    DaVis}
  \item{id}{(int) give the index to frames to be read for calculations}
  \item{vect_scale}{(numeric) scaling factor to scale the vector of 
    centroid-to-tail vector. The end of the scaled vector is where the center 
    of the flux line locate}
  \item{body_width}{(numeric) this control the length of the flux line, should be
    given in pixel to match \code{body_land}}
  \item{body_height}{(numeric) similar to \code{body_width} for \code{side_flux}}
  \item{scal}{scale to be used for width/ height, should be in same unit as 
    displacement used for calculation of vector field}
  \item{plot}{(logical) if \code{TRUE} produce plot with \code{body_land}, 
    flux line and interpolated/projected vectors on flux line}
  \item{relative}{(logical) if \code{TRUE} calculation will be based on relative
    velocity rather than absolute velocity}
  \item{t_stamp}{(numeric) time, should be in same unit as time differences used
    for calculation of vector field}
}

\value{
(numeric) return flux in area/time
}

\details{
It is assumed that body landmarks are unscaled (in unit of pixels) when output 
from \code{TPSDig2} but vector fields from \code{DaVis} are scaled. \code{body_width} / \code{body_height} is expected to be unscaled too.

Relative velocity use moving body of nauplius as frame of reference while
absolute velocity just takes the velocity fields as they are: using earth as 
frame of reference. Relative velocity makes use of \code{body_land} to calculate
body velocity and substract vectors' velocity from this calculated body velocity.

\code{side_flux}'s \code{body_land} also contains 3 landmarks, but the landmarks 
are: 1) posterior margin along long axis; 2) anterior margin along long axis; 
3) a point indicating the direction of ventral side, i.e. a point anywhere on 
the right side of the posterior-anterior vector. The third landmark guides the extension of the flux line in \code{side_flux}.
}