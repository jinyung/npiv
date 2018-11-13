\name{ang_sep}
\alias{ang_sep}
\title{Calculate angular separation}
\description{
Given time series on appendages' positions, calculate their angular positions, 
stroke amplitude and separation between appendages.}
\usage{
ang_sep(body_land, limb_land)
}
\arguments{
  \item{body_land}{(3d-array) body landmarks or contours to calculate
            centroid, 3rd dimension is the frames}
  \item{limb_land}{(3d-array) limb landmarks, 3rd dimension is the frames}
}
\value{
  \item{ang_sep}{(matrix) angular separations between ant1-ant2; ant1-mand; 
    ant2-mand}
  \item{amplitude}{(num) amplitudes for ant1, ant2 and mand}
  \item{langle}{(matrix) angle in degree for ant1, ant2 and mand}
}
\details{
Abbreviations for the appendages:
\itemize{
  \item \code{ant1} = antennule
  \item \code{ant2} = antenna
  \item \code{mand} = mandibles
}
Please refer to the details of \code{\inst\python\kinematics.py} for more details of the underlying functions (\code{vangle}, \code{sangle}, \code{langle}, \code{tsangle}, \code{tlangle}).
}
}