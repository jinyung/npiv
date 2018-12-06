\name{npiv-package}
\alias{npiv-package}
\alias{npiv}
\docType{package}
\title{Kinematics and hydrodynamics analysis on nauplii}
\description{
This is a companion package to the manuscript:

Wong, J. Y., Chan, B. K. K., Chan, K. Y. K. Swimming kinematics of naupliar 
larvae with contrasting feeding modes.

To reproduce the results, please refer to readme at Github repository for 
explanation on how to use the package and updates.
}

\details{
Functions included in the package were developed for analyses of data generated 
from \href{http://life.bio.sunysb.edu/morph/soft-dataacq.html}{\code{TPSDig2}} 
and \href{https://www.lavision.de/en/products/davis-software/}{\code{DaVis}} 
softwares for kinematics and hydrodynamics analyses, respectively. 
\code{TPSDig2} is a landmark registration tool and was used for digitization of 
body and appendages positions. \code{DaVis} is a software used for particle 
image velocimetry, and outputs velocity and vorticity fields. There is no 
utility tools for direct data import from other softwares but all underlying 
calculations use \code{R}'s \code{\link{matrix}} or \code{\link{array}} classes 
of data.

Most of the functions were written specifically for the manuscript, especially 
those for kinematics analyses as they rely on the specific definitions of 
landmarks on nauplius larva. However, some functions such as 
\code{\link{calc_flux}} and \code{\link{calc_r}}, are suitable for general 
use with any velocity field data.
}