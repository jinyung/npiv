\name{npiv-internal}
\alias{npiv-internal}
\alias{angle}
\alias{arrow}
\alias{calc_u}
\alias{ceiling_dec}
\alias{cent}
\alias{dif}
\alias{equi_rad}
\alias{floor_dec}
\alias{get_file}
\alias{lseq}
\alias{makelist}
\alias{max_u}
\alias{perm_ttest}
\alias{plotvf}
\alias{range01}
\alias{se}
\alias{unitv}
\alias{vf2surface}
\alias{xycoords2mat}
\title{Miscellaneous small functions}
\description{
  Small functions lacking their own dedicated documentation:
  
  \describe{
    \item{\code{angle}}{calculate angle between target vector and a vector 
      extends from starting point to positive direction of x-axis}
    \item{\code{arrow}}{a different way to plot arrow}
    \item{\code{calc_u}}{calculate speed from velocity vectors output read
      with \code{\link{read_davis}}}
    \item{\code{ceiling_dec}}{modified \code{ceiling} function, rounded up to 
      decimal point specified}
    \item{\code{cent}}{calculate centroid}
    \item{\code{dif}}{modified \code{diff}, which keep length of vector by 
      inserting NAs}
    \item{\code{equi_rad}}{given area, calculate radius of circle of equivalent 
      area}
    \item{\code{floor_dec}}{modified \code{floor} function, rounded down to 
      decimal point specified}
    \item{\code{get_file}}{utility tool to get i-th frame filename from Davis 
      output}
    \item{\code{lseq}}{modified \code{seq}, to get logarithmic sequence}
    \item{\code{makelist}}{make an empty named list}
    \item{\code{max_u}}{utility tool to find the peak fluid speed (at certain
      quantile) of a frame}
    \item{\code{perm_ttest}}{perform permutational T-test, two tailed}
    \item{\code{plotvf}}{plot vector field using data read from 
      \code{\link{read_davis}}}
    \item{\code{range01}}{rescale the values into range of 0-1, with min value 
      as 0 and max as 1}
    \item{\code{se}}{calculate standard error}
    \item{\code{unitv}}{rescale vector into unit vector, i.e. divide the vector 
      by its magnitude}
    \item{\code{vf2surface}}{turn format by \code{\link{read_davis}} into 
      \code{fields} package \code{surface}'s format}
    \item{\code{xycoords2mat}}{turn list of x and y into a 2 column xy matrix}
  }
}
\details{
Explanations on how to use these functions were documented in the comments of the 
source codes in
\code{file.path(system.file("R", package = "npiv"), "npiv-internal.R")}
}

\references{
\code{floor_dec}/ \code{ceiling_dec} came from a solution provided in 
\url{https://stackoverflow.com/a/39611375}

\code{lseq} came from a solution provided in 
\url{https://stackoverflow.com/a/29963530}
}