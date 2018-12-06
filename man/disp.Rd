\name{disp}
\alias{disp}
\title{Calculate displacement}
\description{
Calculate forward/backward displacement
}
\usage{
disp(body_land, keep_len = TRUE, win_size)
}
\arguments{
  \item{body_land}{(3d array) landmarks on nauplii carapace}
  \item{keep_len}{(logical) if \code{TRUE}, \code{NA} will be inserted to keep 
    the length of displacement to be same as length of time frames}
  \item{win_size}{(int) size of the moving window to compute difference, 
    default to 3}
}

\value{
(num) displacement, with positive value = forward displacement and negative 
  value = reverse displacement
}

\details{
Direction of displacement is determined by dot product between displacement 
vector and vector from centroid to tail processes (pointing in posterior 
direction).

Calculation of displacement depends on computation of lag differences 
(see \code{\link{diff}}) and thus the returned vector of values is shorter than 
the number of input time frames. \code{win_size} is equivalent to value of 
\code{lag} argument of \code{\link{diff}} +1. \code{keep_len} prevent the 
shortening of returned values by inserting \code{NA}(s) padding, which makes 
subsequent indexing covenient. However, presence of \code{NA}s could also 
interfere with subsequent processing such as smoothing (e.g. Savitzky-Golay 
filter in \link[signal]{sgolayfilt}). 
}