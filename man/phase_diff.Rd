\name{phase_diff}
\alias{phase_diff}
\title{Calculate phase of beating}
\description{
Given angular positions, calculate proportion of appendages moving in same directions.
}
\usage{
phase_diff(langle_smooth) 
}
\arguments{
  \item{langle_smooth}{(matrix) angle of appendages}
}

\value{
  \item{two_pairs_result}{pairwise proportion of appendages moving in same phase}
  \item{all_pairs_result}{proportion of all pairs moving in same phase}
}

\details{
This function works with returned value of \code{\link{ang_sep}} (the \code{langle} result). It is recommended to smooth the angle data before use of this function to remove digitization errors.
}
\seealso{
\code{\link{ang_sep}}
}