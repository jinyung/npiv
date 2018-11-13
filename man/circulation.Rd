\name{circulation}
\alias{circulation}
\title{Calculate circulation}
\description{
Calculate circulation from vorticity field.
}
\usage{
circulation(dir, body_land, id, dim = c(1280, 1024), 
            win_size = 16, plot = FALSE, side = c('r', 'l'), 
            scal, threshold)
}
\arguments{
  \item{dir}{(char) path to directory containing the vorticity fields files}
  \item{body_land}{ (3d array) landmarks for separating left and right}
  \item{id}{(int) frames to read}
  \item{dim}{(int) size of the photos in pixels, two values consists of width and 
    height}
  \item{win_size}{(int) size of each vorticity value represent, one value. if 16
    is given, it is assume each value covers area of 16x16 pixels}
  \item{plot}{(logical) if \code{TRUE} plots of vorticity field with separating 
    mid line will be produced}
  \item{side}{(char) choices of '\code{r}' = right side or '\code{l}' = left side 
    for which the side of vorticity field will be used for calculation}
  \item{scal}{(num) scale, default in mm/pix}
  \item{threshold}{(num) vorticity values under this value are not used in 
    calculation and treated as noise}
}
\value{
  \item{cir_pos}{circulation in anti-clockwise direction}
  \item{cir_neg}{circulation in clockwise direction}
}
\details{
Vorticity field is splitted into left/right by mid line computed with the 
provided landmarks on nauplii carapace (\code{body_land}), and calculation is 
proceeded with values from only one side, specified by \code{side}.

Returned values of circulation are separated into rotation in clockwise/ anti-
clockwise directions. Circulation calculated from vorticity values <0 is the 
rotation in clockwise direction (right hand rule) and vice versa. 
}