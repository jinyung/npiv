\name{read_davis}
\alias{read_davis}
\alias{read_davis_2}
\title{Read DaVis files}
\description{Read vector and scalar field outputs from DaVis software
}
\usage{
read_davis(file)

read_davis_2(dir, id, type = c("df", "list"), dim = c(1280, 1024), 
             win_size = 16)
}
\arguments{
  \item{file}{(char) path to file to be read}
  \item{dir}{(char) path to directory containing the file to be read}
  \item{id}{(int) index to frame to be read}
  \item{type}{(char) 'df' or 'list', see details}
  \item{dim}{(int) two values indicating width and height of field in pixel}
  \item{win_size}{size of each value represent, one value. if 16
    is given, it is assume each value covers area of 16x16 pixels}
}

\value{
Vector/ scalar field
}

\details{
\code{read_davis} is suitable to read velocity field (vector field) from DaVis.

\code{read_davis_2} is suitable to read vorticity field (scalar field) from 
DaVis. Argument \code{type} determines the resulted formats, which should be 
chosen based on compatility with functions from other popular packages for 
subsequent analysis: '\code{df}' is in the format used by 
\link[graphics]{image} and \link[fields]{image.plot} while '\code{list}' is in 
the format used by \link[fields]{quilt.plot}.
}