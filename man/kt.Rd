\name{kt}
\alias{kt}
\title{kinematics analyses}
\description{Perform kinematics analyses on swimming appendages of nauplius}
\details{
Functions for kinematics analyses were written in \code{python} and is located in
\code{kinematics.py module}. The module is imported into
\code{R} when the package is loaded via
\code{reticulate} package as \code{kt}
object*, and the functions can be called with \code{kt$<function name>}. Help file
can be accessed with \code{reticulate::py_help} function, e.g.
\code{reticulate::py_help(kt)} will show the help messages for all functions in the
module**.
}
\note{
* Please refrain from naming \code{R} object as \code{kt}.

** May not work inside RStudio. Terminal (Ubuntu)/ cmd (Windows) are ok.
}
