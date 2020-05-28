# npiv

This is a companion package to the manuscript:

> Wong, J. Y., Chan, B. K. K., Chan, K. Y. K. (2020). Evolution of feeding shapes swimming kinematics of barnacle naupliar larvae: a comparison between trophic modes. *Integrative Organismal Biology*, obaa011. [doi:10.1093/iob/obaa011](https://doi.org/10.1093/iob/obaa011)

Functions included in the package were developed for analyses of data generated
from [`TPSDig2`](http://life.bio.sunysb.edu/morph/soft-dataacq.html) and
[`DaVis`](https://www.lavision.de/en/products/davis-software/) softwares for
kinematics and hydrodynamics analyses, respectively. `TPSDig2` is a landmark
registration tool and was used for digitization of body and appendages
positions. `DaVis` is a software used for particle image velocimetry, and
outputs velocity and vorticity fields. There is no utility tools for direct data
import from other softwares but all underlying calculations use `R`'s `matrix`
or `array` classes of data.

Most of the functions were written specifically for the manuscript. However, some functions such as 
[`calc_flux`](R/calc_flux.R) and [`calc_r`](R/calc_r.R), are suitable for 
general use with any velocity field data.

## Reproducing the results
Data are not bundled with the package due to the large file size and are stored
in [OSF repository](https://osf.io/r9abn/). 

___

**Notes on using python module from the package**
Functions for kinematics analysis was written in `python` and is located in
[`kinematics.py` module](inst/python/kinematics.py). The module is imported into
`R` when the package is loaded via
[`reticulate`](https://rstudio.github.io/reticulate/index.html) package as `kt`
object*, and the functions can be called with `kt$<function name>`. Help file
can be accessed with `reticulate::py_help` function, e.g.
`reticulate::py_help(kt)` will show the help messages for all functions in the
module**.

NOTES: 

\* Please refrain from naming `R` object as `kt`.

\*\* May not work inside RStudio. Terminal (Ubuntu)/ cmd (Windows) are ok.
