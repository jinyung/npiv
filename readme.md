# npiv

This is a companion package to the manuscript:

> Wong, J. Y., Chan, B. K. K., Chan, K. Y. K. Swimming kinematics of naupliar 
larvae with contrasting feeding modes.

Functions included in the package were developed for analyses of data generated
from [`TPSDig2`](http://life.bio.sunysb.edu/morph/soft-dataacq.html) and
[`DaVis`](https://www.lavision.de/en/products/davis-software/) softwares for
kinematics and hydrodynamics analyses, respectively. `TPSDig2` is a landmark
registration tool and was used for digitization of body and appendages
positions. `DaVis` is a software used for particle image velocimetry, and
outputs velocity and vorticity fields. There is no utility tools for direct data
import from other softwares but all underlying calculations use `R`'s `matrix`
or `array` classes of data.

## Notes on using python module from the package

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

## Reproducing the results
Data are not bundled with the package due to the large file size and are stored
in [OSF repository](https://osf.io/r9abn/). Landmarks data for kinematics analysis, data exported from `Davis` program (velocity field and vorticity field data) can be accessed and directly downloaded inside `R`*:

```r
# download path and local temp saving path (or change to path you want to save)
tps_url <- 'https://files.osf.io/v1/resources/r9abn/providers/googledrive/kinematics/?zip='
vf_url <- 'https://files.osf.io/v1/resources/r9abn/providers/googledrive/vf_pxscale/?zip='
vort_url <- 'https://files.osf.io/v1/resources/r9abn/providers/googledrive/Vorticity/?zip='
tps_temp <- tempfile(fileext = ".zip")
vf_temp <- tempfile(fileext = ".zip")
vort_temp <- tempfile(fileext = ".zip")

# download and unzip the downloaded data
tps_zip <- download.file(tps_url, tps_temp)
vf_zip <- download.file(vf_url, vf_temp)
vort_zip <- download.file(vort_url, vort_temp)
tps_dir <- unzip(tps_temp, exdir = tempdir())
vf_dir <- unzip(vf_temp, exdir = tempdir())
vort_dir <- unzip(vort_temp, exdir = tempdir())
```

`vf_dir` and `vort_dir` store list of paths for velocity and vorticity fields 
files downloaded into temporary folders. It takes time for `osf` to zip the 
files, so expect waiting time before download commence. 

\* NOTE: method not tested, will only work after the project goes public, now the API will block 
access.