# npiv

This is a companion package to the manuscript:

> Wong, J. Y., Chan, B. K. K., Chan, K. Y. K. (2020). Evolution of feeding shapes 
swimming kinematics of barnacle naupliar larvae: a comparison between trophic modes. 
*Integrative Organismal Biology*, obaa011. 
[doi:10.1093/iob/obaa011](https://doi.org/10.1093/iob/obaa011)

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

These data can also be downloaded directly inside `R`:

```r
# download path 
tps_url <- 'https://files.osf.io/v1/resources/r9abn/providers/googledrive/kinematics/?zip='
vf_url <- 'https://files.osf.io/v1/resources/r9abn/providers/googledrive/vf_pxscale/?zip='
vort_url <- 'https://files.osf.io/v1/resources/r9abn/providers/googledrive/Vorticity/?zip='

# create local temporary saving paths 
# (or change to any path you want to save, temporary files are.. temporary)
tps_temp <- tempfile(fileext = ".zip")
vf_temp <- tempfile(fileext = ".zip")
vort_temp <- tempfile(fileext = ".zip")

# download and unzip the downloaded data
# it takes time for osf to zip the files, so expect waiting time before download 
# commence (and download for vector/vorticity fields may fail sometimes due to 
# long response time resulted from large file sizes). 
tps_zip <- download.file(tps_url, tps_temp, mode = 'wb')
vf_zip <- download.file(vf_url, vf_temp, mode = 'wb')
vort_zip <- download.file(vort_url, vort_temp, mode = 'wb')
tps_dir <- unzip(tps_temp, exdir = tempdir())
vf_dir <- unzip(vf_temp, exdir = tempdir())
vort_dir <- unzip(vort_temp, exdir = tempdir())

# now you can read the downloaded files, for example: 'tps_dir' is now a list
# of file paths of all the .TPS files downloaded, to read the first .TPS file
# using the 'readtps' function:
npiv::kt$readtps(tps_dir[1])
```
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

\* Please refrain from naming `R` object as `kt`.

\*\* May not work inside RStudio. Terminal (Ubuntu)/ cmd (Windows) are ok.
