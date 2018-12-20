# --- normalize a vector into 0-1 range ---
# x: (num)
range01 <- function(x) {
  ((x - min(x)) / (max(x) - min(x)))
}

# --- centroid ---
cent <- function(landmark) {
  apply(landmark, 2, mean)
}

# --- xy coords format convert---
# a small covenient wrapper to turn list of x and y into a 2column xy matrix
xycoords2mat <- function(xycoord) {
  matrix(unlist(xycoord), length(xycoord$x), 2,
         dimnames = list(NULL, c('x', 'y')))
}

# --- diff() with NAs, same length ---
# a small covenient wrapper for diff()
dif <- function(x, lag = 2, differences = 1) {
  out <- rep(NA, length(x))
  idx <- c(floor((1+lag/2)):floor((length(x)-lag/2)))
  out[idx] <- diff(x, lag = lag, differences = differences)
  out
}

# --- create logarithmic sequence ---
# arguments similar to seq, see `?seq`
# https://stackoverflow.com/a/29963530
lseq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}

# --- print filename of davis output for i-th slice of snapshot ---
# slice_num: (int) i-th snapshot
get_file <- function(slice.num) {
  paste0(sprintf("B%05d", slice.num), ".txt")
}

# --- radius of circle of equivalent area ---
# area: (num)
equi_rad <- function(area) {
  sqrt(area / pi)
}

# --- plot vector field ---
# dat: (matrix) object read with read_davis
# width, image: (int) size of image, default 1280 x 1024
# type: (char) data values in either mm or pixel
# fac: (int) scaling factor for arrow length
# scale: (num) scale in mm/pix
plotvf <- function(dat, width = 1280, height = 1024,
                   type = c("mm", "pixel"), fac = 30, scale,
                   add = FALSE) {

  type <- match.arg(type)

  if (type == "mm") {
    lab <- "mm"
  } else if (type == "pixel") {
    # scale
    if (!missing(scale)) {
      width <- width * scale
      height <- height * scale
      dat <- dat * scale
      lab <- "mm"
    } else {
      lab <- "pixel"
    }
  }

  # size of arrows heads vary with speed
  arrowsize <- range01(calc_u(dat)) * 0.2

  # empty plot
  if (add != TRUE)
    plot(NULL, type = "n", xlim = c(0, width), ylim = c(height, 0),
         xlab = lab, ylab = lab, asp = 1)
  # arrows
  shape::Arrows(dat$x, dat$y, dat$x + dat$u * fac,
                dat$y + dat$v * fac, code = 2,
                arr.length = arrowsize)
}

# --- calculate speed from velocity vectors --- dat: (matrix) object read with
# read_davis(), davis output is already velocity no need to divide by time time:
# ()
calc_u <- function(dat) {
  sqrt(dat$u^2 + dat$v^2)
}
# calc_u <- function(dat, time = 1/2000) {
#   sqrt(dat$u^2 + dat$v^2) / time
# }

# angle between target vector and a vector extends from starting point &
# parallel to x axis
# x and y are two points of the target vector, with x = from, y = to
# x = c(x1, y1); y = c(x2, y2)
angle <- function(x, y, degree = TRUE) {
  p1 <- x
  p2 <- y
  p3 <- c(p1[1] + 50, p1[2])  # 50 is arbitrarily set
  a = c(p2[1] - p1[1], p2[2] - p1[2])  # vect 1
  b = c(p3[1] - p1[1], p3[2] - p1[2])  # vect 2
  result <- acos(sum(a * b) / (sqrt(sum(a * a)) * sqrt(sum(b * b))))
  if (degree)
    result <- result * 180 / pi  # convert from rad to deg
  if (p2[2] > p1[2])
    result <- 360 - result
  result
}

# --- simplified arrows plotting function ---
# v1 = c(x1, y1); v2 = c(x2, y2)
arrow <- function(v1, v2, ...) {
  x1 <- v1[1]
  x2 <- v1[2]
  y1 <- v2[1]
  y2 <- v2[2]
  arrows(x1, x2, y1, y2, ...)
}

# --- find the unit vector of a vector ---
# v = c(x, y, ...) vector
unitv <- function(v) {
  v / sqrt(sum(v^2))  # v divide by magnitude of v
}

# --- vector field format convert ---
# turn format returned by read_davis into fields package `surface` format, with
# y axis inverted
# out: to set the output component (u/v component of vf)
# invert: whether to invert y-axis. default is true as i set my vf origin at the
# bottom left corner in DaVis but surface format use topleft as origin
vf2surface <- function(vf, out = c('u', 'v'), invert = TRUE) {
  # arguments
  out <- match.arg(out)
  z <- eval(parse(text = paste0('vf$', out)))

  # extract
  x <- unique(vf$x)
  y <- unique(vf$y)
  z <- matrix(z, length(x), length(y))
  z <- z[, length(y):1]

  # correct change of y axis problem
  if (out == 'v') {
    if (invert)
      z <- -z
  }


  # return
  return(list(x = x, y = y, z = z))
}

# --- calculate standard error ---
se <- function(x) sd(x)/ sqrt(length(x))

# --- make an empty named list ---
# n (int) = length of list
# name (vector of string) = name of items in list
makelist <- function(n, name) {
  setNames(vector('list', n), name)
}

# --- floor/ceiling in decimal points---
# https://stackoverflow.com/a/39611375
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

# permutational t-test, two tailed
perm_ttest <- function(x, y, perm = 999, seed = 8888, plot = FALSE) {
  set.seed(seed)
  xn <- length(x)
  yn <- length(y)
  n <- xn + yn
  pool <- c(x, y)
  obs_diff <- mean(x) - mean(y)
  null_diff <- NULL
  for (i in 1:perm) {
    resampled <- sample(pool)
    nullx <- resampled[1:xn]
    nully <- resampled[(xn+1):n]
    null_diff[i] <- mean(nullx) - mean(nully)
  }
  pval = sum(abs(c(null_diff, obs_diff)) >= abs(obs_diff)) / (perm+1)
  if (plot == TRUE) {
    hist(null_diff, col = 'gray', main = '', freq = FALSE,
         xlab = 'Mean difference')
    abline(v = obs_diff, col = 'blue')
    legend('topleft', bg = 'white', pch = c(22, NA), lty = c(NA, 1),
           col = c('black', 'blue'), pt.bg = c('grey', NA),
           legend = c('Null distribution', 'Observed'),
           inset = 0.01)
    box()
  }
  return(pval)
}

# utility tool to find the peak fluid speed (at certain quantile) of a frame
# peak_id = frame to find the peak fluid speed
# (e.g. middle of power stroke with peak power)
# quan = which quantile, 95% u of that frame (95% excluding the min, where u
# below min is regarded as noise)
max_u <- function(dir, peak_id, quan = 0.99, min = 0.0005) {
  dat <- read_davis(file.path(dir, get_file(peak_id)))
  U <- calc_u(dat)
  quantile(U[U>min], quan)
}
