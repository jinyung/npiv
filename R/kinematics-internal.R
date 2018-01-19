# centroid
cent <- function(landmark) {
  apply(landmark, 2, mean)
}

# angle between target vector and a vector extends from starting point &
# parallel to x axis
# x and y are two points of the target vector
# x = c(x1, y1); y = c(x2, y2)
angle <- function(x, y) {
  p1 <- x
  p2 <- y
  p3 <- c(p1[1] + abs(p2[1] - p1[1]), p1[2])
  a = c(p1[1] - p2[1], p1[2] - p2[2])
  b = c(p1[1] - p3[1], p1[2] - p3[2])
  result <- acos(sum(a * b) / (sqrt(sum(a * a)) * sqrt(sum(b * b))))
  result <- result * 180/ pi
  if (p2[2] > p1[2]) 
    result <- 360 - result
  result
}

# simplified arrows plotting function
# v1 = c(x1, y1); v2 = c(x2, y2)
arrow <- function(v1, v2, ...) {
  x1 <- v1[1]
  x2 <- v1[2]
  y1 <- v2[1]
  y2 <- v2[2]
  arrows(x1, x2, y1, y2, ...)
}

# find the unit vector of a vector
# v = c(x, y, ...) vector
unitv <- function(v) {
  v / sqrt(sum(v^2))  # v divide by magnitude of v
}

# # 
# ild <- function(a, b) {
#   sqrt(sum((b-a)^2))
# }

# turn format returned by read_davis into fields package `surface`` format, with
# y axis inverted
vf2surface <- function(vf, out = c('u', 'v')) {
  # arguments
  out <- match.arg(out)
  z <- eval(parse(text = paste0('vf$', out)))
  
  # extract
  x <- unique(vf$x)
  y <- unique(vf$y)
  z <- matrix(z, length(x), length(y))
  z <- z[, length(y):1]
  
  # correct change of y axis problem
  if (out == 'v')
    z <- -z
  
  # return
  return(list(x = x, y = y, z = z))
}

xycoords2mat <- function(xycoord) {
  matrix(unlist(xycoord), length(xycoord$x), 2, 
         dimnames = list(NULL, c('x', 'y')))
}