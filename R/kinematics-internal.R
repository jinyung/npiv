# taken from otolith project

cent <- function(landmark) {
  apply(landmark, 2, mean)
}

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
