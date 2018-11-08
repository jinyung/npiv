# use dot product to differentiate direction, with dot prooduct >0 in same
# direction, dot product <0 in different direction and dot product = 0 as
# perpendicular

# depends on whether need to keep same time sampling points (same length easier
# for indexing, but will insert NA) or need to smooth the result (Savitzky-Golay
# filter will give undesirable result with NAs), adjust keep_len parameter.
# Generally it is easier to remove NA than to add later

disp <- function(centroid, body_land, keep_len = TRUE, win_size) {
  # size of moving window
  if (missing(win_size))
    win_size <- 3
  pad <- (win_size - 1) / 2
  
  # calculate distance
  if (keep_len == TRUE)
    xy_diff <- apply(centroid, 2, dif, win_size - 1)  # no n-2
  else
    xy_diff <- apply(centroid, 2, diff, win_size - 1)  # n-2
  xy_square <- apply(xy_diff, 2, `^`, 2)
  xy_square_sum <- apply(xy_square, 1, sum)
  xy_dist <- sqrt(xy_square_sum)
  
  # calculate dot product sign to determine forward/reverse
  tail_land <- t(body_land[2, , ])
  ref_vect <- tail_land - centroid
  if (keep_len == FALSE)
    ref_vect <- ref_vect[(1+pad):(dim(ref_vect)[1] - pad), ]  # to n-2
  # ref_vect <- (ref_vect[-dim(ref_vect)[1], ] + ref_vect[-1, ]) / 2  # to n-2
  disp_sign <- NULL
  for (i in 1:dim(ref_vect)[1])
    disp_sign[i] <- sign(xy_diff[i, ] %*% ref_vect[i, ])  # sign(0) = 0
  
  # distance times sign = displacement
  xy_disp <- xy_dist * -disp_sign  # '-' cos rev to tail is forward & vice versa
  
  return(xy_disp)
}