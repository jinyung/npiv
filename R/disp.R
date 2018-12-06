disp <- function(body_land, keep_len = TRUE, win_size) {
  # size of moving window
  if (missing(win_size))
    win_size <- 3
  pad <- (win_size - 1) / 2
  
  # centroid
  centroid <- t(apply(body_land, 3, cent))
  
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
  disp_sign <- NULL
  for (i in 1:dim(ref_vect)[1])
    disp_sign[i] <- sign(xy_diff[i, ] %*% ref_vect[i, ])  # sign(0) = 0
  
  # distance times sign = displacement
  xy_disp <- xy_dist * -disp_sign  # '-' cos rev to tail is forward & vice versa
  
  return(xy_disp)
}