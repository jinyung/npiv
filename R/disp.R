disp <- function(centroid, body_land) {
  # calculate displacement
  xy_diff <- apply(centroid, 2, diff)  # n-1
  xy_square <- apply(xy_diff, 2, `^`, 2)
  xy_square_sum <- apply(xy_square, 1, sum)
  xy_disp <- sqrt(xy_square_sum)
  
  # calculate angle to determine forward/reverse
  tail_land <- t(body_land[2, , ])
  ref_vect <- tail_land - centroid
  ref_vect <- (ref_vect[-dim(ref_vect)[1], ] + ref_vect[-1, ]) / 2  # to n-1
  angle <- NULL
  for (i in 1:dim(ref_vect)[1]) 
    angle[i] <- kt$vangle(xy_diff[i, ], ref_vect[i, ])
  disp_sign <- ifelse(angle > 90, 1, -1)  
  
  # disp times sign
  xy_disp <- xy_disp * disp_sign
  
  return(xy_disp)
}