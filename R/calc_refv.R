
calc_refv <- function(body_land_array, type = c('top', 'side'),
                      scal, t_stamp) {
  type <- match.arg(type)
  if (type == 'top') {
    centroid_mat <- t(apply(body_land_array, 3, cent))
    tail_land_mat <- t(body_land_array[2, , ])
  } else {
    centroid_mat <- t(apply(body_land_array[-3, , ], 3, cent))  # 3rd lm only for direction
    tail_land_mat <- t(body_land_array[1, , ])
  }
  
  cent_diff <- apply(centroid_mat, 2, dif)   # with 2 less values
  # NA-ed for 1st and last one
  # refv defined as velocity of centroid in opposite direction to cent-tail-vec
  cent_tail_vect <- tail_land_mat - centroid_mat
  ref_disp <- NULL
  for (i in 1:dim(cent_diff)[1])
    ref_disp[i] <- cent_diff[i, ] %*% -unitv(cent_tail_vect[i, ])
  refv <- ref_disp * scal / dif(t_stamp)
  return(refv)
}
