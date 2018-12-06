get_fluxln <- function(body_land, vect_scal = 1, ln_width,
                       ln_width_scal = 1, type = c('top', 'side')) {
  # function things
  type <- match.arg(type)

  # centroid, tail, and centroid-to-tail-vector
  # for side, 3rd lm only for direction so not used
  if (type == 'top') {
    centroid <- cent(body_land)
    tail_land <- body_land[2, ]
  }
  else {
    centroid <- cent(body_land[-3, ])
    tail_land <- body_land[1, ]
  }
  cent_tail_vect <- tail_land - centroid

  # scale centroid-to-tail-vector
  cent_tail_vect_scaled <- cent_tail_vect * vect_scal

  # center point of flux line (top)
  # starting point of flux line (side)
  tail_extend <- centroid + cent_tail_vect_scaled

  # unit normal vector to cent_tail_vect
  unormv1 <- unitv(c(cent_tail_vect[2], -cent_tail_vect[1]))  # clockwise
  unormv2 <- unitv(c(-cent_tail_vect[2], cent_tail_vect[1]))  # anti-clockwise

  if (type == 'side') {
    ## determine the unit normal vector in ventral side
    # centroid-to-3rd-landmark vector
    cent_3rd_vect <- body_land[3, ] - centroid
    # use dot product to determine
    if (unormv1 %*% cent_3rd_vect > 0)
      unormv_vent <- unormv1
    else
      unormv_vent <- unormv2
  }

  # scale unormv and put it on the center point / starting point
  if (type == 'top') {
    # resulting flux line is from left to right of nauplius body
    pt1 <- tail_extend + unormv1 * ln_width * ln_width_scal / 2
    pt2 <- tail_extend + unormv2 * ln_width * ln_width_scal / 2
  } else {
    pt1 <- tail_extend
    pt2 <- tail_extend + unormv_vent * ln_width * ln_width_scal
  }

  # return
  return(xy.coords(x = c(pt1[1], pt2[1]), y = c(pt1[2], pt2[2])))
}
