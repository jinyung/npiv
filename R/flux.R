flux <- function(body_tps_file, vf_dir, vect_scale, body_width, scale, 
                 plot = TRUE) {
  
  # read tps
  body_tps <- kt$readtps(body_tps_file)
  id <- as.integer(body_tps$id)
  body_land <- body_tps$coords
  
  # centroid to tail vector
  tail_land <- t(body_land[2, , ])
  centroid <- t(apply(body_land, 3, cent))
  cent_tail_vect <- tail_land - centroid
  
  # scale centroid-to-tail-vector
  cent_tail_vect_scaled <- cent_tail_vect * vect_scale
  tail_extend <- centroid + cent_tail_vect_scaled
  
  # perpendicular line to cent-tail vect
  slope_ct <- cent_tail_vect[, 2] / cent_tail_vect[, 1]  # cent-tail vector slope
  slope_ct_perp <- -1/slope_ct  # slope of perpendicular line
  intercept_ct_perp <- tail_extend[, 2] - (slope_ct_perp * tail_extend[, 1])
  
  ## solution1
  # # get lines from the fh, not used. body width preferred over fh width
  # # find the interception points from these three lines
  # left_fh_land <- t(body_land[3, , ])
  # right_fh_land <- t(body_land[1, , ])
  # intercept_left_fh <- left_fh_land[, 2] - (slope_ct * left_fh_land[, 1])
  # intercept_right_fh <- right_fh_land[, 2] - (slope_ct * right_fh_land[, 1])
  
  ## alternative solution
  # create arbitrarily extended vectors to both sides and turn it into unit vector
  # Note: did not consider the case where y = x (vertical line)
  ext <- 50  # arbitrarily set
  ext1 <- as.matrix(data.frame(x = ext, y = ext*slope_ct_perp))  # arbitraty vector1
  unit_ext1 <- t(apply(ext1, 1, unitv))  # unit vector of ext1
  flux_line_pt1 <- tail_extend + unit_ext1 * body_width/2
  # neg unit vector1 gives opposite direction vector, i.e. unit_ext2 = -unit_ext1
  flux_line_pt2 <- tail_extend - unit_ext1 * body_width/2  
  
  # sample equidistant points along the line for calculation of flux
  npoints <- 17  # arbitrarily set, have to be odd number > 3, produce n-1/2 cells 
  cell.idx <- seq(2, npoints - 1, 2)
  interp_ln <- list()
  for (i in 1:length(id)) {
    approx_result <- approx(x = c(flux_line_pt1[i, 1], flux_line_pt2[i, 1]), 
                            y = c(flux_line_pt1[i, 2], flux_line_pt2[i, 2]), 
                            n = npoints)
    interp_ln[[i]] <- lapply(approx_result, `[`, cell.idx)
  }
  
  # create arbitrary normal vector on interp_ln and scale to unit vector
  ext <- 50  # arbitrarily set
  normv <- as.matrix(data.frame(x = ext, y = ext*slope_ct))  # arbitraty normal vector
  unit_normv <- t(apply(normv, 1, unitv))  # unit normal vector
  
  # determine direction of unit normal vector by dot product with tail vector
  unv_direction <- NULL
  for (i in 1:dim(unit_normv)[1])
    unv_direction[i] <- sign(cent_tail_vect[i, ] %*% unit_normv[i, ])  
  
  # unit_normv should be in same direction as tail vector, so flux going into
  # labrum will be in negative direction (following convention of 'in'=-ve and
  # 'out'=+ve)
  unit_normv <- unit_normv * c(unv_direction)
  
  projected_vf <- interp_vf <- vf_u <- vf_v <- list()
  for (i in seq_along(id)) {
    
    # read the vector field and change the format for interpolation
    vf <- read_davis(file.path(vf_dir, get_file(id[i])))
    vf_u[[i]] <- vf2surface(vf, 'u')
    vf_v[[i]] <- vf2surface(vf, 'v')
    
    # interpolate u and v on the line
    interp_vf_tmp <- matrix(NA, length(cell.idx), 2)
    interp_vf_tmp[, 1] <- fields::interp.surface(vf_u[[i]], 
                                                 xycoords2mat(interp_ln[[i]]))
    interp_vf_tmp[, 2] <- fields::interp.surface(vf_v[[i]], 
                                                 xycoords2mat(interp_ln[[i]]))
    interp_vf[[i]] <- interp_vf_tmp
    
    # projection by dot product with the normal unit vector
    projected_vf[[i]] <- interp_vf[[i]] %*% unit_normv[i,]
  }
  
  flux <- sapply(projected_vf, sum) * body_width/length(cell.idx)
  
  if (!missing(scale))
    flux <- flux * scale * scale
  
  # optional plot check
  if (plot) {
    for (i in seq_along(id)) {
      plot(body_land[, , i], asp = 1, ylim = c(0, 1024), xlim = c(0, 1280), 
           xlab = 'x', ylab = 'y')
      segments(flux_line_pt1[i, 1], flux_line_pt1[i, 2], 
               flux_line_pt2[i, 1], flux_line_pt2[i, 2], 
               col = 2)
      for (j in 1:8)
        arrow(xycoords2mat(interp_ln[[i]])[j, ], 
              xycoords2mat(interp_ln[[i]])[j, ] + interp_vf[[i]][j, ]  * 10000, 
              length = 0.02, col = "gray")
      for (j in 1:8)
        arrow(xycoords2mat(interp_ln[[i]])[j, ], 
              xycoords2mat(interp_ln[[i]])[j, ] + unit_normv[i,] * 
                projected_vf[[i]][j, ]  * 10000, 
              length = 0.02, col = 1)
    }
  }
  
  return(flux)
}
