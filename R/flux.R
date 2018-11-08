# vf_dir: (path) directory containing the vector field output from Davis id:
# (int) give the index to frames to be read for calculations vect_scale:
# (numeric) where the center of the flux line should locate, aend point of body
# centroid-tail vector after scaling with this scaling factor body_width:
# (numeric) this control the length of the flux line, should be given in pixel
# to match body_land scal: scale of body_width, in the same unit of displacement
# of vector field, it is assumed that vector field is properly scaled body_land:
# should be given in pixel relative (boolean): if true use the body as frame of
# reference, i.e. velocity passing thru the flux line is relative to the
# swimming speed

# examples:
# vf is in m/s, scal has to be changed to m/pixel, end result will be in m2/s
# 
  # body_land = unlist(body_lm[[8]])
  # vf_dir = vf_list[8]
  # id = id[[8]]
  # vect_scale = vect_scale[8]
  # body_width = body_width[8]
  # scal = scal[8]
  # t_stamp = tlist

flux <- function(body_land, vf_dir, id, vect_scale, body_width, scal, 
                 plot = FALSE, relative = TRUE, t_stamp) {
  
  # # read tps
  # body_tps <- kt$readtps(body_tps_file)
  # id <- as.integer(body_tps$id)
  # body_land <- body_tps$coords
  
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
  
  # create arbitrarily extended vectors to both sides and turn it into unit vector
  # Note: did not consider the case where y = x (vertical line)
  ext <- 50  # arbitrarily set
  ext1 <- as.matrix(data.frame(x = ext, y = ext * slope_ct_perp))  # arbitraty vector1
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
  
  unit_normv <- t(apply(normv, 1, unitv))  # unit normal vector, with magnitude 
                                           # of 1, doesn't matter px or mm = 1
  
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
    
    # projection by dot prod with the normal unit vector, give magnitude w sign
    projected_vf[[i]] <- interp_vf[[i]] %*% unit_normv[i, ]  # this is in the 
                                                             # unit of the vf 
  }
  
  # change the frame of reference
  # with 2 less value, NA-ed for first and last one
  
  cent_diff <- apply(centroid, 2, dif)
  ref_v <- NULL
  for (i in 1:dim(cent_diff)[1])
    ref_v[i] <- cent_diff[i, ] %*% -unitv(cent_tail_vect[i, ])  # opposite direction to cent-tail
  
  if (relative == TRUE) {
    # eudist <- function(x) sqrt(x[1]^2 + x[2]^2)
    # displ <- apply(tail_diff, 1, eudist) * scal
    # ref_speed <- displ / dif(t_stamp)
    
    # reference velocity (= nauplius body speed)

    # relative velocity 
    projected_vf <- mapply(`+`, projected_vf, ref_v * scal / dif(t_stamp), 
                           SIMPLIFY = FALSE)
  }
  
  flux <- sapply(projected_vf, sum) * (body_width * scal)/length(cell.idx)
  
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
  # return(list(projected_vf = projected_vf, flux=flux, ref_v = ref_v * scal / dif(t_stamp)))
}
