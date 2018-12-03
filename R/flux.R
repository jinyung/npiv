# ... arguments to be passed to FUN, default to get_fluxln
# in case of get_fluxlinm include body_land, vect_scal, ln_width, scal,
# ln_width_scal
calc_flux <- function(FUN = get_fluxln, vf_path,
                      side = c('right', 'left'),
                      refv = 0, width, ...) {
  side <- match.arg(side)
  fluxlnpts <- FUN(...)
  projected_vf <- proj_vf(fluxlnpts, vf_path = vf_path, side = side)
  
  # return sum of vectors * width they represent, after adjustment by relative
  # velocity
  if (as.character(substitute(FUN)) == 'get_fluxln')
    width <- ln_width*scal*ln_width_scal
  return((mean(projected_vf)+refv) * width)
}

tcalc_flux <-function(body_land_array, vf_dir, id, vect_scal,
                      type = c('top', 'side'), relative = TRUE,
                      ln_width, ln_width_scal,
                      side = c('right', 'left'),
                      t_stamp, scal) {
  
  # function things
  if(dim(body_land_array)[3] != length(id))
    stop('body_land_array and id needs to contain same number of frames')
  type = match.arg(type)
  side = match.arg(side)
  
  # relative
  if (relative) {
    if (missing(t_stamp) | missing (scal)) {
      stop('t_stamp/ scal need to be provided if relative = TRUE')
    }
    refv <- calc_refv(body_land_array, type = type, scal = scal,
                      t_stamp = t_stamp)
  } else {
    refv = 0
  }
  
  # calc flux
  result <- NULL
  for (i in 1:length(id)) {
    result[i] <- calc_flux(FUN = get_fluxln,
                           vf_path = file.path(vf_dir, get_file(id[i])),
                           refv = refv[i],
                           vect_scal = vect_scal,
                           ln_width, ln_width_scal = ln_width_scal,
                           side = side,
                           body_land = body_land_array[, , i])
  }
  
  return(result)
}

# given two points, calculate projected_vf
# nseg = 17 is arbitrarily set
# nseg = how many equal segments are the flux line cut into
# xycoord use xy.coords format, with $x and $y list, with x1 and y1 as the
# starting point and x2 and y2 as the end point of flux line vector (which
# determine the direction of the flux)
# default is right, because that is the correct direction for top flux estimation
# when used with get_fluxln
#xycoord <- get_fluxln(body_land, ln_width = 106)

proj_vf <- function(xycoord, vf_path, nseg = 8, side = c('right', 'left')) {

  side <- match.arg(side)

  ## define flux line
  # let's say 8 segments
  # indices of center points of each segments on a line of 8 seg from 17 points
  cell_idx <- seq(2, nseg * 2 + 1, 2)
  # get the xycoords of these 17 points
  approx_out <- approx(xycoord, n = nseg*2+1)  # in/output in xy.coords format
  # from the indices, get the center points of each segment
  flux_ln <- lapply(approx_out, `[`, cell_idx)

  ## interpolate vector field on flux line
  # get the vector field
  vf <- read_davis(vf_path)
  vf_u <- vf2surface(vf, 'u')  # list with x, y, z
  vf_v <- vf2surface(vf, 'v')

  # interpolate u and v on the flux line
  interp_vf<- matrix(NA, length(cell_idx), 2)
  interp_vf[, 1] <- fields::interp.surface(vf_u, xycoords2mat(flux_ln))
  interp_vf[, 2] <- fields::interp.surface(vf_v, xycoords2mat(flux_ln))

  ## projection
  # let's make flux line a vector
  flux_ln_mat <- xycoords2mat(flux_ln)  # just change format
  flux_ln_vec <- unitv(flux_ln_mat[2, ] - flux_ln_mat[1, ])

  # unit normal vector, have two in either direction
  # find normal vector by swapping x and y, and minus to either one
  unormv1 <- c(flux_ln_vec[2], -flux_ln_vec[1])  # clockwise (right)
  unormv2 <- c(-flux_ln_vec[2], flux_ln_vec[1])  # anti-clockwise (left)

  # projection by dot prod with the normal unit vector
  if (side == 'right')
    projected_vf <- interp_vf %*% unormv1
  else
    projected_vf <- interp_vf %*% unormv2

  return(projected_vf)
}
