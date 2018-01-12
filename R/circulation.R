# function to calculate circulation with feature of splitting the vorticity
# field into left and right of body

circulation <- function(dir, body_tps_file, dim = c(1280, 1024), 
                        win_size = 16, plot = FALSE, side = c('r', 'l'), 
                        scale, threshold) {
  
  # argument processing/ get dimension data
  side <- match.arg(side)
  if (!missing(threshold))
    threshold <- abs(threshold)
  xdim <- dim[1]  # width of the field
  ydim <- dim[2]  # height of the field
  gridx <- seq(win_size/2, xdim - win_size/2, win_size)
  gridy <- seq(win_size/2, ydim - win_size/2, win_size)
  gridxn <- xdim / win_size  # number of cell along width
  gridyn <- ydim / win_size  # number of cell along height
  
  # read body tps file, body landmarks used to define left-right
  body_tps <- kt$readtps(body_tps_file)
  id <- as.integer(body_tps$id)
  body_land <- body_tps$coords
  
  # read the vorticity field files into a list
  dat_list = list()
  for (i in seq_along(id)) {
    dat_list[[i]] <- read_davis_2(dir, id[i], type = 'df')
  }
  names(dat_list) <- id
  
  # draw a mid body line, get the slope and intercept of the line
  tail_land <- t(body_land[2, , ])
  centroid <- t(apply(body_land, 3, cent))
  ref_vect <- tail_land - centroid
  slope <- ref_vect[, 2] / ref_vect[, 1]
  intercept <- centroid[, 2] - (slope * centroid[, 1])
  
  # initialize for loop
  left <- right <- list()
  
  for (i in seq_along(id)) {
    # calculate the xy coordinates on grid on mid line
    yy <- slope[i] * gridx + intercept[i]
    xx <- (gridy - intercept[i]) / slope[i]
    
    # determine left or right of mid line, first by calculate angle of mid-line
    ang <- angle(centroid[i, ], tail_land[i, ])  # calculate angle of tail 
    
    # angle determine orientation, orientation determine left side
    if (ang == 0) {
      comparex <- `|`  # give TRUE not matter ==/>/<
      comparey <- `<`
    } else if (ang == 180) {
      comparex <- `|`
      comparey <- `>`
    } else if (ang == 90) {
      comparex <- `<`
      comparey <- `|`
    } else if (ang == 270) {
      comparex <- `>`
      comparey <- `|`
    } else if (ang < 90) {
      comparex <- `<`
      comparey <- `<`
    } else if (ang < 180) {
      comparex <- `<`
      comparey <- `>`
    } else if (ang < 270) {
      comparex <- `>`
      comparey <- `>`
    } else if (ang < 360) {
      comparex <- `>`
      comparey <- `<`
    }
    xx_bool_i <- comparex(dat_list[[i]]$x, rep(xx, each = gridxn))
    yy_bool_i <- comparey(dat_list[[i]]$y, rep(yy, gridyn))
    left[[i]] <- xx_bool_i & yy_bool_i
    right[[i]] <- !left[[i]]
  }
  
  # calculate circulation
  cir_neg <- cir_pos <- NULL
  # which side?
  if (side == 'l')
    side_idx <- left
  else
    side_idx <- right
    
  for (i in seq_along(id)) {
    vor <- dat_list[[i]]$z[side_idx[[i]]]
    vor_pos <- vor[vor > 0]
    vor_neg <- vor[vor < 0]
    
    # threshold is for vorticity and use absolute number
    if (!missing(threshold)) {
      vor_pos <- vor_pos[vor_pos > threshold]
      vor_neg <- vor_neg[abs(vor_neg) > threshold]
    }
    
    cir_pos[i] <- sum(vor_pos * scal^2 * win_size, na.rm = TRUE)
    cir_neg[i] <- sum(vor_neg * scal^2 * win_size, na.rm = TRUE)
  }
  
  # visual check
  if (plot == TRUE) {
    # t_0 <- id[1]  # to calculate time based on id
    
    for (i in seq_along(id)) {
      # tiff(filename = sprintf("%03d.tif", i), width = 1400, height = 1152, 
      #      res = 192, compression = 'lzw')
      fields::quilt.plot(x = dat_list[[i]]$x, y = dat_list[[i]]$y, 
                         z = dat_list[[i]]$z, asp = 1,
                         nx = gridxn, ny = gridyn, add = FALSE, add.legend = TRUE)
      
      # labels
      points(dat_list[[1]]$x[left[[i]]], dat_list[[1]]$y[left[[i]]], 
             col = 'gray40', pch = 22)
      points(dat_list[[1]]$x[right[[i]]], dat_list[[1]]$y[right[[i]]], 
             col = 1, pch = 22)
      abline(a = intercept[i], b = slope[i], lwd = 2)
      # legend('bottomleft', legend = paste((id[i] - t_0)*0.5, 'ms'), bg = 'white', 
      #        inset = 0.05)
      # dev.off()
    }
  }
  
  return(list(cir_pos, cir_neg))
} 