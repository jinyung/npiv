ang_sep <- function(body_tps_file, limb_tps_file, smooth = as.integer(5)) {
  # read
  limb_tps <- kt$readtps(limb_tps_file)
  body_tps <- kt$readtps(body_tps_file)
  limb_land <- limb_tps$coords
  body_land <- body_tps$coords
  
  # limb angle
  langle <- kt$tlangle(body_land, limb_land)
  langle_smooth <- apply(langle, 2, kt$smooth, win_size = smooth)
  
  # define power and recovery stroke by ant2
  t_len <- length(body_tps$id)
  langle_power <- langle_smooth[1:(t_len/2), ]
  langle_recor <- langle_smooth[(t_len - t_len/2 + 1):t_len, ]
  
  # define mid of power stroke and mid of recovery stroke, by ant2
  power_mid_idx <- which.max(langle_power[, 2]) + round((which.min(langle_smooth[, 2]) - which.max(langle_power[, 2]))/2)
  recor_mid_idx <- which.min(langle_smooth[, 2]) + round((t_len/2 + which.max(langle_recor[, 2]) - which.min(langle_smooth[, 2]))/2)
  
  # calculate angular separation
  ang_sep = kt$tsangle(body_land, limb_land)
  ang_sep_smooth = apply(ang_sep, 2, kt$smooth, smooth)
  
  # result
  power_sep <- ang_sep_smooth[power_mid_idx, ]
  recor_sep <- ang_sep_smooth[recor_mid_idx, ]
  edge <- (smooth - 1) / 2
  ang_sep_smooth <- ang_sep_smooth[-c(1:edge, (t_len - edge + 1):t_len), ]
  langle_smooth <- langle_smooth[-c(1:edge, (t_len - edge + 1):t_len), ]
  max_sep <- apply(langle_smooth, 2, max) - apply(langle_smooth, 2, min)
  return(list(power_sep = power_sep, recor_sep = recor_sep, ang_sep = ang_sep, 
              ang_sep_smooth = ang_sep_smooth, max_sep = max_sep, 
              langle_smooth = langle_smooth, langle = langle))
}