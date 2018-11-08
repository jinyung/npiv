# assume id of tps file as number and the difference between id represent time difference
# mid of power and recovery stroke defined as the max speed of power stroke and recovery stroke
# but here assume the first half consist of power stroke and second half with recovery stroke
# and take the steepest change in angular speed
#
#
# return
# max_sep = amplitude
# ang_sep <- function(body_tps_file, limb_tps_file) {
#   # read
#   limb_tps <- kt$readtps(limb_tps_file)
#   body_tps <- kt$readtps(body_tps_file)
#   limb_land <- limb_tps$coords
#   body_land <- body_tps$coords
#
#   # limb angle
#   langle <- kt$tlangle(body_land, limb_land)
#   rownames(langle) <- body_tps$id
#   colnames(langle) <- c('ant1', 'ant2', 'mand')
#   # langle_smooth <- apply(langle, 2, kt$smooth, win_size = smooth)
#   langle_smooth <- apply(langle, 2, signal::sgolayfilt, 3, 19)
#
#   # define power and recovery stroke by ant2
#   # t_len <- length(body_tps$id)
#   # langle_power <- langle_smooth[1:(t_len/2), ]
#   # langle_recor <- langle_smooth[(t_len - t_len/2 + 1):t_len, ]
#
#   # define mid of power stroke and mid of recovery stroke, by ant2
#   # tlist <- as.integer(body_tps$id)
#   # d_t <- diff(tlist, 2)
#   # langle_diff_ant2 <- diff(langle_smooth[, 'ant2'], 2)
#   # langle_dt_ant2 <- langle_diff_ant2 / d_t
#   # power_mid_idx <- which.min(langle_dt_ant2[1:(t_len/2)]) + 1  # idx +1 because diff shorten length
#   # recor_mid_idx <- c((t_len/2 + 1):t_len)[which.max(langle_dt_ant2[(t_len/2 + 1):t_len]) + 1]
#
#   # calculate angular separation
#   ang_sep = kt$tsangle(body_land, limb_land)
#   # ang_sep_smooth = apply(ang_sep, 2, kt$smooth)
#   ang_sep_smooth = apply(ang_sep, 2, signal::sgolayfilt, 3, 19)
#
#
#   # result
#   # power_sep <- ang_sep_smooth[power_mid_idx, ]
#   # recor_sep <- ang_sep_smooth[recor_mid_idx, ]
#   # edge <- (smooth - 1) / 2
#   # ang_sep_smooth <- ang_sep_smooth[-c(1:edge, (t_len - edge + 1):t_len), ]
#   langle_smooth <- langle_smooth[-c(1:edge, (t_len - edge + 1):t_len), ]
#   max_sep <- apply(langle_smooth, 2, max) - apply(langle_smooth, 2, min)
#   return(list(power_sep = power_sep, recor_sep = recor_sep, ang_sep = ang_sep,
#               ang_sep_smooth = ang_sep_smooth, max_sep = max_sep,
#               langle_smooth = langle_smooth, langle = langle,
#               power_mid_idx = power_mid_idx, recor_mid_idx = recor_mid_idx))
# }

ang_sep <- function(body_land, limb_land) {

  # calculate limb angle
  langle <- kt$tlangle(body_land, limb_land)
  colnames(langle) <- limb <- c('ant1', 'ant2', 'mand')
  # langle_smooth <- apply(langle, 2, signal::sgolayfilt, 3, 19)

  # calculate angular separation
  ang_sep = kt$tsangle(body_land, limb_land)
  colnames(ang_sep) <- apply(combn(limb, 2), 2, paste, collapse = '-')
  # ang_sep_smooth = apply(ang_sep, 2, signal::sgolayfilt, 3, 19)

  # result
  amplitude <- apply(langle, 2, max) - apply(langle, 2, min)
  # return(list(ang_sep = ang_sep, ang_sep_smooth = ang_sep_smooth,
  #             amplitude = amplitude, langle_smooth = langle_smooth,
  #             langle = langle))
  return(list(ang_sep = ang_sep,
              amplitude = amplitude,
              langle = langle))
}
