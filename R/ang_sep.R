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
  return(list(ang_sep = ang_sep,
              amplitude = amplitude,
              langle = langle))
}
