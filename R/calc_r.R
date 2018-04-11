# --- calculate r as a function of U* ---
# dir: (path) directory containing the velocity output file
# from: (int) frame at end of power stroke
# to: (int) frame at start of recovery stroke
# num: (int) number of snapshots leading to `from` at equal intervals
# thre: (num) threshold values to adjust U*, two values required
# thre.len: (int) number of U* values to use
# val.area: (num) area each velocity vector value covers in mm^2, ..
# .. by default 16pix*16pix*scale in mm/pix

calc_r <- function(dir, from, to, num = 5, thre = c(0.5, 0.95),
                   thre.len = 15, val.area) {

  # frame sequence defined as specific number of frames leading to power stroke
  frame.seq <- round(seq(from, to, length.out = num + 1))
  U <- list()

  # read the files from sequences of frames and calculate U
  for (i in seq_along(frame.seq)) {
    dat <- read_davis(file.path(dir, get_file(frame.seq[i])))
    U[[i]] <- calc_u(dat)
  }

  # set U* based on distribution of U
  # I define boundary of U* as min and max of designated quantile values of U of
  # ..all frames, default are 50% and 95%, same set of U* for all frames.
  U.star <- lseq(min(sapply(U, quantile, thre[1])),
                 max(sapply(U, quantile, thre[2])),
                 length.out = thre.len)

  # get the area and and r of equivalent circle at diff U*
  # initiate first
  S <- r <- matrix(NA, nrow = length(U.star), ncol = length(frame.seq),
                   dimnames = list(U.star = U.star, frame = frame.seq))
  for (i in seq_along(frame.seq)) {
    for (j in seq_along(U.star)) {
      S[j, i] <- area <- sum(U[[i]] >= U.star[j]) * val.area
      r[j, i] <- equi_rad(area)
    }
  }

  # output
  return(list(r = r, S = S, u = U.star))
}
