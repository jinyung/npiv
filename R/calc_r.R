calc_r <- function(dir, idlist, thre = c(0.5, 0.95), thre.len = 15, val.area) {
  
  # # frame sequence defined as specific number of frames leading to power stroke
  # frame.seq <- round(seq(from, to, length.out = num + 1))
  frame.seq <- idlist
  
  # read the files from sequences of frames and calculate U
  U <- list()
  for (i in seq_along(frame.seq)) {
    dat <- read_davis(file.path(dir, get_file(frame.seq[i])))
    # U[[i]] <- calc_u(dat, time_scal)
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


calc_s <- function(dir, idlist, ustar = 0.0005, val.area) {
  
  frame.seq <- idlist
  
  # read the files from sequences of frames and calculate U
  U <- list()
  for (i in seq_along(frame.seq)) {
    dat <- read_davis(file.path(dir, get_file(frame.seq[i])))
    U[[i]] <- calc_u(dat)
  }
  
  s <- r <-  NULL
  for (i in seq_along(frame.seq)) {
    s[i] <- area <- sum(U[[i]] >= ustar) * val.area
    r[i] <- equi_rad(area)
  }

  # output
  return(list(r = r, s = s))
}
  
