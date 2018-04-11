# --- calculate r as a function of U* ---
# dir: (path) directory containing the vector field output file
# idlist: (int vector) index to frames for calculations
# thre: (num) threshold values to adjust U*, two values required
# thre.len: (int) number of U* values to use
# val.area: (num) area each velocity vector value covers in mm^2, ..
# .. by default 16pix*16pix*(scale in mm/pix)^2
# time_scal: (num) time interval between frames for calculation of U 
# Note: the scale of the displacement used for calculation of U follows that..
# .. of vector field

calc_r <- function(dir, idlist, thre = c(0.5, 0.95),thre.len = 15, val.area, 
                   time_scal = 1/2000) {
  
  # # frame sequence defined as specific number of frames leading to power stroke
  # frame.seq <- round(seq(from, to, length.out = num + 1))
  frame.seq <- idlist
  
  # read the files from sequences of frames and calculate U
  U <- list()
  for (i in seq_along(frame.seq)) {
    dat <- read_davis(file.path(dir, get_file(frame.seq[i])))
    U[[i]] <- calc_u(dat, time_scal)
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
  
  # # scale
  # if (!missing(scal)) {
  #   U.star <- U.star * scal 
  #   S <- S * scal * scal
  #   r <- r * scal
  # }
  
  # output
  return(list(r = r, S = S, u = U.star))
}
