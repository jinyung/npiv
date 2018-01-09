read_davis_2 <- function(dir, slice.num, dim = c(1280, 1024), win_size = 16) {
  # read
  vor <- as.matrix(read.table(file.path(dir, get_file(slice.num)), sep = "\t", 
                              skip = 1, dec = ","))
  
  # set x and y
  xdim <- dim[1]
  ydim <- dim[2]
  y <- as.integer(rep(seq(win_size/2, ydim - win_size/2, win_size), 
                      xdim/win_size))
  x <- as.integer(rep(seq(win_size/2, xdim - win_size/2, win_size), 
                      each = ydim/win_size))
  
  return(data.frame(x = x, y = y, z = c(vor)))
}

read_davis <- function(file) {
  dat <- read.table(file, sep = "\t", skip = 1, dec = ",")
  colnames(dat) <- c("x", "y", "u", "v")
  return(dat)
}

