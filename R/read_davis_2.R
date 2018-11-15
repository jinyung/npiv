read_davis_2 <- function(dir, id, type = c("df", "list"),
                         dim = c(1280, 1024), win_size = 16) {
  
  # arg proc
  type <- match.arg(type)
  
  # read
  vor <- as.matrix(read.table(file.path(dir, get_file(id)), sep = "\t", 
                              skip = 1, dec = ","))
  
  # reverse the order of rows so that y axis range from 0-1000 not 1000-0
  vor <- apply(vor, 2, rev)
  vor[vor == 0] <- NA  # turn mask 0-value into NA
  vor <- t(vor)
  
  # set x and y
  xdim <- dim[1]
  ydim <- dim[2]
  x <- seq(win_size/2, xdim - win_size/2, win_size)
  y <- seq(win_size/2, ydim - win_size/2, win_size)
  xx <- as.integer(rep(x, ydim/win_size))  # rep for x takes ydim, vice versa
  yy <- as.integer(rep(y, each = xdim/win_size))   
  
  # return
  if (type == 'list')
    result <- list(x = x, y = y, z = vor)
  else
    result <- data.frame(x = xx, y = yy, z = c(vor))
    
  return(result)
}


read_davis <- function(file) {
  dat <- read.table(file, sep = "\t", skip = 1, dec = ",")
  colnames(dat) <- c("x", "y", "u", "v")
  return(dat)
}

