phase_diff <- function (langle_smooth) {
  limb_comb <- list(ant1_ant2 = c(1, 2), ant1_mand = c(1, 3), 
                    ant2_mand = c(2, 3))
  # get first derivative
  d_langle_smooth <- apply(langle_smooth, 2, dif)
  # get the sign of the first derivative to determine +/- slopes
  sign_langle_smooth <- apply(d_langle_smooth, 2, sign)
  two_pairs_result <- sapply(limb_comb, function(x) {
                                          sign1 <- sign_langle_smooth[, x[1]];
                                          sign2 <- sign_langle_smooth[, x[2]];
                                          inphase <- sum(sign1 == sign2, na.rm = TRUE);
                                          sum(inphase) / length(sign1)})
  all_pairs_result <- sum(sign_langle_smooth[, 1] == sign_langle_smooth[, 2] & 
                          sign_langle_smooth[, 2] == sign_langle_smooth[, 3], na.rm = TRUE) / 
                        length(sign_langle_smooth[, 1])
  return(c(two_pairs_result, all_pairs_result))
}

# phase_diff2 <- function (langle_smooth) {
#   limb_comb <- list(ant1_ant2 = c(1, 2), ant1_mand = c(1, 3), 
#                     ant2_mand = c(2, 3))
#   # get first derivative
#   d_langle_smooth <- apply(langle_smooth, 2, dif)
#   # get the sign of the first derivative to determine +/- slopes
#   sign_langle_smooth <- apply(d_langle_smooth, 2, sign)
#   two_pairs_result <- sapply(limb_comb, function(x) {
#     sign1 <- sign_langle_smooth[, x[1]];
#     sign2 <- sign_langle_smooth[, x[2]];
#     inphase <- sum(sign1 == sign2);
#     sum(inphase) / length(sign1)})
#   all_pairs_result <- sum(sign_langle_smooth[, 1] == sign_langle_smooth[, 2] & 
#                             sign_langle_smooth[, 2] == sign_langle_smooth[, 3]) / 
#     length(sign_langle_smooth[, 1])
#   return(c(two_pairs_result, all_pairs_result))
# }
# 

