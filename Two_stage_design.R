optimalTwoStageDesign <- function(a, b, p0, p1) {
  # save results
  .ls <- vector("list", 8)
  names(.ls) <- c("n", "n1", "r1", "r", "a", "power", "EN", "PET")
  # probability of rejecting the drug
  probReject <- function(r1, r, p, n1, n2) {
    p_stage_1 <- pbinom(r1, n1, p)
    p_stage_2 <- 0
    for (x in (r1 + 1):min(n1, r)) {
      p_stage_2 <- p_stage_2 + dbinom(x, n1, p) * pbinom(r - x, n2, p)
    }
    return(p_stage_1 + p_stage_2)
  }
  # n range
  p_bar <- (p0 + p1) / 2
  n_start <- floor(p_bar * (1 - p_bar) * ((qnorm(1 - a) + qnorm(1 - b)) / (p1 - p0))**2)
  n_min <- max(n_start - 20, 10)
  n_max <- min(n_start + 30, 90)
  cat("The recommended lower value of n is:", n_start, "\n")
  cat("Search between", n_min, "and", n_max, "\n\n")
  # record the smallest EN and n
  n_optimal <- 1000
  en_min <- 1000
  for (n in n_min:n_max) {
    for (n1 in 1:(n - 1)) {
      # if the sample size for the first stage is larger than the min(EN), stop
      if (n1 >= n_optimal) {break}
      for (r1 in 1:(n1 - 1)) {
        # prob of early stop
        PET <- pbinom(r1, n1, p0)
        # expected total
        EN <- n1 + (1 - PET) * (n - n1)
        r_max <- NULL
        p1_reject_current <- NULL
        for (r in (r1 + 1):(n - 1)) {
          # check type 2 error
          p1_reject <- probReject(r1, r, p1, n1, n - n1)
          if (p1_reject < b) {
            r_max <- r
            p1_reject_current <- p1_reject
          }
        }
        # if found the r met the type II error
        if (!is.null(r_max)) {
          # check type 1 error
          p0_reject <- probReject(r1, r_max, p0, n1, n - n1)
          if (p0_reject >= 1 - a) {
            # record EN and n for optimal design
            if (EN < en_min) {
              en_min <- EN
              n_optimal <- n
            }
            .ls[["n"]] <- c(.ls[["n"]], n)
            .ls[["n1"]] <- c(.ls[["n1"]], n1)
            .ls[["r1"]] <- c(.ls[["r1"]], r1)
            .ls[["r"]] <- c(.ls[["r"]], r_max)
            .ls[["a"]] <- c(.ls[["a"]], 1 - p0_reject)
            .ls[["power"]] <- c(.ls[["power"]], 1 - p1_reject_current)
            .ls[["EN"]] <- c(.ls[["EN"]], EN)
            .ls[["PET"]] <- c(.ls[["PET"]], PET)
          }
        }
      }
    }
  }
  if (is.null(.ls[["n"]])) {
    stop("No eligible results found")
  } else {
    .df <- as.data.frame(.ls)
    rownames(.df) <- NULL
    # print out optimal design
    cat("\nThe Optimal design output is:\n")
    print(round(.df[which.min(.df$EN), ], 4))
    # print out minmax design
    cat("\nThe MiniMax design output is:\n")
    print(round(subset(a, n == min(n) & EN == min(EN[n == min(n)])), 4))
    return(.df)
  }
}
tic()
a <- optimalTwoStageDesign(0.1, 0.1, 0.25, 0.4)
toc()

# library(data.table)

# TwoStageDesign <- function(a, b, p0, p1) {
#   # n range
#   p_bar <- (p0 + p1) / 2
#   n_start <- floor(p_bar * (1 - p_bar) * ((qnorm(1 - a) + qnorm(1 - b)) / (p1 - p0))**2)
#   n_min <- max(n_start - 20, 10)
#   n_max <- min(n_start + 30, 80)
#   # probability of rejecting the drug
#   probReject <- function(r1, r, p, n1, n2) {
#     p_stage_1 <- pbinom(r1, n1, p)
#     p_stage_2 <- 0
#     for (x in (r1 + 1):min(n1, r)) {
#       p_stage_2 <- p_stage_2 + dbinom(x, n1, p) * pbinom(r - x, n2, p)
#     }
#     return(p_stage_1 + p_stage_2)
#   }
#   # a data.table to save results
#   DT <- data.table(
#     n = n_min:n_max,
#     p0 = p0,
#     p1 = p1,
#     a = a,
#     b =  b
#   )
#   DT <- DT[, data.table(n1 = 1:(n - 1)), by = n][
#     , data.table(r1 = 1:(n1 - 1)), .(n, n1)
#   ][
#     , data.table(r = (r1 + 1):(n - 1)), .(n, n1, r1)
#   ]
#   DT[, p1_reject := probReject(r1, r, p1, n1, n - n1), by = seq_len(nrow(DT))]
#   DT <- DT[DT[p1_reject <= b, .I[which.max(r)], .(n, n1, r1)]$V1]
#   DT[, p0_reject := probReject(r1, r, p0, n1, n - n1), by = seq_len(nrow(DT))]
#   DT <- DT[p0_reject >= 1 - a]
#   # prob of early stop
#   DT[, PET := pbinom(r1, n1, p0)]
#   # expected total
#   DT[, EN := n1 + (1 - PET) * (n - n1)]
#   return(DT)
# }
# tic()
# b <- TwoStageDesign(0.05, 0.2, 0.05, 0.15)
# toc()


