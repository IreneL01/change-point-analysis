generate_data <- function(n, seed, left_censored_rate = 0.1, right_censored_rate = 0.1) {
  set.seed(seed)  
  
  N <- n * 500
  X1 <- rbinom(N, 1, 0.5)
  X2 <- runif(N, 0, 1)
  
  W <- runif(N)
  target_value <- -log(W)
  t_seq <- seq(18, 80, by = 0.1) 
  Lambda_t <- 0.6 * (0.04 * seq(0, 62, by = 0.1))^2 
  U <- t_seq[sapply(target_value/exp(gamma[1] * X1 + gamma[2] * X2),
                    function(value) which.min(abs(Lambda_t - value)))]
  
  # error terms
  eta <- rnorm(N, 0, 3)
  
  # follow-up time generation
  followup_start <- sample(seq(18, 80, by = 0.1), N, replace = TRUE)
  visit_counts <- sample(3:6, N, replace = TRUE) 
  total_visits <- sum(visit_counts)
  intervals <- sample(seq(1, 3, by = 0.1), total_visits, replace = TRUE) 
  ind  <- c(1, 1 + cumsum(visit_counts[1:(N-1)])) 
  intervals[ind] <- 0
  subjid <- rep(1:N, times = visit_counts)
  start_long <- rep(followup_start, times = visit_counts)
  
  # assemble follow-up table
  dt <- data.table(subjid = subjid, Interval = intervals,
                   followup_start = start_long)
  dt[, visdyAge := followup_start + cumsum(Interval), by = subjid]
  followup_end <- dt[, .(FU_end = max(visdyAge)), by = subjid][order(subjid)]$FU_end
  
  # censoring
  censoring <- rep("interval", N)
  censoring[U > followup_end] <- "right"
  censoring[U < followup_start] <- "left"
  censoring[U == followup_end] <- "deleted"
  censoring[U == followup_start] <- "deleted"
  #censoring[U > followup_end & (U - followup_end) >= 10] <- "deleted"
  #censoring[U < followup_start & (followup_start - U) >= 10] <- "deleted"
  censoring <- factor(censoring, levels = c("left", "interval", "right", "deleted"))
  
  df <- data.frame(subjid = 1:N, X1, X2, OnsetTime = U, eta,
                   censoring, followup_start, followup_end, visit_count = visit_counts)
  
  # stratified sampling
  n_left  <- round(n * left_censored_rate)
  n_right <- round(n * right_censored_rate)
  n_intv  <- n - n_left - n_right
  
  sampled_ids <- c(
    sample(df[df$censoring == "left", "subjid"], n_left),
    sample(df[df$censoring == "interval", "subjid"], n_intv),
    sample(df[df$censoring == "right", "subjid"], n_right)
  )
  
  final_df <- df[df$subjid %in% sampled_ids, ]
  final_df <- final_df[match(sampled_ids, final_df$subjid), ]
  final_IDs <- final_df$subjid
  
  followup_df <- dt[subjid %in% final_IDs]
  followup_df[, VisitIndex := 1:.N, by = subjid]
  followup_df <- followup_df[, .(subjid, visdyAge, VisitIndex)]
  
  final_df_info <- final_df %>%
    select(subjid, X1, X2, OnsetTime, eta, censoring, followup_start, followup_end, visit_count)
  
  final_long_df <- followup_df %>%
    left_join(final_df_info, by = "subjid") %>%
    rename(Censoring = censoring, VisitCount = visit_count)
  
  # compute U_interval
  U_bounds <- final_long_df %>%
    group_by(subjid) %>%
    summarise(U = first(OnsetTime),
              Times = list(visdyAge))
  get_surv_interval <- function(u, times_vec) {
    times_vec <- sort(times_vec)
    if (u < min(times_vec)) {
      surv1 <- 0
      surv2 <- times_vec[1]
    } else if (u > max(times_vec)) {
      surv1 <- tail(times_vec, 1)
      surv2 <- Inf
    } else {
      ind <- max(which(times_vec <= u))
      if(ind == length(times_vec)){
        left <- times_vec[ind - 1]
        right <- times_vec[ind]
      }else{
        left <- times_vec[ind]
        right <- times_vec[ind + 1]
      }
      surv1 <- left
      surv2 <- right
    }
    return(c(surv1, surv2))
  }
  LRinterval  <- mapply(get_surv_interval, U_bounds$U, U_bounds$Times, SIMPLIFY = FALSE)
  LRintervals <- do.call(rbind, LRinterval)
  U_bounds$survStartFmt <- LRintervals[, 1] # mapply(function(i){LRintervals[[i]][1]}, 1:n) # 
  U_bounds$survEndFmt <- LRintervals[, 2] # mapply(function(i){LRintervals[[i]][2]}, 1:n) # 
  
  final_long_df <- final_long_df %>%
    left_join(U_bounds, by = "subjid")
  final_long_df$eps <- rnorm(nrow(final_long_df), mean = 0, sd = sigma )
  # simulate yij
  final_long_df <- final_long_df %>%
    mutate(
      yij_true = a + b * visdyAge + c * pmax(visdyAge - OnsetTime + Delta, 0) +
        beta[1] * X1 + beta[2] * X2 + eta,
      yij = yij_true + eps
    )
  
  final_long_df$survEndFmt <- ifelse(is.infinite(final_long_df$survEndFmt), NA, final_long_df$survEndFmt)
  final_long_df$survStartFmt <- ifelse(final_long_df$survStartFmt == 0, NA, final_long_df$survStartFmt)
  
  final_long <- final_long_df %>%
    select(subjid, X1, X2, yij, yij_true, OnsetTime, visdyAge, survStartFmt, survEndFmt, eta, eps, Censoring)
  return(final_long_df)
}