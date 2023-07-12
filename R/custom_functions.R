# Adapted from 'as_detect' from the asdetect package to control the MAD threshold
# wl: window length

as_detect_mad <-
  function (ts, dt = NA, lowwl = 5, highwl = "default", mad_thr=3, mad_cst=1.4826){
    l <- length(ts)
    if (highwl == "default") {
      higwl <- floor(l/3)
    } else {
      higwl <- highwl
    }
    tip_distrib <- rep(0, l) # counts how often time points are involved in windows exceeding 3*MAD
    wls <- lowwl:higwl # set of window lengths
    # grad <- list() # added: list of gradients

    for (j in 1:length(wls)) { # For the different window length
      breaks <- l/wls[j]
      if (breaks != floor(breaks)) {
        remainder <- l - floor(breaks) * wls[j]
        ts_temp <- ts[(floor(remainder/2) + 1):(floor(remainder/2) +
                                                  floor(breaks) * wls[j])]
        # ts_temp: longest time series divisible by window length
      }
      if (breaks == floor(breaks)) {
        remainder <- NA
        ts_temp <- ts
      }

      lm_coeffs <- array(NA, dim = c(breaks, 2))
      # lm_coeff: array to store intercept and slope of lm from each window of a given length

      for (i in 1:breaks) { # For the different windows
        lm_mod <- lm(ts_temp[1:wls[j] + ((i - 1) * wls[j])] ~
                       c(1:wls[j])) # lm_mod: result of the linear regression made on a given window
        lm_coeffs[i, ] <- lm_mod$coefficients
      }

      outlier_ind_up <- which((lm_coeffs[, 2] - median(lm_coeffs[, 2]))/mad(lm_coeffs[, 2], constant=mad_cst) > mad_thr)
      outlier_ind_down <- which((lm_coeffs[, 2] - median(lm_coeffs[, 2]))/mad(lm_coeffs[, 2], constant=mad_cst) < -mad_thr)
      # grad[[j]] <- rep((lm_coeffs[, 2] - median(lm_coeffs[, 2]))/mad(lm_coeffs[, 2]),each=wls[j])
      # outlier_ind_up: list of windows with a slope above 3*MAD
      # outlier_ind_down: list of windows with a slope below -3*MAD


      # tip_distrib: counts how often time points are involved in windows exceeding 3*MAD
      if (is.na(remainder)) {
        for (k in outlier_ind_up) {
          tip_distrib[1:wls[j] + ((k - 1) * wls[j])] <- tip_distrib[1:wls[j] +
                                                                      ((k - 1) * wls[j])] + 1
        }
        for (k in outlier_ind_down) {
          tip_distrib[1:wls[j] + ((k - 1) * wls[j])] <- tip_distrib[1:wls[j] +
                                                                      ((k - 1) * wls[j])] - 1
        }
      }
      if (!is.na(remainder)) {
        for (k in outlier_ind_up) {
          tip_distrib[(floor(remainder/2) + 1) + 1:wls[j] +
                        ((k - 1) * wls[j])] <- tip_distrib[(floor(remainder/2) +
                                                              1) + 1:wls[j] + ((k - 1) * wls[j])] + 1
        }
        for (k in outlier_ind_down) {
          tip_distrib[(floor(remainder/2) + 1) + 1:wls[j] +
                        ((k - 1) * wls[j])] <- tip_distrib[(floor(remainder/2) +
                                                              1) + 1:wls[j] + ((k - 1) * wls[j])] - 1
        }
      }
    }

    # dt: time step of the time series, if it is one leave default dt=NA
    if (is.na(dt)) {
      result <- tip_distrib/length(wls)
      # result: cumulative contribution from all window lengths divided by the number of window lengths
    } else {
      t <- rep(NA, length(tip_distrib))
      t[1] <- 0
      for (i in 2:length(t)) {
        t[i] <- t[i - 1] + dt
      }
      result <- data.frame(t = t, detect = tip_distrib/length(wls))
    }
    return(result)
    # return(list("result"=result, "grad"=grad))
  }

environment(as_detect_mad) <- asNamespace('asdetect') # to call hidden functions from asdetect



# Adapted from 'shift_type' from the asdetect package to solution a limit case

custom_shift_type <- function (ts, where_as_pos, dt = FALSE, width = "tenth")
{
  if (width == "tenth") {
    w <- floor(length(ts)/20) # shift_type needs time series of 20 or more
  }
  else {
    w <- floor(width/2)
  }
  if (dt == FALSE) {
    pos <- where_as_pos
  }
  else {
    pos <- where_as_pos/dt
    if (width != "tenth") {
      w <- w/dt
    }
  }
  if (pos <= w) { # include case when pos == w otherwise (pos - w) will equal 0 and cause length problem below
    poss <- 1:(pos + w)
  }
  else if (pos > (length(ts) - w)) {
    poss <- (pos - w):length(ts)
  }
  else {
    poss <- (pos - w):(pos + w)
  }
  as_mod <- lm(ts[poss] ~ c(1:length(poss)))
  full_mod <- lm(ts ~ c(1:length(ts)))
  as_grad <- as_mod$coefficients[2]
  full_grad <- full_mod$coefficients[2]
  if (abs(as_grad) < abs(full_grad)) {
    marker <- 0
  }
  else {
    marker <- 1
  }
  return(marker)
}

environment(custom_shift_type) <- asNamespace('asdetect') # to call hidden functions from asdetect


# Adapted from 'where_as' from the asdetect package to make it silent

where_as_quiet <- function (ts, dt = NA, thresh = 0.7, quiet = FALSE)
{
  inds <- which(abs(ts) > thresh) # inds: position of time points with detection line above threshold

  if (length(inds) > 0) { # if any point above threshold

    encoding <- rle(diff(inds)) # encoding: lengths and values of runs of lag-1 difference

    # numtip <- length(which(encoding$length != 1))
    numtip <- length(which(encoding$value == 1)) # numtip: number of runs above threshold

    # whichtip <- which(encoding$length != 1)
    whichtip <- which(encoding$value == 1) # whichtip: position of runs in encoding

    tip_pos <- rep(NA, numtip) # tip_pos: list of position(s) with min or max detection score
    tip_prob <- rep(NA, numtip) # tip_prob: list of detection score(s) of min or max

    for (k in 1:numtip) {
      if (k == 1) {
        if (ts[inds[1]] > 0) {
          ## if maximum spread on several time point, it returns the first position:
          # tip_pos[k] <- inds[which.max(ts[inds[1:(encoding$length[1])]] + 1)]
          # tip_pos[k] <- inds[which.max( ts[ inds[1:(encoding$length[1] + 1)] ] )]

          ## if maximum spread on several time point, it returns the median position:
          tip_pos[k] <- floor(median(inds[which(ts[ inds[1:(encoding$length[1] + 1)] ] ==
                                                  max( ts[ inds[1:(encoding$length[1] + 1)] ]))]))
        }
        if (ts[inds[1]] < 0) {
          # tip_pos[k] <- inds[which.min(ts[inds[1:(encoding$length[1])]] + 1)]
          # tip_pos[k] <- inds[which.min( ts[ inds[1:(encoding$length[1] + 1)] ] )]

          tip_pos[k] <- floor(median(inds[which(ts[ inds[1:(encoding$length[1] + 1)] ] ==
                                                  min( ts[ inds[1:(encoding$length[1] + 1)] ]))]))
        }
        tip_prob[k] <- ts[tip_pos[k]]
      }
      else {

        inds_temp <- inds[sum(encoding$length[1:(whichtip[k] -
                                                   1)]):sum(encoding$length[1:whichtip[k]]) +
                            1] # inds_temp: for run k, positions of the kth run

        if (ts[inds_temp[1]] > 0) {
          # tip_pos[k] <- inds_temp[which.max(ts[inds_temp])]
          tip_pos[k] <- floor(median(inds_temp[which(ts[inds_temp] == max(ts[inds_temp]))]))
        }
        if (ts[inds_temp[1]] < 0) {
          # tip_pos[k] <- inds_temp[which.min(ts[inds_temp])]
          tip_pos[k] <- floor(median(inds_temp[which(ts[inds_temp] == min(ts[inds_temp]))]))
        }
        tip_prob[k] <- ts[tip_pos[k]]
      }
    }
    if (!is.na(dt)) {
      t <- rep(NA, length(ts))
      t[1] <- 0
      for (i in 2:length(t)) {
        t[i] <- t[i - 1] + dt
      }
      tip_pos <- t[tip_pos]
    }
    results <- list()
    results$as_pos <- tip_pos
    results$dt_val <- tip_prob
    return(results)
  }
  if (length(inds) == 0) { # if no point above threshold, still return the maximum with a warning
    results <- list()
    results$as_pos <- which.max(abs(ts))
    if (max(abs(ts)) == max(ts)) {
      results$dt_val <- max(ts)
    }
    if (-max(abs(ts)) == min(ts)) {
      results$dt_val <- min(ts)
    }
    if (!quiet) print("Threshold not detected, maximum returned instead")

    return(results)
  }
}

environment(where_as_quiet) <- asNamespace('asdetect') # to call hidden functions from asdetect

