
# Data preparation --------------------------------------------------------

#' Trajectory classification (from Rigal et al. 2020)
#'
#' @param df a data frame stemming from run_simu function
#' @param thr unfeasibility threshold
#' @param type a character specifying whether the series come from simulations ("sim") or from RAMLDB ("RAM")
#'
#' @return a list of data frames ready for analyses
#' @export

prep_data <- function(df, thr=0, type="sim", apriori){

  if (type=="sim"){

    if ("iter" %in% names(df)){

      ts_type <- "TB"
      dataset <- df %>%
        {if (apriori==TRUE) dplyr::select(., scen, iter, year, TB,
                                          expected_traj, expected_class) else .} %>%
        {if (apriori==FALSE) dplyr::select(., scen, iter, year, TB) else .} %>%
        dplyr::rename(Y=TB)

      # %>%
      # dplyr::mutate(scen=factor(scen, levels = unique(df$scen)))
      iter <- seq_len(max(df$iter))
      scen_list <- df %>% dplyr::distinct(scen) %>% pull()
    }

  } else if (type=="RAM"){

    ts_type <- names(df)[3]
    iter <- 1
    dataset <- df %>%
      dplyr::rename(Y=ts_type) %>%
      dplyr::mutate(iter = 1)

    scen_list <- df %>% dplyr::distinct(scen) %>% pull()

  } else if (type=="data"){

    ts_type <- names(df)[3]
    iter <- 1
    dataset <- df %>%
      dplyr::rename(Y=ts_type) %>%
      dplyr::mutate(iter = 1)

    scen_list <- df %>% dplyr::distinct(scen) %>% pull()

  }

  sets <- list()
  for (j in 1:length(scen_list)){

    by_scen <- dataset %>%
      dplyr::filter(scen==scen_list[j]) %>%
      split(.$iter)

    dats <- lapply(by_scen, function(dat){

      X <- dat$year
      Y <- dat$Y
      Y_SE <- abs(rnorm(length(X), mean = 0.05*(max(Y)-min(Y)), sd = 0.005*(max(Y)-min(Y))))

      if (apriori==TRUE){

        expected_traj <- dat$expected_traj
        expected_class <- dat$expected_class

        set <- data.frame(X=X, Y=Y, Y_SE=Y_SE, expected_traj=expected_traj, expected_class=expected_class)

      } else if (apriori==FALSE){

        set <- data.frame(X=X, Y=Y, Y_SE=Y_SE)
      }

    })

    sets <- c(sets, dats)
  }
  if (type=="sim"){
    names(sets) <- paste0(rep(scen_list, each=max(df$iter)),"_iter", sprintf("%02d", 1:max(df$iter)))
  } else {
    names(sets) <- scen_list
  }


  # sets <- list()
  # for (j in 1:length(scen_list)){
  #
  #   for (i in iter){
  #
  #     dat <- dataset %>% dplyr::filter(scen==scen_list[j] & iter==i)
  #
  #     X <- dat$year
  #     Y <- dat$Y
  #     Y_SE <- abs(rnorm(length(X), mean = 0.05*(max(Y)-min(Y)), sd = 0.005*(max(Y)-min(Y))))
  #
  #     if (apriori==TRUE){
  #
  #       expected_traj <- dat$expected_traj
  #       expected_class <- dat$expected_class
  #
  #       set <- data.frame(X=X, Y=Y, Y_SE=Y_SE, expected_traj=expected_traj, expected_class=expected_class)
  #
  #     } else if (apriori==FALSE){
  #
  #       set <- data.frame(X=X, Y=Y, Y_SE=Y_SE)
  #
  #     }
  #
  #     sets[[ paste0(scen_list[j],"_iter", sprintf("%02d", i)) ]] <- set
  #   }
  #
  # }

  if (!is.null(thr)){
    sets <- sets %>%
      lapply(function(x) x %>% dplyr::filter(Y>thr))
  }


return(list("ts"=sets, "ts_type"=ts_type))

}


# Scale simulated time series (custom scaling)

scale_simu <- function(simu_list, name){

  scaling <- simu_list[[1]] %>%
    dplyr::group_by(scen, iter) %>%
    dplyr::mutate(magn = max(TB)-min(TB)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(scen) %>%
    dplyr::summarise(magn=mean(magn))

  scl_simu_list <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::left_join(scaling, by="scen") %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::mutate(TB = ifelse(expected_class!="no_change", (TB-mean(TB))/magn*2, TB),
                                          # TB = (TB-mean(TB))/magn*2,
                                          TB = TB+(5-first(TB))
                            ) %>%
                            dplyr::ungroup()
  )

  saveRDS(scl_simu_list, paste0("data/03_simulations/all_simu_", name, "_scl.rds"))

  return(scl_simu_list)

}

# Scale by mean simulated time series

scale_mean_simu <- function(simu_list, name){

  scale_mean_simu_list <- lapply(simu_list,
                           function(x) x %>%
                             dplyr::group_by(scen, iter) %>%
                             dplyr::mutate(TB = TB/mean(TB)) %>%
                             dplyr::ungroup()
  )

  saveRDS(scale_mean_simu_list, paste0("data/03_simulations/all_simu_", name, "_mnscl.rds"))

  return(scale_mean_simu_list)

}


# Normalize simulated time series (rescale min-max)

norm_simu <- function(simu_list, name){

  norm_simu_list <- lapply(simu_list,
                           function(x) x %>%
                             dplyr::group_by(scen, iter) %>%
                             dplyr::mutate(
                               # TB = ifelse(expected_class!="no_change", (TB-min(TB))/(max(TB)-min(TB)), TB)
                               TB = (TB-min(TB))/(max(TB)-min(TB))
                             ) %>%
                             dplyr::ungroup()
                           )

  saveRDS(norm_simu_list, paste0("data/03_simulations/all_simu_", name, "_nrm.rds"))

  return(norm_simu_list)

}


# Standardize simulated time series

strd_simu <- function(simu_list, name){

  strd_simu_list <- lapply(simu_list_nonstd,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::mutate(
                              # TB = ifelse(expected_class!="no_change", (TB-mean(TB))/sd(TB)+min(TB), TB)
                              TB = (TB-mean(TB))/sd(TB),
                              TB = TB - min(TB)
                            ) %>%
                            dplyr::ungroup()
  )

  saveRDS(strd_simu_list, paste0("data/03_simulations/all_simu_", name, "_std.rds"))

  return(strd_simu_list)

}





# List of stocks with time series available for classification:

list_available <- function(ts_type, min_len) {

  ts <- timeseries_values_views %>%
    dplyr::group_by(stockid) %>%
    dplyr::select(stockid, year, tidyselect::all_of(ts_type)) %>%
    na.omit(ts_type) %>%
    dplyr::summarise(first=min(year),
                     last=max(year)) %>%
    dplyr::mutate(length = last-first+1) %>%
    dplyr::mutate(scen = paste(ts_type, stockid, sep="_")) %>%
    dplyr::relocate(scen) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metadata, by="stockid") %>%
    dplyr::select(-c(FisheryType, taxGroup)) %>%
    dplyr::left_join(taxonomy, by="scientificname") %>%
    dplyr::left_join(gaps_info(ts_type) %>% dplyr::select(stockid, gaps, missing), by="stockid") %>%
    dplyr::mutate(gaps = ifelse(is.na(gaps), FALSE, TRUE)) %>%
    dplyr::filter(length>=min_len & !gaps)

  return(ts)
}


# List of stocks with time series available for classification specifying start and end:

list_available_lim <- function(ts_type, min_len, min_y, max_y) {

  ts <- timeseries_values_views %>%
    dplyr::group_by(stockid) %>%
    dplyr::select(stockid, year, tidyselect::all_of(ts_type)) %>%
    na.omit(ts_type) %>%
    dplyr::summarise(first=min(year),
                     last=max(year)) %>%
    dplyr::mutate(length = last-first+1) %>%
    dplyr::mutate(scen = paste(ts_type, stockid, sep="_")) %>%
    dplyr::relocate(scen) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metadata, by="stockid") %>%
    dplyr::select(-c(FisheryType, taxGroup)) %>%
    dplyr::left_join(taxonomy, by="scientificname") %>%
    dplyr::left_join(gaps_info(ts_type) %>% dplyr::select(stockid, gaps, missing), by="stockid") %>%
    dplyr::mutate(gaps = ifelse(is.na(gaps), FALSE, TRUE)) %>%
    dplyr::filter(length>=min_len & !gaps & first<=min_y & last>=max_y)

  return(ts)
}




# Trajectory classification -----------------------------------------------

#' #' Trajectory classification (from Rigal et al. 2020)
#' #'
#' #' @param Y vector of the time series values
#' #' @param X vector of the timesteps
#' #' @param dataset two-column dataframe with timesteps and values of the time series
#' #' @param interval_size
#' #'
#' #' @return a single-row data frame with infos about the trajectory
#' #' @export
#'
#' class_trajectory <- function (Y = NULL, X = NULL, dataset = NULL, interval_size = 0.5)
#' {
#'   if (is.null(X) == TRUE & is.null(Y) == TRUE & is.null(dataset) == TRUE){
#'     stop("either 'dataset' or at least 'Y' and 'X' must be specified")
#'   }
#'   if (is.null(X) == TRUE & is.null(Y) == TRUE) {
#'     Y <- dataset[,2]
#'     X <- dataset[,1]
#'   }else{
#'     if (class(Y) == "character" & class(X) == "character") {
#'       if (is.null(dataset) == TRUE) {
#'         stop("if 'Y' and 'X' are character, 'dataset' must exist")
#'       }else{
#'         Y <- dataset[, Y]
#'         X <- dataset[, X]
#'       }
#'     }else{
#'       if (!(class(Y) %in% c("numeric","integer")) == TRUE & !(class(X) %in% c("numeric","integer")) == TRUE) {stop("'Y' and 'X' must be either characters or vector but 'class' must be similar")}
#'     }
#'   }
#'
#'   data <- data.frame(cbind(Y, X))
#'   data <- data[order(data$X),] # ordering the X values
#'
#'   if (length(X)<4){
#'     stop("time series length must be at least 4")
#'   }
#'
#'   Y <- data$Y
#'   X <- data$X
#'
#'   linear.model <- lm(Y~X)
#'   aic_lin <- AIC(linear.model)
#'
#'   orthogonal_polynomial <- lm(Y~poly(X,2, raw=F))
#'   aic_pol <- AIC(orthogonal_polynomial)
#'
#'   # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial
#'   # we have to perform a variable change to obtain relevant values in the X interval
#'   # for first_order_coefficient, second_order_coefficient and intercept,
#'   # knowing that X'= alpha*X + beta
#'   # and chi = eta*X'^2 + theta
#'
#'   gammab  <-  orthogonal_polynomial$coefficients[3]
#'   delta  <-  orthogonal_polynomial$coefficients[2]
#'   epsilon  <-  orthogonal_polynomial$coefficients[1]
#'
#'   alpha  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[2]
#'   beta  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[1]
#'
#'   eta  <-  1/lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[2]
#'   theta  <-  (-lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[1])*eta
#'
#'   Y2<-Y*(max(X)-min(X))/(max(Y)-min(Y)) # p2 and p3 are relevant when Y and X amplitudes are equivalent,
#'   # in particular when studying scaled-to-1 indices, Y and X amplitudes
#'   # may be very different, so we scaled the amplitudes to calculate p2 and p3
#'   polynomial_orthonormal_basis<-lm(Y2~poly(X,2, raw=T))$coefficients
#'
#'   if(summary(orthogonal_polynomial)$coefficients[3, 4] <= 0.05){  # non linear case
#'     classification <- data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
#'                                  first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
#'                                  second_order_coefficient = (alpha^2)*gammab*eta,
#'                                  second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
#'                                  strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
#'                                  intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
#'                                  x_m = (X[length(X)]-X[1])/2+X[1],
#'                                  p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta), # points of interest
#'                                  p2 = (-polynomial_orthonormal_basis[2]+1)/(2*polynomial_orthonormal_basis[3]),
#'                                  p3 = (-polynomial_orthonormal_basis[2]-1)/(2*polynomial_orthonormal_basis[3]))
#'
#'     aic_max_shape <- aic_pol
#'
#'   }else{ # linear case
#'     classification <- data.frame(first_order_coefficient = delta*alpha,
#'                                  first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
#'                                  second_order_coefficient = 0,
#'                                  second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
#'                                  strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
#'                                  intercept = epsilon+delta*beta,
#'                                  x_m = (X[length(X)]-X[1])/2+X[1],
#'                                  p1 = NA,
#'                                  p2 = NA,
#'                                  p3 = NA)
#'     aic_max_shape <- aic_lin
#'
#'   }
#'
#'   classification$r.sq <- summary(orthogonal_polynomial)$adj.r.squared # retrieve the adjusted coefficient of determination
#'
#'   # compute the derivative at xm-delta and at xm + delta with delta being half of the input interval size
#'   derivative  <-  2*(classification$x_m-(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
#'   derivative2  <-  2*(classification$x_m+(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
#'
#'
#'   if(sign(derivative) != sign(derivative2)){ # non consistent direction around x_m
#'     classification$derivative  <-  NA
#'     classification$intercept_derivative  <-  NA
#'   }else{ # consistent direction around x_m
#'     classification$derivative  <-  mean(c(derivative, derivative2))
#'     classification$intercept_derivative  <-  (classification$second_order_coefficient*classification$x_m^2+classification$first_order_coefficient*classification$x_m+classification$intercept)-classification$x_m*classification$derivative
#'   }
#'
#'   # compute the derivative of the curvature function
#'   classification$derivated_curvature  <-  -12*(classification$second_order_coefficient^2)*(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)*(classification$second_order_coefficient/abs(classification$second_order_coefficient))/
#'     ((1+(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)^2)^(2.5))
#'
#'   if(classification$second_order_pvalue>0.05){classification$derivated_curvature <- NA}
#'
#'   classification$direction <- NA # classify the direction
#'   classification$direction[which(classification$derivative > 0)] <- "increase"
#'   classification$direction[which(classification$derivative < 0)] <- "decrease"
#'   classification$direction[which(is.na(classification$derivative))] <- "stable"
#'   classification$direction[which(as.numeric(classification$first_order_pvalue)>0.05 & as.numeric(classification$second_order_pvalue)>0.05)] <- "stable"
#'
#'   classification$acceleration <- NA # classify the acceleration
#'   classification$acceleration[which(classification$derivated_curvature < 0)] <- "accelerated"
#'   classification$acceleration[which(classification$derivated_curvature > 0)] <- "decelerated"
#'   classification$acceleration[which(classification$direction == "stable" &
#'                                       classification$second_order_coefficient < 0)] <- "concave"
#'   classification$acceleration[which(classification$direction == "stable" &
#'                                       classification$second_order_coefficient > 0)] <- "convex"
#'   classification$acceleration[which(is.na(classification$derivated_curvature))] <- "constant"
#'
#'   classification$shape_class <- paste(classification$direction, # give the final classification combining direction and acceleration
#'                                       classification$acceleration,
#'                                       sep="_")
#'
#'   linear.model.summary <- summary(linear.model) # provide the linear approach results for comparison
#'
#'   classification$linear_slope <- linear.model.summary$coefficients[2, 1]
#'   classification$linear_slope_pvalue <- linear.model.summary$coefficients[2, 4]
#'   classification$linear_intercept <- linear.model.summary$coefficients[1, 1]
#'
#'   classification$first_X_value <- X[1]
#'   classification$last_X_value <- X[length(X)]
#'   classification$aic_lin <- aic_lin
#'   classification$aic_pol <- aic_pol
#'   classification$aic_max_shape <- aic_max_shape
#'
#'   row.names(classification) <- "Y"
#'
#'   return(classification)
#'
#' }


#' Trajectory classification (from Rigal et al. 2020)
#'
#' @param Y vector of the time series values
#' @param X vector of the timesteps
#' @param dataset two-column dataframe with timesteps and values of the time series
#' @param interval_size
#'
#' @return a two-row data frame with infos about the trajectory (polynomial and linear)
#' @export

class_trajectory_both <- function (Y = NULL, X = NULL, dataset = NULL, interval_size = 0.5){

  if (is.null(X) == TRUE & is.null(Y) == TRUE & is.null(dataset) == TRUE){
    stop("either 'dataset' or at least 'Y' and 'X' must be specified")
  }
  if (is.null(X) == TRUE & is.null(Y) == TRUE) {
    Y <- dataset[,2]
    X <- dataset[,1]
  }else{
    if (class(Y) == "character" & class(X) == "character") {
      if (is.null(dataset) == TRUE) {
        stop("if 'Y' and 'X' are character, 'dataset' must exist")
      }else{
        Y <- dataset[, Y]
        X <- dataset[, X]
      }
    }else{
      if (!(class(Y) %in% c("numeric","integer")) == TRUE & !(class(X) %in% c("numeric","integer")) == TRUE) {stop("'Y' and 'X' must be either characters or vector but 'class' must be similar")}
    }
  }

  data <- data.frame(cbind(Y, X))
  data <- data[order(data$X),] # ordering the X values

  if (length(X)<4){
    stop("time series length must be at least 4")
  }

  Y <- data$Y
  X <- data$X

  null.model <- lm(Y~1)
  # aic_nch <- AIC(null.model) # stats::AIC.default
  # rsq_nch <- summary(null.model)$adj.r.squared
  nrmse_nch <- sqrt(sum(summary(null.model)$residuals^2)/length(Y))/sd(Y)
  aic_nch <- MuMIn::AICc(null.model)
  # stats4::logLik(null.model) # (df=2)

  linear.model <- lm(Y~X)
  # aic_lin <- AIC(linear.model)
  # rsq_lin <- summary(linear.model)$adj.r.squared
  nrmse_lin <- sqrt(sum(summary(linear.model)$residuals^2)/length(Y))/sd(Y)
  aic_lin <- MuMIn::AICc(linear.model)
  # stats4::logLik(linear.model) # (df=3)
  # aic_lin <- extractAIC(linear.model)[2]

  orthogonal_polynomial <- lm(Y~poly(X,2, raw=F))
  # aic_pol <- AIC(orthogonal_polynomial)
  # rsq_pol <- summary(orthogonal_polynomial)$adj.r.squared
  nrmse_pol <- sqrt(sum(summary(orthogonal_polynomial)$residuals^2)/length(Y))/sd(Y)
  aic_pol <- MuMIn::AICc(orthogonal_polynomial)
  # stats4::logLik(orthogonal_polynomial) # (df=4)
  # aic_pol <- extractAIC(orthogonal_polynomial)[2]

  # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial
  # we have to perform a variable change to obtain relevant values in the X interval
  # for first_order_coefficient, second_order_coefficient and intercept,
  # knowing that X'= alpha*X + beta
  # and chi = eta*X'^2 + theta

  gammab  <-  orthogonal_polynomial$coefficients[3]
  delta  <-  orthogonal_polynomial$coefficients[2]
  epsilon  <-  orthogonal_polynomial$coefficients[1]

  alpha  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[2]
  beta  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[1]

  eta  <-  1/lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[2]
  theta  <-  (-lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[1])*eta

  Y2<-Y*(max(X)-min(X))/(max(Y)-min(Y)) # p2 and p3 are relevant when Y and X amplitudes are equivalent,
  # in particular when studying scaled-to-1 indices, Y and X amplitudes
  # may be very different, so we scaled the amplitudes to calculate p2 and p3
  polynomial_orthonormal_basis<-lm(Y2~poly(X,2, raw=T))$coefficients

  # non linear case
  classification <- data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
                               first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                               second_order_coefficient = (alpha^2)*gammab*eta,
                               second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                               strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                               intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
                               x_m = (X[length(X)]-X[1])/2+X[1],
                               p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta), # points of interest
                               p2 = (-polynomial_orthonormal_basis[2]+1)/(2*polynomial_orthonormal_basis[3]),
                               p3 = (-polynomial_orthonormal_basis[2]-1)/(2*polynomial_orthonormal_basis[3]),
                               aic = aic_pol,
                               nrmse = nrmse_pol)

  # linear case
  classification[2,] <- data.frame(first_order_coefficient = delta*alpha,
                                   first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                   second_order_coefficient = 0,
                                   second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                   strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                   intercept = epsilon+delta*beta,
                                   x_m = (X[length(X)]-X[1])/2+X[1],
                                   p1 = NA,
                                   p2 = NA,
                                   p3 = NA,
                                   aic = aic_lin,
                                   nrmse = nrmse_lin)

  # no change case
  classification[3,] <- data.frame(first_order_coefficient = 0,
                                   first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                   second_order_coefficient = 0,
                                   second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                   strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                   intercept = null.model$coefficients,
                                   x_m = (X[length(X)]-X[1])/2+X[1],
                                   p1 = NA,
                                   p2 = NA,
                                   p3 = NA,
                                   aic = aic_nch,
                                   nrmse = nrmse_nch)


  # Significance of the coefficients:
  if(summary(orthogonal_polynomial)$coefficients[3, 4] <= 0.05){  # non linear case
    classification$best_model <- c("best", NA, NA)
  } else if(summary(orthogonal_polynomial)$coefficients[2, 4]<=0.05){ # linear case
    classification$best_model <- c(NA, "best", NA)
  } else { # no change case
    classification$best_model <- c(NA, NA, "best")
  }

  # classification$r.sq <- summary(orthogonal_polynomial)$adj.r.squared # retrieve the adjusted coefficient of determination

  for (i in 1:3){

    # compute the derivative at xm-delta and at xm + delta with delta being half of the input interval size
    derivative  <-  2*(classification$x_m[i]-(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient[i]+classification$first_order_coefficient[i]
    derivative2  <-  2*(classification$x_m[i]+(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient[i]+classification$first_order_coefficient[i]


    if(sign(derivative) != sign(derivative2) | i==3){ # non consistent direction around x_m
      classification$derivative[i]  <-  NA
      classification$intercept_derivative[i]  <-  NA
    }else{ # consistent direction around x_m
      classification$derivative[i]  <-  mean(c(derivative, derivative2))
      classification$intercept_derivative[i]  <-  (classification$second_order_coefficient[i]*classification$x_m[i]^2+classification$first_order_coefficient[i]*classification$x_m[i]+classification$intercept[i])-classification$x_m[i]*classification$derivative[i]
    }

    # compute the derivative of the curvature function
    classification$derivated_curvature[i]  <-  -12*(classification$second_order_coefficient[i]^2)*(2*classification$second_order_coefficient[i]*classification$x_m[i]+classification$first_order_coefficient[i])*(classification$second_order_coefficient[i]/abs(classification$second_order_coefficient[i]))/
      ((1+(2*classification$second_order_coefficient[i]*classification$x_m[i]+classification$first_order_coefficient[i])^2)^(2.5))

    if(classification$second_order_pvalue[i]>0.05 & i != 1){ # keep derivated curvature even if not significant for polynomial fit
      classification$derivated_curvature[i] <- NA
    }

    classification$direction[i] <- NA # classify the direction
    classification$direction[i][which(classification$derivative[i] > 0)] <- "increase"
    classification$direction[i][which(classification$derivative[i] < 0)] <- "decrease"
    classification$direction[i][which(is.na(classification$derivative[i]))] <- "stable"
    classification$direction[i][which(as.numeric(classification$first_order_pvalue[i])>0.05 & as.numeric(classification$second_order_pvalue[i])>0.05)] <- "stable"

    classification$acceleration[i] <- NA # classify the acceleration
    classification$acceleration[i][which(classification$derivated_curvature[i] < 0)] <- "accelerated"
    classification$acceleration[i][which(classification$derivated_curvature[i] > 0)] <- "decelerated"
    classification$acceleration[i][which(classification$direction[i] == "stable" &
                                           classification$second_order_coefficient[i] < 0)] <- "concave"
    classification$acceleration[i][which(classification$direction[i] == "stable" &
                                           classification$second_order_coefficient[i] > 0)] <- "convex"
    classification$acceleration[i][which(is.na(classification$derivated_curvature[i]))] <- "constant"

    classification$shape_class[i] <- paste(classification$direction[i], # give the final classification combining direction and acceleration
                                           classification$acceleration[i],
                                           sep="_")
  }

  if(classification$shape_class[2] == "stable_constant"){
    aic_lin <- classification[2,]$aic <- AIC(linear.model)-2 # To remove what corresponds to one useless parameter
  }

  linear.model.summary <- summary(linear.model) # provide the linear approach results for comparison

  classification$linear_slope <- linear.model.summary$coefficients[2, 1]
  classification$linear_slope_pvalue <- linear.model.summary$coefficients[2, 4]
  classification$linear_intercept <- linear.model.summary$coefficients[1, 1]

  classification$first_X_value <- X[1]
  classification$last_X_value <- X[length(X)]

  row.names(classification) <- c("Y_pol","Y_lin", "Y_nch")

  return(classification)

}



#' Markov-chain simulations of classification
#'
#' @param dataset a data frame ready for analyses
#' @param niter a number of MC simulations (no resampling done if niter = 1)
#' @param ref_year reference year, by default equal to the mid year of the interval
#' @param correction set the reference value to 100 and correct values below 0 before logtransformation
#' @param fit a character to specify the type of fit ("pol" or "lin") default is best fit
#'
#' @return a data frame with as many rows as iterations
#' @export

mc_trend <- function(dataset,         # data
                     niter,           # number of MC simulations
                     ref_year=NULL,   # reference year, by default equal to the mid year of the interval
                     correction=FALSE, # set the reference value to 100 and correct values below 0 before logtransformation
                     fit=NULL) # specify the type of fit ("pol" or "lin") default is best fit
{

  b <- data.frame(t(rep(NA, 15)))
  attributes(b)$names <- c("second_order_coefficient",
                           "first_order_coefficient",
                           "strd_error",
                           "shape_class",
                           "intercept",
                           "p_1",
                           "p_2",
                           "p_3",
                           "slope_p_value",
                           "slope",
                           "ref_year",
                           "aic",
                           "best_model",
                           "second_order_pvalue",
                           "nrmse"
                           # "aic_lin",
                           # "aic_pol",
                           # "aic_max_shape"
  )

  if(is.null(ref_year)){
    ref_year <- dataset$X[round(nrow(dataset)/2)+1]
  }

  if(correction == TRUE){
    ref_value <- dataset$Y[dataset$X == ref_year]
    if(ref_value == 1){
      dataset$Y <- 100*dataset$Y # set reference year value to 100
      dataset$Y[which(dataset$Y <= 1)] <- 1 # set values < 1 to 1
      dataset$Y_SE[which(dataset$Y <= 1)] <- 0 # and their SE to 0
      dataset$Y_SE <- 100*dataset$Y_SE
    }
    if(ref_value!=1 & ref_value!=100){
      if(ref_value>1){
        dataset$Y <- dataset$Y/ref_value
        dataset$Y_SE <- dataset$Y_SE/ref_value
      }
      if(ref_value<1){
        stop("use 'correction = FALSE' when value of the reference year is strictly below 1")
      }
      dataset$Y <- 100*dataset$Y
      if(length(which(dataset$Y <= 1))>0){print("caution, low values corrected, if strongly decreasing or increasing trajectory, use respectively  first or last year as referential")}
      dataset$Y[which(dataset$Y <= 1)] <- 1
      dataset$Y_SE[which(dataset$Y <= 1)] <- 0
      dataset$Y_SE <- 100*dataset$Y_SE
    }
    if(ref_value == 100){
      dataset$Y[which(dataset$Y <= 1)] <- 1
      dataset$Y_SE[which(dataset$Y <= 1)] <- 0
    }

    dataset$sd <- dataset$Y_SE/dataset$Y # SE of log transformed data
    dataset$log <- log(dataset$Y) # log transforme Y
    for(j in 1:nrow(dataset)){
      if(dataset$sd[j]>(dataset$log[j])){dataset$sd[j] <- (dataset$log[j])} # set SE to value amplitude if SE > value (if not, it leads to huge values when resampling in the next loop)
    }
  }

  if(correction == FALSE){
    if(min(dataset$Y)<0){
      min_value <- abs(min(dataset$Y))
      dataset$Y <- dataset$Y+min_value
    }else{min_value <- 0}
    dataset$sd <- dataset$Y_SE
    dataset$log <- dataset$Y
  }

  for(i in 1:niter){

    if (niter==1){ # no resampling if only one iteration
      a <- dataset$Y

    } else { # simulate Y values from normal distribution (mean= original Y value, sd = original SE)
      a <- rnorm(nrow(dataset), mean=dataset$log, sd=dataset$sd)
    }

    if(correction == TRUE){
      a <- exp(a)/exp(a[which(dataset$X == ref_year)])*100 # set reference year value to 100 and retransform values if logtranformed
    }else{a <- a-min_value}

    a <- class_trajectory_both(a, dataset$X)
    best_model <- a %>% dplyr::filter(best_model=="best") %>% rownames() %>% sub("Y_","",.)
    if(is.null(fit)){
      a <- a %>% dplyr::filter(best_model=="best")
    }else if(fit == "nch"){
      a <- a %>% dplyr::filter(rownames(.)=="Y_nch")
    }else if(fit == "lin"){
      a <- a %>% dplyr::filter(rownames(.)=="Y_lin")
    }else if(fit == "pol"){
      a <- a %>% dplyr::filter(rownames(.)=="Y_pol")
    }


    b[i, 1] <- a$second_order_coefficient
    b[i, 2] <- a$first_order_coefficient
    b[i, 3] <- a$strd_error
    b[i, 4] <- a$shape_class
    b[i, 5] <- a$intercept
    if(a$second_order_coefficient!=0){
      if(findInterval(a$p1,  c(min(dataset$X), max(dataset$X))) == 1){ # record changing point inside time series
        b[i, 6] <- a$p1}else{b[i, 6] <- NA}
      if(findInterval(a$p2,  c(min(dataset$X), max(dataset$X))) == 1){
        b[i, 7] <- a$p2}else{b[i, 7] <- NA}
      if(findInterval(a$p3,  c(min(dataset$X), max(dataset$X))) == 1){
        b[i, 8] <- a$p3}else{b[i, 8] <- NA}
    }else{
      b[i, 6] <- NA
      b[i, 7] <- NA
      b[i, 8] <- NA
    }
    b[i, 9] <- a$linear_slope_pvalue
    b[i, 10] <- a$linear_slope
    b[i, 12] <- a$aic
    b[i, 13] <- best_model
    b[i, 14] <- a$second_order_pvalue
    b[i, 15] <- a$nrmse
    # b[i, 12] <- a$aic_lin
    # b[i, 13] <- a$aic_pol
    # b[i, 14] <- a$aic_max_shape
  }
  b[, 4] <- as.factor(b[, 4])
  b[, 11] <- rep(ref_year, nrow(b))

  return(b)
}


#' Summary of MC simulations of classification
#'
#' @param sets a list of data frame ready for analyses
#' @param niter a number of MC simulations
#' @param ref_year reference year, by default equal to the mid year of the interval
#' @param correction set the reference value to 100 and correct values below 0 before logtransformation
#' @param fit a character to specify the type of fit ("pol" or "lin") default is best fit
#'
#' @return a data frame with as many rows as simulations
#' @export

res_trend <- function(sets,
                      niter,
                      ref_year=NULL,
                      correction=FALSE,
                      fit=NULL){

  res <- data.frame()

  for (i in seq_len(length(sets))){

    if(nrow(sets[[i]])>3 & anyNA(sets[[i]]$Y_SE) == FALSE){

      simulated <- mc_trend(sets[[i]], niter, ref_year=NULL, correction=FALSE, fit)
      simulated <- simulated %>% dplyr::mutate(best_model=as.factor(best_model))

      if(length(levels(simulated$shape_class))>1){ # test the significance of the most numerous class
        test <- RVAideMemoire::multinomial.theo.multcomp(simulated$shape_class, p = rep(1/length(levels(simulated$shape_class)),
                                                                                        length(levels(simulated$shape_class))), prop=TRUE)
        # Not working if there's no linear alternative among the best shapes
        #   if(min(test$p.value2[test$observed>test$expected])<0.05){
        #     max_shape <- row.names(test$p.value)[which(test$observed == max(test$observed[test$observed>test$expected]))]
        #   }else{ # if non significant class, the linear one is chosen
        #     max_shape <-c("increase_constant", "decrease_constant", "stable_constant")[which.max(c(length(grep("increase", simulated$shape_class)), length(grep("decrease", simulated$shape_class)), length(grep("stable", simulated$shape_class))))]
        #   }

        # So temporarily not relying on the test
        max_shape <- row.names(test$p.value)[which(test$observed == max(test$observed[test$observed>test$expected]))]
        if (grepl("constant", max_shape)){
          best_model <- "lin"
        }else{best_model <- "pol"}
      }

      if(length(levels(simulated$shape_class)) == 1){max_shape <- levels(simulated$shape_class)
      best_model <- levels(simulated$best_model)}

      alpha2 <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 1]))
      sd_alpha2 <- sd(as.numeric(simulated[simulated$shape_class == max_shape, 1]))
      alpha1 <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 2]))
      sd_alpha1 <- sd(as.numeric(simulated[simulated$shape_class == max_shape, 2]))
      inter <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 5]))
      strd <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 3]))
      p_1 <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 6]), na.rm=T)
      sd_p_1 <- sd(as.numeric(simulated[simulated$shape_class == max_shape, 6]), na.rm=T)
      if( !is.na(p_1) && findInterval(p_1, c(min(sets[[i]]$X), max(sets[[i]]$X))) != 1){p_1 <- sd_p_1 <- as.numeric(NA)}
      p_2 <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 7]), na.rm=T)
      sd_p_2 <- sd(as.numeric(simulated[simulated$shape_class == max_shape, 7]), na.rm=T)
      if( !is.na(p_2) && findInterval(p_2, c(min(sets[[i]]$X), max(sets[[i]]$X))) != 1){p_2 <- sd_p_2 <- as.numeric(NA)}
      p_3 <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 8]), na.rm=T)
      sd_p_3 <- sd(as.numeric(simulated[simulated$shape_class == max_shape, 8]), na.rm=T)
      if( !is.na(p_3) && findInterval(p_3, c(min(sets[[i]]$X), max(sets[[i]]$X))) != 1){p_3 <- sd_p_3 <- as.numeric(NA)}
      slope_p_value <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 9]), na.rm=T)
      slope <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 10]), na.rm=T)
      slope_sd <- sd(as.numeric(simulated[simulated$shape_class == max_shape, 10]), na.rm=T)
      aic <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 12]))
      second_order_pvalue <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 14]), na.rm=T)
      nrmse <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 15]))
      # aic_lin <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 12]))
      # aic_pol <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 13]))
      # aic_max_shape <- mean(as.numeric(simulated[simulated$shape_class == max_shape, 14]))

    }else{alpha2 <- alpha1 <- sd_alpha1 <- inter <- strd <- p_1 <- sd_p_1 <- p_2 <- sd_p_2 <- p_3 <- sd_p_3 <- max_shape <- slope_p_value <- slope <- slope_sd <- aic <- best_model <- second_order_pvalue <- nrmse <- NA}

    ref_year <- simulated[1,11]

    res <- res %>% rbind(data.frame(alpha2, alpha1, sd_alpha1, inter, strd, p_1, sd_p_1, p_2, sd_p_2, p_3, sd_p_3, second_order_pvalue, slope_p_value, slope, slope_sd, nrmse, ref_year, max_shape=as.factor(max_shape), aic=aic, best_model=best_model, best_class=NA, trend=NA))

    res <- res %>%
      dplyr::mutate(
        best_class = ifelse(!is.na(max_shape),
                            dplyr::case_when(
                              stringr::str_detect(max_shape,"stable_constant") ~ "no_change",
                              stringr::str_detect(max_shape, "constant") ~ "linear",
                              stringr::str_detect(max_shape, "acc|dec|conc|conv") ~ "quadratic"
                            ), NA),
        trend = ifelse(!is.na(max_shape),
                       sub("_.*", "", max_shape), NA)
      )

  }


  return(res)
}



#' Summary of breakpoints analyses
#'
#' @param sets a list of data frame ready for analyses
#' @param abr_mtd a vector with abbreviation(s) corresponding to the breakpoints method(s) to use
#' ("cp", "cp2", "str_m", "str_t", "ecpt_m", "ecpt_t", "asd", "mcp", "chg")
#' @param asd_thr a numeric threshold for as_detect method
#' @param asd_chk a logical paramater for check_true_shift in asd_fct2
#'
#' @return a data frame with results of breakpoints analyses
#' @export

abrupt_classif <- function(sets, abr_mtd, asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst){

  ts <- sets %>%
    lapply(function(x) dplyr::select(x, c(X,Y)) %>% dplyr::rename(year=X))

  shifts_res <- ts %>%
    lapply(function(x) shifts(x, abr_mtd, asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst))

  abt_res <- list()
  for(k in abr_mtd){

    abt_res_mtd <- shifts_res %>%
      lapply(function(x) x$res_table %>%
               dplyr::filter(abbr==k)) %>%
      do.call(what=rbind.data.frame)
    abt_res[[k]] <- abt_res_mtd
  }


  return(list("shifts_res" = shifts_res, "abt_res" = abt_res))
}



#' Makes trajectory fitting
#'
#' @param sets a list of data frame ready for analyses
#' @param abr_mtd a list of abbreviation with breakpoints methods to use (see function abrupt_classif for the list)
#' @param type a character specifying whether the series come from simulations ("sim") or RAMLDB ("RAM")
#' @param asd_thr a numeric threshold for as_detect method
#' @param asd_chk a logical paramater for check_true_shift in asd_fct2
#' @param apriori logical to state whether expected class and trajectories are indicated
#'
#' @return a data frame with the results of different trajectory fitting,
#' @export

fit_models <- function(sets, abr_mtd, type, asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst, apriori){

  # Running the different model fitting:
  res_nch <- res_trend(sets, niter=1, correction=FALSE, fit="nch")
  res_lin <- res_trend(sets, niter=1, correction=FALSE, fit="lin")
  res_pol <- res_trend(sets, niter=1, correction=FALSE, fit="pol")
  res_abt <- abrupt_classif(sets, abr_mtd=abr_mtd, asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst)
  lengths <- data.frame(first = sapply(sets, function(x) min(x$X), simplify="array"),
                        last = sapply(sets, function(x) max(x$X), simplify="array"),
                        length = sapply(sets, function(x) nrow(x), simplify="array"))


  summ_res <- cbind(
    res_nch %>% dplyr::select(contains("max_shape")|contains("aic")|contains("trend")|contains("nrmse")) %>% dplyr::rename_with(~str_c(., "_nch")),
    res_lin %>% dplyr::select(contains("max_shape")|contains("aic")|contains("trend")|contains("nrmse")) %>% dplyr::rename_with(~str_c(., "_lin")),
    res_pol %>% dplyr::select(contains("max_shape")|contains("aic")|contains("best_model")|contains("trend")|contains("nrmse")) %>% dplyr::rename_with(~str_c(., "_pol")),
    do.call("cbind",
            lapply(abr_mtd, function(x) res_abt$abt_res[[x]] %>%
                     dplyr::select(contains("brk")|aic|trend|nrmse) %>%
                     dplyr::mutate(max_shape=paste0(n_brk,"_breakpoint")) %>%
                     dplyr::rename_with(~str_c(., paste0("_", x)))
                   )
            ),
    res_abt$abt_res[["chg"]] %>% dplyr::select(mag, SDbef, SDaft, step_size),
    lengths
  ) %>% dplyr::rename(signif_model = best_model_pol)

  # Add expectations if known
  if(type=="sim" & apriori){

    expected <- do.call(rbind, (lapply(sets, function(x) x[1,]))) %>%
      dplyr::select(expected_traj, expected_class)
    summ_res <- cbind(summ_res, expected)

  }

  return(list("summ_res"=summ_res,
              "res_fit"=list("res_nch"=res_nch,
                             "res_lin"=res_lin,
                             "res_pol"=res_pol,
                             "res_abt"=res_abt)
              )
         )
}



#' Comparison of trajectory fitting
#'
#' @param sets a list of data frame ready for analyses
#' @param str a character specifying the type of structure to use for classification ("aic","2_step","aic_asd")
#' @param abr_mtd a list of abbreviation with breakpoints methods to use (see function abrupt_classif for the list)
#' @param type a character specifying whether the series come from simulations ("sim") or RAMLDB ("RAM")
#' @param asd_thr a numeric threshold for as_detect method
#' @param asd_chk a logical paramater for check_true_shift in asd_fct2
#' @param noise_comb
#' @param showplots logical that indicate whether or not to show plots (slows down classification)
#' @param save_plot logical, whether to save plots or not
#' @param save_plot_bis logical, whether to save plot_bis or not
#' @param apriori logical to state whether expected class and trajectories are indicated
#' @param run_loo logical to state whether to perform leave-one-out process
#' @param two_bkps (experimental) logical to state if you look for more than one breakpoint (not fully tested)
#' @param smooth_signif logical, should significance of the coefficient be taken into account in the selection of best model?
#' @param outplot logical, last plot in return value?
#' @param ind_plot character to return the last plot for a given trajectory, either "nch", "lin", "pol", or "abt"
#' @param detection_plot logical to plot detection plot when abrupt
#' @param dirname directory name where to save plots
#'
#' @return a data frame with the results of different trajectory fitting,
#' a three-column data frame with the output of the best fitting model, trajectory, and class,
#' also plots the different fits
#' @export

traj_class <- function(sets, str, abr_mtd, asd_chk, type="sim", asd_thr,
                       noise_comb=NULL, showplots=FALSE, apriori=FALSE,
                       run_loo, two_bkps, smooth_signif, outplot=FALSE,
                       lowwl=5, highwl="default", mad_thr=3, mad_cst=1.4826,
                       ind_plot=NULL, save_plot=TRUE, edge_lim=5, congr_brk=5,
                       plot_one_in=10, detection_plot=TRUE, dirname=NULL,
                       save_plot_bis=TRUE){

  ### Original time series

  # Make trajectory fitting:
  res <- fit_models(sets$ts, abr_mtd, type, asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst, apriori)

  # Define best trajectory according to decision tree:
  best_traj <- best_traj_aic(class_res=res, type=type, apriori=apriori, aic_selec=str, smooth_signif=smooth_signif,
                             edge_lim=edge_lim, congr_brk=congr_brk)

  # If looking for two breakpoints:
  if(two_bkps){

    # Keep time series with breakpoints:
    list_brk_ts <- best_traj %>%
      dplyr::filter(traj=="1_breakpoint") %>%
      dplyr::select(simu_id, loc_brk_chg)

    # If some time series have breakpoints:
    if (nrow(list_brk_ts) != 0){

    # Split time series (keeping the breakpoint in each):
    sets2 <- list()
    for (i in 1:nrow(list_brk_ts)){

      sets2[[paste0(list_brk_ts$simu_id[i],"_1st")]] <- sets$ts[[list_brk_ts$simu_id[i]]] %>%
        dplyr::filter(X <= list_brk_ts$loc_brk_chg[i])

      sets2[[paste0(list_brk_ts$simu_id[i],"_2nd")]] <- sets$ts[[list_brk_ts$simu_id[i]]] %>%
        dplyr::filter(X >= list_brk_ts$loc_brk_chg[i])
      # %>% dplyr::mutate(X = 1:nrow(.))
    }

    # Remove too short sub time series:
    sets2 <- sets2[lapply(sets2, nrow)>=15]

    # If any time series eligible:
    if(length(sets2)>0){

      # Fit models of splitted time series:
      res2 <- fit_models(sets2, abr_mtd, type, asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst, apriori)
      best_traj_split <- best_traj_aic(class_res=res2, type=type, apriori=apriori, aic_selec=str, smooth_signif=smooth_signif,
                                       edge_lim=edge_lim, congr_brk=congr_brk)

      # Locate auxiliary breakpoints:
      list_brk_ts2 <- best_traj_split %>%
        dplyr::mutate(part = factor(gsub("^.*_", "", simu_id), levels=c("1st", "2nd")),
                      simu_id = gsub("_1st|_2nd", "", simu_id)) %>%
        dplyr::filter(traj == "1_breakpoint") %>%
        dplyr::select(simu_id, loc_brk_chg, part) %>%
        # tidyr::pivot_wider(names_from = part, values_from = loc_brk_chg) %>% # not able to disable drop levels
        tidyr::spread(key = part, value = loc_brk_chg, drop=FALSE) %>%
        dplyr::rename(loc_aux1_chg = "1st",
                      loc_aux2_chg = "2nd")

      best_traj <- dplyr::left_join(best_traj, list_brk_ts2, by="simu_id") %>%
        dplyr::mutate(nb_brk = dplyr::case_when(
          class != "abrupt" ~ 0,
          is.na(loc_aux1_chg) & is.na(loc_aux2_chg) ~ 1,
          is.na(loc_aux1_chg) | is.na(loc_aux2_chg) ~ 2,
          !(is.na(loc_aux1_chg) & is.na(loc_aux2_chg)) ~ 3
        ))
    } else {

      best_traj <- best_traj %>%
        dplyr::mutate(loc_aux1_chg = NA,
                      loc_aux2_chg = NA,
                      nb_brk = ifelse(class == "abrupt", 1, 0)
        )
    }

    } else {

      best_traj <- best_traj %>%
        dplyr::mutate(loc_aux1_chg = NA,
                      loc_aux2_chg = NA,
                      nb_brk = ifelse(class == "abrupt", 1, 0)
        )
    }
  }


  ### Leave-One-Out time series
  if (run_loo){

    # Make all leave-one-out time series and name them accordingly:
    loo <- lapply(sets$ts, function(x) lapply(1:nrow(x), function(k) x[-k,]))
    for (i in 1:length(loo)) names(loo[[i]]) <- paste0(names(loo)[i], "_", paste0("loo",1:length(loo[[i]])))

    # Make trajectory fitting for loo:
    res_loo <- lapply(loo, function(x) fit_models(x, abr_mtd, type, asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst, apriori))

    # Define best trajectory according to decision tree for loo:
    best_traj_loo <- lapply(res_loo, function(x) best_traj_aic(class_res=x, type=type, apriori=apriori, aic_selec=str, smooth_signif=smooth_signif,
                                                               edge_lim=edge_lim, congr_brk=congr_brk))

    best_traj_loo <- lapply(1:length(sets$ts), function(i) best_traj_loo[[i]] %>%
                              dplyr::mutate(X=sets$ts[[i]]$X))

    # LOO frequencies for class:
    loo_class <- do.call("rbind",  lapply(1:length(sets$ts), function(i) {

      class_freq <- best_traj_loo[[i]] %>%
        dplyr::group_by(class) %>%
        dplyr::summarise(class_freq_loo = n()/nrow(best_traj_loo[[i]])) %>%
        tidyr::pivot_wider(names_from = class, values_from = class_freq_loo)

      missing <- setdiff(levels(best_traj_loo[[i]]$class), names(class_freq))

      class_freq[missing] <- 0

      class_freq <- class_freq[levels(best_traj_loo[[i]]$class)] %>%
        dplyr::select_all(list(~ paste0("loo_", .))) %>%
        dplyr::mutate(simu_id = names(sets$ts)[i])


    }))

    best_traj <- best_traj %>%
      dplyr::left_join(loo_class, by="simu_id")

  }

  # store asd parameters
  if (str == "aic_asd") best_traj <- best_traj %>%
      dplyr::mutate(asd_thr=asd_thr,
                    lowwl=lowwl,
                    highwl=highwl,
                    mad_thr=mad_thr)

  # Plot the 4 different fits with best one colored:
  if(showplots == TRUE){

    if (!run_loo) best_traj_loo <- NULL

    ## Make plots
    plots_traj_nch <- plot_traj_multi_abt(sets, rslt=res$res_fit$res_nch, best_traj, plot_class="no_change", best_traj_loo=best_traj_loo)
    plots_traj_lin <- plot_traj_multi_abt(sets, rslt=res$res_fit$res_lin, best_traj, plot_class="linear", best_traj_loo=best_traj_loo)
    plots_traj_pol <- plot_traj_multi_abt(sets, rslt=res$res_fit$res_pol, best_traj, plot_class="quadratic", best_traj_loo=best_traj_loo)
    plots_traj_chg <- plots_traj_abt <- plot_traj_multi_abt(sets, rslt=res$res_fit$res_abt, best_traj, plot_class="abrupt", asd_thr, best_traj_loo=best_traj_loo, detection_plot=detection_plot)

    if(!is.null(ind_plot)){
      if(ind_plot=="best"){

        ind_plots <- get(paste0("plots_traj_",
                                sub("max_shape_","",best_traj$best_model)))
      } else {
        ind_plots <- get(paste0("plots_traj_",ind_plot))
      }
    }


    for (i in 1:nrow(res$summ_res)){


      if(str=="aic_asd") class_plot <- cowplot::plot_grid(plots_traj_nch[[i]], plots_traj_lin[[i]],
                                                          plots_traj_pol[[i]], plots_traj_abt$plots[[i]],
                                                          nrow=2, ncol=2, align = 'h')
      if(str=="aic") class_plot <- cowplot::plot_grid(plots_traj_nch[[i]], plots_traj_lin[[i]],
                                                      plots_traj_pol[[i]], plots_traj_abt[[i]],
                                                      nrow=2, ncol=2, align = 'h')

      if(str=="aic_asd") endname <- paste0("_",str,"_coef",smooth_signif,"asd_thr",asd_thr,"l",lowwl, "h", highwl,"mad", mad_thr,".png")
      if(str=="aic") endname <- paste0("_",str,"_coef",smooth_signif,".png")

      if(type=="sim" & save_plot & i%%plot_one_in==0){

        png(filename = paste0(dirname,
                              sets$ts[[i]]$expected_class[1],"/", paste(head(stringr::str_split(names(sets$ts)[i],"_")[[1]],-1), collapse="_"),"_",
                              noise_comb,"_", gsub("^.*_", "", names(sets$ts)[i]),"__", endname), width=8, height=6, units="in", res=300)
        print(class_plot)
        dev.off()


        if(str=="aic_asd" & save_plot_bis){

          png(filename = paste0(dirname,
                                sets$ts[[i]]$expected_class[1],"/abr_", paste(head(stringr::str_split(names(sets$ts)[i],"_")[[1]],-1), collapse="_"),"_",
                                noise_comb,"_", gsub("^.*_", "", names(sets$ts)[i]), endname), width=6, height=6, units="in", res=300)
          print(plots_traj_abt$plot_bis[[i]])
          dev.off()
        }

      }

      if (type=="RAM" & save_plot){

        class_plot <- class_plot +
          patchwork::plot_annotation(sub("_iter01","",names(sets$ts)[i]) %>% sub("_"," ",.),
                                     theme = theme(plot.title = element_text(hjust = 0.5))
                                     )

        png(filename = paste0(dirname,
                              sub("_iter01","",names(sets$ts)[i]), endname), width=8, height=6, units="in", res=300)
        print(class_plot)
        dev.off()


        if(str=="aic_asd" & save_plot_bis){

          png(filename = paste0(dirname,"/abr_",
                                sub("_iter01","",names(sets$ts)[i]), endname), width=6, height=6, units="in", res=300)
          print(plots_traj_abt$plot_bis[[i]])
          dev.off()
          }
        }

      if (type=="data" & save_plot){

        class_plot <- class_plot +
          patchwork::plot_annotation(sub("_iter01","",names(sets$ts)[i]) %>% sub("_"," ",.),
                                     theme = theme(plot.title = element_text(hjust = 0.5))
          )

        png(filename = paste0(dirname,sub("_iter01","",names(sets$ts)[i]), endname),
            width=8, height=6, units="in", res=300)
        print(class_plot)
        dev.off()

        if(str=="aic_asd" & save_plot_bis){

          png(filename = paste0(dirname,"/abr_",sub("_iter01","",names(sets$ts)[i]), endname),
              width=6, height=6, units="in", res=300)
          print(plots_traj_abt$plot_bis[[i]])
          dev.off()

          }
      }

    }

  }

  if (outplot){ # For last plot only

    if(!is.null(ind_plot)){
      return(list("res"=res$res, "best_traj"=best_traj,
                  "res_detail"=res$res_fit, "class_plot"=ind_plots))

    } else {
      return(list("res"=res$res, "best_traj"=best_traj,
                  "res_detail"=res$res_fit, "class_plot"=class_plot))
    }
  }

  return(list("res"=res$res, "best_traj"=best_traj,
              "res_detail"=res$res_fit))

}




#' Determine best trajectory based on AIC comparisons
#'
#' @param class_res a dataframe with results from the different classifications
#' @param type a character specifying whether the series come from simulations ("sim") or RAMLDB ("RAM")
#' @param apriori logical to state whether expected class and trajectories are indicated
#' @param aic_selec a character specifying whether the best trajectory is the lowest ("min") or consider simplest model with delta AIC < 2 ("2units")
#' @param smooth_signif logical, should significance of the coefficient be taken into account in the selection of best model?
#' @param edge_lim numeric, minimal breakpoint distance to start or end dates
#' @param congr_brk numeric, maximal acceptable distance between chngpt and as_detect breaks
#'
#'
#' @return a data frame with the output of the best fitting model, trajectory, and class
#' @export

best_traj_aic <- function(class_res, type, apriori, aic_selec, smooth_signif, edge_lim, congr_brk){

  traj_lvl <- c("stable_constant","increase_constant","decrease_constant",
                "stable_concave","stable_convex","decrease_accelerated",
                "decrease_decelerated","increase_accelerated",
                "increase_decelerated","1_breakpoint")

  class_lvl <- c("no_change","linear","quadratic","abrupt")

  best_traj <- class_res$summ_res %>%

    # Find minimal AIC and subtract to all AIC:
    tibble::rownames_to_column(var="simu_id") %>%
    dplyr::rowwise() %>% # necessary?

    dplyr::mutate(min_aic = min(dplyr::c_across(dplyr::contains("aic")),
                                na.rm=TRUE),
                  dplyr::across(dplyr::contains("aic") & !"min_aic", ~ .x - min_aic),
                  dplyr::across(.names = "weight_{.col}", .cols = dplyr::contains("aic") & !"min_aic", .fns = ~ exp(-0.5*.x)),
                  dplyr::across(.cols = dplyr::contains("weight_aic"), .fns = ~ round(.x / sum(dplyr::c_across(dplyr::contains("weight_aic")), na.rm=TRUE), digits=2))) %>%
    dplyr::ungroup()

  # Divide multiple asd breakpoints into several columns &
  # Compute distance between chg and asd breakpoints:
  if (aic_selec=="aic_asd"){

    best_traj <- best_traj %>%
      dplyr::mutate(loc_brk_asd_sep = strsplit(loc_brk_asd, ";")) %>%
      tidyr::unnest(loc_brk_asd_sep) %>%
      dplyr::mutate(loc_brk_asd_sep = as.numeric(loc_brk_asd_sep),
                    loc_brk_chg = as.numeric(loc_brk_chg)) %>%
      dplyr::group_by(simu_id) %>%
      dplyr::mutate(row = row_number()) %>%
      tidyr::pivot_wider(names_from=row, names_prefix="loc_brk_asd_", values_from=loc_brk_asd_sep) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::contains('loc_brk_asd_'), .names = "dup_{.col}")) %>%
      dplyr::rename_with(~sub('dup_loc_brk_asd_',"diff_loc_",.), dplyr::contains('dup_loc_brk_asd_')) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::contains("diff_loc_"), ~abs(loc_brk_chg-.x))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(loc_brk_asd = as.numeric(loc_brk_asd),
                    diff_loc_min = min(dplyr::c_across(dplyr::contains("diff_loc_")), na.rm=TRUE),
                    loc_brk_asd_min = min(dplyr::c_across(dplyr::contains("loc_brk_asd_")), na.rm=TRUE),
                    loc_brk_asd_max = max(dplyr::c_across(dplyr::contains("loc_brk_asd_")), na.rm=TRUE),
                    dplyr::across(.cols = dplyr::contains("loc"), ~ dplyr::na_if(., Inf)),
                    loc_brk_asd = as.character(loc_brk_asd)) %>%
    suppressWarnings()
  }

  # Expand in length to select by aic below:
  best_traj <- best_traj %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(contains("aic") & !contains("weight_aic"), names_to = "best_aic", values_to = "vals")

  ## for delta AIC < 2 keep model with lowest complexity according to the ranking:
  if (aic_selec=="2units"){
    best_traj <- best_traj %>%
      dplyr::filter(best_aic != "min_aic") %>%
      dplyr::mutate(best_aic = factor(best_aic, levels=c("aic_nch", "aic_lin", "aic_pol", "aic_chg"))) %>%
      dplyr::filter(vals < 2) %>%
      dplyr::group_by(simu_id) %>%
      dplyr::mutate(n = n()) %>%
      dplyr::filter(n==1 | rank(best_aic) == 1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(best_aic = sub("aic_","max_shape_",best_aic)) %>%
      dplyr::select(-c(n,vals))
  }

  ## [experimental] keep only model with lowest AIC regardless of delta AIC:
  if (aic_selec=="aic"){

    # order models according to their AIC:
    rank_aic <- best_traj %>%
      dplyr::group_by(simu_id) %>%
      dplyr::filter(!best_aic %in% c("aic_asd", "min_aic")) %>%
      dplyr::mutate(aic_rank = dplyr::dense_rank(vals)) %>%
      dplyr::select(contains("aic") & !contains("weight_aic") | "simu_id") %>%
      tidyr::pivot_wider(names_from = aic_rank, names_prefix = "aic_rank", values_from = best_aic)

    best_traj <- best_traj %>%
      dplyr::filter(vals == 0) %>%
      dplyr::select(-vals) %>%
      dplyr::left_join(rank_aic, by="simu_id") %>%

      dplyr::mutate(
        # If the step size is too low:
        best_aic = ifelse(best_aic=="aic_chg" & abs(step_size) <1,
                          aic_rank2, best_aic)) %>%

      dplyr::mutate(best_aic = sub("aic_","max_shape_", best_aic))
  }


  ## [original] keep only model with lowest AIC regardless of delta AIC:
  # if (aic_selec=="min"){
  #
  #   best_traj <- best_traj %>%
  #     dplyr::filter(vals == 0) %>%
  #     dplyr::mutate(best_aic = sub("aic_","max_shape_",best_aic)) %>%
  #     dplyr::select(-vals)
  # }


  ## lowest AIC but confirm with asdetect if abrupt:
  if (aic_selec=="aic_asd"){

    # order models according to their AIC:
    rank_aic <- best_traj %>%
      dplyr::group_by(simu_id) %>%
      dplyr::filter(!best_aic %in% c("aic_asd", "min_aic")) %>%
      dplyr::mutate(aic_rank = dplyr::dense_rank(vals)) %>%
      dplyr::select(contains("aic") & !contains("weight_aic") | "simu_id") %>%
      tidyr::pivot_wider(names_from = aic_rank, names_prefix = "aic_rank", values_from = best_aic)

    # keep the model with the lowest AIC:
    best_traj <- best_traj %>%
      dplyr::filter(vals == 0) %>%
      dplyr::select(-vals)

    # if best AIC model is abrupt but do not fill in the following criteria
    # take 2nd best model in the ranking:
    best_traj <- best_traj %>%
      dplyr::left_join(rank_aic, by="simu_id") %>%

      dplyr::mutate(

        not_abr = "not_lowest_aic",

        # If asdetect doesn't find any breakpoints:
        not_abr = ifelse(best_aic=="aic_chg" &
                           max_shape_asd=="0_breakpoint",
                         "no_asd_brk", not_abr),

        best_aic = ifelse(best_aic=="aic_chg" &
                            max_shape_asd=="0_breakpoint",
                          aic_rank2, best_aic),

        # If the breakpoints found are not congruent:
        # best_aic = ifelse(best_aic=="aic_chg" &
        #                     diff_loc > 3,
        #                   aic_rank2, best_aic),
        not_abr = ifelse(best_aic=="aic_chg" & diff_loc_min > congr_brk,
                         "no_congr_brk", not_abr),

        best_aic = ifelse(best_aic=="aic_chg" & diff_loc_min > congr_brk,
                                        aic_rank2, best_aic),

        # If the step size is too low:
        # best_aic = ifelse(best_aic=="aic_chg" & abs(step_size) <1,
        #                   aic_rank2, best_aic),

        # If the breakpoints are too close to the edges:
        # only for chngpt
        not_abr = ifelse(best_aic=="aic_chg" &
                           (loc_brk_chg-first<edge_lim |
                              last-loc_brk_chg<edge_lim),
                         "edge_lim", not_abr),

        best_aic = ifelse(best_aic=="aic_chg" &
                            (loc_brk_chg-first<edge_lim |
                               last-loc_brk_chg<edge_lim),
                          aic_rank2, best_aic),

        # best_aic = ifelse(best_aic=="aic_chg" &
        #                     ((loc_brk_chg-first<edge_lim & loc_brk_asd_min-first<edge_lim) |
        #                        (last-loc_brk_chg<edge_lim & last-loc_brk_asd_max<edge_lim)),
        #                   aic_rank2, best_aic)


        not_abr = ifelse(best_aic=="aic_chg", NA, not_abr)

      ) %>%

      dplyr::mutate(best_aic = sub("aic_","max_shape_", best_aic))

  }


  # Display best model and trajectory for each simulation
  best_traj <- best_traj %>%

    dplyr::select(simu_id|contains("max_shape")|contains("expected")|
                  best_aic|signif_model|contains("loc")|contains("weight")|
                  mag|contains("SD")|contains("nrmse")|step_size|contains("not_abr")) %>%

    tidyr::pivot_longer(contains("max_shape"), names_to = "best_model", values_to = "traj") %>%

    # Keep the trends for each fit (and then only keep the best):
    dplyr::left_join(best_traj %>%
                       # dplyr::select(simu_id | (contains("trend") & !"trend_asd")) %>%
                       dplyr::select(simu_id | contains("trend")) %>%
                       {if(aic_selec == "aic_asd") dplyr::select(., -"trend_asd") else .} %>%
                       tidyr::pivot_longer(contains("trend"), names_to = "trend_model", values_to = "trend") %>%
                       dplyr::rename(best_model = trend_model) %>%
                       dplyr::mutate(best_model = sub("trend_","max_shape_",best_model)),
                     by = c("simu_id", "best_model"))


  # Check significance of coefficients for smooth (non-abrupt) trajectories:
  if(smooth_signif){

    best_traj <- best_traj %>%
      dplyr::mutate(signif_model = paste0("max_shape_", signif_model),
                    best_aic = ifelse(best_aic != "max_shape_chg" & signif_model != best_aic,
                                      signif_model, best_aic))
  }

  # Final reshaping of the results:
  best_traj <- best_traj %>%
    dplyr::filter(best_aic==best_model) %>%
    {if(aic_selec == "aic_asd") dplyr::select(., -best_aic & -weight_aic_asd) else dplyr::select(., -best_aic)} %>%
    dplyr::mutate(simu_id = factor(simu_id),
                  traj = traj %>% factor(levels = traj_lvl),
                  class = dplyr::case_when(
                    stringr::str_detect(traj,"stable_constant") ~ "no_change", # acts as smooth significant for 1st order
                    stringr::str_detect(traj, "constant") ~ "linear",
                    stringr::str_detect(traj, "breakpoint") ~ "abrupt",
                    stringr::str_detect(traj, "acc|dec|conc|conv") ~ "quadratic"
                  ) %>% factor(levels = class_lvl),
                  trend = trend %>% factor(levels = c("stable", "decrease", "increase"))) %>%
    dplyr::rename(weight_aic_no_change = weight_aic_nch,
                  weight_aic_linear = weight_aic_lin,
                  weight_aic_quadratic = weight_aic_pol,
                  weight_aic_abrupt = weight_aic_chg,
                  nrmse_no_change = nrmse_nch,
                  nrmse_linear = nrmse_lin,
                  nrmse_quadratic = nrmse_pol,
                  nrmse_abrupt = nrmse_chg) %>%
    dplyr::relocate(where(is.factor)) %>%
    dplyr::relocate(c(simu_id, best_model))

  # Make expected trajectories factors:
  if(type=="sim" & apriori){
    best_traj <- best_traj %>%
      dplyr::mutate(expected_traj = expected_traj %>% factor(levels = traj_lvl),
                    expected_class = expected_class %>% factor(levels = class_lvl))
  }

  return(best_traj)

}



# Wrapper functions -------------------------------------------------------

#' Run classification for all time series of a given type
#'
#' @param ts_type a type of empirical time series (TCbest, TBbest, SProd...)
#' @param min_len integer, the minimal length of time series to include (default 20 time points)
#' @param str
#' @param normalize logical, whether or not to normalize by setting maximal value to one
#' @param showplots logical that indicate whether or not to show plots (slows down classification)
#' @param run_loo logical to state whether to perform leave-one-out process
#'
#' @return a data frame with the results of different trajectory fitting,
#' a three-column data frame with the output of the best fitting model, trajectory, and class,
#' also plots the different fits
#' @export

run_classif_RAM <- function(ts_type, min_len=20, str, asd_thr, normalize, showplots=TRUE, run_loo, two_bkps, smooth_signif){

  # List of stocks with time series available for classification:
  # stocks_avail <- list_available(ts_type=ts_type, min_len=min_len)
  stocks_avail <- list_available(ts_type=ts_type, min_len=min_len)

  # Classify for all time series of this type:
  traj_ts <- data.frame()
  outlist <- list()

  thr <- NULL
  if(ts_type %in% c("TCbest","TBbest","R")) thr <- 0

  for (i in 1:nrow(stocks_avail)){

    set <- extract_RAM(stocks_avail$stockid[i], ts_type) %>%
      prep_data(thr=thr, type="RAM", apriori=FALSE)

    # scale by average
    if(ts_type == "SProd"){

      set_TB <- extract_RAM(stocks_avail$stockid[i], "TBbest") %>%
        prep_data(thr=thr, type="RAM", apriori=FALSE)

      mean_TB <- set_TB$ts[[1]]$Y %>% mean()

      set$ts[[1]] <- set$ts[[1]] %>%
        dplyr::mutate(Y_SE = Y_SE/mean_TB, # Y_SE not used
                      Y = Y/mean_TB)

    } else {

      set$ts[[1]] <- set$ts[[1]] %>%
        dplyr::mutate(Y_SE = Y_SE/mean(Y), # Y_SE not used
                      Y = Y/mean(Y))
    }

    set_length <- nrow(set$ts[[1]])

    if (normalize){

      # Normalization
      # set$ts[[1]] <- set$ts[[1]] %>%
      #   dplyr::mutate(Y_SE = Y_SE/max(Y)*100,
      #                 Y = Y/max(Y)*100)

      # Standardization
      set$ts[[1]] <- set$ts[[1]] %>%
        dplyr::mutate(Y_SE = Y_SE/mean(Y),
                      Y = (Y-mean(Y))/sd(Y))
    }

    if(str == "aic") abr_mtd <- c("chg")
    if(str == "aic_asd") abr_mtd <- c("chg", "asd")

    trajs <- traj_class(sets=set, str=str, abr_mtd=abr_mtd, asd_thr=asd_thr, asd_chk=FALSE,
                        type="RAM", showplots=showplots, apriori=FALSE, run_loo=run_loo,
                        two_bkps=two_bkps, smooth_signif=smooth_signif,
                        lowwl=5, highwl="default", mad_thr=3,
                        edge_lim=0, congr_brk=5)


    # Add linear slope, breakpoints location, magnitude:
    trajs$best_traj <- trajs$best_traj %>%
      dplyr::mutate(slope = trajs$res_detail$res_lin$slope
                    # ,
                    # length_diff = set_length-stocks_avail[i,]$length
                    # loc_brk_chg = ifelse(trajs$best_traj$class=="abrupt",
                    #                      trajs$res_detail$res_abt$abt_res$chg$loc_brk, NA),
                    # loc_brk_asd = ifelse(trajs$best_traj$class=="abrupt" & str == "aic_asd",
                    #                      trajs$res_detail$res_abt$abt_res$asd$loc_brk, NA),
                    # diff_loc = ifelse(is.na(loc_brk_chg)|is.na(loc_brk_asd), NA, diff_loc),
                    # mag_brk_chg = ifelse(trajs$best_traj$class=="abrupt",
                    #                      trajs$res_detail$res_abt$abt_res$chg$mag, NA),
                    # stepsize_brk_chg = ifelse(trajs$best_traj$class=="abrupt",
                    #                           trajs$res_detail$res_abt$abt_res$chg$step_size, NA)
      )

    traj_ts <- traj_ts %>% dplyr::bind_rows(trajs$best_traj)
    outlist[[stocks_avail$stockid[i]]] <- trajs

    if (i%%10 == 0) print(paste0(i,"/",nrow(stocks_avail)))

  }

  traj_ts_full <- traj_ts %>%
    dplyr::mutate(stockid = simu_id %>%
                    sub("_iter01","", .) %>%
                    sub(paste0(ts_type,"_"),"", .),
                  ts_type = sub("_.*", "", simu_id)) %>%
    dplyr::left_join(stocks_avail, by="stockid")

  saveRDS(outlist, paste0("res/04_RAM_class/outlist_",ts_type,"_",str,"_coeff",smooth_signif,".rds"))
  saveRDS(traj_ts_full, paste0("res/04_RAM_class/traj_",ts_type,"_",str,"_coeff",smooth_signif,".rds"))

  return(return(list("traj_ts_full"=traj_ts_full, "outlist"=outlist)))

}


#' Run classification for all time series of a given type cutting at start and end
#'
#' @param ts_type a type of empirical time series (TCbest, TBbest, SProd...)
#' @param min_len integer, the minimal length of time series to include (default 20 time points)
#' @param str
#' @param normalize logical, whether or not to normalize by setting maximal value to one
#' @param showplots logical that indicate whether or not to show plots (slows down classification)
#' @param run_loo logical to state whether to perform leave-one-out process
#' @param min_y start year
#' @param max_y end year
#'
#' @return a data frame with the results of different trajectory fitting,
#' a three-column data frame with the output of the best fitting model, trajectory, and class,
#' also plots the different fits
#' @export

run_classif_RAM_lim <- function(ts_type, min_len=20, str, asd_thr, normalize, showplots=TRUE, run_loo, two_bkps, smooth_signif, min_y, max_y){

  # List of stocks with time series available for classification:
  # stocks_avail <- list_available(ts_type=ts_type, min_len=min_len)
  stocks_avail <- list_available_lim(ts_type=ts_type, min_len=min_len, min_y=min_y, max_y=max_y)

  # Classify for all time series of this type:
  traj_ts <- data.frame()
  outlist <- list()

  thr <- NULL
  if(ts_type %in% c("TCbest","TBbest")) thr <- 0

  for (i in 1:nrow(stocks_avail)){

    set <- extract_RAM_lim(stocks_avail$stockid[i], ts_type, min_y=min_y, max_y=max_y) %>%
      prep_data(thr=thr, type="RAM", apriori=FALSE)

    # scale by average
    if(ts_type == "SProd"){

      set_TB <- extract_RAM_lim(stocks_avail$stockid[i], "TBbest", min_y=min_y, max_y=max_y) %>%
        prep_data(thr=thr, type="RAM", apriori=FALSE)

      mean_TB <- set_TB$ts[[1]]$Y %>% mean()

      set$ts[[1]] <- set$ts[[1]] %>%
        dplyr::mutate(Y_SE = Y_SE/mean_TB, # Y_SE not used
                      Y = Y/mean_TB)

    } else {

      set$ts[[1]] <- set$ts[[1]] %>%
        dplyr::mutate(Y_SE = Y_SE/mean(Y), # Y_SE not used
                      Y = Y/mean(Y))
    }

    set_length <- nrow(set$ts[[1]])

    if (normalize){

      # Normalization
      # set$ts[[1]] <- set$ts[[1]] %>%
      #   dplyr::mutate(Y_SE = Y_SE/max(Y)*100,
      #                 Y = Y/max(Y)*100)

      # Standardization
      set$ts[[1]] <- set$ts[[1]] %>%
        dplyr::mutate(Y_SE = Y_SE/mean(Y),
                      Y = (Y-mean(Y))/sd(Y))
    }

    if(str == "aic") abr_mtd <- c("chg")
    if(str == "aic_asd") abr_mtd <- c("chg", "asd")

    trajs <- traj_class(sets=set, str=str, abr_mtd=abr_mtd, asd_thr=asd_thr, asd_chk=FALSE,
                        type="RAM", showplots=showplots, apriori=FALSE, run_loo=run_loo,
                        two_bkps=two_bkps, smooth_signif=smooth_signif,
                        lowwl=5, highwl="default", mad_thr=3,
                        edge_lim=0, congr_brk=5)


    # Add linear slope, breakpoints location, magnitude:
    trajs$best_traj <- trajs$best_traj %>%
      dplyr::mutate(slope = trajs$res_detail$res_lin$slope
                    # ,
                    # length_diff = set_length-stocks_avail[i,]$length
                    # loc_brk_chg = ifelse(trajs$best_traj$class=="abrupt",
                    #                      trajs$res_detail$res_abt$abt_res$chg$loc_brk, NA),
                    # loc_brk_asd = ifelse(trajs$best_traj$class=="abrupt" & str == "aic_asd",
                    #                      trajs$res_detail$res_abt$abt_res$asd$loc_brk, NA),
                    # diff_loc = ifelse(is.na(loc_brk_chg)|is.na(loc_brk_asd), NA, diff_loc),
                    # mag_brk_chg = ifelse(trajs$best_traj$class=="abrupt",
                    #                      trajs$res_detail$res_abt$abt_res$chg$mag, NA),
                    # stepsize_brk_chg = ifelse(trajs$best_traj$class=="abrupt",
                    #                           trajs$res_detail$res_abt$abt_res$chg$step_size, NA)
      )

    traj_ts <- traj_ts %>% dplyr::bind_rows(trajs$best_traj)
    outlist[[stocks_avail$stockid[i]]] <- trajs

    if (i%%10 == 0) print(paste0(i,"/",nrow(stocks_avail)))

  }

  traj_ts_full <- traj_ts %>%
    dplyr::mutate(stockid = simu_id %>%
                    sub("_iter01","", .) %>%
                    sub(paste0(ts_type,"_"),"", .),
                  ts_type = sub("_.*", "", simu_id)) %>%
    dplyr::left_join(stocks_avail, by="stockid")

  saveRDS(outlist, paste0("res/04_RAM_class/outlist_",ts_type,"_",str,"_coeff",smooth_signif,"_",min_y,"-",max_y,".rds"))
  saveRDS(traj_ts_full, paste0("res/04_RAM_class/traj_",ts_type,"_",str,"_coeff",smooth_signif,"_",min_y,"-",max_y,".rds"))

  return(return(list("traj_ts_full"=traj_ts_full, "outlist"=outlist)))

}




run_classif_data <- function(df_list, min_len=20, str, normalize, showplots=TRUE, run_loo, two_bkps, smooth_signif,
                             group, time, variable, outplot=FALSE, ind_plot=NULL, dirname=NULL, save_plot=TRUE){

  # List of time series meeting the length criterion:
  df_list <- df_list[df_list %>% lapply(nrow)>min_len]

  # Classify for all time series of this type:
  traj_ts <- data.frame()
  outlist <- list()

  thr <- NULL
  # if(ts_type %in% c("TCbest","TBbest")) thr <- 0

  for (i in 1:length(df_list)){

    set <- df_list[[i]] %>%
      dplyr::rename(scen = .data[[group]],
                    year = .data[[time]]) %>%
      dplyr::select(scen, year, .data[[variable]]) %>%
      prep_data(thr=thr, type="RAM", apriori=FALSE)

    set_length <- nrow(set$ts[[1]])

    if (normalize){

      # Normalization
      # set$ts[[1]] <- set$ts[[1]] %>%
      #   dplyr::mutate(Y_SE = Y_SE/max(Y)*100,
      #                 Y = Y/max(Y)*100)

      # Standardization
      set$ts[[1]] <- set$ts[[1]] %>%
        dplyr::mutate(Y_SE = Y_SE/mean(Y),
                      Y = (Y-mean(Y))/sd(Y))
    }

    if(str == "aic") abr_mtd <- c("chg")
    if(str == "aic_asd") abr_mtd <- c("chg", "asd")

    trajs <- traj_class(sets=set, str=str, abr_mtd=abr_mtd, asd_thr=0.15, asd_chk=FALSE,
                        type="data", showplots=showplots, apriori=FALSE, run_loo=run_loo,
                        two_bkps=two_bkps, smooth_signif=smooth_signif, outplot=outplot,
                        ind_plot=ind_plot, lowwl=5, highwl="default", mad_thr=3,
                        edge_lim=5, congr_brk=5, dirname=dirname, save_plot=save_plot)

    # Add linear slope, breakpoints location, magnitude:
    trajs$best_traj <- trajs$best_traj %>%
      dplyr::mutate(slope = trajs$res_detail$res_lin$slope)

    traj_ts <- traj_ts %>% dplyr::bind_rows(trajs$best_traj)
    outlist[[names(df_list)[i]]] <- trajs

    if (i%%10 == 0) print(paste0(i,"/",length(df_list)))
  }

  traj_ts_full <- traj_ts %>%
    dplyr::mutate(species = simu_id %>%
                    sub("_iter01","", .))

  # saveRDS(outlist, paste0(dirname, "outlist_",str,"_coeff",smooth_signif,"_data.rds"))
  # saveRDS(traj_ts_full, paste0(dirname, "traj_",str,"_coeff",smooth_signif,"_data.rds"))

  return(list("traj_ts_full"=traj_ts_full, "outlist"=outlist))

}




#' Compute one column of the confusion matrix
#'
#' @param best_traj a data frame with the output of the best fitting model, trajectory, class, LOO proportion
#' @param expected a character specifying the reference trajectory or class used to generate the simulations
#' @param type a character specifying the type of classification ("traj" or "class")
#'
#' @return a single column dataframe corresponding to a column of the confusion matrix
#' @export

make_conf_mat <- function(best_traj){

  conf_mat_traj <- caret::confusionMatrix(data=best_traj$traj, reference=best_traj$expected_traj)
  conf_mat_class <- caret::confusionMatrix(data=best_traj$class, reference=best_traj$expected_class)

  return(list("conf_mat_traj" = conf_mat_traj, "conf_mat_class"=conf_mat_class))
}





# [Server] Run classifications in parallel
# i: integer indicating noise combination

classif_noise_comb <- function(simu_list, str, run_loo, asd_thr, i){

  print(paste0("noise combination: ",i,"/",length(simu_list)))
  noise_comb <- names(simu_list)[i]

  sets <- prep_data(df=simu_list[[i]], thr=NULL, type="sim", apriori=TRUE)
  # sets <- c(sets[1])

  if(str == "aic") abr_mtd <- c("chg")
  if(str == "aic_asd") abr_mtd <- c("chg", "asd")

  ## workflow 3: AIC + asdetect:
  trajs <- traj_class(sets, str=str, abr_mtd=abr_mtd, asd_thr=asd_thr, asd_chk=FALSE,
                      type="sim", noise_comb=noise_comb, smooth_signif=TRUE,
                      showplots=FALSE, apriori=TRUE, run_loo=run_loo, two_bkps = TRUE,
                      lowwl=5, highwl="default", mad_thr=3, mad_cst=1.4826,
                      edge_lim=5, congr_brk=5, plot_one_in=10)

  # traj_list[[noise_comb]] <- trajs

  conf_mat <- make_conf_mat(trajs$best_traj)
  # mat_list[[noise_comb]] <- conf_mat


  return(list("trajs" = trajs,
              "conf_mat" = conf_mat))
}



# Make plots --------------------------------------------------------------


#' Summary of MC simulations of classification
#'
#' @param sets a list of data frame ready for analyses
#' @param rslt data frame returned by res_trend function
#' @param plot_class class of trajectory corresponding the rslt
#' @param best_traj a data frame with the output of the best fitting model, trajectory, class, LOO proportion
#' @param best_traj_loo output from LOO process
#' @param detection_plot logical to plot the detection plot
#'
#' @return a list of plots with best fit
#' @export

plot_traj_multi_abt <- function(sets, rslt, best_traj, plot_class, asd_thr=NULL, best_traj_loo=NULL, detection_plot=TRUE){

  if (plot_class != "abrupt") plots <- vector(mode = "list", length = nrow(rslt))
  if (plot_class == "abrupt") plots <- plot_bis <- vector(mode = "list", length = length(rslt$shifts_res))

  for (i in 1:length(sets$ts)){

    # Define y-axis name:
    if (sets$ts_type %in% c("TB", "TBbest")) ts_type <- "Biomass"
    # if (sets$ts_type %in% c("TB", "TBbest")) ts_type <- "State variable (Biomass)"
    if (sets$ts_type %in% c("TC", "TCbest")) ts_type <- "Catch"
    if (sets$ts_type %in% c("SProd")) ts_type <- "Surplus production"
    if (sets$ts_type %in% c("Index", "index")) ts_type <- "Index"
    if (sets$ts_type %in% c("R")) ts_type <- "Recruitment"

    # Plot time series and model fit [smooth]:

    if (plot_class != "abrupt"){

      p <- local({
        i <- i
        ggplot(sets$ts[[i]], aes(x = X, y = Y))+
          geom_line()+
          theme_light(base_size = 7)+
          labs(x = "time unit (year)", y = ts_type)+
          # expand_limits(y=c(0,10))+
          expand_limits(y=0)+
          stat_function(fun=function(x){rslt[i,]$alpha2*x^2+rslt[i,]$alpha1*x+rslt[i,]$inter}, color="blue")+
          stat_function(fun=function(x){rslt[i,]$alpha2*x^2+rslt[i,]$alpha1*x+rslt[i,]$inter-rslt[i,]$strd}, linetype = "dashed", color="blue")+
          stat_function(fun=function(x){rslt[i,]$alpha2*x^2+rslt[i,]$alpha1*x+rslt[i,]$inter+rslt[i,]$strd}, linetype = "dashed", color="blue")
      })

      # Plot time series and model fit [abrupt]:
    } else {

      table_chg <- rslt$abt_res$chg %>% dplyr::slice(i) %>%
        dplyr::mutate(loc_brk = as.numeric(loc_brk))

      # as_detect score
      if(!is.null(rslt$abt_res$asd)){

        asd_ts <- rslt$shifts_res[[i]]$asd_detect %>%
          tibble::as_tibble() %>%
          tibble::rownames_to_column(var="year") %>%
          dplyr::mutate(year=as.numeric(year)+sets$ts[[i]]$X[1]-1)

        p_asd <- asd_ts %>%
          ggplot(aes(x=year, y=value))+
          geom_hline(yintercept = c(-asd_thr, asd_thr), col="red", linetype="dotted")+
          geom_line()+
          theme_light(base_size = 7)+
          labs(y="Detection")+
          theme(plot.title = element_text(hjust = 0.5))+
          expand_limits(y=c(-1,1))+
          ggtitle("as_detect detection score")

        p_loc <- rslt$abt_res$asd %>% dplyr::slice(i) %>%
          dplyr::pull(loc_brk) %>% strsplit(";") %>% `[[`(1) %>%
          as.numeric() %>% suppressWarnings()

        if (!is.na(p_loc[1])) p_asd <- p_asd +
          geom_vline(xintercept = p_loc, col="red", linetype="dashed")

        if(rslt$abt_res$asd[i,"loc_brk"]!="NA"){

          asd_loc <- rslt$abt_res$asd %>% dplyr::slice(i) %>%
            dplyr::pull(loc_brk) %>% strsplit(";") %>% `[[`(1) %>%
            as.numeric()

        } else {asd_loc <- NA}

      } else {
        # table_asd <- table_chg ; table_asd[1,] <- NA
        asd_loc <- NA
      }

      p <- local({
        i <- i
        ggplot(sets$ts[[i]], aes(x = X, y = Y))+
          geom_line()+
          theme_light(base_size = 7)+
          labs(x = "time unit (year)", y = ts_type)+
          # expand_limits(y=c(0,10))+
          expand_limits(y=0)+

          geom_vline(xintercept = table_chg$loc_brk, col="blue", linetype="dashed")+
          geom_line(data = rslt$shifts_res[[i]]$chg_outlist$pred_chg,
                    aes(x=year, y=bp), col = "blue", alpha=0.7)+
          scale_colour_manual(values = rep("red", table_chg$n_brk+1))+
          theme(legend.position = "none") %>%
          suppressWarnings()
      })

      # Plot additional breaktimes [chg]
      if ("loc_aux1_chg" %in% names(best_traj)) {

        loc_aux_chg <- c(best_traj[["loc_aux1_chg"]][i], best_traj[["loc_aux2_chg"]][i])
        p <- p +
          geom_vline(xintercept = loc_aux_chg[!is.na(loc_aux_chg)],
                     col="dodgerblue4", linetype="dotted", alpha=0.5)
      }

      # Plot as_detect breaktimes [asd]
      if (!is.na(asd_loc[1])){

        p <- p +
          geom_vline(xintercept = asd_loc, col="red", linetype="dashed")
      }
    }

    # Define title:
    if(plot_class == "no_change"){
      title_part <- paste0("<br>Intercept = ", format(rslt[i,]$inter, digits=2, scientific=TRUE))
    }

    if(plot_class == "linear"){
      if(rslt[i,]$slope_p_value<0.001){
        pval <- " p_val<0.001"
      } else {
        pval <- paste0(" p_val = ", round(rslt[i,]$slope_p_value, digits=3))
      }
      title_part <- paste0("<br>Slope = ", format(rslt[i,]$alpha1, digits=2, scientific=TRUE), pval)
    }

    if(plot_class == "quadratic"){
      if(rslt[i,]$second_order_pvalue<0.001){
        pval <- " p_val<0.001"
      } else {
        pval <- paste0(" p_val = ", round(rslt[i,]$second_order_pvalue, digits=3))
      }
      title_part <- paste0("<br>2nd_order_coeff = ", format(rslt[i,]$alpha2, digits=2, scientific=TRUE), pval)
    }

    if(plot_class == "abrupt"){
      title_part <- paste0("<br>Breaktime(s): <span style='color:blue'>", table_chg$loc_brk, "</span>",
                          if(!(is.na(best_traj$loc_aux1_chg[i]) & is.na(best_traj$loc_aux2_chg[i]))) {
                            paste0(" (<span style='color:dodgerblue4'>",best_traj$loc_aux1_chg[i],"</span>",
                                     ", <span style='color:dodgerblue4'>",best_traj$loc_aux2_chg[i],"</span>)")
                                   # if(!is.na(best_traj$loc_aux2_chg[i])){
                                   #   paste0(", <span style='color:dodgerblue2'>",best_traj$loc_aux2_chg[i],"</span>")
                                   # },")")
                          },
                          if(!is.na(asd_loc[1])) {paste0("; <span style='color:red'>",paste(asd_loc,collapse=","),"</span>")},
                          "  Step size = ", round(table_chg$step_size, digits=2))
      }

    # Define title [LOO]:
    if (!is.null(best_traj_loo)){

      # Add title and LOO points [LOO]:
      if (plot_class != "abrupt"){

        p <- p +
          ggtitle(paste0("<b>",stringr::str_to_title(sub("_"," ",rslt[i,]$best_class))," ", rslt[i,]$trend,"</b>",
                        "<br>AICc = ", round(rslt[i,]$aic, digits=2),
                        "  wAICc = ", best_traj[[paste0("weight_aic_", plot_class)]][i],
                        "  LOO = ", round(best_traj[[paste0("loo_", plot_class)]][i], digits=2),
                        "  NRMSE = ", round(best_traj[[paste0("nrmse_", plot_class)]][i], digits=2),
                        title_part)) +
          theme(plot.title = ggtext::element_markdown(size=10, lineheight = 1.1))+
          geom_point(data=best_traj_loo[[i]] %>%
                       dplyr::left_join(sets$ts[[i]], by="X") %>%
                       dplyr::filter(class==plot_class), aes(x=X, y=Y), col="orange", alpha=0.7)
      } else {

        p <- p +
          ggtitle(paste0("<b>Abrupt ", rslt$abt_res$chg[i,]$trend,"</b>",
                        "<br>AICc = ", round(table_chg$aic, digits=2),
                        "  wAICc = ", best_traj[[paste0("weight_aic_", plot_class)]][i],
                        "  LOO = ", round(best_traj[[paste0("loo_", plot_class)]][i], digits=2),
                        "  NRMSE = ", round(best_traj[[paste0("nrmse_", plot_class)]][i], digits=2),
                        title_part)) +
          theme(plot.title = ggtext::element_markdown(size=10, lineheight = 1.1))+

          # Add histograms for LOO breaktimes:
          geom_histogram(data=best_traj_loo[[i]] %>%
                           dplyr::filter(class=="abrupt"),
                         aes(x=loc_brk_chg, y=after_stat(ncount)*max(sets$ts[[i]][[names(sets$ts[[i]])[2]]])/2),
                         fill="blue", alpha=.3, colour=NA, binwidth = 1)

          if(!is.null(rslt$abt_res$asd)){
            p <- p +
              geom_histogram(data=best_traj_loo[[i]] %>%
                               dplyr::filter(class=="abrupt") %>%
                               dplyr::mutate(loc_brk_asd = strsplit(loc_brk_asd, ";")) %>%
                               tidyr::unnest(loc_brk_asd) %>%
                               dplyr::mutate(loc_brk_asd = as.numeric(loc_brk_asd)),
                             aes(x=loc_brk_asd, y=after_stat(ncount)*max(sets$ts[[i]][[names(sets$ts[[i]])[2]]])/2),
                             fill="red", alpha=.3, colour=NA, binwidth = 1)
          }

        p <- p +
          geom_point(data=best_traj_loo[[i]] %>%
                       dplyr::left_join(sets$ts[[i]], by="X") %>%
                       dplyr::filter(class==plot_class), aes(x=X, y=Y), col="orange", alpha=0.7)

      }


    } else {
      # Add title [no LOO]:
      if (plot_class != "abrupt"){

      p <- p +
        ggtitle(paste0("<b>",stringr::str_to_title(sub("_"," ",rslt[i,]$best_class))," ", rslt[i,]$trend,"</b>",
                      "<br>AICc = ", round(rslt[i,]$aic, digits=2),
                      "  wAICc = ", best_traj[[paste0("weight_aic_", plot_class)]][i],
                      "  NRMSE = ", round(best_traj[[paste0("nrmse_", plot_class)]][i], digits=2),
                      title_part))+
        theme(plot.title = ggtext::element_markdown(size=10, lineheight = 1.1))

      } else {

        p <- p +
          ggtitle(paste0("<b>Abrupt ", rslt$abt_res$chg[i,]$trend,"</b>",
                        "<br>AICc = ", round(table_chg$aic, digits=2),
                        "  wAICc = ", best_traj[[paste0("weight_aic_", plot_class)]][i],
                        "  NRMSE = ", round(best_traj[[paste0("nrmse_", plot_class)]][i], digits=2),
                        title_part)) +
          theme(plot.title = ggtext::element_markdown(size=10, lineheight = 1.1))

      }
    }


    # Color the background accordingly
    if (best_traj$class[i] == plot_class){

      bkg_col <- "lightskyblue"
      # if("expected_class" %in% names(best_traj)){
      #
      #   if (best_traj$class[i] == best_traj$expected_class[i]) { bkg_col <- "#9af394"
      #   } else { bkg_col <- "lightcoral" }
      #
      # # } else { bkg_col <- "lightblue" }
      # } else { bkg_col <- "lightskyblue" }

    } else { bkg_col <- "grey80" }

    p <- p + theme(plot.background = element_rect(fill = bkg_col))


    # Add as_detect plot
    if (plot_class == "abrupt" & !is.null(rslt$abt_res$asd)){

      p_bis <- cowplot::plot_grid(p,
                               p_asd + theme(plot.background = element_rect(fill = bkg_col)),
                               align = 'hv', ncol=1, rel_heights = c(3, 2))


    # Weight/LOO radar plot
      wrad <- best_traj %>%
        dplyr::select(simu_id | dplyr::contains("weight")) %>%
        `colnames<-`(c("simu_id","nch","lin","qdr","abt")) %>%
        dplyr::slice(i)

      wrad_plot <-
        ggradar::ggradar(wrad, axis.label.size = 1.5,
                         grid.label.size = 0, group.point.size = 1, group.line.width = .2,
                         group.colours = "red",
                         background.circle.transparency=0, centre.y=0, gridline.mid.colour="grey20",
                         gridline.min.colour="grey20", gridline.max.colour="grey20",
                         axis.line.colour="grey20", grid.line.width=0.25)+
        theme(
          plot.background = element_blank(),
          panel.background = element_blank(),
          plot.caption = element_text(hjust = 0.5, vjust = 0, size = 6))+
        labs(caption = "wAICc")

      p_bis <- cowplot::ggdraw() +
        cowplot::draw_plot(p_bis, x = 0, y = 0, width = 1, height = 1) +
        cowplot::draw_plot(wrad_plot, x = 0.75, y = .86, width = .15, height = .15)

      if (!is.null(best_traj_loo)){

        loorad <- best_traj %>%
          dplyr::select(simu_id | dplyr::contains("loo")) %>%
          `colnames<-`(c("simu_id","nch","lin","qdr","abt")) %>%
          dplyr::slice(i)

        loorad_plot <-
          ggradar::ggradar(loorad, axis.label.size = 1.5,
                           grid.label.size = 0, group.point.size = 1, group.line.width = .2,
                           group.colours = "red",
                           background.circle.transparency=0, centre.y=0, gridline.mid.colour="grey20",
                           gridline.min.colour="grey20", gridline.max.colour="grey20",
                           axis.line.colour="grey20", grid.line.width=0.25)+
          theme(
            plot.background = element_blank(),
            panel.background = element_blank(),
            plot.caption = element_text(hjust = 0.5, vjust = 0, size = 6))+
          labs(caption = "LOO")

        p_bis <- p_bis +
          cowplot::draw_plot(loorad_plot, x = 0.85, y = .86, width = .15, height = .15)
      }

      plot_bis[[i]] <- p_bis

    }

    plots[[i]] <- p

  }

  if (plot_class == "abrupt" & !is.null(rslt$abt_res$asd) & detection_plot){

    return(list("plots"=plots, "plot_bis"=plot_bis))

  } else {

    return(plots)
  }
}




# Plot all time series for a given noise combination

plot_given_noise <- function(path, noise_comb, save=TRUE){

  simu_list_plot <- readRDS(path)
  scen_exp <- unique(simu_list_plot[[1]] %>% pull(scen))

  gl <- lapply(1:length(scen_exp), function(i){

    plot_simu(simu_list_plot[[noise_comb]] %>%
                dplyr::filter(scen==scen_exp[i]))[[2]]

  })


  if (save){
    png(paste0("res/03_ricker_library/simu_library/scenarios_",noise_comb,".png"), width = 800, height = 1200)
    facet <- gridExtra::grid.arrange(grobs=gl, ncol=2, clip=TRUE)
    dev.off()
  }

  return(gl)
}


# LOO by class and classification correctness
weight_boxplot_class_fun <- function(traj, save, workflow, title){

  traj_TCbest_wAIC <- traj %>%
    tidyr::pivot_longer(cols=contains("weight_"),
                        names_to = "weight_class",
                        values_to = "weight") %>%
    dplyr::mutate(weight_class = sub("weight_aic_","",weight_class)) %>%
    dplyr::filter(weight_class==class)

  traj_TCbest_wAIC %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(mean_weight = mean(weight))

  wAIC_boxplot <- traj_TCbest_wAIC %>%
    mutate(correct = ifelse(class==expected_class, TRUE, FALSE),
           length = factor(length, levels = c("l100","l50","l20"))) %>%
    ggplot(aes(x = class, y = weight, fill = correct))+
    # geom_jitter(aes(colour=correct), width=0.1, height = 0, alpha=0.5)+
    geom_point(aes(colour=correct), alpha=0.02,
               position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.1))+
    # geom_violin(alpha = 0.7)+
    geom_boxplot(alpha = 0.3, outlier.shape = NA)+
    stat_summary(fun.data = function(x){
      return(c(y = 0.1, label = length(x))) },
      geom = "text", fun = median, position = position_dodge(width = 0.8))+
    theme_light()+
    expand_limits(y=0)+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
    guides(fill="none", colour="none")+
    facet_grid(vars(length))+
    labs(x="Trajectory class", y="AICc weight")+
    ggtitle(title)
  wAIC_boxplot

  # if("expected_class" %in% names(traj)) type <- "sim" else type <- "RAM"
  # ts_type <- traj$simu_id[1] %>% stringr::str_split("_") %>% `[[`(1) %>% `[`(1)

  if(save) ggsave(plot = wAIC_boxplot, width=8, height=10,
                  paste0("res/03_ricker_library/simu_library/wAIC_boxplot_class_",workflow,".png"))

  return(wAIC_boxplot)
}



# Post hoc analyses -------------------------------------------------------



#' Pool confusion matrices into one (by length)
#'
#' @param path the path to the .rds file output from traj function
#' @param title a character to display as title
#'
#' @return a pheatmap
#' @export

conf_mat_pool <- function(path, title, show_legend=FALSE) {

  outlist <- readRDS(path)
  pool_mat <- outlist$mat_list[[1]]$conf_mat_class$table

  for(i in 2:length(outlist$mat_list)){
    pool_mat <- pool_mat + outlist$mat_list[[i]]$conf_mat_class$table
  }
  colnames(pool_mat) <- rownames(pool_mat) <- sub("_"," ",colnames(pool_mat))

  if(show_legend){
    ftsz <- 15 ; ftsz_nb <- 20
  } else { ftsz <- 150 ; ftsz_nb <- 200 }

  byclass <- sum(pool_mat[,1])
  p <- pheatmap::pheatmap(pool_mat/byclass,
                          color = grDevices::colorRampPalette(
                            rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),
                          breaks = seq(0, 1, .01),
                          cluster_rows = FALSE, cluster_cols = FALSE,
                          legend = show_legend, display_numbers = TRUE,
                          number_format = "%.2f", fontsize = ftsz,
                          number_color="black", fontsize_number = ftsz_nb,
                          # main=stringr::str_split(path, "/")[[1]][4] %>% sub(x=., "_1",""),
                          main=title)
  return(p)
}


