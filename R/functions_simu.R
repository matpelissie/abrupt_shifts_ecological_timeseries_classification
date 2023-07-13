###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###
#
# 22/03/2022 mathieu.pelissie@ens-lyon.fr
#
# Functions for stock dynamics simulations
#
# functions_simu.R
#
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###


# I) Simulated data -----------------------------------------------------

#' Make scenario
#'
#' @param shape shape of the time series to generate (cst, lin, qdr, abt)
#' @param l length of the time series
#' @param trend (lin, qdr, abt) "pos" or "neg"
#' @param min (lin, qdr, abt) minimal value
#' @param max (lin, qdr, abt) maximal value
#' @param velocity (qdr) "acc" or "dec"
#' @param breaktime (abt) time step at which the step occurred
#'
#' @return a vector with the desired features and a string with the parameters used
#'
#' @export

make_scen <- function(shape, l=100, trend=NULL,
                      min=NULL, max=NULL, velocity=NULL,
                      breaktime=floor(l/2)){

  if(shape=="cst") {

    scen <- rep(min, l)
    name <- sprintf("l%g_%s_X%g", l, shape, min)


  } else if(shape=="lin") { # Shape linear

    scen <- seq(min, max, length.out = l)
    name <- sprintf("l%g_%s_%s_X%g-%g", l, shape, trend, min, max)
    if(trend=="neg"){
      scen <- rev(scen)
      name <- sprintf("l%g_%s_%s_X%g-%g", l, shape, trend, max, min)
    }

  } else if(shape=="ccv") { # Shape concave

    scen <- c(seq(max, min, length.out = breaktime),
              seq(min, max, length.out = l-breaktime))
    name <- sprintf("l%g_%s_X%g-%g-%g", l, shape, max, min, max)

  } else if(shape=="cvx") { # Shape convex

    scen <- c(seq(min, max, length.out = breaktime),
              seq(max, min, length.out = l-breaktime))
    name <- sprintf("l%g_%s_X%g-%g-%g", l, shape, min, max, min)


  } else if(shape=="qdr") { # Shape quadratic

    if(velocity=="dec"){

      if(trend=="pos"){
        a <- (min-max)/l^2
        b <- 2*(max-min)/l
        c <- min
        name <- sprintf("l%g_%s_%s_%s_X%g-%g", l, shape, trend, velocity, min, max)
      }

      if(trend=="neg"){
        a <- (max-min)/l^2
        b <- 2*(min-max)/l
        c <- max
        name <- sprintf("l%g_%s_%s_%s_X%g-%g", l, shape, trend, velocity, max, min)
      }
    }

    if(velocity=="acc"){

      if(trend=="pos"){
        a <- (max-min)/l^2
        b <- 0
        c <- min
        name <- sprintf("l%g_%s_%s_%s_X%g-%g", l, shape, trend, velocity, min, max)
      }

      if(trend=="neg"){
        a <- (min-max)/l^2
        b <- 0
        c <- max
        name <- sprintf("l%g_%s_%s_%s_X%g-%g", l, shape, trend, velocity, max, min)
      }

    }

    pol2deg <- function(x) a*x^2+b*x+c
    scen <- pol2deg(1:l)


  } else if(shape=="abt") { # Shape abrupt

    if(breaktime>l) stop("Break point should occur within the time series")

    if(trend=="pos"){
      scen <- c(rep(min, breaktime-1), rep(max, l-breaktime+1))
      name <- sprintf("l%g_%s_%s_brk%g_X%g-%g", l, shape, trend, breaktime, min, max)
    }

    if(trend=="neg"){
      scen <- c(rep(max, breaktime-1), rep(min, l-breaktime+1))
      name <- sprintf("l%g_%s_%s_brk%g_X%g-%g", l, shape, trend, breaktime, max, min)
    }
  }


  return(list("scen"=scen, "name"=name))

}



#' Make multiple scenario with multiple breakpoints
#'
#' @param l length of the time series
#' @param trend "ccv", "cvx", "pos", or "neg"
#' @param vals list of plateau values
#' @param brks list of time steps at which the breaks occurred
#'
#' @return a vector with the desired features and a string with the parameters used
#'
#' @export

make_scen_multibrk <- function(l=100, trend=NULL,
                      val1=NULL, val2=NULL, val3=NULL,
                      brk1=NULL, brk2=NULL){

  if(brk1>l | brk2>l) stop("Breakpoints should occur within the time series")

  shape <- "abt"
  scen <- c(rep(val1, brk1), rep(val2, brk2-brk1), rep(val3, l-brk2))
  name <- sprintf("l%g_%s_%s_brk%g-%g_X%g-%g-%g", l, shape, trend, brk1, brk2, val1, val2, val3)

  return(list(scen, name))
}





#' Get csv data
#'
#' @param file relative path and file name with noise scenarios
#' @return a csv with noise or simulation paramaters for scenarios in row
#' @export

get_data <- function(file){
  readr::read_csv(file)
}


#' Make simulations and plots
#'
#' @param name character for scenario name
#' @param Ts length of time series
#' @param r vector of growth rates (same length as time series)
#' @param F vector of harvesting rates (same length as time series)
#' @param K carrying capacity
#' @param P exponent for sigmoid functional response
#' @param H half-saturation rate
#' @param sr process (i.e. demographic) error
#' @param se environmental noise
#' @param su observational noise (i.e. measurement error)
#' @param iter number of repetitions of time series
#' @param thr unfeasibility threshold
#' @param init percentage of carrying capacity as initial size
#' @param trs length of the transient period (to be removed, default = 50 steps)
#' @param expected_class expectation trajectory class
#' @param jfr frequency of jump events
#' @param jsz base size of the jump event (buffered by process error and population size)
#'
#' @return a list of four elements:
#' - a list of matrix with time varying inputs and outputs for all iterations
#' - a data frame with time varying inputs and outputs for all iterations
#' - patchwork plot with TB, TC, ER, SProd time series
#' - a character indicating the percentage of iterations that went extinct
#' @export

run_simu <- function(name, Ts, r, F, K, P, H, sr, se, su, iter, thr, init, trs=50, expected_class, jfr, jsz){


  # if and only if F is variable
  full_name <- sprintf("%s_r%#.2g_H%g", sub("X", "F", name), r, H)

  if (length(r)==1) r <- rep(r,Ts)
  if (length(H)==1) H <- rep(H,Ts)

  if (length(r)!=length(F) | length(F)!=length(H) | length(F)!= Ts){
    stop("Parameters should be vectors the same length as time series length Ts")
  }

  TBts <- matrix(NA, nrow=Ts, ncol=iter)
  TCts <- matrix(NA, nrow=Ts, ncol=iter)
  rts <- matrix(NA, nrow=Ts, ncol=iter)
  Fts <- matrix(NA, nrow=Ts, ncol=iter)
  Hts <- matrix(NA, nrow=Ts, ncol=iter)
  ERts <- matrix(NA, nrow=Ts, ncol=iter)
  SProdts <- matrix(NA, nrow=Ts, ncol=iter)
  EIGts <- matrix(NA, nrow=Ts, ncol=iter)
  srts <- matrix(NA, nrow=Ts, ncol=iter)
  sets <- matrix(NA, nrow=Ts, ncol=iter)
  suts <- matrix(NA, nrow=Ts, ncol=iter)

  b <- r/K

  i <- 1
  extra <- 0
  while (i <= iter){

    # Begin at a given % of carrying capacity:
    x <- K*init

    # Let the equilibrium be reached:
    for (t in 1:trs){

      N <- x*exp(r[1]-b[1]*x)
      TC <- F[1]*x^P/(H[1]^P+x^P)
      x <- N - TC

      if (x<thr){
        if (x<(-10^1)) x <- -10^1 # to avoid reaching -Inf
        print("Warning: the systems fails to stabilize, check the parameters")
      }

    }

    # First time point with stochasticity:
    TBts[1, i] <- x*exp(su*rnorm(1))

    for (t in 2:Ts){

      r1 <- r[t]*exp(rnorm(n=1, sd=sr)) # process error (demographic stochasticity)
      N <- x*exp(r1-b[t]*x+rnorm(n=1, sd=se))
      TC <- F[t]*x^P/(H[t]^P+x^P)

      if (runif(1) < jfr){ # In case of jump is present

        jump_ok <- FALSE
        jump_count <- 0

        while(!jump_ok & jump_count<50){ # Avoid jumping to extinction

          J <- rnorm(n=1, sd=se)*jsz*x
          if (N - TC + J > 0) jump_ok <- TRUE
          if (N - TC + J < 0) jump_count <- jump_count + 1
        }

      } else {
        J <- 0
      }

      x <- N - TC + J

      # Compute stability:
      eig <- exp(r[t]-b[t]*K*init+se*rnorm(1))*(1-b[t]*K*init) - F[t]*((2*K*init/((K*init)^2+H[t]^2))-(2*(K*init)^3/((K*init)^2+H[t]^2)^2))


      if (x<thr){

        extra <- extra+1
        break
        # if (x<(-10^1)) x <- -10^1 # to avoid reaching -Inf

      } else {

        TBts[t, i] <- x*exp(rnorm(n=1, sd=su)) # observation error
        TCts[t, i] <- TC
        rts[t, i] <- r1
        Fts[t, i] <- F[t]
        Hts[t, i] <- H[t]
        EIGts[t, i] <- eig
        srts[t, i] <- sr
        sets[t, i] <- se
        suts[t, i] <- su

      }

    }

    # Is constant time series not collapsed:
    # if (stringr::str_detect(name,"cst") & x>0.1*TBts[1, i]) i <- i+1
    if (!expected_class %in% c("abrupt") & x>0.1*TBts[1, i] & TBts[90, i]>0.1*TBts[1, i]) i <- i+1

    # Is time series complete:
    # else if (!stringr::str_detect(name,"cst") & x>thr) i <- i+1
    else if (expected_class %in% c("abrupt") & x>thr) i <- i+1


    # To escape if not possible to complete the time series:
    if (extra>2*max(iter,50)){
      return(warning("Warning: too many simulations fail to avoid extinction, check parameters and noise"))
    }
  }

  for (k in 2:Ts){

    ERts[k,] <- TCts[k-1,]/TBts[k,]
    SProdts[k-1,] <- TBts[k,] - TBts[k-1,] + TCts[k-1,]

  }

  ERts <- apply(ERts, 2, function(x) replace(x, is.infinite(x) | is.nan(x), NA))


  lib_ts <- data.frame(scen = as.vector(full_name),
                       iter = rep(1:iter, each=Ts),
                       year = rep(1:Ts, iter),
                       TB = as.vector(TBts),
                       r = as.vector(rts),
                       F = as.vector(Fts),
                       H = as.vector(Hts),
                       P = rep(P, Ts*iter),
                       EIG = as.vector(EIGts),
                       TC = as.vector(TCts),
                       ER = as.vector(ERts),
                       SProd = as.vector(SProdts),
                       sr = as.vector(srts),
                       se = as.vector(sets),
                       su = as.vector(suts)
                       )

  scen_mat <- list("TB"=TBts, "TC"=TCts, "r"=rts, "F"=Fts, "ER"=ERts, "SProd"=SProdts)


  return(list(
    "scen_mat"=scen_mat,
    "df"=lib_ts
    ))

}



#' Is stock collapsed?
#'
#' @param m total biomass matrix
#' @param thr threshold of collapsed state (default is 0.2)
#'
#' @return a list of two objects:
#' - a list the same length as iterations with boolean indicating if stocks
#' reached long-lasting collapsed state (i.e. > 5 consecutive years below threshold)
#' - a matrix of the same dimensions as the input matrix with logical info
#' regarding instantaneous collapsed state of stock
#' @export

is_collapsed <- function(m, thr = 0.2){

  popmax <- apply(m, 2, max)
  punc_collapsed <- apply(m, 2, function(x) ifelse(x<popmax*thr, TRUE, FALSE))

  lcollapsed <- punc_collapsed %>%
    as.data.frame %>%
    lapply(function(x) x)

  collapsed <- lcollapsed %>%
    lapply(function(x) max( c(with(rle(x), lengths[values]), 0)) >5 )

  return(list(collapsed, punc_collapsed))
}



# Wrapper functions -------------------------------------------------------

#' Make and save simulated data
#'
#' @param param_df a data frame with scenarios parameters
#' @param scen_fct a character to specify the scenario-making function to use
#'
#' @return a data frame with simulated data for all scenarios
#' @export

make_store_simu <- function(param_df, scen_fct, se, sr, su, jfr, jsz){

  ## Repetitive?
  # for (k in 1:length(noise_names)){
  #
  #   assign(noise_names[k], noise_df %>%
  #            dplyr::pull(noise_names[k]) %>%
  #            `[`(n))
  # }

  param_names <- colnames(param_df)
  all_simu <- data.frame(df=matrix(ncol = 17))[-1,]

  for (j in 1:nrow(param_df)){ # For each scenario

    for (i in 1:length(param_names)){ # import parameter values

      assign(param_names[i], param_df %>%
               dplyr::pull(param_names[i]) %>%
               `[`(j))
    }

    if (scen_fct == "make_scen"){
      Fseq <- make_scen(shape=shape, trend=trend, min=min, max=max,
                        velocity=velocity, l=l, breaktime = breaktime)

    } else if (scen_fct == "make_scen_multibrk"){
      Fseq <- make_scen_multibrk(trend=trend, l=l,
                                 val1=val1, val2=val2, val3=val3,
                                 brk1=brk1, brk2=brk2)
    }

    sim <- run_simu(name=Fseq[[2]], r=r, F=Fseq[[1]], Ts=length(Fseq[[1]]),
                    K=10, P=2, H=H, sr=sr, se=se, su=su, iter=iter,
                    thr=0, init=0.9, expected_class=expected_class,
                    jfr=jfr, jsz=jsz)

    if(class(sim) != "character"){ # if no warnings

      sim$df <- sim$df %>%
        dplyr::mutate(expected_traj=expected_traj,
                      expected_class=expected_class)
      all_simu <- all_simu %>% rbind(sim$df)
    }

  }

  # saveRDS(all_simu,
  #         sprintf("data/03_simulations/all_simu_se%g_sr%g_su%g_l%g_iter%g.rds",
  #                 se, sr, su, l, iter))
  # write.csv(all_simu,
  # sprintf("data/03_simulations/all_simu_se%g_sr%g_su%g_l%g_iter%g.rds",
  # se, sr, su, l, iter))

  return(all_simu)
}




# Bifurcation plot --------------------------------------------------------

#' Plot stable states and evidence potential alternative stable states
#'
#' @param rval a value of growth rate
#' @param P exponent for sigmoid functional response
#' @param H half-saturation rate
#' @param Fmax the maximal harvesting rate
#' @param Fstep the interval between two values of harvesting rate
#' @param TBhigh the starting value for high biomass
#'
#' @return a ggplot with values of stable states along a range of F
#' @export

ASS_row <- function(Ts, rval, P, H, K, Fmax, Fstep, TBhigh){

  Fseq <- seq(0, Fmax, Fstep)

  row_inc <- row_dec <- states_inc <- states_dec <- stab_inc <- stab_dec <-
    data.frame(matrix(NA, length(rval), length(Fseq))) %>%
    `rownames<-`(paste0("r", rval)) %>%
    `colnames<-`(paste0("F", Fseq))

  # Increasing F & high starting biomass
  init <- TBhigh
  for (j in Fseq){

    F <- rep(j, Ts)
    r <- rep(rval, Ts)
    simu <- run_simu(name="test", Ts, r, F, K, P, H=H, sr=0, se=0, su=0, iter=1, thr=0, init,
                     expected_class="abrupt", jfr=0, jsz=0)[[2]] %>% try()
    # expected_class abrupt to avoid removing collapsed time series

    if(class(simu) !="try-error"){

      states_inc[paste0("r",rval), paste0("F",j)] <- eq_state(simu, thr=0)[[2]]

      # init <- simu[nrow(simu), "TB"]/K # take last value as stable state value
      init <- mean(simu[(nrow(simu)*0.9):nrow(simu), "TB"])/K # take mean of 100 last values for cycles
      init <- max(init, 10^(-5)) # to avoid reaching zero
      row_inc[paste0("r", rval), paste0("F", j)] <- init*K
      stab_inc[paste0("r", rval), paste0("F", j)] <- simu %>% dplyr::pull("EIG") %>% tail(1)

    } else {
      states_inc[paste0("r",rval), paste0("F",j)] <- "unfeasible"
      row_inc[paste0("r", rval), paste0("F", j)] <- NA
    }

    # print(j)
  }

  states_inc_long <- states_inc %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "states_inc") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  row_inc_long <- row_inc %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "value_inc") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  stab_inc_long <- stab_inc %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "stab_inc") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  # Decreasing F & low starting biomass
  # init <- TBlow
  # start from stable equilibrium of the increasing F, no need for TBlow
  for (j in rev(Fseq)){

    F <- rep(j, Ts)
    r <- rep(rval, Ts)
    simu <- run_simu(name=NA, Ts, r, F, K, P, H=H, sr=0, se=0, su=0, iter=1, thr=0, init,
                     expected_class="abrupt", jfr=0, jsz=0)[[2]] %>% try()

    if(class(simu) !="try-error"){

    states_dec[paste0("r",rval), paste0("F",j)] <- eq_state(simu, thr=0)[[2]]

    # init <- simu[nrow(simu), "TB"]/K # take last value as stable state value
    init <- mean(simu[(nrow(simu)*0.9):nrow(simu), "TB"])/K # take mean of 100 last values for cycles
    init <- max(init, 10^(-5)) # to avoid reaching zero
    row_dec[paste0("r", rval), paste0("F", j)] <- init*K
    stab_dec[paste0("r", rval), paste0("F", j)] <- simu %>% dplyr::pull("EIG") %>% tail(1)

    } else {
      states_dec[paste0("r",rval), paste0("F",j)] <- "unfeasible"
      row_dec[paste0("r", rval), paste0("F", j)] <- NA
    }

    # print(j)
  }

  states_dec_long <- states_dec %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "states_dec") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  row_dec_long <- row_dec %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "value_dec") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  stab_dec_long <- stab_dec %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "stab_dec") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  bif_wide <- row_inc_long %>%
    dplyr::left_join(row_dec_long, by=c("r", "F")) %>%
    dplyr::left_join(states_inc_long, by=c("r", "F")) %>%
    dplyr::left_join(states_dec_long, by=c("r", "F")) %>%
    dplyr::left_join(stab_inc_long, by=c("r", "F")) %>%
    dplyr::left_join(stab_dec_long, by=c("r", "F")) %>%
    # arbitrary threshold of 10^(-1)
    dplyr::mutate(ASS = ifelse(abs(value_inc - value_dec)>10^(-1), TRUE, FALSE)) %>%
    dplyr::mutate(diffstate = ifelse(states_inc != states_dec, TRUE, FALSE))


  bif_long <- bif_wide %>%
    tidyr::pivot_longer(cols = contains("_"),
                        names_to = c(".value","sense"),
                        names_pattern = "(\\w+)_(\\w+)")


  # bifColor <- c("#E64B35B2", "#4DBBD5B2")
  coeff <- max(abs(bif_long %>% pull(stab))/10)

  bif_plot <- ggplot(bif_long, aes(x=F, y=value, color=as.factor(sense)))+
    geom_line()+
    geom_point(size=0.5)+
    # geom_line(aes(x=F, y=stab/coeff, color=as.factor(sense)), lty=3)+
    # scale_y_continuous(
    #   name = "TB at equilibrium",
    #   sec.axis = sec_axis(~.*coeff, name = "Stability"))+
    labs(y="TB at equilibrium", color="Sense")+
    theme_light()+
    expand_limits(y = 0)+
    theme(legend.position = c(0.9,0.8))
    # labs(caption = paste0("r = ", rval, " ; H = ", H, " ; P = ", P))

  stab_plot <- ggplot(bif_long, aes(x=F, y=stab, color=as.factor(sense)))+
    geom_hline(yintercept=0, lty=3)+
    geom_line(lty=1)+
    labs(y="Stability")+
    theme_light()+
    expand_limits(y = 0)+
    theme(legend.position = "none")+
    labs(caption = paste0("r = ", rval, " ; H = ", H, " ; P = ", P))

  bif_stab <- bif_plot / stab_plot +
  patchwork::plot_layout(heights = unit(c(4, 1), c('null')))

  return(list(bif_wide, bif_stab))

}

#' Save the plot of (alternative) stable states
#'
#' @param H half-saturation rate
#' @param P exponent for sigmoid functional response
#' @param rval a value of growth rate
#' @param Fmax the maximal harvesting rate
#' @param Fstep the interval between two values of harvesting rate
#' @param TBhigh the starting value for high biomass
#'
#' @return No return value
#' @export

ASS_save <- function(Ts, rval, P, H, Fmax, Fstep, TBhigh){

  ASS_plot <- ASS_row(Ts=Ts, H=H, P=P, rval=rval,
                      Fmax=Fmax, Fstep=Fstep,
                      TBhigh=TBhigh)[[2]]
  filename <- sprintf("res/03_ricker_library/h%g/h%g_p%g/03_ASS_h%g_p%g_r%#.2g_tbhigh%g.png",
                      H, H, P, H, P, rval, TBhigh)
  ggsave(ASS_plot, filename=filename,
         height=1900, width=2100, unit="px")

  return(invisible(NULL))

}



#' Give the state at the equilibrium for a given time series
#'
#' @param lib_ts data frame with time varying inputs and outputs for all iterations
#' @param thr unfeasibility threshold (default 0)
#'
#' @return a character indicating state at equilibrium (unfeasible, stable, cycles, chaos)
#' @export

eq_state <- function(lib_ts, thr=0){

  sample_length <- nrow(lib_ts)/10

  unique_diff <- lib_ts %>%
    tail(n=sample_length) %>%
    dplyr::mutate(diff_lag1=TB-dplyr::lag(TB)) %>%
    tidyr::drop_na(diff_lag1) %>%
    dplyr::distinct(diff_lag1) %>%
    dplyr::pull()

  # Temporary ###-###-###-##-##
  plot_diff <- ggplot(lib_ts %>%
         tail(n=100) %>%
         dplyr::mutate(diff_lag1=TB-dplyr::lag(TB)) %>%
         tidyr::drop_na(diff_lag1),
         aes(x=year, y=diff_lag1))+
    geom_line()+
    geom_point()
  ###-###-###-###-#

  state <- NA

  if (lib_ts %>% dplyr::pull(TB) %>% tail(1) < thr # last value below unfeasibility threshold
      ) {
    state <- "unfeasible"

    } else if (unique_diff %>% length() == 1 | # fully stabilized
               max(abs(unique_diff)) < 10^(-5) | # if oscillating very little
               unique_diff %>% sign() %>% unique() %>% length() == 1 # if not yet stabilized

               ) {
      state <- "stable"

    } else if (unique_diff %>% length() < (sample_length-1) | # regular cycles
               sign(abs(unique_diff[-1]) - abs(stats::lag(unique_diff, k=2)[-1])) %>% unique() %>% length() == 1 # if lag-2 cycles of decreasing amplitude

               ) {
      state <- "cycle"

    } else if (unique_diff %>% length() == sample_length-1) { # all lag-1 diff values different
      state <- "chaos"

    }

  return(list(unique_diff, state, plot_diff))

}


#' Plot system stable states under a range of constant r and F
#'
#' @param Ts length of time series
#' @param K carrying capacity
#' @param rmax the maximal growth rate value (default 3)
#' @param rstep the interval between two values of growth rate (default 0.1)
#' @param Fmax the maximal harvesting rate value (default 3)
#' @param Fstep the interval between two values of harvesting rate (default 0.1)
#' @param H half-saturation rate
#' @param P exponent for sigmoid functional response
#' @param thr extinction threshold
#' @param init percentage of carrying capacity as initial size
#'
#' @return a list of two objects:
#' a ggplot with stable states along ranges of r and F
#' the data frame used for the ggplot
#' @export

bifurc_analysis <- function(Ts, K, rmax=3, rstep=0.1,
                            Fmax=3, Fstep=0.1, P, H, thr, init,
                            save, dirname="res/03_ricker_library"){

  rseq <- seq(0, rmax, rstep)
  Fseq <- seq(0, Fmax, Fstep)

  ASS <- data.frame()

  for (i in rseq){

    ASS <- ASS %>%
      rbind(ASS_row(Ts=Ts, rval = i, Fmax=Fmax, K=K, Fstep=Fstep, H=H, P=P, TBhigh=init)[[1]])

    print(paste0(i, "/", rmax))
  }

  # Plot states

  bif_ana <- ASS %>%
    dplyr::mutate(ASS = ifelse(ASS, TRUE, NA))

  bif_ana$states_inc <- factor(bif_ana$states_inc,
                          levels = c("unfeasible","stable","cycle","chaos"))

  bif_ana$states_dec <- factor(bif_ana$states_dec,
                               levels = c("unfeasible","stable","cycle","chaos"))
  state.col <- c(unfeasible="#c44601", stable="#5ba300", cycle="#0073e3", chaos="#e79b01")

  # Types of stable equilibrium and ASS when increasing F
  plot_state_inc <- ggplot(bif_ana, aes(x=F, y=r))+
    geom_tile(aes(fill=states_inc, alpha=value_inc))+
    scale_fill_manual(values = state.col)+
    scale_alpha(range=c(0.3,1), limits=c(0,10), na.value=1)+
    geom_point(data=subset(bif_ana, !is.na(ASS) & !states_inc %in% c("unfeasible", "chaos")),
               # Display ASS only for stable equilibrium
               aes(x=F, y=r), size=1)+
    scale_color_manual(values = c("TRUE"="black"))+
    # geom_hline(yintercept=c(1, 1.6, 2.1))+
    theme_minimal()+
    labs(caption = paste0("h = ", H, " ;  init TB inc = ", init*100, "% K ; increasing F"))


  # Types of stable equilibrium and ASS when decreasing F
  plot_state_dec <- ggplot(bif_ana, aes(x=F, y=r))+
    geom_tile(aes(fill=states_dec, alpha=value_dec))+
    scale_fill_manual(values = state.col)+
    scale_alpha(range=c(0.3,1), limits=c(0,10))+
    geom_point(data=subset(bif_ana, !is.na(ASS) & !states_dec %in% c("unfeasible", "chaos")),
               # Display ASS only for stable equilibrium
               aes(x=F, y=r), size=1)+
    scale_color_manual(values = c("TRUE"="black"))+
    # geom_hline(yintercept=c(1, 1.6, 2.1))+
    theme_minimal()+
    labs(caption = paste0("h = ", H, " ;  init TB inc = ", init*100, "% K ; decreasing F"))


  # Plot values

  bif_ana_val <- bif_ana %>%
    dplyr::mutate(value_inc = ifelse(!states_inc %in% c("unfeasible", "chaos"), value_inc, NA),
                  value_dec = ifelse(!states_dec %in% c("unfeasible", "chaos"), value_dec, NA))

  # Equilibrium values when increasing F
  plot_val_inc <- ggplot(bif_ana_val, aes(x=F, y=r))+
    geom_tile(aes(fill=value_inc))+
    geom_point(data=subset(bif_ana, !is.na(ASS) & !states_inc %in% c("unfeasible", "chaos")),
               # Display ASS only for stable equilibrium
               aes(x=F, y=r), size=1)+
    scale_color_manual(values = c("TRUE"="black"))+
    # geom_hline(yintercept=c(1, 1.6, 2.1))+
    theme_minimal()+
    labs(caption = paste0("h = ", H, " ;  init TB inc = ", init*100, "% K ; increasing F"))

  # Equilibrium values when decreasing F
  plot_val_dec <- ggplot(bif_ana_val, aes(x=F, y=r))+
    geom_tile(aes(fill=value_dec))+
    geom_point(data=subset(bif_ana, !is.na(ASS) & !states_dec %in% c("unfeasible", "chaos")),
               # Display ASS only for stable equilibrium
               aes(x=F, y=r), size=1)+
    scale_color_manual(values = c("TRUE"="black"))+
    # geom_hline(yintercept=c(1, 1.6, 2.1))+
    theme_minimal()+
    labs(caption = paste0("h = ", H, " ;  init TB inc = ", init*100, "% K ; decreasing F"))

  plots_bif_ana <- list(plot_state_inc, plot_state_dec, plot_val_inc, plot_val_dec)
  print(plots_bif_ana)


  # save <- ifelse(readline(prompt="Do you want to save these plots? \n 1: Yes\n 0: No ")=="1", TRUE, FALSE)

  if (save){

    sprintf("%s/h%g", dirname, H) %>% dir.create(showWarnings = FALSE)
    sprintf("%s/h%g/h%g_p%g", dirname, H, H, P) %>% dir.create(showWarnings = FALSE)
    filename <- sprintf("%s/h%g/h%g_p%g/03_bif_ana_h%g_p%g_rmax%g_Fmax%g_init%g",
                        dirname, H, H, P, H, P, rmax, Fmax, init)

    saveRDS(bif_ana, paste0(filename,".rds"))

    plot_exts <- c("stateinc","statedec","valueinc","valuedec")
    lapply(plot_exts,
           function(x){
             ggsave(plots_bif_ana[[which(plot_exts==x)]],
                    filename=paste0(filename,"_", x, ".pdf"),
                    height=1900, width=2100,
                    unit="px")
           }
    )

  }

  return(list(plot_state_inc, plot_state_dec, plot_val_inc, plot_val_dec, bif_ana))

}



# Breakpoint algorithms ---------------------------------------------------


#' Breakpoints analysis using 'asdetect' package (based on change in gradient)
#'
#' @param stock_ts time series to analyse as 'ts' object
#' @param asd_thr numerical, value for detection threshold
#'
#' @return a one-row data frame with info about potentially detected breakpoints
#' @export

asd_fct <- function(stock_ts, asd_thr, check_true_shift, lowwl, highwl, mad_thr, mad_cst){

  len <- length(stock_ts)
  detect <- as_detect_mad(stock_ts, lowwl=lowwl, highwl=highwl, mad_thr=mad_thr, mad_cst=mad_cst)
  # returns position of shift (or of the maximum detected value):
  where <- try(where_as_quiet(detect, thresh = asd_thr, quiet=TRUE), silent=TRUE)

  if(class(where)[1]=="try-error"){

    where_low <- try(where_as_quiet(detect, thresh = 1, quiet=TRUE))

    where <- list()
    where$as_pos <- where_low$as_pos
    where$dt_val <- where_low$dt_val

  }

  # One line multiple as_detect breakpoints

  asd_out <- data.frame(abbr="asd",
                        mtd="asdetect",
                        n_brk=0,
                        loc_brk=NA,
                        aic=NA,
                        trend=NA,
                        mag=NA,
                        SDbef=NA,
                        SDaft=NA,
                        step_size=NA,
                        nrmse=NA)

  if (length(where$as_pos)>0){

    for (i in 1:length(where$as_pos)){

      if(check_true_shift){

        if(abs(where$dt_val[i])>asd_thr){

          # returns 1 if true abrupt shift detected:
          # true_shift <- asdetect::shift_type(stock_ts, where, width=7)
          if (len<20) true_shift <- custom_shift_type(stock_ts, where$as_pos[i], width = 2)
          else true_shift <- custom_shift_type(stock_ts, where$as_pos[i], width = 7)

        } else {true_shift <- 0}

        if (abs(where$dt_val[i])>asd_thr & true_shift == 1){

          asd_out["loc_brk"] <- paste0(asd_out["loc_brk"],";",
                                       c(where$as_pos[i], NA)[1] + start(stock_ts)[1] - 1)
          asd_out["n_brk"] <- asd_out["n_brk"] + 1
        }


      } else {

        if (abs(where$dt_val[i])>asd_thr){

          asd_out["loc_brk"] <- paste0(asd_out["loc_brk"],";",
                                       c(where$as_pos[i], NA)[1] + start(stock_ts)[1] - 1)
          asd_out["n_brk"] <- asd_out["n_brk"] + 1
        }

      }
    }
  }

  asd_out["loc_brk"] <- sub("NA;", "", asd_out["loc_brk"])

  return(list("df" = asd_out, "detect" = detect))

}

#' Breakpoints analysis using 'chngpt' package
#'
#' @param ts time series to analyse as 'ts' object
#'
#' @return a one-row data frame with info about potentially detected breakpoints
#' @export

chg_fct <- function(ts){

  Y <- tail(names(ts),1)

  chg <- chngpt::chngptm(formula.1 = as.formula(paste(Y,"~1")),formula.2 = ~year,
                         type="step", family="gaussian", data=ts)

  # chg$logliks
  pred_chg <- data.frame(year = ts$year,
                         bp = chg$best.fit$fitted.values
                         )
  # k <- 3 # number of variables in the model
  nrmse <- sqrt(sum(residuals(chg)^2)/length(ts$Y))/sd(ts$Y)

  chg_out <- data.frame(abbr = "chg",
                        mtd = "chgnpt",
                        n_brk = 1,
                        loc_brk = chg$chngpt,
                        aic = MuMIn::AICc(chg),
                        trend = ifelse(pred_chg$bp[1] > pred_chg$bp[length(pred_chg$bp)],
                                       "decrease", "increase"),
                        mag = pred_chg$bp[length(pred_chg$bp)] - pred_chg$bp[1],
                        SDbef = ts %>%
                          dplyr::filter(year<=chg$chngpt) %>%
                          dplyr::pull(Y) %>%
                          sd(),
                        SDaft = ts %>%
                          dplyr::filter(year>=chg$chngpt) %>%
                          dplyr::pull(Y) %>%
                          sd(),
                        nrmse = nrmse
                        )
  chg_out <- chg_out %>% mutate(step_size = mag/((SDbef+SDaft)/2))

  # stats4::logLik(chg) # (df=4)

  chg_outlist <- list("chg_out" = chg_out, "pred_chg" = pred_chg)

  return(chg_outlist)

}


#' Multiple breakpoints analysis
#'
#' @param ts time series to analyses as data frame (3 columns: year, Y, Y_SE)
#' @param abr_mtd a vector with abbreviation(s) corresponding to the breakpoints method(s) to use
#' @param asd_thr a numeric threshold for as_detect method
#' @param asd_chk a logical paramater for check_true_shift in asd_fct
#'
#' @return a list of three objects:
#' one data frame with info about potentially detected breakpoints
#' two lists of three data frames for bpm and bpt
#' one list of EnvCpt output
#' @export

shifts <- function(ts, abr_mtd, asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst){

  stock_ts <- ts %>%
    dplyr::pull(ncol(ts)) %>%
    ts(start=ts$year[1], end=tail(ts$year,1))

  # colnames(ts)[1] <- "year" # no more relevant

  res_table <- data.frame()


  if ("asd" %in% abr_mtd){

    asd_out <- asd_fct(stock_ts, asd_thr, check_true_shift=asd_chk, lowwl, highwl=highwl, mad_thr, mad_cst)
    res_table <- rbind(res_table, asd_out$df)
  }

  if ("chg" %in% abr_mtd){

    chg_outlist <- chg_fct(ts)
    chg_out <- chg_outlist$chg_out

    res_table <- rbind(res_table, chg_out)

  }

  if ("str_m" %in% abr_mtd & "chg" %in% abr_mtd){

    return(list("res_table" = res_table, "chg_outlist" = chg_outlist,
                "bpm_outlist" = bpm_outlist))

  } else if ("asd" %in% abr_mtd & "chg" %in% abr_mtd) {

    return(list("res_table" = res_table, "chg_outlist" = chg_outlist,
                "asd_detect" = asd_out$detect))

  } else if ("chg" %in% abr_mtd) {

    return(list("res_table" = res_table, "chg_outlist" = chg_outlist))

  } else {

    return(list("res_table" = res_table))

  }

}


# Make plots --------------------------------------------------------


#' #' Plot simulations
#' #'
#' #' @param lib_ts dataframe with simulated time series
#' #'
#' #' @return a patchwork plot with TB, TC, ER, SProd time series
#' #' @export
#'
#' plot_simu <- function (lib_ts){
#'
#'   colourCount <- length(unique(lib_ts$iter))
#'   getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))
#'
#'   # coeff <- max(lib_ts %>% drop_na() %>% pull(F)) / max(lib_ts %>% pull(TB))
#'   coeff <- 2.5
#'   tbColor <- "red"
#'     FColor <- "blue"
#'
#'     TB_plot <- ggplot(lib_ts, aes(x=year, y=TB, color=as.factor(iter), group=iter))+
#'       geom_hline(yintercept=0, lty=3)+
#'       geom_line(col=tbColor, alpha=0.1)+
#'       # geom_point(aes(size=punc_collapse), col="black", alpha=0.2, na.rm=TRUE)+
#'       scale_size_manual(values = c("TRUE"=2, "FALSE"=NA))+
#'       geom_line(data=lib_ts %>% drop_na(F), aes(x=year, y=F/coeff, color=as.factor(iter)),
#'                 col=FColor, lty=2)+
#'       scale_y_continuous(
#'         name = "Total biomass",
#'         sec.axis = sec_axis(~.*coeff, name = "Harvesting rate", breaks=seq(0,25,5)))+
#'       # expand_limits(y = 0, x = max(500))+
#'       # expand_limits(y = c(0, 8), x = max(lib_ts$year))+
#'       expand_limits(y = c(0, 12.5), x = max(lib_ts$year))+
#'       # labs(y="Total Biomass")+
#'       # hrbrthemes::theme_ipsum()+
#'       theme_light()+
#'       theme(legend.position = "none",
#'             axis.title.x = element_blank(),
#'             axis.title.y = element_text(color = tbColor, size=13),
#'             axis.title.y.right = element_text(color = FColor, size=13))+
#'       scale_color_manual(values = getPalette(colourCount))+
#'       labs(caption = paste0("r = ", lib_ts %>% pull(r) %>% tail(1) %>% round(2),
#'                             " ; H = ", lib_ts %>% pull(H) %>% tail(1),
#'                             " ; P = ", lib_ts %>% pull(P) %>% tail(1),
#'                             " ; F range = [", min(lib_ts %>% drop_na() %>% pull(F)) %>% round(2),
#'                             ", ", max(lib_ts %>% drop_na() %>% pull(F)) %>% round(2), "]",
#'                             ", sr = ", lib_ts %>% pull(sr) %>% tail(1),
#'                             ", se = ", lib_ts %>% pull(se) %>% tail(1)))
#'
#'     # locations of the first breakpoint in biomass
#'     # geom_vline(aes(xintercept=loc_cp, color=as.factor(iter)), linetype="solid", alpha=0.5)+
#'     # geom_vline(aes(xintercept=loc_cp2, color=as.factor(iter)), linetype="longdash", alpha=0.5)+
#'     # geom_vline(aes(xintercept=loc_bpm, color=as.factor(iter)), linetype="dashed", alpha=0.5)
#'     # geom_vline(aes(xintercept=loc_bpt, color=as.factor(iter)), linetype="dotted", alpha=0.5)
#'
#'
#'     TC_plot <- ggplot(lib_ts, aes(x=year, y=TC, color=as.factor(iter)))+
#'       geom_line(na.rm=TRUE)+
#'       expand_limits(y = 0, x = max(lib_ts$year))+
#'       labs(y="Total Catch")+
#'       theme_light()+
#'       theme(legend.position = "none",
#'             axis.title.x = element_blank())+
#'       scale_color_manual(values = getPalette(colourCount))
#'
#'     # locations of the first breakpoint in biomass
#'     # geom_vline(aes(xintercept=loc_cp, colour=as.factor(iter)), linetype="solid", alpha=0.5)+
#'     # geom_vline(aes(xintercept=loc_cp2, colour=as.factor(iter)), linetype="longdash", alpha=0.5)+
#'     # geom_vline(aes(xintercept=loc_bpm, colour=as.factor(iter)), linetype="dashed", alpha=0.5)
#'     # geom_vline(aes(xintercept=loc_bpt, colour=as.factor(iter)), linetype="dotted", alpha=0.5)
#'
#'
#'     ER_plot <- ggplot(lib_ts, aes(x=year, y=ER, color=as.factor(iter)))+
#'       geom_line(na.rm=TRUE)+
#'       expand_limits(y = 0, x = max(lib_ts$year))+
#'       labs(y="Exploitation Rate")+
#'       theme_light()+
#'       theme(legend.position = "none",
#'             axis.title.x = element_blank())+
#'       scale_color_manual(values = getPalette(colourCount))
#'
#'     # locations of the first breakpoint in biomass
#'     # geom_vline(aes(xintercept=loc_cp, colour=as.factor(iter)), linetype="solid", alpha=0.5)+
#'     # geom_vline(aes(xintercept=loc_cp2, colour=as.factor(iter)), linetype="longdash", alpha=0.5)+
#'     # geom_vline(aes(xintercept=loc_bpm, colour=as.factor(iter)), linetype="dashed", alpha=0.5)
#'     # geom_vline(aes(xintercept=loc_bpt, colour=as.factor(iter)), linetype="dotted", alpha=0.5)
#'
#'     SProd_plot <- ggplot(lib_ts, aes(x=year, y=SProd, color=as.factor(iter)))+
#'       geom_line(na.rm=TRUE)+
#'       expand_limits(y = 0, x = max(lib_ts$year))+
#'       labs(y="Surplus Production")+
#'       theme_light()+
#'       theme(legend.position = "none")+
#'       scale_color_manual(values = getPalette(colourCount))
#'
#'     # locations of the first breakpoint in biomass
#'     # geom_vline(aes(xintercept=loc_cp, colour=as.factor(iter)), linetype="solid", alpha=0.5)+
#'     # geom_vline(aes(xintercept=loc_cp2, colour=as.factor(iter)), linetype="longdash", alpha=0.5)+
#'     # geom_vline(aes(xintercept=loc_bpm, colour=as.factor(iter)), linetype="dashed", alpha=0.5)
#'     # geom_vline(aes(xintercept=loc_bpt, colour=as.factor(iter)), linetype="dotted", alpha=0.5)
#'
#'     patch <- TB_plot / TC_plot / ER_plot / SProd_plot
#'     patch <- patch + patchwork::plot_annotation(title = paste0("Scenario: ", "   iterations: ", colourCount))
#'
#'     return(list(patch, TB_plot))
#' }


#' Plot simulations without any axis
#'
#' @param lib_ts dataframe with simulated time series
#'
#' @return a minimal plot with TB curve
#' @export

plot_simu_light <- function (lib_ts){

  tbColor <- "black"

  TB_plot <- ggplot(lib_ts, aes(x=year, y=TB, color=as.factor(iter), group=iter))+
    # geom_hline(yintercept=0, lty=3)+
    geom_line(col=tbColor)+
    scale_y_continuous(
      name = "Total biomass")+
    expand_limits(y = c(0, 10), x = max(lib_ts$year))+
    theme_light()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()
    )

  return(TB_plot)
}



#' Plot simulations in a simple way
#'
#' @param lib_ts dataframe with simulated time series
#' @param tbcolor line color
#' @param alpha line transparency
#'
#' @return a minimal plot with TB curve
#' @export

plot_simu_simple <- function (lib_ts, var, tbcolor="black", alpha=1, xname="", yname="", ylim=12.5){

  TB_plot <- ggplot(lib_ts, aes(x=year, y=get(var), color=as.factor(iter), group=iter))+
    geom_line(col=tbcolor, alpha=alpha)+
    scale_x_continuous(name = xname)+
    scale_y_continuous(name = yname)+
    expand_limits(y = c(0, ylim))+
    theme_light()+
    theme(legend.position = "none")

  return(TB_plot)
}



#' Plot a time series with breakpoints
#'
#' @param ts time series dataframe
#' @param shifts_res output from shift function
#' @param abr_mtd a list of abbreviation with breakpoints methods to use (see function abrupt_classif for the list)
#'
#' @return a list of three plots
#'
#' @export

plot_bp <- function(ts, shifts_res, abr_mtd, best_traj=NULL, best_traj_loo=NULL){

  plots_out <- list()

  pop <- ggplot(ts, aes(x=get(grep("X|year",names(ts), value=TRUE)),
                        # y=get(names(ts)[2])/max(.data[[names(ts)[2]]])))+
                        y=get(names(ts)[2])))+

    geom_line()+
    labs(x=grep("X|year",names(ts), value=TRUE),
         y=as.character(names(ts)[2]))+
    theme_light(base_size = 7)+
    expand_limits(y=0)

  if (!is.null(best_traj_loo)){

    coeff <- max(ts[[names(ts)[2]]])

    pop <- pop +

      geom_histogram(data=best_traj_loo %>%
                       dplyr::filter(class=="abrupt"),
                     aes(x=loc_brk_chg, y=after_stat(ncount)*max(ts[[names(ts)[2]]])/2),
                     fill="blue", alpha=.3, colour=NA, binwidth = 1)+

      geom_histogram(data=best_traj_loo %>%
                       dplyr::filter(class=="abrupt"),
                     aes(x=loc_brk_asd, y=after_stat(ncount)*max(ts[[names(ts)[2]]])/2),
                     fill="red", alpha=.3, colour=NA, binwidth = 1)+
      # scale_y_continuous(
      #   sec.axis = sec_axis(~.*400/max(ts[[names(ts)[2]]])/2, name = "Brks freq LOO", breaks=seq(0,100,25)))+


    geom_point(data=best_traj_loo %>%
                 dplyr::left_join(ts, by="X") %>%
                 # dplyr::filter(class=="abrupt"), aes(x=X, y=Y, col=diff_loc), alpha=.7)
                 dplyr::filter(class=="abrupt"), aes(x=X, y=Y), col="orange", alpha=.7)
    pop

  }

  plots_out[["pop"]] <- pop

  # chgnpt chg
  if ("chg" %in% abr_mtd){

    table_chg <- shifts_res$res_table %>% dplyr::filter(abbr=="chg")
    loc_brk_chg <- table_chg %>% dplyr::pull(loc_brk) %>% stringr::str_split(";") %>% `[[`(1) %>% as.numeric()

    if ("loc_aux1_chg" %in% names(best_traj)) {
      loc_aux_chg <- c(best_traj$loc_aux1_chg, best_traj$loc_aux2_chg)
    }

    bp_chg <- pop +
      geom_vline(xintercept = loc_brk_chg, col="blue", linetype="dashed")+
      geom_line(data = shifts_res$chg_outlist$pred_chg, aes(x=year, y=bp), col = "blue", alpha=0.7)+
      scale_colour_manual(values = rep("red", table_chg$n_brk+1))+
      theme(legend.position = "none") %>%
      suppressWarnings()

    if ("loc_aux1_chg" %in% names(best_traj)) {

      bp_chg <- bp_chg +
        geom_vline(xintercept = loc_aux_chg[!is.na(loc_aux_chg)], col="blue", linetype="dotted", alpha=0.5)
    }

    # asdetect asd
    if ("asd" %in% abr_mtd){

      table_asd <- shifts_res$res_table %>% filter(abbr=="asd")
      loc_brk_asd <- table_asd %>% pull(loc_brk) %>% str_split(";") %>% `[[`(1) %>% as.numeric()

      if (!is.na(loc_brk_asd)){

        bp_chg <- bp_chg +
          geom_vline(xintercept = loc_brk_asd, col="red", linetype="dashed")

      }
    }
  }

  # Write title:
  if (!is.null(best_traj_loo)){
    bp_chg <- bp_chg +
      ggtitle(paste0("<b>Abrupt ", shifts_res$res_table$trend[2],"</b>",
                     "<br>AIC = ", round(table_chg$aic, digits=2)," ",
                     " AIC weight = ", best_traj[["weight_aic_abrupt"]]," ",
                     " LOO = ", round(best_traj[["loo_abrupt"]], digits=3),
                     "<br>Breaktime(s): <span style='color:blue'>",table_chg$loc_brk,"</span>",
                     if(!is.null(best_traj$loc_aux1_chg)) {
                       paste0("(<span style='color:lightblue'>",best_traj$loc_aux1_chg,"</span>",
                              if(!is.null(best_traj$loc_aux2_chg)){
                                paste0(", <span style='color:lightblue'>",best_traj$loc_aux2_chg,"</span>")

                              },")")
                     },
                     if(!is.na(loc_brk_asd)) {paste0(", <span style='color:red'>",loc_brk_asd,"</span>")},
                     "  Magnitude = ", round(table_chg$mag, digits=3)))+
      theme(plot.title = ggtext::element_markdown(size=10, lineheight = 1.1))

  } else {
    bp_chg <- bp_chg +
      ggtitle(paste0("<b>Abrupt ", shifts_res$res_table$trend[2],"</b>",
                     "<br>AIC = ", round(table_chg$aic, digits=2)," ",
                     " AIC weight = ", best_traj[["weight_aic_abrupt"]],
                     "<br>Breaktime(s): <span style='color:blue'>", table_chg$loc_brk, "</span> ",
                     if(!is.null(best_traj$loc_aux1_chg)) {
                       paste0("(<span style='color:lightblue'>",best_traj$loc_aux1_chg,"</span>",
                              if(!is.null(best_traj$loc_aux2_chg)){
                                paste0(", <span style='color:lightblue'>",best_traj$loc_aux2_chg,"</span>")
                              },")")
                     },
                     if(!is.na(loc_brk_asd)) {
                       paste0(", <span style='color:red'>",loc_brk_asd,"</span>")
                     },
                     "  Magnitude = ", round(table_chg$mag, digits=3)))+
      theme(plot.title = ggtext::element_markdown(size=10, lineheight = 1.1))
    bp_chg
  }

  # Set background colour:
  if (best_traj$class == "abrupt"){

    if("expected_class" %in% names(best_traj)){

      if (best_traj$class == best_traj$expected_class) { bkg_col <- "lightgreen"
      } else { bkg_col <- "lightcoral" }

    } else { bkg_col <- "lightblue" }

  } else { bkg_col <- "grey80" }

  bp_chg <- bp_chg + theme(plot.background = element_rect(fill = bkg_col))

  plots_out[["bp_chg_asd"]] <- bp_chg


  # # asdetect asd integrated before
  # if ("asd" %in% abr_mtd){
  #
  #   table_asd <- shifts_res$res_table %>% filter(abbr=="asd")
  #   loc_brk_asd <- table_asd %>% pull(loc_brk) %>% str_split(";") %>% `[[`(1) %>% as.numeric()
  #
  #   if (!is.na(loc_brk_asd)){
  #
  #     bp_asd <- bp_chg +
  #       geom_vline(xintercept = loc_brk_asd, col="red", linetype="dashed")
  #
  #     plots_out[["bp_chg_asd"]] <- bp_asd
  #
  #   }
  # }


  return(plots_out)
}


# Plot reproducible examples of time series in all noise combinations in noise_df and
# for all scenarios in param_df

plot_ts_exp <- function(param_df, noise_df, length, seed){

  noise_names <- colnames(noise_df)
  param_names <- colnames(param_df)
  scen_fct <- "make_scen"
  jfr <- 0.2

  # for all scenarios:
  for (j in 1:nrow(param_df)){

    gl <- list()
    # import parameter values:
    for (i in 1:length(param_names)){

      assign(param_names[i], param_df %>%
               dplyr::pull(param_names[i]) %>%
               `[`(j))
    }

    # make mortality scenario:
    if (scen_fct == "make_scen"){
      Fseq <- make_scen(shape=shape, trend=trend, min=min, max=max,
                        velocity=velocity, l=l, breaktime = breaktime)
    }

    for (n in 1:nrow(noise_df)){

      # set noise parameters:
      for (m in 1:length(noise_names)){

        assign(noise_names[m], noise_df %>%
                 dplyr::pull(noise_names[m]) %>%
                 `[`(n))
      }

      set.seed(seed)
      sim <- run_simu(name=Fseq[[2]], r=r, F=Fseq[[1]], Ts=length(Fseq[[1]]),
                      K=10, P=2, H=H, sr=sr, se=se, su=su, iter=1,
                      thr=0, init=0.9, expected_class=expected_class,
                      jfr=jfr, jsz=jsz)

      if(length=="l50") {
        sim$df <- sim$df %>%
          dplyr::group_by(scen, iter) %>%
          dplyr::filter(year %% 2 == 1) %>%
          dplyr::mutate(year=seq(1:50),
                        scen = sub("l100", "l50", scen),
                        scen = sub("brk50", "brk25", scen)) %>%
          dplyr::ungroup()
        Fseq[[2]] <- sub("l100", "l50", Fseq[[2]])
      }

      if(length=="l20"){
        sim$df <- sim$df %>%
          dplyr::group_by(scen, iter) %>%
          dplyr::filter(year %% 5 == 1) %>%
          dplyr::mutate(year=seq(1:20),
                        scen = sub("l100", "l20", scen),
                        scen = sub("brk50", "brk10", scen)) %>%
          dplyr::ungroup()
        Fseq[[2]] <- sub("l100", "l20", Fseq[[2]])
      }

      TB_plot <- ggplot(sim$df, aes(x=year, y=TB, color=as.factor(iter), group=iter))+
        geom_line(col="red", alpha=1)+
        scale_size_manual(values = c("TRUE"=2, "FALSE"=NA))+
        # scale_y_continuous(
        #   name = "Total biomass")+
        expand_limits(y = c(0, 12.5))+
        theme_light()+
        theme(legend.position = "none",
              axis.title.x = element_blank()
              # ,axis.title.y = element_text(color = "red", size=13)
        )+
        labs(caption = paste0("r = ", r,
                              " ; H = ", H,
                              " ; P = ", 2,
                              " ; F range = [", head(Fseq[[1]],1),
                              ", ", tail(Fseq[[1]],1), "]",
                              ", sr = ", sr,
                              ", se = ", se,
                              ", jfr = ", jfr,
                              ", jsz = ", jsz
        ))

      if(n==1) TB_plot <- TB_plot + ggtitle("se = 0.001 & sr = 0.001")
      if(n==2) TB_plot <- TB_plot + ggtitle("sr = 0.001 & se = 0.025") +
        scale_y_continuous(name = "jump = 0*se")
      if(n %in% c(3,4)) TB_plot <- TB_plot + scale_y_continuous(name = paste0("jump = ",jsz,"*se"))
      if(n %in% c(5,8,11)) TB_plot <- TB_plot + ggtitle(paste0("sr = ", sr, " & se = 0.025"))

      gl <- c(gl, list(TB_plot))
    }

    p <- cowplot::plot_grid(plotlist = gl[c(1,2,5,8,11)], nrow=1) / cowplot::plot_grid(plotlist = c("",gl[c(3,6,9,12)]), nrow=1) / cowplot::plot_grid(plotlist = c("",gl[c(4,7,10,13)]), nrow=1)

    cowplot::save_plot(filename = paste0("res/03_ricker_library/simu_library/ts_exp/",Fseq[[2]],"_noise_lvl_seed",seed,".png"),
                       plot=p, base_height = 8, base_width = 22)

    print(paste0("scenario: ", j, "/", nrow(param_df)))
  }

  return(invisible(NULL))
}


