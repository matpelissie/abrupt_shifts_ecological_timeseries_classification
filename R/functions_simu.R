###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###
#
# 22/03/2022 mathieu.pelissie@ens-lyon.fr
#
# Functions for stock dynamics simulations
#
# functions_simu.R
#
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###


# I) Simulate data -----------------------------------------------------

#' Get csv data
#'
#' @param file relative path and file name with noise scenarios
#'
#' @return csv with noise or simulation paramaters for scenarios in row
#'
#' @export

get_data <- function(file){
  readr::read_csv(file)
}



#' Make mortality scenario
#'
#' @param shape shape of the timeseries to generate, either "cst" (constant),
#' "lin" (linear), "qdr" (quadratic monotonous), "ccv" (concave),
#' "cvx" (convex), or "abt" (abrupt).
#' @param l length of the timeseries
#' @param trend (lin, qdr, abt) character "pos" (positive) or "neg" (negative)
#' @param min (lin, qdr, abt) minimal value
#' @param max (lin, qdr, abt) maximal value
#' @param velocity (qdr) "acc" (accelerating) or "dec" (decelerating)
#' @param breaktime (abt) timestep at which the break occurred
#'
#' @return vector with the desired features and string with the parameters
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
        name <- sprintf("l%g_%s_%s_%s_X%g-%g",
                        l, shape, trend, velocity, min, max)
      }

      if(trend=="neg"){
        a <- (max-min)/l^2
        b <- 2*(min-max)/l
        c <- max
        name <- sprintf("l%g_%s_%s_%s_X%g-%g",
                        l, shape, trend, velocity, max, min)
      }
    }

    if(velocity=="acc"){

      if(trend=="pos"){
        a <- (max-min)/l^2
        b <- 0
        c <- min
        name <- sprintf("l%g_%s_%s_%s_X%g-%g",
                        l, shape, trend, velocity, min, max)
      }

      if(trend=="neg"){
        a <- (min-max)/l^2
        b <- 0
        c <- max
        name <- sprintf("l%g_%s_%s_%s_X%g-%g",
                        l, shape, trend, velocity, max, min)
      }

    }

    pol2deg <- function(x) a*x^2+b*x+c
    scen <- pol2deg(1:l)


  } else if(shape=="abt") { # Shape abrupt

    if(breaktime>l) stop("Break point should occur within the timeseries")

    if(trend=="pos"){
      scen <- c(rep(min, breaktime-1), rep(max, l-breaktime+1))
      name <- sprintf("l%g_%s_%s_brk%g_X%g-%g",
                      l, shape, trend, breaktime, min, max)
    }

    if(trend=="neg"){
      scen <- c(rep(max, breaktime-1), rep(min, l-breaktime+1))
      name <- sprintf("l%g_%s_%s_brk%g_X%g-%g",
                      l, shape, trend, breaktime, max, min)
    }
  }

  return(list("scen"=scen, "name"=name))

}



#' Make mortality scenario with multiple breakpoints
#'
#' @param l length of the timeseries
#' @param trend character "ccv" (concave), "cvx" (convex),
#' "pos" (positive), or "neg" (negative)
#' @param val* value of mortality for the first (val1), middle (val2),
#' and last (val3) steps.
#' @param brk* timestep at which the first (brk1) and second (brk2) breaks
#' occurred
#'
#' @return vector with the desired features and string with the parameters
#'
#' @export

make_scen_multibrk <- function(l=100, trend=NULL,
                      val1=NULL, val2=NULL, val3=NULL,
                      brk1=NULL, brk2=NULL){

  if(brk1>l | brk2>l) stop("Breakpoints should occur within the timeseries")

  shape <- "abt"
  scen <- c(rep(val1, brk1), rep(val2, brk2-brk1), rep(val3, l-brk2))
  name <- sprintf("l%g_%s_%s_brk%g-%g_X%g-%g-%g",
                  l, shape, trend, brk1, brk2, val1, val2, val3)

  return(list("scen"=scen, "name"=name))
}



#' Make simulations and plots
#'
#' @param name character for scenario name
#' @param Ts length of timeseries
#' @param r vector of growth rates (same length as timeseries)
#' @param F vector of mortality rates (same length as timeseries)
#' @param K carrying capacity
#' @param P exponent for sigmoid functional response
#' @param H half-saturation rate
#' @param sr demographic stochasticity
#' @param se environmental stochasticity
#' @param su observational noise (i.e. measurement error), not used
#' @param iter number of repetitions of timeseries
#' @param thr unfeasibility threshold
#' @param init percentage of carrying capacity as initial size
#' @param trs length of the transient period (to be removed, default = 50 steps)
#' @param expected_class expected trajectory class
#' @param jfr frequency of jump events
#' @param jsz base size of the jump event (buffered by demographic stochasticity
#' and population size)
#'
#' @return list of two elements:
#' - list of matrix with time varying inputs and outputs for all iterations
#' - data frame with time varying inputs and outputs for all iterations
#'
#' @export

run_simu <- function(name, Ts, r, F, K, P, H, sr, se, su,
                     iter, thr, init, trs=50, expected_class, jfr, jsz){

  # if and only if F is variable:
  full_name <- sprintf("%s_r%#.2g_H%g", sub("X", "F", name), r, H)

  if (length(r)==1) r <- rep(r,Ts)
  if (length(H)==1) H <- rep(H,Ts)

  if (length(r)!=length(F) | length(F)!=length(H) | length(F)!= Ts){
    stop("Parameters should be vectors the same length as timeseries length Ts")
  }

  TBts <- matrix(NA, nrow=Ts, ncol=iter)
  TCts <- matrix(NA, nrow=Ts, ncol=iter)
  rts <- matrix(NA, nrow=Ts, ncol=iter)
  Fts <- matrix(NA, nrow=Ts, ncol=iter)
  Hts <- matrix(NA, nrow=Ts, ncol=iter)
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

    # First timepoint with stochasticity:
    TBts[1, i] <- x*exp(su*rnorm(1))

    for (t in 2:Ts){

      r1 <- r[t]*exp(rnorm(n=1, sd=sr)) # demographic stochasticity
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


      if (x<thr){

        extra <- extra+1
        break

      } else {

        TBts[t, i] <- x*exp(rnorm(n=1, sd=su)) # observation error (not used)
        TCts[t, i] <- TC
        rts[t, i] <- r1
        Fts[t, i] <- F[t]
        Hts[t, i] <- H[t]
        srts[t, i] <- sr
        sets[t, i] <- se
        suts[t, i] <- su

      }

    }

    # Is constant timeseries not collapsed:
    if (!expected_class %in% c("abrupt") &
        x>0.1*TBts[1, i] & TBts[90, i]>0.1*TBts[1, i]) i <- i+1

    # Is timeseries complete:
    else if (expected_class %in% c("abrupt") & x>thr) i <- i+1

    # To escape if not possible to complete the timeseries:
    if (extra>2*max(iter,50)){
      return(warning("Warning: too many simulations fail to avoid extinction,
                     check parameters and noise"))
    }
  }

  lib_ts <- data.frame(scen = as.vector(full_name),
                       iter = rep(1:iter, each=Ts),
                       year = rep(1:Ts, iter),
                       TB = as.vector(TBts),
                       r = as.vector(rts),
                       F = as.vector(Fts),
                       H = as.vector(Hts),
                       P = rep(P, Ts*iter),
                       TC = as.vector(TCts),
                       sr = as.vector(srts),
                       se = as.vector(sets),
                       su = as.vector(suts)
                       )

  scen_mat <- list("TB"=TBts, "TC"=TCts, "r"=rts, "F"=Fts)

  return(list(
    "scen_mat"=scen_mat,
    "df"=lib_ts
    ))

}



#' Make and save simulated data
#'
#' @param param_df data frame with scenarios parameters
#' @param scen_fct character to specify the scenario-making function to use
#' @param sr demographic stochasticity
#' @param se environmental stochasticity
#' @param su observational noise (i.e. measurement error), not used
#' @param jfr frequency of jump events
#' @param jsz base size of the jump event (buffered by demographic stochasticity
#' and population size)
#'
#' @return data frame with simulated data for all scenarios
#'
#' @export

make_store_simu <- function(param_df, scen_fct, sr, se, su, jfr, jsz){

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

  return(all_simu)

}



# II) Bifurcation plot --------------------------------------------------------

#' Give the state at the equilibrium for a given timeseries
#'
#' @param lib_ts data frame with inputs and outputs for all iterations
#' @param thr unfeasibility threshold (default 0)
#'
#' @return character indicating state at equilibrium (unfeasible, stable,
#' cycles, chaos)
#'
#' @export

eq_state <- function(lib_ts, thr=0){

  sample_length <- nrow(lib_ts)/10

  # Unique lag-1 difference of the last 10% timepoints:
  unique_diff <- lib_ts %>%
    tail(n=sample_length) %>%
    dplyr::mutate(diff_lag1=TB-dplyr::lag(TB)) %>%
    tidyr::drop_na(diff_lag1) %>%
    dplyr::distinct(diff_lag1) %>%
    dplyr::pull()

  # Determine equilibrium state based on a series of criteria:
  state <- NA

  if (lib_ts %>%
      dplyr::pull(TB) %>%
      tail(1) < thr) { # last value below unfeasibility threshold
    state <- "unfeasible"

  } else if (unique_diff %>% length() == 1 | # fully stabilized
             max(abs(unique_diff)) < 10^(-5) | # if oscillating very little
             unique_diff %>% sign() %>% unique() %>%
             length() == 1) { # if not yet stabilized
    state <- "stable"

  } else if (unique_diff %>% length() < (sample_length-1) | # regular cycles
             sign(abs(unique_diff[-1]) -
                  abs(stats::lag(unique_diff, k=2)[-1])) %>% unique() %>%
             length() == 1) { # if lag-2 cycles of decreasing amplitude
    state <- "cycle"

  } else if (unique_diff %>%
             length() == sample_length-1) { # all lag-1 diff values different
    state <- "chaos"
  }

  return(state)

}



#' Plot stable states and evidence potential alternative stable states
#'
#' @param Ts length of timeseries
#' @param rval value of growth rate
#' @param P exponent for sigmoid functional response
#' @param H half-saturation rate
#' @param K carrying capacity
#' @param Fmax the maximal mortality rate
#' @param Fstep the interval between two values of mortality rate
#' @param TBhigh the starting value for high biomass
#'
#' @return table and ggplot with values of stable states along a range of F
#'
#' @export

ASS_row <- function(Ts, rval, P, H, K, Fmax, Fstep, TBhigh){

  Fseq <- seq(0, Fmax, Fstep)

  row_inc <- row_dec <- states_inc <- states_dec <- stab_inc <- stab_dec <-
    data.frame(matrix(NA, length(rval), length(Fseq))) %>%
    `rownames<-`(paste0("r", rval)) %>%
    `colnames<-`(paste0("F", Fseq))

  # Find equilibrium for increasing F & high starting biomass:
  init <- TBhigh
  for (j in Fseq){

    F <- rep(j, Ts)
    r <- rep(rval, Ts)
    simu <- run_simu(name="test", Ts, r, F, K, P, H=H, sr=0, se=0, su=0, iter=1,
                     thr=0, init,expected_class="abrupt", jfr=0, jsz=0)[[2]] %>%
      try() # expected_class abrupt to avoid removing collapsed timeseries

    if(class(simu) !="try-error"){

      states_inc[paste0("r",rval), paste0("F",j)] <- eq_state(simu, thr=0)

      # take mean of last 10% timepoints:
      init <- mean(simu[(nrow(simu)*0.9):nrow(simu), "TB"])/K
      init <- max(init, 10^(-5)) # to avoid reaching zero

      row_inc[paste0("r", rval), paste0("F", j)] <- init*K
      stab_inc[paste0("r", rval), paste0("F", j)] <- simu %>%
        dplyr::pull("EIG") %>%
        tail(1)

    } else {
      states_inc[paste0("r",rval), paste0("F",j)] <- "unfeasible"
      row_inc[paste0("r", rval), paste0("F", j)] <- NA
    }

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

  # Find equilibrium for decreasing F & low starting biomass:
  for (j in rev(Fseq)){

    F <- rep(j, Ts)
    r <- rep(rval, Ts)
    simu <- run_simu(name=NA, Ts, r, F, K, P, H=H, sr=0, se=0, su=0, iter=1,
                     thr=0, init,expected_class="abrupt", jfr=0, jsz=0)[[2]] %>%
      try()

    if(class(simu) !="try-error"){

    states_dec[paste0("r",rval), paste0("F",j)] <- eq_state(simu, thr=0)

    # take mean of last 10% timepoints:
    init <- mean(simu[(nrow(simu)*0.9):nrow(simu), "TB"])/K
    init <- max(init, 10^(-5)) # to avoid reaching zero

    row_dec[paste0("r", rval), paste0("F", j)] <- init*K
    stab_dec[paste0("r", rval), paste0("F", j)] <- simu %>%
      dplyr::pull("EIG") %>% tail(1)

    } else {
      states_dec[paste0("r",rval), paste0("F",j)] <- "unfeasible"
      row_dec[paste0("r", rval), paste0("F", j)] <- NA
    }

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

  # Compare equilibria:
  bif_wide <- row_inc_long %>%
    dplyr::left_join(row_dec_long, by=c("r", "F")) %>%
    dplyr::left_join(states_inc_long, by=c("r", "F")) %>%
    dplyr::left_join(states_dec_long, by=c("r", "F")) %>%
    dplyr::left_join(stab_inc_long, by=c("r", "F")) %>%
    dplyr::left_join(stab_dec_long, by=c("r", "F")) %>%
    # arbitrary threshold of 10^(-1) for equal equilibrium value:
    dplyr::mutate(ASS = ifelse(abs(value_inc - value_dec)>10^(-1),
                               TRUE, FALSE)) %>%
    dplyr::mutate(diffstate = ifelse(states_inc != states_dec,
                                     TRUE, FALSE))


  bif_long <- bif_wide %>%
    tidyr::pivot_longer(cols = contains("_"),
                        names_to = c(".value","sense"),
                        names_pattern = "(\\w+)_(\\w+)")


  bif_plot <- ggplot(bif_long, aes(x=F, y=value, color=as.factor(sense)))+
    geom_line()+
    geom_point(size=0.5)+
    labs(y="TB at equilibrium", color="Sense")+
    theme_light()+
    expand_limits(y = 0)+
    theme(legend.position = c(0.9,0.8))


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

  return(list("table" = bif_wide, "plot" = bif_stab))

}



#' Plot system stable states under a range of constant r and F
#'
#' @param Ts length of timeseries
#' @param K carrying capacity
#' @param rmax the maximal growth rate value (default 3)
#' @param rstep the interval between two values of growth rate (default 0.1)
#' @param Fmax the maximal mortality rate value (default 3)
#' @param Fstep the interval between two values of mortality rate (default 0.1)
#' @param P exponent for sigmoid functional response
#' @param H half-saturation rate
#' @param thr extinction threshold
#' @param init percentage of carrying capacity as initial size
#'
#' @return list of plots with stable states and along ranges of r and F and
#' the data frame used for the plots
#'
#' @export

bifurc_analysis <- function(Ts, K, rmax=3, rstep=0.1,
                            Fmax=3, Fstep=0.1, P, H, thr, init){

  rseq <- seq(0, rmax, rstep)
  Fseq <- seq(0, Fmax, Fstep)
  ASS <- data.frame()

  # Find equilibrium states for every growth rate value:
  for (i in rseq){

    ASS <- ASS %>%
      rbind(ASS_row(Ts=Ts, rval = i, Fmax=Fmax, K=K,
                    Fstep=Fstep, H=H, P=P, TBhigh=init)$table)

    print(paste0(i, "/", rmax))
  }

  # Plot types and values of equilibrium state:
  bif_ana <- ASS %>%
    dplyr::mutate(ASS = ifelse(ASS, TRUE, NA))

  bif_ana$states_inc <- factor(bif_ana$states_inc,
                          levels = c("unfeasible","stable","cycle","chaos"))

  bif_ana$states_dec <- factor(bif_ana$states_dec,
                               levels=c("unfeasible","stable","cycle","chaos"))
  state.col <- c(unfeasible="#c44601", stable="#5ba300",
                 cycle="#0073e3", chaos="#e79b01")

  # Types of equilibrium and ASS when increasing F:
  plot_inc <- ggplot(bif_ana, aes(x=F, y=r))+
    geom_tile(aes(fill=states_inc, alpha=value_inc))+
    scale_fill_manual(values = state.col)+
    scale_alpha(range=c(0.3,1), limits=c(0,10), na.value=1)+
    geom_point(data=subset(bif_ana, !is.na(ASS) &
                             !states_inc %in% c("unfeasible", "chaos")),
               # display ASS only for stable equilibrium
               aes(x=F, y=r), size=1)+
    scale_color_manual(values = c("TRUE"="black"))+
    theme_minimal()+
    labs(caption = paste0("h = ", H, " ;  init TB inc = ",
                          init*100, "% K ; increasing F"))


  # Types of equilibrium and ASS when decreasing F:
  plot_dec <- ggplot(bif_ana, aes(x=F, y=r))+
    geom_tile(aes(fill=states_dec, alpha=value_dec))+
    scale_fill_manual(values = state.col)+
    scale_alpha(range=c(0.3,1), limits=c(0,10))+
    geom_point(data=subset(bif_ana, !is.na(ASS) &
                             !states_dec %in% c("unfeasible", "chaos")),
               # display ASS only for stable equilibrium
               aes(x=F, y=r), size=1)+
    scale_color_manual(values = c("TRUE"="black"))+
    theme_minimal()+
    labs(caption = paste0("h = ", H, " ;  init TB inc = ",
                          init*100, "% K ; decreasing F"))


  return(list("plot_inc" = plot_inc,
              "plot_dec" = plot_dec,
              "bif_ana" = bif_ana))
}



# III) Make plots --------------------------------------------------------



#' Plot simulations in a simple way
#'
#' @param lib_ts dataframe with simulated timeseries
#' @param tbcolor line color
#' @param alpha line transparency
#'
#' @return minimal plot with TB curve
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


