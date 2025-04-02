# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library("logger"))

f_log <- file.path("log.txt")
log_appender(appender_file(f_log))
# message(Sys.time(), " Log file initialised ", f_log)
log_info("*** START UP ***")

suppressWarnings(suppressPackageStartupMessages(library(cmdstanr)))
library(ggplot2)
library(parallel)
library(mcprogress)
suppressWarnings(suppressPackageStartupMessages(library(tictoc)))
library(poisson)
library(data.table)
library(survival)
library(flexsurv)
library(mvnfast)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_02"
} else {
  log_info("Run method ", args[1])
}

tic()

set.seed(4682467)

mc_cores <- 50
mod_01 <- cmdstanr::cmdstan_model("./stan/log-logistic-aft-01.stan")
output_dir_mcmc <- paste0(getwd(), "/tmp")
output_dir_res <- paste0(getwd(), "/out")


# function to build a default return object in the case of exceptions or
# warnings which might otherwise stop the simulation
#
get_default_return_obj <- function(N){

  d_pars <- data.table(
    id_analy = seq_along(N),
    N_enrol = N,
    p0_sim = NA_real_,
    p1_sim = NA_real_,
    p_sup = NA_real_,
    p_fut = NA_real_,
    mu_rmst0 = NA_real_,
    mu_rmst1 = NA_real_,
    mu_drmst = NA_real_,
    # median tte
    scale0 = NA_real_,
    scale1 = NA_real_,
    # linear predictor parameters
    gamma1 = NA_real_,
    gamma2 = NA_real_,
    shape = NA_real_,
    p0_360 = NA_real_,
    p1_360 = NA_real_,
    rd_360 = NA_real_,
    # probability thresholds
    pr_drmst_gt_0 = NA_real_,
    pr_gamma2_gt_0 = NA_real_
  )

  d_pars
}


integrand_surv <- function(x, mu, shape) {
  a = exp(mu)
  S = 1 / (1 + (x/a)^shape)
  S
}

# id_sim = 1
# id_scen = 1
# N = c(1500, 2000, 2500, 3000)
# # p_0, p_1
# p = c(0.1, 0.07)

do_trial <- function(
    id_sim = 1,
    N = c(1500, 2000, 2500, 3000),
    # probability of event by day 360
    p = c(0.1, 0.07),
    p_thresh_sup = 0.95,
    p_thresh_fut = 0.30,
    # 1 = flexsurvreg, 2 = pathfinder, 3 = mcmc
    estimation_method = 1
){

  log_info(paste0(match.call()[[1]]), " sim ", id_sim)

  # Number of analyses (will always be 2, but anyway)
  K <- length(N)
  N_tot <- max(N)
  # Cohorts
  N_c <- c(N[1], diff(N))
  # Start and end indexes
  is <- c(0, cumsum(N_c)[-length(N_c)]) + 1
  ie <- cumsum(N_c)

  # For simplicity just start the follow up clock at the time of birth.
  t_0 <- c(0, poisson::nhpp.event.times(
    rate = 5.5,
    N_tot - 1,
    prob.func = function(t){pmin(t/360, 1)}))

  # All data
  d <- data.table(
    id = 1:N_tot,
    trt = rep(0:1, length = N_tot),
    t_0 = t_0
  )

  d_pars <- get_default_return_obj(N)

  # Construct log-logistic event distributions assuming a
  # cumulative incidence of p[1] \approx 10% at day 360 for the control RSVpreF
  # group.

  # The hazard function of the log-logistic is
  # hump-shaped and is a bit like a log normal but with longer tails.
  # It initially increases, reaches a maximum and then decreases
  # toward 0 as lifetimes become larger and larger.
  # see flexsurv::dllogis
  # This isn't exactly right, but it will probably do ok over the first 360 days
  # of life.
  # Unlike lognormal, loglistic has the benefit of a closed form hazard function.

  # Under the following parameterisation, the interpretation of the parameters
  # is the same as in the standard Weibull distribution (i.e. not the ph version).
  # We normally model on the scale parameter, i.e. \alpha in the following:

  # f = \frac{(b/a)(t/a)^{b-1}}{(1 + (t/a)^b)^2}
  # S = \frac{1}{1 + (t/a)^b}

  # Further details here:
  # https://maj-biostat.github.io/notebooks/log-logistic-aft-in-stan.html

  # We want a 10% event rate at 360 days,
  # S(360) = 1 - 0.1 = 0.9, so choose a value for b = 2 and solve for a

  # S(360) = \frac{1}{1 + (360/a)^b} = 0.9
  # 1/0.9      = 1 + (360/a)^b
  # (1/0.9) - 1 = (360/a)^b
  # (0.1/0.9)^{1/b}  =  360/a
  # a = 360 / (0.1/0.9)^{1/b}
  # a = 1080

  # Generally
  # a = \frac{t}{ ( 1/(1-p)  - 1 )^{1/b} }
  # where p is the desired event rate.

  # Where alpha is the shape parameter controling whether the hazard will
  # increase or fall over time and where alpha > 1 implies increasing hazard.
  # Maintain the same shape parameter for both groups.
  b <- 2

  a1 <- 360 / (  (1/(1-p[1])) - 1  )^(1/b)
  a2 <- 360 / (  (1/(1-p[2])) - 1  )^(1/b)

  # Generate event times
  d[, t_evt := fifelse(
    trt == 0,
    rllogis(.N, shape = b, scale = a1),
    rllogis(.N, shape = b, scale = a2))]

  # log_info("median0 ", a1)
  # log_info("median1 ", a2)
  # surv 360
  # 1 / (1 + (360/a1)^b)
  # 1 / (1 + (360/a2)^b)
  # log_info("s0_360 ", 1 / (1 + (360/a1)^b))
  # log_info("s1_360 ", 1 / (1 + (360/a2)^b))

  ## Trial loop --------
  continue_trial <- TRUE
  ii <- 1

  while(continue_trial){

    # Pick up the current cohort
    d_c <- d[1:N[ii], ]

    # Determine how much follow up time for each of the units in
    # the cohort so that we can work out what events have happened and whether
    # we observe them or not (remember that censoring occurs at 360 days)

    # If we are at the final analysis then we assume that we have waited the
    # 360 days to gather the full followup, i.e. the last person enrolled will have
    # a follow up of 360 days.
    if(ii == length(N)){
      d_c[, t_fu := max(d_c$t_0) + 360 - t_0]
    } else {
      d_c[, t_fu := max(d_c$t_0) - t_0]
    }

    # In order to see the event occur we need to have followed them for longer
    # than the time to event but if the event occurs after 360 days then we
    # will not see it.
    d_c[, evt := as.numeric(t_fu >= t_evt & t_evt <= 360)]

    # Ensure that the time we use in the analysis stops at the 360 days
    d_c[, t_evt_obs := copy(t_evt)]
    # Censored events are set to 360 days
    d_c[evt == 0, t_evt_obs := 360]

    # Survival approach - log-logistic regression

    # Prepare data for Stan
    ld <- list(
      N = nrow(d_c),
      P = 2,
      X = cbind(1, d_c$trt),
      y = d_c$t_evt_obs,
      event = d_c$evt,
      N_pred = 361,
      t_surv = 0:360,
      # hyperparameters
      mu0_gamma = c(5, 0),
      sd0_gamma = c(2, 2),
      rho_shape = 0.5
    )

    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-seq-sim", id_sim, "-intrm", ii)

    # MCMC fit - sink to remove the noise
    # snk <- capture.output(

    if(estimation_method == 1){

      fs1 <- flexsurvreg(Surv(t_evt_obs, evt) ~ trt, data = dd, dist = "llogis")

      d_post <- data.table(rmvn(1000, fs1$opt$par, fs1$cov))
      names(d_post) <- names(fs1$opt$par)
      setcolorder(d_post, c("scale", "trt1", "shape"))
      d_post[, shape := exp(shape)]
      names(d_post) <- c(paste0("gamma", 1:2), "shape")
      d_post[, scale0 := exp(gamma1)]
      d_post[, scale1 := exp(gamma1 + gamma2)]

      d_post[, p0_360 := pllogis(360, shape = shape, scale = scale0)]
      d_post[, p1_360 := pllogis(360, shape = shape, scale = scale1)]

    } else if(estimation_method %in% 2:3){

      if(estimation_method == 2){
        # For resolving high pareto values see
        # https://users.aalto.fi/~ave/casestudies/Birthdays/birthdays.html
        m1 <- mod_01$pathfinder(
          ld,
          init = function() {list(
            gamma = c(runif(1, 5, 10), runif(1, -1, 1)),
            shape = runif(1, 0, 4)
          )},
          num_paths=4, single_path_draws=250,
          history_size=50, max_lbfgs_iters=100,
          refresh = 0, draws = 1000,
          output_dir = output_dir_mcmc,
          output_basename = foutname)

      } else if(estimation_method == 3){

        m1 <- mod_01$sample(
          ld,
          init = function() {list(
            gamma = c(runif(1, 5, 10), runif(1, -1, 1)),
            shape = runif(1, 0, 4)
          )},
          iter_warmup = 1000, iter_sampling = 1000,
          parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = T,
          max_treedepth = 10, adapt_delta = 0.95,
          output_dir = output_dir_mcmc,
          output_basename = foutname)
      }

      d_post <- data.table(
        m1$draws(
          variables = c("gamma", "shape",
                        "scale0", "scale1",
                        "p0_360", "p1_360"),
          format = "matrix"))

      names(d_post) <- c(paste0("gamma", 1:2),
                         "shape", "scale0", "scale1",
                         "p0_360", "p1_360")
    }


    m_rmst <- matrix(NA, ncol = 2, nrow = nrow(d_post))
    for(vv in 1:nrow(d_post)){
      m_rmst[vv, 1] <- integrate(
        integrand_surv, lower = 0, upper = 360,
        mu = d_post$gamma1[vv],
        shape = d_post$shape[vv])$value
      m_rmst[vv, 2] <- integrate(
        integrand_surv, lower = 0, upper = 360,
        mu = d_post$gamma1[vv] + d_post$gamma2[vv],
        shape = d_post$shape[vv])$value
    }

    d_rmst <- data.table(m_rmst)
    names(d_rmst) <- paste0("trt", 0:1)

    d_pars[id_analy == ii, `:=`(
      N_enrol = N[ii],
      p0_sim = p[1],
      p1_sim = p[2],
      p_sup = p_thresh_sup,
      p_fut = p_thresh_fut,
      mu_rmst0 = mean(d_rmst$trt0),
      mu_rmst1 = mean(d_rmst$trt1),
      mu_drmst = mean(d_rmst$trt1 - d_rmst$trt0),
      scale0 = mean(d_post$scale0),
      scale1 = mean(d_post$scale1),
      gamma1 = mean(d_post$gamma1),
      gamma2 = mean(d_post$gamma2),
      shape = mean(d_post$shape),
      p0_360 = mean(d_post$p0_360),
      p1_360 = mean(d_post$p1_360),
      rd_360 = mean(d_post$p1_360 - d_post$p0_360),
      pr_drmst_gt_0 = mean((d_rmst$trt1 - d_rmst$trt0)>0),
      pr_gamma2_gt_0 = mean((d_post$gamma2)>0)
    )]

    # test for superiority and futility conditionss
    if(d_pars[id_analy == ii, pr_drmst_gt_0 > p_thresh_sup] |
       d_pars[id_analy == ii, pr_drmst_gt_0 < p_thresh_fut]
    ){
      continue_trial = F
    } else {
      ii <- ii + 1
    }

    if(ii > length(N)){
      continue_trial = F
    }


    # m_surv0 <- matrix(NA, ncol = 360, nrow = nrow(d_post))
    # m_surv1 <- matrix(NA, ncol = 360, nrow = nrow(d_post))
    # for(vv in 1:nrow(d_post)){
    #
    #   m_surv0[vv, ] <- 1/ (1 + (1:360/d_post$scale0[vv])^d_post$shape[vv])
    #   m_surv1[vv, ] <- 1/ (1 + (1:360/d_post$scale1[vv])^d_post$shape[vv])
    #
    # }
    # d_fig <- rbind(
    #   cbind(trt = 0, data.table(m_surv0)),
    #   cbind(trt = 1, data.table(m_surv1))
    # )
    # d_fig <- melt(d_fig, id.vars = "trt")
    # d_fig[, day := as.numeric(gsub("V", "", variable))]
    #
    # d_fig <- d_fig[, .(
    #   mu_surv = mean(value),
    #   q_025 = quantile(value, prob = 0.025),
    #   q_975 = quantile(value, prob = 0.975)
    # ), keyby = .(trt, day)]
    # d_fig[, trt := factor(trt)]
    #
    # ggplot(d_fig, aes(x = day, y = mu_surv, group = trt))  +
    #   geom_line(lwd = 0.2) +
    #   theme_minimal()


  }


  return(
    list(
      d_pars = d_pars,
      d = d
    )
  )

}

# n_sim = 1000
# id_sim = i
# N = c(2000, 2500, 3000)
# p0 = 0.1
# p1_lwr = 0.07
# p1_upr = 0.07
# p_thresh_sup = 0.9
# p_thresh_fut = 0.4
# p1 <- seq(p1_lwr, p1_upr, length = n_sim)
# p = c(p0, p1[i])

n_sim <- 1000
N <- c(2000, 2500, 3000)
p0 <- 0.1
p1 <- seq(0.07, 0.1, by = 0.01)
p_thresh_sup <- 0.95
p_thresh_fut <- 0.30
id_sim <- 1

run_sim_02 <- function(
    n_sim = 1000,
    N = c(2000, 2500, 3000),
    p0 = 0.1,
    p1 = seq(0.07, 0.1, by = 0.01),
    p_thresh_sup = 0.95,
    p_thresh_fut = 0.30
){

  log_info(paste0(match.call()[[1]]))

  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    mc_cores <- 5
  }

  l_res <- list()

  # repeat each of the elements in p1 n_sim times
  # we will smooth over the true effect size, but we need some way to
  # stratify by interim analysis
  p1 <- rep(p1, each = n_sim)
  i <- 1
  p = c(p0, p1[i])

  log_info("Total simulations to run ", length(p1))

  l_res <- pmclapply(1:length(p1), function(i){

    do_trial(
      id_sim = i, N = N,
      p = c(p0, p1[i]),
      p_thresh_sup = 0.96, p_thresh_fut = 0.35,
      estimation_method = 1
    )

  }, mc.cores = mc_cores, title = paste0("Sequential trial simulation"))

  foutname <- paste0(
    format(Sys.time(), format = "%Y%m%d%H%M%S"), "-sim-res.qs")
  qs::qsave(l_res, file = file.path(output_dir_res, foutname))



  toc(log = TRUE, quiet = TRUE)

  log_info("run_none_sim_05: .", unlist(tic.log(format = TRUE)))

  foutname
}


run_none_sim_02 <- function(){
  log_info("run_none_sim_05: Nothing doing here bud.")
}

main_sim_02 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim_02()
