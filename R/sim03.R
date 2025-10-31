


source("./R/init.R")
source("./R/data.R")

args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library("logger"))

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim03"
  args[2] = "./sim03/cfg-sim03-v01.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}

# Log setup
f_log_sim <- file.path("./logs", "log-sim.txt")
log_appender(appender_file(f_log_sim))
log_info("*** START UP ***")

f_cfgsc <- file.path("./etc", args[2])
g_cfgsc <- config::get(file = f_cfgsc)
stopifnot("Config is null" = !is.null(g_cfgsc))

ix <- 1

output_dir_mcmc <- paste0(getwd(), "/tmp")
output_dir_res <- paste0(getwd(), "/out")





# Main trial loop.
run_trial03_v1 <- function(
    ix,
    l_spec,
    return_posterior = F
){

  log_info("Entered  run_trial for trial ", ix)

  # Get enrolment times for arbitrary large set of patients
  # Simpler to produce enrol times all in one hit rather than try to do them
  # incrementally with a non-hom pp.

  # events per day
  lambda = l_spec$pt_per_day
  # ramp up over x months
  rho = function(t) pmin(t/l_spec$ramp_up_days, 1)


  # rr <- unlist(pblapply(1:100, cl = 4, FUN=function(ii){
  #   ttt <- get_enrol_time(sum(l_spec$N), lambda, rho)
  #   max(ttt)
  # }))
  # mean(rr) / 365

  # day of enrolment
  loc_t0 <- get_enrol_time(sum(l_spec$N_max), lambda, rho)

  # d_fig <- data.table(id = 1:length(loc_t0), loc_t0)
  # d_fig[, yr := loc_t0 / 365]
  # ggplot(d_fig, aes(x = yr, y = id)) +
  #   scale_x_continuous("Year", lim = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
  #   geom_step() +
  #   geom_vline(xintercept = 2.5, col = "red")

  # loop controls
  N_current <- 0
  evt_obs <- 0
  d_all <- data.table()

  ## LOOP -------
  while(evt_obs < l_spec$evt_target & N_current < l_spec$N_max){

    trt <- rep(c(0, 1), len = l_spec$N_batch)
    probs <- fifelse(trt == 1, l_spec$pr_ii, l_spec$pr_mv)
    y <- rbinom(l_spec$N_batch, 1, probs)

    evts_new <- sum(y)

    if(evt_obs + evts_new < l_spec$evt_target){
      d_all <- rbind(
        d_all,
        data.table(trt = trt, y = y)
      )
      evt_obs <- evt_obs + evts_new
      N_current <- N_current + l_spec$N_batch

    } else {
      evt_reqd <- l_spec$evt_target - evt_obs

      evt_cum <- cumsum(y)
      ix_stop <- which(evt_cum >= evt_reqd)[1L]
      d_all <- rbind(
        d_all,
        data.table(trt = trt[1:ix_stop], y = y[1:ix_stop])
      )

      evt_obs <- evt_obs + sum(y[1:ix_stop])
      N_current <- N_current + ix_stop
    }

  }

  d_all[, t_0 := loc_t0[1:.N]]

  if (N_current >= l_spec$N_max) {
    msg <- paste0("Reached ", l_spec$N_max, ". Reconfigure event probabilities or max sample size.")
    stop(msg)
  }


  y0 <- rbeta(l_spec$N_ptcl, 1 + d_all[trt == 0, sum(y)], 1 + d_all[trt == 0, .N] - d_all[trt == 0, sum(y)])
  y1 <- rbeta(l_spec$N_ptcl, 1 + d_all[trt == 1, sum(y)], 1 + d_all[trt == 1, .N] - d_all[trt == 1, sum(y)])

  win <- as.numeric(mean(y1 - y0 < l_spec$dec_delta_sup) > l_spec$dec_thresh_sup)


  return(
    list(
      d_all = d_all,
      win = win

    )
  )

}



run_sim03 <- function(){

  log_info(paste0(match.call()[[1]]))

  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 3
  }

  l_spec <- list()

  l_spec$nex <- g_cfgsc$nex
  l_spec$pt_per_day <- g_cfgsc$pt_per_day
  l_spec$ramp_up_days <- g_cfgsc$ramp_up_days

  l_spec$N_ptcl <- g_cfgsc$N_ptcl

  l_spec$N_max <- g_cfgsc$N_max
  l_spec$N_batch <- g_cfgsc$N_batch

  # Think in terms of a nominal sample size.
  # If you were thinking of a sample size that is around 4000 and rand was 1:1
  # then if pr_mv = 0.1 and pr_ii = 0.07, the expected event count would be
  # 4000 * 0.5 * 0.1 + 4000 * 0.5 * 0.07 = 340
  l_spec$pr_mv <- g_cfgsc$pr_mv
  l_spec$pr_ii <- g_cfgsc$pr_ii
  # number of events targeted for mv vs ii
  l_spec$evt_target <- g_cfgsc$evt_target


  l_spec$prior$pri_a <- g_cfgsc$pri_a
  l_spec$prior$pri_b <- g_cfgsc$pri_b

  l_spec$dec_delta_sup <- g_cfgsc$dec_delta_sup
  l_spec$dec_thresh_sup <- g_cfgsc$dec_thresh_sup

  if(l_spec$nex > 0){
    log_info("Creating ", l_spec$nex, " example trials with full posterior")
    l_spec$ex_trial_ix <- sort(sample(1:g_cfgsc$nsim, size = l_spec$nex, replace = F))
  }
  return_posterior <- F
  str(l_spec)
  e = NULL
  ix <- 1


  r <- pbapply::pblapply(
    X=1:g_cfgsc$nsim, cl = g_cfgsc$mc_cores, FUN=function(ix) {
      log_info("Simulation ", ix);

      if(ix %in% l_spec$ex_trial_ix){
        return_posterior = T
      } else {
        return_posterior = F
      }

      ll <- tryCatch({
        run_trial03_v1(
          ix,
          l_spec,
          return_posterior = return_posterior
        )
      },
      error=function(e) {
        log_info("ERROR in pblapply LOOP (see terminal output):")
        message(" ERROR in pblapply LOOP " , e);
        log_info("Traceback (see terminal output):")
        message(traceback())
        stop(paste0("Stopping with error ", e))
      })

      ll
    })



  d_pr_dec <- data.table()
  for(i in 1:length(r)){

    log_info("Appending d_pr_dec for result ", i)

    if(is.recursive(r[[i]])){
      d_pr_dec <- rbind(
        d_pr_dec,
        cbind(
          sim = i, win = r[[i]]$win
        )
      )
    } else {
      log_info("Value for r at this index is not recursive ", i)
      message("r[[i]] ", r[[i]])
      message(traceback())
      stop(paste0("Stopping due to non-recursive element "))
    }

  }


  # data from each simulated trial
  d_all <- rbindlist(lapply(1:length(r), function(i){
    r[[i]]$d_all
  } ), idcol = "sim")



  d_pr_dec[, mean(win)]

  d_all[, .N, keyby = sim][, mean(N)]


}


run_none_sim03 <- function(){
  log_info("run_none_sim03: Nothing doing here bud.")
}

main_sim03 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim03()
