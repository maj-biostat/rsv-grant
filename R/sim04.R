


source("./R/init.R")
source("./R/data.R")

args = commandArgs(trailingOnly=TRUE)


# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim04"
  args[2] = "./sim04/cfg-sim04-v01.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}

# Log setup
f_log_sim <- file.path("./logs", "log-sim04.txt")
log_appender(appender_file(f_log_sim))
log_info("*** START UP ***")

f_cfgsc <- file.path("./etc", args[2])
g_cfgsc <- config::get(file = f_cfgsc)
stopifnot("Config is null" = !is.null(g_cfgsc))

ix <- 1

output_dir_mcmc <- paste0(getwd(), "/tmp")
output_dir_res <- paste0(getwd(), "/out")





# Main trial loop.
run_trial04_v1 <- function(
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
  loc_t0 <- get_enrol_time(sum(l_spec$N), lambda, rho)

  # d_fig <- data.table(id = 1:length(loc_t0), loc_t0)
  # d_fig[, yr := loc_t0 / 365]
  # ggplot(d_fig, aes(x = yr, y = id)) +
  #   scale_x_continuous("Year", lim = c(0, 4), breaks = seq(0, 4, by = 0.5)) +
  #   geom_step() +
  #   geom_vline(xintercept = 2.5, col = "red")

  # loop controls
  stop_enrol <- FALSE
  N_analys <- length(l_spec$N)
  # set interim number
  l_spec$ic <- 1
  d_all <- data.table()

  # decisions
  g_rule_type <- c("sup",  "fut")

  d_pr_dec <- CJ(
    ic = 1:N_analys,
    rule = factor(g_rule_type, levels = g_rule_type),
    par = factor(c("rd")),
    p = NA_real_,
    dec = NA_integer_
  )

  d_anlys <- CJ(
    ic = 1:N_analys,
    N_enrol = NA_real_,
    N_anlys = NA_real_,
    t_anlys = NA_real_,
    stopped = NA_integer_
  )

  ## LOOP -------
  while(!stop_enrol){

    log_info("Trial ", ix, " cohort ", l_spec$ic)

    # next chunk of data on pts.
    if(l_spec$ic == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
    } else {
      l_spec$is <- l_spec$ie + 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
    }

    trt <- rep(c(0, 1), len = l_spec$N[l_spec$ic])
    probs <- fifelse(trt == 1, l_spec$pr_ii, l_spec$pr_mv)
    y <- rbinom(l_spec$N[l_spec$ic], 1, probs)
    t_0 <- loc_t0[l_spec$is:l_spec$ie]
    t_fu <- t_0 + 365


    # combine the existing and new data
    d_all <- rbind(
      d_all,
      data.table(
        id = l_spec$is:l_spec$ie,
        t_0 = t_0, t_fu = t_fu,
        trt = trt,
        y = y
      )
      )


    if(l_spec$ie == sum(l_spec$N)){
      log_info("Trial ", ix, " final analysis, using all pt")
      t_now <- d_all[, max(t_fu)]
      d_mod <- copy(d_all)

    } else {
      # t0 is the entry time and is repeated
      t_now <- d_all[, max(t_0)]

      # include all those who have reached the first follow up at the time of
      # the interim analysis
      incl_ids <- d_all[t_fu <= t_now, id]
      d_mod <- d_all[id %in% incl_ids]

    }

    d_post <- data.table(
      p_0 = rbeta(
        l_spec$N_ptcl,
        1 + d_mod[trt == 0, sum(y)],
        1 + d_mod[trt == 0, .N] - d_mod[trt == 0, sum(y)]),
      p_1 = rbeta(
        l_spec$N_ptcl,
        1 + d_mod[trt == 1, sum(y)],
        1 + d_mod[trt == 1, .N] - d_mod[trt == 1, sum(y)])
    )
    d_post[, rd := p_1 - p_0]

    d_post_long <- melt(d_post, measure.vars = names(d_post),
                        variable.name = "par")



    d_pr_dec[

      rbind(
        d_post_long[par %in% c("rd"), .(
          ic = l_spec$ic,
          rule = factor("sup", levels = c("sup", "fut")),
          p = mean(value < l_spec$delta$sup),
          dec = as.integer(mean(value < l_spec$delta$sup) > l_spec$thresh$sup)
        ), keyby = par],
        d_post_long[par %in% c("rd"), .(
          ic = l_spec$ic,
          rule = factor("fut", levels = c("sup", "fut")),
          p = mean(value > l_spec$delta$fut),
          dec = as.integer(mean(value > l_spec$delta$fut) > l_spec$thresh$fut)
        ), keyby = par]
      ),

      on = .(ic, rule, par), `:=`(
        p = i.p, dec = i.dec
      )
    ]

    d_anlys[ic == l_spec$ic, N_enrol := nrow(d_all)]
    d_anlys[ic == l_spec$ic, N_anlys := nrow(d_mod)]
    d_anlys[ic == l_spec$ic, t_anlys := t_now]
    d_anlys[ic == l_spec$ic, stopped := d_pr_dec[ic == l_spec$ic, as.integer(any(dec))]  ]
    d_anlys[, stopped := cumsum(stopped)]
    d_anlys[stopped > 1, stopped := 1]

    d_stop <- d_pr_dec[
      ic <= l_spec$ic & par %in% c("rd"),
      .(resolved = as.integer(sum(dec) > 0)), keyby = .(par)]

    # if(d_stop$resolved[1] == 1){
    #   log_info("Stop trial as all questions addressed ", ix)
    #   stop_enrol <- T
    # }

    # next interim
    l_spec$ic <- l_spec$ic + 1

    if(l_spec$ic > N_analys){
      stop_enrol <- T
    }


  }


  # did we stop (for any reason) prior to the final interim?
  if(any(d_pr_dec$dec == 1)){
    stop_at <- d_pr_dec[dec == 1, ic[1]]
  } else {
    stop_at <- l_spec$ic - 1
  }


  l_ret <- list(
    # data collected in the trial
    d_all = d_all,

    d_pr_dec = d_pr_dec,
    d_anlys = d_anlys,

    stop_at = stop_at
  )

  # if(return_posterior){
  #   l_ret$d_post_all <- copy(d_post_all)
  # }
  #
  log_info("Return trial ", ix)
  return(l_ret)

}



run_sim04 <- function(){

  log_info(paste0(match.call()[[1]]))

  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 3
  }

  l_spec <- list()
  l_spec$desc <- g_cfgsc$desc
  l_spec$nsim <- g_cfgsc$nsim
  l_spec$nex <- g_cfgsc$nex

  l_spec$nex <- g_cfgsc$nex
  l_spec$pt_per_day <- g_cfgsc$pt_per_day
  l_spec$ramp_up_days <- g_cfgsc$ramp_up_days

  l_spec$N_ptcl <- g_cfgsc$N_ptcl

  l_spec$N <- g_cfgsc$N_pt

  l_spec$pr_mv <- g_cfgsc$pr_mv
  l_spec$pr_ii <- g_cfgsc$pr_ii

  l_spec$prior$pri_a <- g_cfgsc$pri_a
  l_spec$prior$pri_b <- g_cfgsc$pri_b

  l_spec$delta$sup <- g_cfgsc$dec_delta_sup
  l_spec$thresh$sup <- g_cfgsc$dec_thresh_sup

  l_spec$delta$fut <- g_cfgsc$dec_delta_fut
  l_spec$thresh$fut <- g_cfgsc$dec_thresh_fut

  l_spec$print_summary <- g_cfgsc$print_summary

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
        run_trial04_v1(
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





  # data from each simulated trial
  d_pr_dec <- rbindlist(lapply(1:length(r), function(i){
    r[[i]]$d_pr_dec
  } ), idcol = "sim")

  # data from each simulated trial
  d_all <- rbindlist(lapply(1:length(r), function(i){
    r[[i]]$d_all
  } ), idcol = "sim")

  # data from each simulated trial
  d_anlys <- rbindlist(lapply(1:length(r), function(i){
    r[[i]]$d_anlys
  } ), idcol = "sim")



  if(l_spec$print_summary){


    # cumulative probability of superiority / futility.

    d_dec_1 <- copy(d_pr_dec[par %in% c("rd")])

    # number of enrolments at each interim (interim sample size sequence)
    d_N <- data.table(ic = seq_along(l_spec$N), N = cumsum(l_spec$N))

    # compute the cumulative instances of a decision being made by sim, each
    # decision type and by parameter
    # The decision rule has to be both active at this interim and triggered
    d_dec_1[, cdec := as.integer(cumsum(dec)>0), keyby = .(sim, rule)]
    d_dec_1[, cdec := as.logical(cdec)]
    d_dec_1 <- merge(d_dec_1, d_N, by = "ic")

    # cumulative proportion for which each decision quantity has been met by
    # analysis and domain
    d_dec_cumprob <- d_dec_1[, .(pr_val = mean(cdec)), keyby = .(ic, N, rule)]

    d_cprob_dec <- data.table(
      desc = l_spec$desc,
      rd_trt = l_spec$pr_ii - l_spec$pr_mv,
      d_dec_cumprob)

    d_tbl <- dcast(d_cprob_dec, desc + rd_trt + rule ~ N, value.var = "pr_val")
    d_tbl[rule == "sup", ref := as.numeric(l_spec$delta$sup)]
    d_tbl[rule == "sup", pr_thresh := l_spec$thresh$sup]

    d_tbl[rule == "fut", ref := as.numeric(l_spec$delta$fut)]
    d_tbl[rule == "fut", pr_thresh := l_spec$thresh$fut]

    setcolorder(d_tbl, c("desc", "rd_trt", "rule", "ref", "pr_thresh"))

    pander(d_tbl, split.table = 120)


    # average sample size
    # based on the enrolment rate and primary endpoint assumed to be at 12 months

    # d_anlys[, .(N = )]
    d_tbl <- d_anlys[stopped == 1, .SD[1], keyby = sim][
      , .(desc = l_spec$desc,
          note = "Avg sample size",
          E_N_enrol = mean(N_enrol),
          sd_N_enrol = sd(N_enrol),
          E_N_analys = mean(N_anlys),
          E_t = mean(t_anlys))]

    pander(d_tbl, split.table = 120)


    # average timing of first analysis

    # d_anlys[, .(N = )]
    d_tbl <- d_anlys[ic == 1, .SD[1], keyby = sim][
      , .(desc = l_spec$desc, note = "Days to analys 1", E_t = mean(t_anlys))]

    pander(d_tbl, split.table = 120)


    # expected duration sample size





  }




}


run_none_sim04 <- function(){
  log_info("run_none_sim04: Nothing doing here bud.")
}

main_sim04 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim04()
