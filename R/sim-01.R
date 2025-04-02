library(cmdstanr)
library(data.table)
library(parallel)
library(ggplot2)

m1 <- cmdstanr::cmdstan_model("stan/logistic.stan")
n_sim = 100
p0 = 0.1
p1 = 0.1
delta_ni = 0.015

run_01 <- function(n_sim = 100, N = 2000, p0 = 0.10, p1 = 0.10, delta_ni = 0.03){

  y0 <- rbinom(n_sim, N, p0)
  y1 <- rbinom(n_sim, N, p1)

  o0 <- exp(qlogis(p0))
  o1 <- exp(qlogis(p1))
  or <- o1/o0

  ii <- 1

  r <- parallel::mclapply(
    X=1:n_sim, mc.cores = 5, FUN=function(ii) {

      w <- tryCatch({

        ld <- list(n = rep(N, 2), y = c(y0[ii], y1[ii]))

        # f1 <- m1$sample(
        #   ld, iter_warmup = 1000, iter_sampling = 2000,
        #   parallel_chains = 2, chains = 2, refresh = 0, show_exceptions = F,
        #   max_treedepth = 13)

        snk <- capture.output(
          f1 <- m1$pathfinder(ld, num_paths=20, single_path_draws=200,
                              history_size=50, max_lbfgs_iters=100,
                              refresh = 0, draws = 2000)
        )

        post_1 <- data.table(f1$draws(variables = c("p1", "p0"), format = "matrix"))

        mean(post_1$p1 - post_1$p0 < delta_ni) > 0.975

      },
      error=function(e) {
        stop(paste0("Stopping with error ", e))
      })

      w
    })

  r <- unlist(r)
  pwr <- mean(r)

  message("Num. sims ", n_sim)
  message("p0    ", p0)
  message("p1    ", p1)
  message("delta_ni ", delta_ni)
  message("or    ", or)
  message("Power ", pwr)
  pwr
}

N_opts <- seq(1000, 2000, by = 200)
pwr <- numeric(length(N_opts))

for(i in seq_along(N_opts)){
  pwr[i] <- run_01(n_sim = 1000, N = N_opts[i], p0 = 0.10, p1 = 0.10, delta_ni = 0.03)
}

d_fig <- data.table(N_opts, pwr)

ggplot(d_fig, aes(x = N_opts, y = pwr)) +
  geom_point()



