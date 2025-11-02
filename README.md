# RSV proposal simulations

Simulation for grant submission CLINICAL TRIALS AND COHORT STUDIES GRANT 2025 Application ID: 2054492 (Snelling).

Based on evaluation of risk by treatment group using sequential design.

## Run

For example, using the `sim04` implementation (the one for the grant):

```
Rscript --vanilla ./R/sim04.R run_sim04 ../etc/sim04/cfg-sim04-v01.yml
```

```
 $ desc         : chr "Null"
 $ nsim         : int 1000
 $ nex          : int 10
 $ pt_per_day   : num 4.6
 $ ramp_up_days : int 90
 $ N_ptcl       : int 10000
 $ N            : int [1:5] 2000 500 500 500 500
 $ pr_mv        : num 0.1
 $ pr_ii        : num 0.1
 $ prior        :List of 2
  ..$ pri_a: int 2
  ..$ pri_b: int 10
 $ delta        :List of 2
  ..$ sup: int 0
  ..$ fut: num -0.01
 $ thresh       :List of 2
  ..$ sup: num 0.99
  ..$ fut: num 0.7
 $ print_summary: int 1
```



