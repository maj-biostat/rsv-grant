# RSV proposal simulations

Current reference sim is `sim02` 




## Run

E.g.

```
Rscript --vanilla ./R/sim-02.R run_sim_02 
```


## Grant text

The primary analysis will be based on a sequential design that evaluates the time to occurrence of RSV by treatment group. 
Analyses will start after 2100 enrolments and be run after every 300 subsequent participants. 
The maximum sample size targets 3000 enrolments, which will be analysed after the 360 day follow up is complete. 
Inference will be based on a Bayesian accelerated failure time model under weakly regularising priors to estimate and compare restricted mean survival time to day 360. 
Decision rules are specified to stop enrolment for superiority and futility. 



Initial simulation work has calibrated the design to use a 99% probability threshold applied to the posterior difference in RMST at the first analysis, falling to 97.25% at the final analysis. Futility is calibrated to use a decision threshold of 35% throughout. Under this setup with a baseline cumulative incidence of 10% at 360 days and 7% in the treatment group, the cumulative probability of declaring superiority is greater than 80% by the maximum sample size. Under the null scenario, the type-i assertion probability is controlled at less than 5% and the probability of declaring futility at the first analysis is greater than 35%. Missingness and the existence of competing risks are expected to be negligible and will therefore be handled through censoring. Further simulation will explore variations to the design and sensitivity options with a focus on efficiency and robustness to misspecification. 
