# plotsurv
Function to plot survival curves

``` r
n <- 1000
library(tidyverse)
tte.table <- tibble::tibble(
  group = rbinom(n,1,0.5),
  ev1 = rpois(n,15) - group*2,
  ev2 = rpois(n,16) - group*7,
  cens = rpois(n,10) - group*3,
  time = pmin(ev1, ev2, cens),
  eventtype = case_when(
    time == ev1 ~ 1,
    time == ev2 ~ 2,
    time == cens ~ 0
  )
)

library(survival)
fit <- survfit(Surv(time, factor(eventtype)) ~ group, data = tte.table)

library(plotsurv)
plotsurv(fit,
         risk.table=T,
         x.breaks=7,
         ticks = T,
         ticksize = 6,
         tickalpha = 0.5,
         include_surv=F,
         conf.int=T,
         group_labels = c("Males, Survival", "Males, ev1", "Females, survival", "Females, ev1"),
         color_lab = "Outcome",
         fill_lab = "Outcome", 
         line_lab = "Outcome") 


```


