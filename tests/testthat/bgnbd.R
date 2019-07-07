
library(dplyr)
library(testthat)

data(cdnowSample)

getMf <- function(set){
    set$week <- as.numeric(set$date - min(set$date)) / 7
    maxWeek <- max(set$week)
    set %>%
        group_by(s_id) %>% ## for each customer ID
        summarise(
            ## total time being a client
            Tx = maxWeek - min(week),
            ## time between last and first transactions
            tx = max(week) - min(week),
            ## weeks with purchases (ignoring first purchase week)
            n = length(unique(week)) - 1,
            )
}

splitDate <- as.Date("1997-09-30")
calibrateSet <- getMf(filter(cdnowSample, date <= splitDate))

testthat::context("bgnbd")

testthat::test_that("bgnbd estimation",{

    ## model estimation
    set.seed(123)
    model <- shopping::bgnbd(Tx = calibrateSet$Tx,
                   tx = calibrateSet$tx,
                   n = calibrateSet$n,
                   start = c(0.243, 4.414, 0.793,  2.426))
    truePar <- c(0.242603727054404, 4.41380777966353, 0.79289875692843, 2.42583420527876)
    expect_equal(model$par, truePar, tolerance=1e-3)

    ##expected sales
    pars <- c(0.24259395230803,
              4.41359091416604,
              0.792919955839573,
              2.42589404751842)
    n <- 2
    tx <- 30.43
    Tx <- 38.8571428571429
    t <- 39
    ss <- shopping::predict.bgnbd(list(par = pars), Tx = Tx, tx=tx, n=n, t=t)
    expect_equal(ss, 1.2260, tolerance=1e-3)

    ## expected sales vectorized
    n <- rep(n, 3)
    tx <- rep(tx, 3)
    Tx <- rep(Tx, 3)
    t <- rep(t, 3)
    ss <- shopping::predict.bgnbd(list(pars = pars), Tx = Tx, tx=tx, n=n, t=t)

    ##palive
    pa <- shopping::palive(model, Tx = Tx, tx=tx, n=n)
    expect_equal(pa, rep(0.7266, 3), tolerance=1e-3)
    
})



