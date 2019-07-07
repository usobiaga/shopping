#' CDnow dataset sample
#' 
#' Dataset pubished by Bruce Hardie see \url{http://www.brucehardie.com/datasets/}.
#' 
#' \enumerate{
#'   \item \code{m_id} is the client id in the master dataset.
#'   \item \code{s_id} is the client id.
#'   \item \code{date} is the date of the purchase
#'   \item \code{n} is the number of purchased CDs
#'   \item \code{value} is the dollars payed.
#' }
#' 
#' @name cdnowSample
#' @docType data
#' 
NULL

#' Estimates the bgnbd model base on “Counting Your Customers” the Easy Way: An Alternative
#' to the Pareto/NBD Model (2005) by Peter S FaderPeter S FaderBruce G. S. HardieKa Lok Lee
#'
#' @param Tx lenght of time since the consumers first purchase. 
#' @param tx recency of the last purchase (time of the last purchase - time of the first).
#' @param n count fo repeat purchases.
#' @param start initial values for the bgnbd model.
#' @param lower lower argument for optimization defaults to \code{c(0.01, 0.01, 0.01, 0.01)}.
#'
#' @export
#'
#'@examples
#'
#' data(cdnowSample)
#'
#'## helper function to get the model variables
#' 
#'getMf <- function(set){
#'    set$week <- as.numeric(set$date - min(set$date)) / 7
#'    maxWeek <- max(set$week)
#'    set %>%
#'        group_by(s_id) %>% ## for each customer ID
#'        summarise(
#'            ## total time being a client
#'            Tx = maxWeek - min(week),
#'            ## time between last and first transactions
#'            tx = max(week) - min(week),
#'            ## weeks with purchases (ignoring first purchase week)
#'            n = length(unique(week)) - 1,
#'            )
#'}
#'
#' splitDate <- as.Date("1997-09-30")
#' calibrateSet <- getMf(filter(cdnowSample, date <= splitDate))
#' set.seed(345)
#' model <- bgnbd(Tx = calibrateSet$Tx, tx = calibrateSet$tx,
#' n = calibrateSet$n, start = c(1, 3, 1, 2.4))
#' pred <- predict(model, Tx = calibrateSet$Tx, tx = calibrateSet$tx, n = calibrateSet$n, 39)
#' 
#'
bgnbd <- function(Tx, tx, n, start = rep(1.5, 4), lower = rep(0.01, 4), ...){

    stopifnot(
        all(lower > 0),
        length(Tx) == length(tx),
        length(Tx) == length(n),
        all(Tx >= tx),
        all(Tx >= 0),
        all(tx >= 0),
        !anyNA(Tx),
        !anyNA(tx),
        !anyNA(n),
        all(n >= 0),
        sum(n) > 1)

    opt_bgnbd <- function(pars, Tx, tx, n){
        r <- tryCatch({
            -sum(bgnbd_ll(pars, Tx, tx, n))
        }, warning=function(cond) {
            1
        })
        if (is.na(r)) r <- 1
        r
    }
    
    opt <- stats::nlminb(
        objective = opt_bgnbd,
        Tx = Tx,
        tx = tx,
        n = n,
        start = start,
        lower = lower,
        ...)
    class(opt) <- 'bgnbd'
    opt
}


#' Predict expected sales
#'
#' Expected conditional sales in a \code{t} time period. \eqn{E(Y(t)| X=x, t_x, T, r, \alpha, a, b)}
#'
#' @param object a object of class \code{bgnbd}.
#' @param Tx lenght of time since the consumers first purchase. 
#' @param tx recency of the last purchase (time of the last purchase - time of the first).
#' @param n count fo repeat purchases.
#' @param t time period in which the expeted sales take place. 
#'
#' @import BMS
#' @export
#'
predict.bgnbd <- function(object, Tx, tx, n, t){

    stopifnot(
        length(Tx) == length(tx),
        length(Tx) == length(n),
        all(Tx >= tx),
        all(Tx >= 0),
        all(tx >= 0),
        all(t >= 0),
        !anyNA(Tx),
        !anyNA(tx),
        !anyNA(n),
        all(n >= 0),
        sum(n) > 1)
    
    pars <- object$par
    A1 <- (pars[3]+pars[4]+n-1)/(pars[3]-1)
    A2 <- ((pars[2]+Tx)/(pars[2]+Tx+t))^(pars[1]+n)
    A3 <- Vectorize(BMS::f21hyper)(pars[1]+n, pars[4]+n, pars[3]+pars[4]+n-1, t/(pars[2]+Tx+t))
    A4 <- 1 + (n>0)*(pars[3]/(pars[4]+n-1)*((pars[2]+Tx)/(pars[2]+tx))^(pars[1]+n))
    expt <- A1*(1-A2*A3)/A4
    expt
}

## palives

#' Probability of being alive.
#'
#' generic function to compute the probability of a customer making a purchase in the future.
#'
#' @param object estimation
#' @param ... parameters passed to other methods
#' @export
palive <- function(object, ...) UseMethod('palive')



#' Probability of being alive.
#'
#' Palive for bgnbd
#'
#' @param object object of class bgnbd
#' @param Tx lenght of time since the consumers first purchase. 
#' @param tx recency of the last purchase (time of the last purchase - time of the first).
#' @param n count fo repeat purchases.
#' @param ... unused
#' 
#' @export
palive.bgnbd <- function(object, Tx, tx, n){

    stopifnot(
        length(Tx) == length(tx),
        length(Tx) == length(n),
        all(Tx >= tx),
        all(Tx >= 0),
        all(tx >= 0),
        !anyNA(Tx),
        !anyNA(tx),
        !anyNA(n),
        all(n >= 0),
        sum(n) > 1)
    
    pars <- object$par
    term1 <- (pars[3]/(pars[4] + n - 1)) * ((pars[2] + Tx)/(pars[2] + tx))^(pars[1] + n)
    1/(1 + as.numeric(n > 0) * term1)
}


#' loglike of BGNBD
#' 
#' @param pars (r, alpha, a, b) parameters of the loglike.
#' @param Tx lenght of time since the consumers first purchase.
#' @param tx recency of the last purchase (time of the last purchase - time of the first)
#' @param n count fo repeat purchases
#' 
bgnbd_ll <- function(pars, Tx, tx, n){
    A1 <- lgamma(pars[1]+n)+pars[1]*log(pars[2]) - lgamma(pars[1])
    A2 <- lgamma(pars[3]+pars[4]) + lgamma(pars[4]+n)-
        lgamma(pars[4])-lgamma(pars[3]+pars[4]+n)
    A3 <- (pars[1]+n)*(-log(pars[2]+Tx))
    A4 <- log(pars[3]) - log(pars[4]+n-1) + (pars[1]+n)*(-log(pars[2]+tx))
    A1 + A2 + log(exp(A3)+(n>0)*exp(A4))
}


