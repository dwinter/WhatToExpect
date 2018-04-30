#' @importFrom Rcpp sourceCpp
#' @useDynLib WhatToExpect
#' @export

fq_AT <- function(path, gz=TRUE, nreads=-1){
    fastq_AT(path.expand(path), gz, nreads)
}

#' @importFrom mixtools normalmixEM
fit_AT_mixture <- function(path){
    x <- fastq_AT(path)
    normalmixEM(x)
}

plot_AT_mix <- function(mix){
    h <- hist(mix$x, breaks="FD", probability=TRUE)
    lines(h$breaks, dnorm(h$breaks, mix$mu[1], mix$sigma[1]) * mix$lambda[1], lwd=3, col="red")
    lines(h$breaks, dnorm(h$breaks, mix$mu[2], mix$sigma[2]) * mix$lambda[2], lwd=3, col="blue")
}

summarize_AT_mix <- function(mix){
    data.frame( component = 1:2,est_prop = mix$lambda, mean_AT = mix$mu, sd_AT = mix$sigma)
}

posterior <- function(x, fitted_model){
    A <- dnorm(x, fitted_model$mu[1], fitted_model$sigma[1])
    B <- dnorm(x, fitted_model$mu[2], fitted_model$sigma[2])
    A/(A+B)
}
