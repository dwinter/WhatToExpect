#' @importFrom Rcpp sourceCpp
#' @useDynLib WhatToExpect
#' @export

fq_AT <- function(path){
    fastq_AT(path.expand(path))
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

