#' @importFrom Rcpp sourceCpp
#' @useDynLib WhatToExpect
#' @export

fq_AT <- function(path, gz=TRUE, nreads=-1){
    fastq_AT(path.expand(path), gz, nreads)
}

#' @importFrom mixtools normalmixEM
fit_AT_mixture <- function(path, gz=TRUE, nreads=-1){
    x <- fastq_AT(path,gz, nreads)
    if(any(is.na(x))){
        x <- x[!is.na(x)]
    }
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

estimate_isochores <- function(R1, R2, insert_size, gz=TRUE, nreads=1e4){
    M <-  cbind( fastq_AT(R1, gz, nreads), 
                fastq_AT(R2, gz, nreads))
    if( any(is.na(M)) ){
        warning("Excluding NAs in AT-richnes  data")
        to_drop <- unique(which(is.na(M), arr.ind=TRUE)[,1])
        M <- M[-to_drop,]
    }
    mixture <- mixtools::normalmixEM(as.vector(M))
    p_comp_A <- apply(M, 1, posterior, mixture)
    classified <- apply(p_comp_A, 2, cut, c(0, 0.01, 0.99, 1), labels=c("low", "inter", "high"))
    resolved <- classified[,apply(classified, 2, function(x) !any(x=="inter"))]
    pi_A <- mean(unlist(resolved) == "high")
    pAA <- mean(apply(resolved,2, function(x) all(x=="high")))
    s <- 1 - (pAA/pi_A)^(1/insert_size)
    list(gc=mixture$mu[1], prop_genome=mixture$lambda[1], mean_len = 1/s, transition_prob=s)
}
    




