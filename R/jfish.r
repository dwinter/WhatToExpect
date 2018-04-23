jfish_count_kmers <- function(read_files, k, starting_size, out_file="mer_counts.jf", gz=FALSE, exe="jellyfish"){
    reads <- paste(read_files, collapse=" ")
    args <- paste ("count -C -m", k, "-s", starting_size)
    if(gz){
        zcat_args <- paste(reads, "| ", exe, args, "/dev/fd/0")
        ret_code <- system2("zcat", zcat_args)
    } else{
        ret_code <- system2(exe, paste(args, reads))
    }
    if(ret_code != 0){
        stop("Jellyfish return error code ", ret_code)
    }
    if( file.exists(out_file) ){
        message("kmer counts written to ",out_file)
        return(out_file)
    }
    stop("jellyfish ran, but outfile doesn't seem to exist?")
}

read_jfish <- function(hist_file, exe="jellyfish"){
    tf <- tempfile()
    ret_code <- system2(exe, paste("histo", hist_file), stdout=tf)
    if (ret_code != 0){
        stop("Jellyfish return error code", ret_code)
    }
    read.table(tf, col.names=c("n", "freq"))
}


plot_kmer_hist <- function(hist, x_cutoff=0.95){
    x <- which(cumsum(kmer_counts$freq)/sum(kmer_counts$freq) > x_cutoff)[1]
    plot(hist[1:x,], type='l', xlab="Number of kmers")
}

find_optima <- function(x, maxima=TRUE){
    diff_of_diffs <- diff(sign(diff(x)))
    if(maxima){
        return(which(diff_of_diffs == -2) + 1)
    }
    which(diff_of_diffs == 2) + 1
}

est_genome_size <- function(kmer_hist, unit="Mb", plot=TRUE, x_cutoff=0.95){
    denom <- switch(tolower(unit), 
             "b" =  1e0,
             "kb" = 1e3,       
             "mb" = 1e6,
             "gb" = 1e9,
              stop("Units should be in 'b', 'Kb', 'Mb', 'Gb', got", unit)
    )
    error_edge <- find_optima(kmer_hist[,2], maxima=FALSE)[1]
    NR <- nrow(kmer_hist)
    peak <- which.max(kmer_hist$freq[(error_edge+1):NR]) + error_edge
    if(plot){
        top <- max(kmer_hist$freq)
        peak_y <- kmer_hist$freq[peak]
        err_colour <- "#C40018"
        peak_colour <- "#292725"
        plot_kmer_hist(kmer_hist)
        rect(xleft=1, xright=error_edge, ytop=top, ybottom=0, col=paste0(err_colour,"99"), border="red")
        text(x=error_edge, y=0.9 * max(kmer_hist$freq), "Excluded kmers", pos=4, col=err_colour)
        points(x=peak, peak_y, col=peak_colour, pch=19)
        lab = paste0("peak coverage (", peak, ")")
        text(x=peak, y=kmer_hist$freq[peak], lab, col=peak_colour, pos=3)
    }
    
    sum(apply(kmer_hist[9:nrow(kmer_hist),], 1, prod)) / (peak*denom)
}





 