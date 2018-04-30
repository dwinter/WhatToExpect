#include <string>
#include <Rcpp.h>
#include <fstream>
#include "gzstream.h"

double AT(std::string sequence){
    uint32_t n = 0;
    uint32_t AT = 0;
    for( auto b : sequence) {
        switch(b){
            case 'A':   n++;
                        AT++;
                        break;
            case 'a':   n++;
                        AT++;
                        break;
            case 'T':   n++;
                        AT++;
                        break;
            case 't':   n++;
                        AT++;
                        break;
            case 'C':   n++;
                        break;
            case 'c':   n++;
                        break;
            case 'G':   n++;
                        break;
            case 'g':   n++;
                        break;
            default:    break;
                      
        }
    }    
    return(AT / double(n));
}

template < typename stream_type >
Rcpp::NumericVector fastq_AT_generic(stream_type &fq_stream, size_t nreads){
    Rcpp::NumericVector res;
    size_t i = 0;
    std::string line;
    while(getline(fq_stream,line)){
        if (i % 4 == 1){
            res.push_back( AT(line) );
        }
        i++;
        if( i / 4 >= nreads) {
            break;
        }
    }
    if (i == 0){
        Rcpp::stop("Couldn't find any reads in fastq file");
    }
    return(res);
}

// [[Rcpp::export]]
Rcpp::NumericVector fastq_AT(std::string fq_path, bool gz, size_t nreads){
    Rcpp::NumericVector res;
    if (gz) {
        igzstream fq(fq_path.c_str());
        res = fastq_AT_generic( fq, nreads ) ;
    }
    else {
        std::ifstream fq(fq_path);
        res = fastq_AT_generic( fq, nreads );
    }
    return(res);
}
