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



// [[Rcpp::export]]
Rcpp::NumericVector fastq_AT_generic(std::string fq_path, bool gz, size_t nreads){
    Rcpp::NumericVector res;
    std::ifstream* fq;
    if (gz) {
        fq = igzstream(fq_path.c_str());
    }
    else {
        fq = std::ifstream(fq_path);
    }
    size_t i = 0;
    std::string line;
    while(getline(fq,line)){
        if (i % 4 == 1){
            res.push_back( AT(line) );
        }
        i++;
        if( i / 4 >= nreads) {
            break;
        }
    }
    return(res);
}
