#include <string>
#include <Rcpp.h>
#include <fstream>

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
Rcpp::NumericVector fastq_AT(std::string fq_path){
    Rcpp::NumericVector res;
    std::ifstream fq(fq_path);
    std::string line;
    size_t i = 0;
    while(getline(fq,line)){
        if (i % 4 == 1){
            res.push_back( AT(line) );
        }
        i++;
    }
    return(res);
}
