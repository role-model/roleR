#include <RcppArmadillo.h>
#include <random>
using namespace Rcpp;
using namespace arma;

// example class where random number generators are members
class randExample {
public:
    std::mt19937 rng; // generator for rand unif nums
    std::uniform_real_distribution<double> dist; // unif dist object
    mat nums; // holds random numbers
    
    // constructor
    randExample(mat nums_) :
        rng(std::mt19937(std::random_device{}())), // this is where a seed could be set
        dist(0, 1),
        nums(nums_) {}
};


// function that fills in random numbers by modifying `randExample x.nums` 
// in place
void updateRand(randExample &x) {
    int n = x.nums.n_rows;
    
    // initialize double to hold random number
    double r;
    for (int i = 0; i < n; i++) {
        r = x.dist(x.rng);
        x.nums(i, 0) = r;
    }
}


// export a wrapper function for testing
//' @name testRand
//' @title testRand
//' @param m m
//' @return outMat
// [[Rcpp::export]]

NumericMatrix testRand(NumericMatrix m) {
    mat numMat = as<mat>(m);
    randExample z = randExample(numMat);
    
    updateRand(z);
    
    NumericMatrix outMat = as<NumericMatrix>(wrap(z.nums));
    
    return outMat;
}
