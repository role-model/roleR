#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;

//' @name rolePhyloCpp
//' @title a C++ class to specify the model phylogeny
//' @param n an \code{integer} containing the number of tips in the phylogeny
//' @param e a two-column \code{numeric matrix} containing all ancestor-child pairs
//' @param l a \code{numeric vector} of edge lengths (in units of time steps = 1/J' generations)
//' @param alive a \code{logical vector} indicating whether tips are extant or not
//' @param tipNames a \code{string vector} of tip names
//' @param scale a \code{numeric} containing time scale translation to years

class rolePhyloCpp {
    public:
        int n;
        NumericMatrix e;
        NumericVector l;
        LogicalVector alive;
        StringVector tipNames;
        double scale;

        //constructor
        rolePhyloCpp(int n_, NumericMatrix e_,NumericVector l_,
                  LogicalVector alive_,StringVector tipNames_,double scale_)
                : n(n_), e(e_), l(l_), alive(alive_), tipNames(tipNames_), scale(scale_)
        {
        }
        
        rolePhyloCpp(){
        }
        
        void birth()
        {
        }

        void death(int i)
        {
            // set tip to dead
            alive[i] = false;
        }
        
        void speciation(int i) {
          
          // nrows of the edge matrix
          int eMax = e.nrow();
          
          // find index of where unrealized edges in edge matrix start
          // eNew <- min(which(e[, 1] == -1))
          int eNew = -1;
          
          for (int k = 0; k < eMax; k++) {
            if (e(k, 0) == -1) {
              eNew = k;
              break;
            }
          }
          
          // return eNew;
          
          // index of the edge matrix of where to add new edge
          // j <- which(e[, 2] == i)
          int j = -1;
          for (int k = 0; k < eMax; k++) {
            if (e(k, 1) == i) {
              j = k;
              break;
            }
          }
          
          // return j;
          
          
          // add one to internal node indices
          //e[e > n] <- e[e > n] + 1
          
          for (int r = 0; r < eNew; r++) {
            for (int c = 0; c < 2; c++){
              if (e(r, c) > n) {
                e(r, c) ++;
              }
            }
          }
          
          // return e(0, 0);
          
          // add new internal node
          int newNode = 2 * n + 1; // index of new node
          e(eNew, 0) = newNode;
          e(1 + eNew, 0) = newNode; // do this more elegantly
          
          // add tips
          e(eNew, 1) = e(j, 1); // add old tip
          e(eNew + 1, 1) = n + 1; // add new tip
          
          // update ancestry of internal nodes
          e(j, 1) = newNode;
          
          // augment edge lengths
          l[eNew] = 0;
          l[1 + eNew] = 0;
          
          // increase all tip edge lengths by 1 time step
          for (int r = 0; r <= eNew + 1; r++) {
            if (e(r, 1) <= n + 1) {
              l(r) ++;
            }
          }
          
          
          // update alive vector
          alive(n) = TRUE;
          
          // update n
          n++;
          
        }

        void speciationOld(int i)
        {
            bool print = true; 

            // find index of where unrealized edges in edge matrix start
            // eNew <- min(which(e[, 1] == -1))
            int eNew = -1;
            for (int k = 0; k < n - 1; k++) {
                if (e(k,1) == -1)
                {
                    eNew = k;
                    break;
                }
            }
            
            if(print){Rcout << "found unrealized edge \n";}
            
            //index of where to add new edge
            // j <- which(e[, 2] == i)
            int j = -1;
            for (int k = 0; k < n; k++) {
                if (e(k,2) == i)
                {
                    j = k;
                    break;
                }
            }
            
            if(print){Rcout << "found index of where to add new edge \n";}

            // add one to internal node indices
            //e[e > n] <- e[e > n] + 1

            for (int r = 0; r < e.nrow() - 1; r++) {
              for(int c = 0; c < e.ncol() - 1; c++){
                if(e(r,c) > n)
                {
                  e(r,c) += 1;
                }
              }
            }
            
            if(print){Rcout << "added one to internal node indices \n";}
            
            // add new node ///takes index of new internal node //may just be 2n 
            int newNode = 2 * n + 1; // index of new node
            e(1 + eNew, 1) = newNode;
          
            //e[(0:1) + eNew, 1] = newNode; // add internal node
            
            if(print){Rcout << "added new node \n";}
            
            e(eNew, 2) = e(j, 2); // add old tip
            
            if(print){Rcout << "added old tip \n";}
            
            e(eNew + 1, 2) = n + 1; // add new tip

            if(print){Rcout << "added new tip \n";}
            
            e(j, 2) = newNode;
            
            if(print){Rcout << "added something \n";}
            
            // augment edge lengths
            l[0 + eNew] = 0; // add new edges (adding 2 edges?)
            l[1 + eNew] = 0;
            //x@l[x@e[, 2] <= n + 1] <- x@l[x@e[, 2] <= n + 1] + 1 // increase tip length
            
            if(print){Rcout << "augmented edge lengths \n";}
            
            // over-write other slots of `x` with new info
            alive.push_back(true);
            tipNames.push_back("new");
     
            n += 1;
            
            if(print){Rcout << "overwrote other slots of x with new info \n";}
        }

        void immigration()
        {

        }
};

RCPP_EXPOSED_CLASS(rolePhyloCpp)
