#pragma once

#include <Rcpp.h>

using namespace Rcpp;

class rolePhyloCpp {
    public:
        int n;
        NumericMatrix e;
        NumericVector l;
        LogicalVector alive;
        StringVector tipNames;
        long scale;

        //constructor
        rolePhyloCpp(int n_, NumericMatrix e_,NumericVector l_,
                  LogicalVector alive_,StringVector tipNames_,long scale_)
                : n(n_), e(e_), l(l_), alive(alive_), tipNames(tipNames_), scale(scale_)
        {
        }

        void birth()
        {
        }

        void death(int i)
        {
            // set tip to dead
            alive[i] = false;
        }

        void speciation(int i)
        {
            //number of tips
            //n <- x@n

            //index of where unrealized edges in edge matrix start
            //eNew <- min(which(x@e[, 1] == -1))

            int eNew = -1;
            for (int k = 0; k < n; k++) {
                if (e(k,1) == -1)
                {
                    eNew = k;
                    break;
                }
            }

            //index of where to add new edge
            //j <- which(x@e[, 2] == i)
            int j = -1;
            for (int k = 0; k < n; k++) {
                if (e(k,2) == i)
                {
                    j = k;
                    break;
                }
            }

            // add one to internal node indices
            //e[e > n] <- e[e > n] + 1

            for (int k = 0; k < n; k++) {
                if (k > n)
                {
                    e(k,1) += 1;
                }
            }

            // add new node
            int newNode = 2 * n + 1; // index of new node
            e(1 + eNew, 1) = newNode;
            //e[(0:1) + eNew, 1] = newNode; // add internal node

            e(eNew, 2) = e(j, 2); // add old tip
            e(eNew + 1, 2) = n + 1; // add new tip


            e(j, 2) = newNode;

            // augment edge lengths
            l[0 + eNew] = 0; // add new edges (adding 2 edges?)
            l[1 + eNew] = 0;
            //x@l[x@e[, 2] <= n + 1] <- x@l[x@e[, 2] <= n + 1] + 1 // increase tip length

            // over-write other slots of `x` with new info
            alive[n + 1] = TRUE;
            tipNames[n + 1] = "new";
            n += 1;
        }

        void immigration()
        {

        }
};

RCPP_EXPOSED_CLASS(rolePhyloCpp)
