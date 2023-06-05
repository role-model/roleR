#include <Rcpp.h>

using namespace Rcpp;

// i is the species index of the parent
// n is the number of species that have ever lived (i.e. index after which to 
//   add new species)
// e is the edge matrix
// l is the vector of edge lengths
// alive is vector of whether or not each tip is extant or extinct 

// [[Rcpp::export]]
// i species index of parent
// n is the index at which you add new species (num species that have ever lived)
void specPhyloRCpp(int i, int n, NumericMatrix e, 
                            NumericVector l, LogicalVector alive) {

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
  
  // update n
  n ++;
  
  // update alive vector
  alive(n) = TRUE;
}

