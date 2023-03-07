// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <string>

using namespace Rcpp;

class roleDataCpp {
public:
    NumericVector indTraitL; 
    NumericVector indSpeciesL; 
    
    NumericVector spAbundL;
    NumericVector spTraitL;
    
    NumericVector spAbundHarmMeanL;
    NumericVector spLastOriginStepL;
    NumericVector spExtinctionStepL; 
    NumericVector spReciprSumL;
    
    NumericVector founderFlagL; 
    
    NumericVector spAbundM;
    NumericVector spTraitM;
    
    NumericVector nTipsP;
    //NumericMatrix edgesP;
    arma::imat edgesP;
    arma::vec lengthsP;
    LogicalVector aliveP;
    CharacterVector tipNamesP;
    NumericVector scaleP;    
    
    arma::mat traitDiffsSq;
    NumericVector traitPow;
    NumericVector envFilterProbs;
    NumericVector compProbs;
    double prevEnvSigma;
    
    std::mt19937 rng; // generator for rand unif nums
    std::uniform_real_distribution<double> unif; // unif distr object that generates between 0 and 1
    std::normal_distribution<double> norm;
    std::normal_distribution<double> tnorm; // norm distr object for trait change
    std::normal_distribution<double> sptnorm; // norm distr object for trait change in speciation 
    
    //roleDataCpp(RObject data){
    //    RObject local = data.slot("localComm");
    //    RObject meta = data.slot("metaComm");
        //indSpTrtL = data.slot("localComm").slot(
    //    indSpTrtL = Rcpp::as<NumericMatrix>(local.slot("indSppTrt"));
        //later parse indSpTrt to spTrt
    //    spAbundTrtM = Rcpp::as<NumericMatrix>(local.slot("sppAbundTrt"));
    //}
    roleDataCpp(RObject local, RObject meta, RObject phylo){

        // parse slots
        //indSpTrtL = Rcpp::as<NumericMatrix>(local.slot("indSppTrt"));
        indTraitL = Rcpp::clone(Rcpp::as<NumericVector>(local.slot("indTrait")));
        indSpeciesL = Rcpp::clone(Rcpp::as<NumericVector>(local.slot("indSpecies")));
        
        //spAbundTrtL = Rcpp::as<NumericMatrix>(local.slot("spAbundTrt"));
        spAbundL = Rcpp::clone(Rcpp::as<NumericVector>(local.slot("spAbund")));
        spTraitL = Rcpp::clone(Rcpp::as<NumericVector>(local.slot("spTrait")));
        
        spAbundHarmMeanL = Rcpp::clone(Rcpp::as<NumericVector>(local.slot("spAbundHarmMean")));
        spLastOriginStepL = Rcpp::clone(Rcpp::as<NumericVector>(local.slot("spLastOriginStep")));
        spExtinctionStepL = Rcpp::clone(Rcpp::as<NumericVector>(local.slot("spExtinctionStep")));
        
        spReciprSumL = NumericVector(10000); // given a 10000 len 
        
        //spAbundTrtM = Rcpp::as<NumericMatrix>(meta.slot("sppAbundTrt"));
        spAbundM = Rcpp::clone(Rcpp::as<NumericVector>(meta.slot("spAbund")));
        spTraitM = Rcpp::clone(Rcpp::as<NumericVector>(meta.slot("spTrait")));
        
        
        //NumericVector nTipsP_vect = Rcpp::as<NumericVector>(phylo.slot("n"));
        //nTipsP = nTipsP_vect
        nTipsP = Rcpp::clone(Rcpp::as<NumericVector>(phylo.slot("n")));
        edgesP = Rcpp::as<arma::imat>(phylo.slot("e"));
        edgesP = edgesP - 1; //subtract one from edges to make them start at 0 
        // fix -1s 
        //for (int e = 0; e < edgesP.nrow(); e++) {
        //    if (edgesP(e, 0) == -2) {
        //        edgesP(e, 0) = -1;
        //        edgesP(e, 1) = -1;
        //    }
        //}
        lengthsP = Rcpp::as<arma::vec>(phylo.slot("l"));
        aliveP = Rcpp::clone(Rcpp::as<LogicalVector>(phylo.slot("alive")));
        tipNamesP = Rcpp::clone(Rcpp::as<CharacterVector>(phylo.slot("tipNames")));
        scaleP = Rcpp::clone(Rcpp::as<NumericVector>(phylo.slot("scale")));
        
        // subtract one from each species to start species at 0 
        indSpeciesL = indSpeciesL - 1;
        
        // n_indv is the number of individuals
        int n_indv = indSpeciesL.length();
        
        // initialize as arma matrix
        traitDiffsSq = arma::mat(n_indv,n_indv);
        
        // add all current squared trait diffs
        for(int i = 0; i < n_indv; i++)
        {
            for(int j = 0; j < n_indv; j++)
            {
                traitDiffsSq(i,j) = pow(indTraitL(i) - indTraitL(j),2);
            }
        }
        
        // initialize founder flag for equilib escape as a vector of 0s of length n_indv 
        founderFlagL = NumericVector(n_indv);
        
        // initialize random number generators
        rng = std::mt19937(std::random_device{}());
        std::uniform_real_distribution<> unif(0, 1);
        std::normal_distribution<double> norm(0, 1);
    }
    
    // copy constructor for roleDataCpp
    roleDataCpp(const roleDataCpp &rhs) {
        
        this->indTraitL = rhs.indTraitL;
        this->indSpeciesL = rhs.indSpeciesL;
        this->spAbundL = rhs.spAbundL;
        this->spTraitL = rhs.spTraitL;
        
        this->spAbundHarmMeanL = rhs.spAbundHarmMeanL;
        this->spLastOriginStepL = rhs.spLastOriginStepL;
        
        this->spExtinctionStepL = rhs.spExtinctionStepL;
        this->spReciprSumL = rhs.spReciprSumL;
        
        this->spAbundM = rhs.spAbundM;
        this->spTraitM = rhs.spTraitM;
        
        this->nTipsP = rhs.nTipsP;
        this->edgesP = rhs.edgesP;
        
        this->lengthsP = rhs.lengthsP;
        this->aliveP = rhs.aliveP;
        
        this->tipNamesP = rhs.tipNamesP;
        this->scaleP = rhs.scaleP;
        
        this->traitDiffsSq = rhs.traitDiffsSq;
        this->envFilterProbs = rhs.envFilterProbs;
        this->prevEnvSigma = rhs.prevEnvSigma;
    }
    
    // copy assignment operator
    roleDataCpp& operator=(const roleDataCpp &rhs) {
        
        this->indTraitL = rhs.indTraitL;
        this->indSpeciesL = rhs.indSpeciesL;
        this->spAbundL = rhs.spAbundL;
        this->spTraitL = rhs.spTraitL;
        
        this->spAbundHarmMeanL = rhs.spAbundHarmMeanL;
        this->spLastOriginStepL = rhs.spLastOriginStepL;
        
        this->spExtinctionStepL = rhs.spExtinctionStepL;
        this->spReciprSumL = rhs.spReciprSumL;
        
        this->spAbundM = rhs.spAbundM;
        this->spTraitM = rhs.spTraitM;
        
        this->nTipsP = rhs.nTipsP;
        this->edgesP = rhs.edgesP;
        
        this->lengthsP = rhs.lengthsP;
        this->aliveP = rhs.aliveP;
        
        this->tipNamesP = rhs.tipNamesP;
        this->scaleP = rhs.scaleP;
        
        this->traitDiffsSq = rhs.traitDiffsSq;
        this->envFilterProbs = rhs.envFilterProbs;
        this->prevEnvSigma = rhs.prevEnvSigma;
    };
};

RCPP_EXPOSED_CLASS(roleDataCpp)