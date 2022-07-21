// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <string>

using namespace Rcpp;

class roleDataCpp {
public:
    //NumericMatrix indSpTrtL;
    NumericVector indTraitL; 
    NumericVector indSpeciesL; 
    
    //NumericMatrix spAbundTrtL; 
    NumericVector spAbundL;
    NumericVector spTraitL;

    //NumericMatrix spAbundTrtM;
    NumericVector spAbundM;
    NumericVector spTraitM;
    
    int nTipsP;
    NumericMatrix edgesP;
    NumericVector lengthsP;
    LogicalVector aliveP;
    CharacterVector tipNamesP;
    double scaleP;    
    
    NumericVector envFilterProbs;
    arma::mat traitDiffsSq;
    
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
        indTraitL = Rcpp::as<NumericVector>(local.slot("indTrait"));
        indSpeciesL =  Rcpp::as<NumericVector>(local.slot("indSpecies"));
        
        //spAbundTrtL = Rcpp::as<NumericMatrix>(local.slot("spAbundTrt"));
        spAbundL = Rcpp::as<NumericVector>(local.slot("spAbund"));
        spTraitL = Rcpp::as<NumericVector>(local.slot("spTrait"));
        
        //spAbundTrtM = Rcpp::as<NumericMatrix>(meta.slot("sppAbundTrt"));
        spAbundM = Rcpp::as<NumericVector>(meta.slot("spAbund"));
        spTraitM = Rcpp::as<NumericVector>(meta.slot("spTrait"));
        
        
        nTipsP = Rcpp::as<int>(phylo.slot("n"));
        edgesP = Rcpp::as<NumericMatrix>(phylo.slot("e"));
        lengthsP = Rcpp::as<NumericVector>(phylo.slot("l"));
        aliveP = Rcpp::as<LogicalVector>(phylo.slot("alive"));
        tipNamesP = Rcpp::as<CharacterVector>(phylo.slot("tipNames"));
        
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
        
        // below code which creates spAbundTrtL from indSpTrtL is deprecated because it doesnt work for seemingly no reason 
        // (object exists in the constructor,
        // but when accessed in iterSim is empty, absolutely no idea why)
        //Rcout << indSpTrtL << "\n";
        //Rcout << spAbundTrtM << "\n";
        
        // // parse local indv indexed matrix to local species indexed matrix 
        // NumericVector species = Rcpp::unique(indSpTrtL(_,0));
        // NumericMatrix spAbundTrtL(100, 2);
        // std::fill(spAbundTrtL.begin(), spAbundTrtL.end(), NumericVector::get_na()); //NumericVector::get_na()
        // 
        // //Rcout << species.length() << "\n";
        // for(int s  = 0; s < species.length(); s++){
        //     Rcout << s << "\n";
        //     int sp = species[s];
        //     Rcout << sp << "\n";
        //     int abund = 0;
        //     NumericVector traits; 
        //     for(int r = 0; r < indSpTrtL.nrow(); r++){
        //         if(indSpTrtL(r,0) == sp){
        //             abund++; 
        //             traits.push_back(indSpTrtL(r,1));
        //         }
        //     }
        //     //Rcout << abund << "\n";
        //     //Rcout << mean(traits) << "\n";
        //     spAbundTrtL(sp,0) = abund;
        //     spAbundTrtL(sp,1) = mean(traits); 
        // }
        // 
        // //spAbundTrtL = spAbundTrtL;
        // //Rcout << spAbundTrtL(0,0) << "\n";
        // //Rcout << spAbundTrtL(0,1) << "\n";
        // //Rcout << "interior" << spAbundTrtL << "\n";
        // 
        // //nTipsP = Rcpp::as<int>(phylo.slot("n"));
        // //edgesP = Rcpp::as<NumericMatrix>(phylo.slot("e"));
        // //lengthsP = Rcpp::as<NumericVector>(phylo.slot("l"));
        // //aliveP = Rcpp::as<LogicalVector>(phylo.slot("alive"));
        // //tipNamesP = Rcpp::as<CharacterVector>(phylo.slot("tipNames"));
    }
};

RCPP_EXPOSED_CLASS(roleDataCpp)