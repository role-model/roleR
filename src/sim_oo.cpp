#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <random>
// #include <fstream>

using namespace Rcpp;
using namespace arma;
using namespace std;


// function to calculate comp matrix
mat compMatCalc(mat x, double sigC) {
    // make a matrix to hold distances
    int n = x.n_rows;
    mat D(n, n);
    D.fill(0.0);

    // init a double to hold one dist
    double trtDist;

    // loop over rows and cols to get distances with `norm` function
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            trtDist = norm(x.row(i) - x.row(j));
            D(i, j) = exp(-pow(trtDist / sigC, 2)); // Gaussian comp kernel
            D(j, i) = D(i, j);
        }
    }
    
    return D;
}


// exposing `compMatCalc` to R for testing purposes
// [[Rcpp::export]]
NumericMatrix compMatCalcTest(NumericMatrix x, double sigC) {
    mat X = as<mat>(x);
    mat m = compMatCalc(X, sigC);

    return wrap(m);
}


// calculate environmental distances
vec envDistCalc(mat x, mat envOptim, double sigE) {
    int n = x.n_rows;
    vec D(n);

    double eDist;

    for (int i = 0; i < n; i++) {
        eDist = norm(x.row(i) - envOptim);
        D(i) = exp(-pow(eDist / sigE, 2));
    }

    return D;
}

// expose `envDistCalc` to R for testing
// [[Rcpp::export]]
NumericVector envDistCalcTest(NumericMatrix x, NumericMatrix envOptim,
                              double sigE) {
    mat X = as<mat>(x);
    mat eo = as<mat>(envOptim);
    // eo = reshape(eo, 1, X.n_cols);

    vec m = envDistCalc(X, eo, sigE);

    return wrap(m);
}


// extract environmental optimum
mat getEnvOptim(S4 x) {
    NumericMatrix mm = as<NumericMatrix>(x.slot("env_optim"));
    mat m = as<mat>(mm);

    return m;
}

// get a param when we know it will be of class `double`
double getParamDouble(S4 p, String s) {
    double x = as<double>(p.slot(s));

    return x;
}

// get a param vector from a function input
// [[Rcpp::export]]
NumericVector getParamFun(S4 p, String s) {
    int n = p.slot("niter");
    n ++; // want to run loop (i=1; i<=`niter`; i++) but for original `niter`
    IntegerVector ii = seq(1, n);
    Function f = p.slot(s);
    NumericVector x = f(ii);
    
    return x;
}


// function to update phylo objects
// maybe we don't need to return anything, maybe pointers would work?
List updatePhylo(int i, int sMax, double scale, imat edge, vec edgeLength,
                 std::vector<bool> alive, std::vector<std::string> tipNames) {
    // index of where unrealized edges in edge matrix start
    int eNew = 2 * sMax - 2;
    
    // check if there is room in the objects for new edges/nodes/tips,
    // if not, then make room
    if (eNew >= edgeLength.size()) {
        // add edges to all relevant objects
        int edgeAdd = 100;
        int startRow = edge.n_rows;
        
        edge.insert_rows(startRow, edgeAdd);
        edge.rows(startRow, startRow + edgeAdd - 1).fill(-1);
        
        edgeLength.resize(startRow + edgeAdd);
        edgeLength.subvec(startRow, edgeLength.n_elem - 1).fill(-1.0);
        
        alive.resize(alive.size() + edgeAdd, false);
        tipNames.resize(tipNames.size() + edgeAdd, "");
    }

    // index of the edge matrix of where to add new edge
    uvec inds = find(edge.col(1) == i);
    int j = inds(0);
    
    // add one to internal nodes
    uvec internalNode = find(edge > sMax); 
    edge.elem(internalNode) += 1;
    
    // add new internal node
    int newNode = 2 * sMax + 1; // index of new node
    edge(eNew, 0) = newNode;
    edge(1 + eNew, 0) = newNode;
    
    // add tips
    edge(eNew, 1) = edge(j, 1); // add old tip
    edge(eNew + 1, 1) = sMax + 1; // add new tip
    
    // update ancestry of internal nodes
    edge(j, 1) = newNode;

    // augment edge lengths
    edgeLength[eNew] = 0;
    edgeLength[1 + eNew] = 0;
    
    // update alive vector
    alive[sMax] = true;
    
    // increase all tip edge lengths by 1 time step
    IntegerVector x = as<IntegerVector>(wrap(edge.col(1)));
    IntegerVector y = seq_len(alive.size());
    LogicalVector alive_rcpp = Rcpp::wrap(alive);
    IntegerVector z = y[alive_rcpp];
    
    LogicalVector ind = in(x, z);
    
    uvec tipi = find((edge.col(1) <= eNew) && as<uvec>(ind));
    
    edgeLength(tipi) += 1 * scale;
    
    // update names
    tipNames[sMax] = "s" + std::to_string(sMax + 1);
    
    // update sMax
    sMax++;
    
    List out = List::create(Named("edge") = edge,
                            Named("edgeLength") = edgeLength,
                            Named("alive") = alive,
                            Named("sMax") = sMax, 
                            Named("tipNames") = tipNames);
    
    return out;
}


// class to hold all objects of a role model
// with methods defined for updating those objects
class roleComm {
private:
    std::mt19937 rng; // generator for rand unif nums
    std::uniform_real_distribution<double> dist; // unif dist object
    IntegerVector localSpp; // passed
    mat localTrt; // passed
    NumericVector spAbund; // passed
    NumericVector harmMean; // passed
    NumericVector pHarmMean; // passed
    IntegerVector lastOriginStep; // passed
    NumericVector lastOriginGen; // passed
    IntegerVector lastUpdate; // passed
    NumericVector invSum; // passed
    NumericVector invSumP; // passed
    NumericVector metaAbund; // passed
    NumericMatrix metaTrt; // passed
    imat edge; // passed
    vec edgeLength; // passed
    std::vector<string> tipNames; // passed
    std::vector<bool> alive; // passed
    S4 params; // passed
    int sMax; // passed
    double sigC; // from params.slot("sigC")
    double sigE; // from params.slot("sigE")
    double sig; // from params.slot("trait_sigma")
    double delta; // passed
    double gamma; // passed
    NumericVector immProb; // from params.slot("imm")
    NumericVector specProb; // from params.slot("speciation_local")
    NumericVector J; // from params.slot("individuals_local")
    mat envOptim; // from params.slot("envOptim")
    mat compMat; // from localTrt_ and sigC
    vec envDist; // from localTrt_ and sigE
    double gens; // passed

public:
    roleComm(IntegerVector localSpp_,
             mat localTrt_,
             NumericVector spAbund_,
             NumericVector harmMean_,
             NumericVector pHarmMean_,
             IntegerVector lastOriginStep_,
             NumericVector lastOriginGen_,
             IntegerVector lastUpdate_,
             NumericVector invSum_,
             NumericVector invSumP_,
             NumericVector metaAbund_,
             NumericMatrix metaTrt_,
             imat edge_,
             vec edgeLength_,
             std::vector<string> tipNames_,
             std::vector<bool> alive_,
             int sMax_,
             double gens_,
             S4 params_) :
    rng(std::mt19937(std::random_device{}())),
    dist(0, 1),
    localSpp(localSpp_),
    localTrt(localTrt_),
    spAbund(spAbund_),
    harmMean(harmMean_),
    pHarmMean(pHarmMean_),
    lastOriginStep(lastOriginStep_),
    lastOriginGen(lastOriginGen_),
    lastUpdate(lastUpdate_),
    invSum(invSum_),
    invSumP(invSumP_),
    metaAbund(metaAbund_),
    metaTrt(metaTrt_),
    edge(edge_),
    edgeLength(edgeLength_),
    tipNames(tipNames_),
    alive(alive_),
    params(params_),
    sMax(sMax_),
    sigC(getParamDouble(params_, "comp_sigma")),
    sigE(getParamDouble(params_, "env_sigma")),
    sig(getParamDouble(params_, "trait_sigma")),
    delta(getParamDouble(params_, "neut_delta")),
    gamma(getParamDouble(params_, "env_comp_delta")),
    immProb(getParamFun(params_, "dispersal_prob")),
    specProb(getParamFun(params_, "speciation_local")),
    J(getParamFun(params_, "individuals_local")),
    envOptim(getEnvOptim(params_)),
    compMat(compMatCalc(localTrt_, sigC)),
    envDist(envDistCalc(localTrt_, envOptim, sigE)), 
    gens(gens_) {}

    // `get` methods
    int getSMax() {
        return sMax;
    }
    
    NumericVector getSpAbund() {
        return spAbund;
    }
    
    List getLocal() {
        IntegerVector spp = as<IntegerVector>(wrap(localSpp));
        NumericMatrix trt = as<NumericMatrix>(wrap(localTrt));

        return List::create(Named("spp") = spp,
                            Named("trt") = trt);
    }

    // method to extract params
    S4 getParams() {
        return params;
    }

    // get all the data
    List getData() {
        // get all the data, trimmed to remove unused elements
        
        // Find rows where edge has -1 values and keep only valid edges
        uvec validEdgeRows = find(edge.col(0) != -1);
        imat trimmedEdge = edge.rows(validEdgeRows);
        vec trimmedEdgeLength = edgeLength.elem(validEdgeRows);
        
        // use sMax to trim vectors relating to tips 
        std::vector<bool> trimmedAlive;
        auto starta = alive.begin();
        auto enda = alive.begin() + sMax;
        trimmedAlive.assign(starta, enda);
        
        std::vector<string> trimmedTipNames;
        auto startt = tipNames.begin();
        auto endt = tipNames.begin() + sMax;
        trimmedTipNames.assign(startt, endt);
        
        // create the output lists with trimmed data
        List locs = List::create(Named("indSpecies") = localSpp,
                                 Named("indTrait") = wrap(localTrt), 
                                 Named("spAbund") = spAbund,
                                 Named("harmMean") = harmMean, 
                                 Named("pHarmMean") = pHarmMean, 
                                 Named("lastOriginStep") = lastOriginStep, 
                                 Named("lastOriginGen") = lastOriginGen,
                                 Named("gens") = gens);
        
        
        List meta = List::create(Named("spAbund") = metaAbund,
                                 Named("spTrait") = metaTrt);
        
        List phylo = List::create(Named("n") = sMax,
                                  Named("e") = wrap(trimmedEdge),
                                  Named("l") = wrap(trimmedEdgeLength),
                                  Named("tipNames") = wrap(trimmedTipNames),
                                  Named("alive") = trimmedAlive);
        
        List out = List::create(Named("localComm") = locs,
                                Named("metaComm") = meta,
                                Named("phylo") = phylo);
        
        return out;
        
    }

    // process methods
    int death(int step) {
        // index of the dead individual
        int idead;

        if (delta == 1) {
            // fully neutral
            idead = sample(J[step], 1)[0] - 1;
        } else {
            // set up vectors to hold death probs
            vec compD;
            vec probs;

            // gamma determines amount of env filtering v. comp
            if (gamma < 1) {
                // competition calcs
                // the `compMat` is symmetric, so could sum either by row or col
                // (i.e. could pass either 0 or 1 for dim); we pass 1 because 
                // it returns a vec (= colvec) which is what we need elsewhere
                compD = sum(compMat, 1);
            } else {
                // set competition term to 0
                compD = zeros(compMat.n_cols);
            }
            
            // non-neutral death probabilities
            // note: death due to env is captured by `1 - envDist` because
            //       the further you are from the optim, the worse your chances
            // note: delta controls how much neutral drift there is, even when
            //       deterministic forces are in play
            probs = delta + (1 - delta) * (
                gamma * (1 - envDist) + (1 - gamma) * compD
            );

            // sample index of dead individual
            idead = sample(J[step], 1, false, wrap(probs))[0] - 1;
        }
        
        return idead;
    }

    int birthImm(int i, int step) {
        // set up indexes
        int inew; // species level index
        int iborn; // individual level index

        // random number to determine which event occurs
        double r = dist(rng);

        // row vec to hold new trait, initialized with random noise
        rowvec newTrt = randn<rowvec>(localTrt.n_cols) * sig *
            2 / J[step]; // re-scale by generation time

        if (r < immProb[step]) { // immigration
            // sample the spp ID of the immigrating individual
            inew = sample(metaAbund.size(), 1, false, metaAbund)[0] - 1;

            // update traits from meta comm
            newTrt += metaTrt.row(inew);
            
            // update last origin
            // check if already present
            bool newImm = spAbund[inew] == 0;
            // bool newImm = std::find(localSpp.begin(), 
            //                         localSpp.end(), 
            //                         inew) == localSpp.end();
            if (newImm) { // if not present, update orig time
                lastOriginStep[inew] = step;
                lastOriginGen[inew] = 1.0 * gens;
                // Rcout << "inew is: " << inew << std::endl;
                // Rcout << "lastOriginGen is: " << lastOriginGen << std::endl;
                // Rcout << "lastOriginStep is: " << lastOriginStep << std::endl;
                // Rcout << "gens is: " << gens << std::endl;
            }
        } else { // local birth
            // sample the individual that gives birth
            iborn = sample(J[step], 1)[0] - 1;
            
            // update traits from local comm
            newTrt += localTrt.row(iborn);
            
            // spp ID of individual that gave birth
            inew = localSpp[iborn];
        }

        // record sp ID of dead individual
        int idead = localSpp[i];
        
        // update local comm spp ID
        localSpp[i] = inew;

        // update traits
        localTrt.row(i) = newTrt;
        
        // return species-level index of sp ID of who died
        return idead;
    }

    int speciation(int i, int step) {
        // random number to determine if speciation happens
        double r = dist(rng);

        // place holder for ID of new species
        int inew;
        
        if (r < specProb[step]) { // speciation occurs
            // determine parent ID from individual ID
            int iparent = localSpp[i];
            
            // scale factor converting iterations to generations
            double scale = 2.0 / J[step];
            
            // updatePhylo assumes R-style indexing starting at 1, so need
            // to add 1 to `iparent` which has C-style indexing starting at 0
            List newPhyInfo = updatePhylo(iparent + 1, sMax, scale, edge, 
                                          edgeLength, alive, tipNames);
            
            // not ideal that we have to cast these things with as<type>
            // *** consider updating
            edge = as<imat>(newPhyInfo["edge"]); 
            edgeLength = as<vec>(newPhyInfo["edgeLength"]);
            tipNames = as<std::vector<string>>(newPhyInfo["tipNames"]);
            alive = as<std::vector<bool>>(newPhyInfo["alive"]);
            
            // update total number of spp
            sMax = newPhyInfo["sMax"];

            // update ID of local individual
            // sMax is now new_sMax = old_sMax + 1; the new species has
            // C++ 0-indexed ID = sMax - 1 (its R 1-indexed tip index = sMax)
            localSpp[i] = sMax - 1;

            // update traits
            rowvec newTrt = localTrt.row(i) +
                randn<rowvec>(localTrt.n_cols) * sig;

            // add record for last origin
            lastOriginStep[sMax - 1] = step + 1; // +1 cause will pass `step-1`

            inew = sMax - 1;
        } else {
            inew = localSpp[i];
        }
        
        return inew;
    }

    
    // book-keeping methods
    void updateDist(int i) {
        // only need to update distances if we're in a non-neutral sim
        if (delta < 1) {
            // update comp distances
            // note: `envDistCalc` can be used because we're calculating 
            //        distances between 1 ind (= the "optimum") and everybody 
            //        else
            vec newComp = envDistCalc(localTrt, localTrt.row(i), sigC);
            compMat.col(i) = newComp;
            compMat.row(i) = newComp.t();
            
            // update env dist
            envDist.row(i) = envDistCalc(localTrt.row(i), envOptim, sigE);
        }
    }
    
    void updateJ(int step) {
        // change size of local comm if necessary 
        int target_size = J[step];
        // int current_size = localSpp.size();
        // 
        // if (target_size > current_size) { // increasing local comm
        //     for (int fill = current_size; fill < target_size; fill++) {
        //         localSpp.push_back(0); // add empty element
        //         birthImm(fill, step); // fill with birthImm
        //         
        //         // always run `speciation` after `birthImm`
        //         speciation(fill, step); 
        //     }
        // } else if (target_size < current_size) { // decreasing local comm
        //     for (int rm = 0; rm < (current_size - target_size); rm++) {
        //         int death_idx = death(step);
        //         localSpp.erase(localSpp.begin() + death_idx);
        //     }
        // }
        
        // record number of generations
        gens += 2.0 / target_size;
    }
    
    void updateSpInfo(int jdead, int jborn, int step) {
        // update vectors pertaining to sp-level info
        // jdead and jborn are species-level indeces
        // Rcout << "step is: " << step
        //       << "; (pre-update) invSum is: " << invSum
        //       << std::endl;
        
        if (jdead == jborn) {
            // death replaced by birth all in same sp, so no need update `a`
            // this method is also used with `jdead == jborn` for *all* spp
            // on iteration steps where we take a snapshot
            
            // make sure all abundances > 0 (only relevant when we are updating
            // all spp on a snapshot step)
            if (spAbund[jborn] > 0) {
                // update from previous steps
                // if step == (last update step - 1), don't do anything
                // we achieve this by `* 0`
                // *** think carefully about lastUpdate and step indeces
                invSum[jborn] += (1 / spAbund[jborn]) * 
                    (step - lastUpdate[jborn] - 1);
                
                // update inverse sum for the current step
                invSum[jborn] += (1 / spAbund[jborn]);
            }
            
            // update `lastUpdate` to current step
            lastUpdate[jborn] = step;
        } else {
            // update dead one from previous steps
            invSum[jdead] += (1 / spAbund[jdead]) * 
                (step - lastUpdate[jdead] - 1);
            
            // if born one had positive abund prior, 
            // update it from previous steps
            if (spAbund[jborn] > 0) {
                invSum[jborn] += (1 / spAbund[jborn]) * 
                    (step - lastUpdate[jborn] - 1);
            }
            
            // update abundances
            spAbund[jborn] ++;
            spAbund[jdead] --;
            
            // update inv sum for the current step
            invSum[jborn] += 1 / spAbund[jborn];
            
            if (spAbund[jdead] > 0) { // only increment if not extirpated
                invSum[jdead] += 1 / spAbund[jdead];
            } else { // if extirpated, set inv sum to 0
                invSum[jdead] = 0;
            }
            
            // update `lastUpdate` to current step
            lastUpdate[jborn] = step; 
            lastUpdate[jdead] = step;
        }
        
        // Rcout << "step is: " << step
        //       << "; (post-update) invSum is: " << invSum
        //       << std::endl;
    }
    
    void updateHMean(int step) {
        // this method only to be called when writing out state in a 
        // snapshot step
        
        // make this a NumericVector so division works 
        NumericVector dur;
        dur = step - 1.0 * lastOriginStep + 1.0;
        dur[dur < 1] = 1;  // species born at this exact snapshot step: dur=0, use 1

        NumericVector h = dur / invSum;
        h[invSum == 0] = 0;
        harmMean = h;
    }
};


// function that takes a `roleData` S4 object and a `roleParams` S4 object, 
// and creates a roleComm rcpp object
roleComm roleCommFromS4(S4 x, S4 p) {
    // local comm stuff
    S4 locs = x.slot("localComm");
    IntegerVector localSpp_ = locs.slot("indSpecies");
    mat localTrt_ = as<mat>(locs.slot("indTrait"));
    NumericVector spAbund_ = locs.slot("spAbund");
    NumericVector harmMean_ = locs.slot("spAbundHarmMean");
    NumericVector pHarmMean_ = locs.slot("spPropHarmMean");
    IntegerVector lastOriginStep_ = locs.slot("spLastOriginStep");
    NumericVector lastOriginGen_ = locs.slot("spLastOriginGen");
    IntegerVector lastUpdate_ = locs.slot("spLastOriginStep");
    NumericVector invSum_ = 1 / harmMean_; 
    NumericVector invSumP_ = 1 / pHarmMean_; 
    double gens_ = locs.slot("gens");

    // meta comm stuff
    S4 meta = x.slot("metaComm");
    NumericVector metaAbund_ = meta.slot("spAbund");
    NumericMatrix metaTrt_ = meta.slot("spTrait");

    // phylo stuff
    S4 phy = x.slot("phylo");
    imat edge_ = as<imat>(phy.slot("e"));
    vec edgeLength_ = as<vec>(phy.slot("l"));
    std::vector<string> tipNames_ = phy.slot("tipNames");
    std::vector<bool> alive_ = phy.slot("alive");
    int sMax_ = as<int>(phy.slot("n"));
    
    // decrement species indeces (so they start at 0)
    localSpp_ = localSpp_ - 1;
    
    // make sure invSum has 0 instead of Inf
    invSum_[harmMean_ == 0] = 0;

    roleComm out = roleComm(localSpp_,
                            localTrt_,
                            spAbund_,
                            harmMean_, 
                            pHarmMean_, 
                            lastOriginStep_, 
                            lastOriginGen_,
                            lastUpdate_, 
                            invSum_,
                            invSumP_,
                            metaAbund_,
                            metaTrt_,
                            edge_,
                            edgeLength_,
                            tipNames_,
                            alive_,
                            sMax_,
                            gens_,
                            p); // recall: p is passed to `roleCommFromS4`

    return out;
}


// expose roleComm to R for testing
// [[Rcpp::export]]
List roleCommTester(S4 x, S4 p) {
    roleComm wow = roleCommFromS4(x, p);

    List l = List::create(Named("dat") = wow.getData(),
                          Named("pzz") = wow.getParams());
    return l;
}

// function to export data from a `roleComm` object back to S4 class of 
// `roleData`. List argument `x` is assumed to be output from 
// `roleComm::.getData()`
// note: we can export data only because we don't need to export params back 
//       out, they're already stored in other R objects that were supplied to 
//       Rcpp
// [[Rcpp::export]]
S4 s4FromRcpp(List x) {
    S4 out("roleData");
    
    // local comm
    S4 locs("localComm");
    List locList = x["localComm"];
    
    // increment by 1 because we have been in cpp 0-based indexing up until now
    IntegerVector rIndexIndSpp = locList["indSpecies"];
    locs.slot("indSpecies") = rIndexIndSpp + 1;
    
    locs.slot("indTrait") = locList["indTrait"];
    locs.slot("indSeqs") = "A"; // place-holder
    locs.slot("spGenDiv") = 1; // place-holder
    locs.slot("spTajD") = 1; // place-holder
    locs.slot("spAbund") = locList["spAbund"];
    locs.slot("spAbundHarmMean") = locList["harmMean"];
    locs.slot("spPropHarmMean") = locList["pHarmMean"];
    locs.slot("spLastOriginStep") = locList["lastOriginStep"];
    locs.slot("spLastOriginGen") = locList["lastOriginGen"];
    locs.slot("gens") = locList["gens"];
    locs.slot("equilibProp") = 1; // *** need to add to simulation
    
    out.slot("localComm") = locs;
    
    // meta comm
    S4 meta("metaComm");
    List metaList = x["metaComm"];
    
    meta.slot("spAbund") = metaList["spAbund"];
    meta.slot("spTrait") = metaList["spTrait"];
    
    out.slot("metaComm") = meta;
    
    // phylo
    S4 phy("rolePhylo");
    List phyList = x["phylo"];
    
    phy.slot("n") = phyList["n"];
    phy.slot("e") = phyList["e"];
    phy.slot("l") = phyList["l"];
    phy.slot("alive") = phyList["alive"];
    phy.slot("tipNames") = phyList["tipNames"];
    
    out.slot("phylo") = phy;
    
    return out;
}


// tester function wrapping the updatePhylo fun
// [[Rcpp::export]]
S4 testUpdatePhylo(S4 tre, int i, double scale) {
    List newTre = updatePhylo(i, tre.slot("n"), scale, tre.slot("e"), 
                              tre.slot("l"), tre.slot("alive"), 
                              tre.slot("tipNames"));
    
    
    // create S4 output
    S4 phy("rolePhylo");
    
    phy.slot("n") = newTre["sMax"];
    phy.slot("e") = newTre["edge"];
    phy.slot("l") = newTre["edgeLength"];
    phy.slot("alive") = newTre["alive"];
    phy.slot("tipNames") = newTre["tipNames"]; 
    
    return phy;
}




// OO version of simulation function
// `x` is a `roleData` object
// `p` is a `roleParams` object
// [[Rcpp::export]]
List simRole(S4 x, S4 p) {
    NumericVector invSumsInit = 1.3; // initial sum of inverse abunds
    
    roleComm wow = roleCommFromS4(x, p); // now need to pass `invSumsInit`

    // number of iterations and output timesteps
    int niter = as<int>(p.slot("niter"));
    int niterTimestep = as<int>(p.slot("niterTimestep"));
    int n = niter / niterTimestep + 1; // number of output values
    
    // Rcout << "n for model steps is: " << n << std::endl;
    
    // niter++; // so we can loop until i < niter
    
    int k; // index for filling in output list `l`
    
    // object to track max species index
    int thisSMax = wow.getSMax();

    // list to hold output, each element will be of class `roleComm`
    List l(n);

    // record initial state
    wow.updateHMean(0);
    l[0] = clone(s4FromRcpp(wow.getData()));
    
    // main loop of sim---starts at 1 because we already recorded the
    // initial state
    for (int i = 1; i <= niter; i++) {
        // check if local comm size has changed and update accordingly
        wow.updateJ(i);
        
        // death
        int idead = wow.death(i); 

        // immigration or local birth
        // returns species-level index of ind that died
        int oldSp = wow.birthImm(idead, i); 
        
        // speciation or not (returns species index for spp-level update)
        // returns species-level index of ind replacing death
        int inew = wow.speciation(idead, i); 

        // update distances
        wow.updateDist(idead);
        
        // update species-level data of local comm
        wow.updateSpInfo(oldSp, inew, i);

        // every `niterTimestep`, record state
        if ((i) % niterTimestep == 0) {
            // update all sp-level data at this time
            // but only for spp with non-0 abund
            
            thisSMax = wow.getSMax();

            // Rcout << "updateSpInfo called from snapshot write-out" << std::endl;
            for (int j = 0; j < thisSMax; j++) {
                wow.updateSpInfo(j, j, i);
            }
            
            // Rcout << "END call from snapshot write-out" << std::endl;
            
            wow.updateHMean(i);
            
            k = (i) / niterTimestep;
            l[k] = clone(s4FromRcpp(wow.getData())); 
        }
        
        // Rcout << "step is: " << i 
        //       << "; spAbund is: " << wow.getSpAbund()
        //       << std::endl;
    }

    return l;
}

