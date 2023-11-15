#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <random>
#include <fstream>
using namespace Rcpp;
using namespace arma;
using namespace std;


void writeToFile(ostream& myfile, std::string text) {
    myfile << text << "\n";
}


void foo(std::string s, int n,  std::string file = "help.txt") {
    ofstream myfile;
    myfile.open(file.c_str());
    
    std::string text = std::to_string(n);
    
    writeToFile(myfile, text+": "+s);
}


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
    IntegerVector ii = seq(1, n);
    Function f = p.slot(s);
    NumericVector x = f(ii);
    
    return x;
}

// function to update phylo objects
// maybe we don't need to return anything, maybe pointers would work?
List updatePhylo(int i, int sMax, double scale, imat edge, vec edgeLength,
                 LogicalVector alive, StringVector tipNames) {

    // index of where unrealized edges in edge matrix start
    int eNew = 2 * sMax - 2;

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
    alive(sMax) = true;
    
    // increase all tip edge lengths by 1 time step
    IntegerVector x = as<IntegerVector>(wrap(edge.col(1)));
    IntegerVector y = seq_len(alive.length());
    IntegerVector z = y[alive];
    
    LogicalVector ind = in(x, z);
    
    uvec tipi = find((edge.col(1) <= eNew) && as<uvec>(ind));
    
    edgeLength(tipi) += 1 * scale;
    
    // update names
    tipNames(sMax) = "s" + std::to_string(sMax + 1);
    
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
    NumericVector metaAbund; // passed
    NumericMatrix metaTrt; // passed
    imat edge; // passed
    vec edgeLength; // passed
    StringVector tipNames; // passed
    LogicalVector alive; // passed
    S4 params; // passed
    int sMax; // passed
    double sigC; // from params.slot("sigC")
    double sigE; // from params.slot("sigE")
    double sig; // from params.slot("trait_sigma")
    double delta; //
    double gamma; //
    NumericVector immProb; // from params.slot("imm")
    NumericVector specProb; // from params.slot("speciation_local")
    mat envOptim; // from params.slot("envOptim")
    mat compMat; // from localTrt_ and sigC
    vec envDist; // from localTrt_ and sigE

public:
    roleComm(IntegerVector localSpp_,
             mat localTrt_,
             NumericVector metaAbund_,
             NumericMatrix metaTrt_,
             imat edge_,
             vec edgeLength_,
             StringVector tipNames_,
             LogicalVector alive_,
             int sMax_,
             S4 params_) :
    rng(std::mt19937(std::random_device{}())),
    dist(0, 1),
    localSpp(localSpp_),
    localTrt(localTrt_),
    // locs.slot("spAbundHarmMean") = 1; // *** need to add to simulation
    // locs.slot("spLastOriginStep") = 1; // *** need to add to simulation
    // locs.slot("spExtinctionStep") = 1; // *** need to add to simulation
    // locs.slot("equilibProp") = 1; // *** need to add to simulation
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
    envOptim(getEnvOptim(params_)),
    compMat(compMatCalc(localTrt_, sigC)),
    envDist(envDistCalc(localTrt_, envOptim, sigE)) {}

    // `get` methods
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
        List locs = List::create(Named("indSpecies") = localSpp,
                                 Named("indTrait") = wrap(localTrt));
        // should be more stuff in above ^

        List meta = List::create(Named("spAbund") = metaAbund,
                                 Named("spTrait") = metaTrt);

        List phylo = List::create(Named("n") = sMax,
                                  Named("e") = wrap(edge),
                                  Named("l") = wrap(edgeLength),
                                  Named("tipNames") = wrap(tipNames),
                                  Named("alive") = alive);

        List out = List::create(Named("localComm") = locs,
                                Named("metaComm") = meta,
                                Named("phylo") = phylo);

        return out;
    }

    // process methods
    int death() {
        // index of the dead individual
        int idead;

        if (delta == 1) {
            // fully neutral
            idead = sample(localSpp.size(), 1)[0] - 1;
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
            probs = gamma * (1 - envDist) + (1 - gamma) * compD;

            // sample index of dead individual
            idead = sample(localSpp.size(), 1, false, wrap(probs))[0] - 1;
        }

        return idead;
    }

    void birthImm(int i, int step) {
        // set up indexes
        int inew;
        int iborn;

        // random number to determine which event occurs
        double r = dist(rng);

        // row vec to hold new trait, initialized with random noise
        rowvec newTrt = randn<rowvec>(localTrt.n_cols) * sig *
            2 / localSpp.size(); // re-scale by generation time

        if (r < immProb[step]) { // immigration
            // sample the spp ID of the immigrating individual
            inew = sample(metaAbund.size(), 1, false, metaAbund)[0] - 1;

            // update traits from meta comm
            newTrt += metaTrt.row(inew);
        } else { // local birth
            // sample the individual that gives birth
            iborn = sample(localSpp.size(), 1)[0] - 1;

            // spp ID of individual that gave birth
            inew = localSpp[iborn];

            // update traits from local comm
            newTrt += localTrt.row(inew);
        }

        // update local comm spp ID
        localSpp[i] = inew;

        // update traits
        localTrt.row(i) = newTrt;
    }

    void speciation(int i, int step) {
        // random number to determine if speciation happens
        double r = dist(rng);

        if (r < specProb[step]) {
            // update phylo
            // Rcout << "iteration is " << step << std::endl;
            // Rcout << "yes speciation; r = " << r << std::endl;
            // Rcout << "specProb = " << specProb[step] << std::endl;
            // Rcout << "specProb size is " << specProb.size() << std::endl;
            
            // scale factor converting time steps in the local comm to generations
            double scale = 2 / localSpp.length();
            
            // run the method to update the phylo
            List newPhyInfo = updatePhylo(i, sMax, scale, edge, edgeLength, alive, 
                                          tipNames);
            
            // not ideal that we have to cast these things with as<type>
            // *** consider updating
            edge = as<imat>(newPhyInfo["edge"]); 
            edgeLength = as<vec>(newPhyInfo["edgeLength"]);
            tipNames = newPhyInfo["tipNames"];
            alive = newPhyInfo["alive"];

            // update total number of spp
            sMax = newPhyInfo["sMax"];

            // update ID of local individual
            localSpp[i] = sMax;

            // update traits
            rowvec newTrt = localTrt.row(i) +
                randn<rowvec>(localTrt.n_cols) * sig; 
                // could re-scale lineage duration
        }
    }

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

            // Rcout << "dist calcs" << std::endl;
            
            // update env dist
            envDist.row(i) = envDistCalc(localTrt.row(i), envOptim, sigE);
            // Rcout << "dist update" << std::endl;
        }
    }
};


// function that takes a `roleData` S4 object and a `roleParams` S4 object, 
// and creates a roleComm rcpp object
roleComm roleCommFromS4(S4 x, S4 p) {
    // local comm stuff
    S4 locs = x.slot("localComm");
    IntegerVector localSpp_ = locs.slot("indSpecies");
    mat localTrt_ = as<mat>(locs.slot("indTrait"));

    // meta comm stuff
    S4 meta = x.slot("metaComm");
    NumericVector metaAbund_ = meta.slot("spAbund");
    NumericMatrix metaTrt_ = meta.slot("spTrait");

    // phylo stuff
    S4 phy = x.slot("phylo");
    imat edge_ = as<imat>(phy.slot("e"));
    vec edgeLength_ = as<vec>(phy.slot("l"));
    StringVector tipNames_ = phy.slot("tipNames");
    LogicalVector alive_ = phy.slot("alive");
    int sMax_ = as<int>(phy.slot("n"));
    
    // decrement species indeces (so they start at 0)
    localSpp_ = localSpp_ - 1;
    // edge_ // need to for phylo stuff????
    // edgeLength_

    // params
    // S4 params_ = x.slot("params");

    roleComm out = roleComm(localSpp_,
                            localTrt_,
                            metaAbund_,
                            metaTrt_,
                            edge_,
                            edgeLength_,
                            tipNames_,
                            alive_,
                            sMax_,
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
    
    locs.slot("indSpecies") = locList["indSpecies"];
    
    locs.slot("indTrait") = locList["indTrait"];
    locs.slot("indSeqs") = "A"; // what to do? make NULL?
    locs.slot("spGenDiv") = 1; // what to do? make NULL?
    // locs.slot("spTrait") = 1; // remove?
    // locs.slot("spAbund") = 1; // remove?
    locs.slot("spAbundHarmMean") = 1; // *** need to add to simulation
    locs.slot("spLastOriginStep") = 1; // *** need to add to simulation
    locs.slot("spExtinctionStep") = 1; // *** need to add to simulation
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
    
    phy.slot("n") = phyList["n"]; // might be sMax, not n
    
    phy.slot("e") = phyList["e"];
    phy.slot("l") = phyList["l"];
    phy.slot("alive") = phyList["alive"];
    phy.slot("tipNames") = "A"; // what to do? remove? or *make intentional?*
    
    out.slot("phylo") = phy;
    
    return out;
}


// tester function wrapping the updatePhylo fun
// [[Rcpp::export]]
S4 testUpdatePhylo(List tre, int i, double scale) {
    // *** don't need to make this list, could just pass right to `updatePhylo`
    List x = List::create(Named("edge") = tre["edge"], 
                          Named("tipNames") = tre["tip.label"],
                          Named("n") = tre["n"],
                          Named("alive") = tre["alive"],
                          Named("edgeLength") = tre["edge.length"],
                          Named("tipNames") = tre["tip.label"]);
    

    List newTre = updatePhylo(i, x["n"], scale, x["edge"], x["edgeLength"], 
                                 x["alive"], x["tipNames"]);
    
    
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
    // consider alternatives to clone????
    x = clone(x);
    roleComm wow = roleCommFromS4(x, p);

    // get params
    // S4 p = wow.getParams();

    // number of iterations and output timesteps
    int niter = as<int>(p.slot("niter"));
    int niterTimestep = as<int>(p.slot("niterTimestep"));
    int n = niter / niterTimestep + 1; // number of output values
    int k; // index for filling in output list `l`

    // list to hold output, each element will be of class `roleComm`
    List l(n);

    // record initial state
    // do we want to output a list or the s4 `roleData` object?
    // probably should output `roleData`
    l[0] = clone(s4FromRcpp(wow.getData()));

    // main loop of sim---starts at 1 because we already recorded the
    // initial state
    for (int i = 1; i <= niter; i++) {
        // death
        int idead = wow.death();
        // foo("dead", i);

        // immigration or local birth
        wow.birthImm(idead, i - 1); // pass i - 1 because loop starts at 1
        // foo("birth", i);
        
        // speciation or not
        wow.speciation(idead, i - 1); // pass i - 1 because loop starts at 1
        // foo("spec", i);

        // update distances
        wow.updateDist(idead);
        // foo("dist", i);

        // every `niterTimestep`, record state
        if (i % niterTimestep == 0) {
            k = i / niterTimestep;
            
            l[k] = clone(s4FromRcpp(wow.getData()));
            // foo("write-out", i);
        }
    }

    return l;
}



void fun(double num) {
    Rcpp::stop("Exception occured!");
}


double takeLog(double val) {
    try {
        fun(val);
    } catch(std::exception &ex) {
        // throw std::range_error("fuuuuuuuuck");
        Rcout << "The value is \n" << val << std::endl;
        return 10;
        // forward_exception_to_r(ex);
    } catch(...) { 
        ::Rf_error("c++ exception (unknown reason)"); 
    }
    return NA_REAL;             // not reached
}

// [[Rcpp::export]]
NumericVector wtf(NumericVector x) {
    x = x - 1;
    return x;
}
