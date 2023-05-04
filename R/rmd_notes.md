RMD notes

------------------------------------------------------------------------

On install, warnings:

- "non-void function does not return a value"
- "no definition for class" roleData, roleExperiment, phylo
- "possible unused arguments" at roleData.20, .40

```
* Updated R/RcppExports.R
There were 18 warnings (use warnings() to see them)

==> R CMD INSTALL --preclean --no-multiarch --with-keep.source roleR

* installing to library ‘/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library’
* installing *source* package ‘roleR’ ...
** using staged installation
clang++ -arch arm64 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c RcppExports.cpp -o RcppExports.o
** libs
clang++ -arch arm64 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c iterModelCpp.cpp -o iterModelCpp.o
In file included from iterModelCpp.cpp:4:
./roleDataCpp.cpp:180:5: warning: non-void function does not return a value [-Wreturn-type]
    };
    ^
iterModelCpp.cpp:688:49: warning: implicit conversion of NULL constant to 'int' [-Wnull-conversion]
                NumericVector probs=NULL, int x=NULL) {
                                               ~^~~~
                                                0
iterModelCpp.cpp:698:1: warning: non-void function does not return a value in all control paths [-Wreturn-type]
}
^
iterModelCpp.cpp:702:55: warning: implicit conversion of NULL constant to 'int' [-Wnull-conversion]
                       RObject params=NULL, int niter=NULL, int i=NULL, //used universally
                                                     ~^~~~
                                                      0
iterModelCpp.cpp:702:67: warning: implicit conversion of NULL constant to 'int' [-Wnull-conversion]
                       RObject params=NULL, int niter=NULL, int i=NULL, //used universally
                                                                 ~^~~~
                                                                  0
iterModelCpp.cpp:703:39: warning: implicit conversion of NULL constant to 'int' [-Wnull-conversion]
                       int dead_index=NULL, // used universally
                                     ~^~~~
                                      0
iterModelCpp.cpp:704:40: warning: implicit conversion of NULL constant to 'int' [-Wnull-conversion]
                       int parent_indv=NULL, // used by call_birth and call_dispersal
                                      ~^~~~
                                       0
iterModelCpp.cpp:705:49: warning: implicit conversion of NULL constant to 'bool' [-Wnull-conversion]
                       bool dispersed_this_iter=NULL, // used by call_speciation and update_speciation_local_meta
                                               ~^~~~
                                                false
iterModelCpp.cpp:706:42: warning: implicit conversion of NULL constant to 'int' [-Wnull-conversion]
                       int speciation_sp=NULL) { // used in update_speciation_local_meta
                                        ~^~~~
                                         0
iterModelCpp.cpp:738:1: warning: non-void function does not return a value in all control paths [-Wreturn-type]
}
^
iterModelCpp.cpp:743:57: warning: implicit conversion of NULL constant to 'int' [-Wnull-conversion]
                         RObject params=NULL, int niter=NULL, int i = NULL){ // used universally 
                                                       ~^~~~
                                                        0
iterModelCpp.cpp:743:71: warning: implicit conversion of NULL constant to 'int' [-Wnull-conversion]
                         RObject params=NULL, int niter=NULL, int i = NULL){ // used universally 
                                                                    ~ ^~~~
                                                                      0
iterModelCpp.cpp:759:1: warning: non-void function does not return a value in all control paths [-Wreturn-type]
}
^
clang++ -arch arm64 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c randExample.cpp -o randExample.o
13 warnings generated.
clang++ -arch arm64 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c roleDataCpp.cpp -o roleDataCpp.o
roleDataCpp.cpp:180:5: warning: non-void function does not return a value [-Wreturn-type]
    };
    ^
clang++ -arch arm64 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c roleParamsCpp.cpp -o roleParamsCpp.o
1 warning generated.
clang++ -arch arm64 -std=gnu++11 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -o roleR.so RcppExports.o iterModelCpp.o randExample.o roleDataCpp.o roleParamsCpp.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/opt/R/arm64/gfortran/lib/gcc/aarch64-apple-darwin20.6.0/12.0.1 -L/opt/R/arm64/gfortran/lib -lgfortran -lemutls_w -lquadmath -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
installing to /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/00LOCK-roleR/00new/roleR/libs
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
in method for ‘getFinalState’ with signature ‘"roleModel"’: no definition for class “roleModel”
in method for ‘getFinalState’ with signature ‘"roleExperiment"’: no definition for class “roleExperiment”
in method for ‘getSumStats’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘getSumStats’ with signature ‘"roleModel"’: no definition for class “roleModel”
in method for ‘getSumStats’ with signature ‘"roleExperiment"’: no definition for class “roleExperiment”
in method for ‘hillAbund’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘hillGenetic’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘hillTrait’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘hillPhylo’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘richness’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘rawAbundance’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘rawSpAbundance’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘rawSppID’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘rawTraits’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘rawGenDiv’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘rawSeqs’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘rawBranchLengths’ with signature ‘"roleData"’: no definition for class “roleData”
in method for ‘rawApePhylo’ with signature ‘"roleData"’: no definition for class “roleData”
Warning: undefined slot classes in definition of "roleData": phylo(class "rolePhylo")
in method for ‘[’ with signature ‘"roleExperiment"’: no definition for class “roleExperiment”
in method for ‘$’ with signature ‘"roleExperiment"’: no definition for class “roleExperiment”
in method for ‘show’ with signature ‘object="roleExperiment"’: no definition for class “roleExperiment”
in method for ‘coerce’ with signature ‘"roleModel","roleExperiment"’: no definition for class “roleModel”
in method for ‘show’ with signature ‘object="roleModel"’: no definition for class “roleModel”
Warning: undefined slot classes in definition of "roleModel": params(class "roleParams")
in method for ‘coerce’ with signature ‘"phylo","rolePhylo"’: no definition for class “phylo”
in method for ‘coerce’ with signature ‘"rolePhylo","phylo"’: no definition for class “phylo”
Note: possible error in 'localComm(data$local$abundance_indv, ': unused arguments (data$local$traits_sp, data$local$pi_sp) at roleData.R:29 
Note: possible error in 'localComm(data$local$abundance_indv, ': unused arguments (data$local$traits_sp, data$local$pi_sp) at roleData.R:40 
** help
*** installing help indices
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (roleR)

```