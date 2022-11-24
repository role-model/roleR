
library(reticulate)
out <- tryCatch(conda_install("r-reticulate",
                              packages="newick",
                              channel="bioconda"),error=function(e) {})
conda_install("r-reticulate",
              packages="msprime",
              channel="conda-forge")
conda_install("r-reticulate",
              packages=c("python-newick","newick_utils"),
              channel="bioconda")

conda_install("r-reticulate",
              packages="newick_utils",
              channel="bioconda")