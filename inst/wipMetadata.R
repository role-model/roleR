# wip metadata code
values <- c(sim$params$values$dispersal_prob, sim$params$values$extinction_meta, 
            sim$params$values$individuals_local, sim$params$values$individuals_meta,
            sim$params$values$speciation_local, sim$params$values$speciation_meta,
            sim$params$values$species_meta, sim$params$values$trait_sigma, 1000) 
names <- c("dispersal prob param","metacomm extinction param",
           "num individuals local param", "num individuals meta param", 
           "localcomm speciation rate param", "metacomm speciation rate param",
           "num species meta param", "trait sigma param","num iterations")
metadata <- data.frame(values)
colnames(metadata) <- "Value"
rownames(metadata) <- names
View(metadata)

# paramter value = inverse gamma 
# expression, text, parse etc to
# maybe executable code? 
# how to capture priors in metadata
# optional field - user specified description of project 
# keyword field - based on description 
# version information R, roleR 