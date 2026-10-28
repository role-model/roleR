# Run a roleExperiment of many models on Koa (Slurm).
#
# Usage: Rscript run_role_koa.R [n_models] [chunk_size] [outdir]
#
# Models are run in chunks. Each finished chunk is saved to `outdir`, so if
# the job is killed (e.g., time limit) you can resubmit the same command and
# it picks up at the first chunk that isn't saved yet.

library(roleR)

args <- commandArgs(trailingOnly = TRUE)
n_prior_samp <- if(length(args) >= 1) as.integer(args[1]) else 10000L
chunk_size <- if(length(args) >= 2) as.integer(args[2]) else 480L
outdir <- if(length(args) >= 3) args[3] else paste0("role_koa_out_", n_prior_samp)

# take the number of workers from Slurm so the script and allocation agree
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

# when run as a Slurm array (indices 1 to n), each task runs every n-th chunk;
# outside an array this is task 1 of 1 and all chunks are run
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
n_tasks <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# priors ----

# draw once and save, so a resumed job uses exactly the same parameters;
# row i of `param_samples` is model_id i
prior_file <- file.path(outdir, "param_samples.rds")

if(file.exists(prior_file)) {
    param_samples <- readRDS(prior_file)
    stopifnot(nrow(param_samples) == n_prior_samp)
} else {
    set.seed(1)
    
    param_samples <- data.frame(
        individuals_local = 1000, 
        individuals_meta = 100000, 
        species_meta = 300, 
        speciation_local = 0, 
        speciation_meta = 0.100,
        extinction_meta = 0.090, 
        dispersal_prob = runif(n_prior_samp, 0.001, 0.1), 
        env_optim = 0, 
        trait_sigma = runif(n_prior_samp, 0.5, 2), 
        env_sigma = runif(n_prior_samp, 2, 10), 
        comp_sigma = runif(n_prior_samp, 2, 10), 
        neut_delta = rep(c(1, 0.15), ceiling(n_prior_samp / 3) * c(1, 2))[1:n_prior_samp], 
        env_comp_delta = rep(0:1, ceiling(n_prior_samp / 2))[1:n_prior_samp], 
        # neut_delta = runif(n_prior_samp), 
        # env_comp_delta = runif(n_prior_samp), 
        mutation_rate = 1e-09, 
        equilib_escape = 1, 
        alpha = 1e+04, 
        num_basepairs = 350, 
        init_type = "oceanic_island",
        niter = 50000,
        niterTimestep = 50000
    )
    
    # write via a temp file and rename, so tasks starting at the same time
    # never see a half-written file
    tmpFile <- paste0(prior_file, ".tmp", task_id)
    saveRDS(param_samples, tmpFile)
    file.rename(tmpFile, prior_file)
}

# functions ----

# list of summary stat functions to be used
funList <- list(
    hill_abund = hillAbund, 
    hill_gen = hillGenetic,
    hill_trait = hillTrait,
    richness = richness,
    mean_taj_D = meanTajD, 
    hill_taj_D = hillTajD, 
    cor_abund_pi = corAbundPi,
    cor_taj_D_abund = corTajDAbund,
    cor_taj_D_pi = corTajDPi
)

# run in chunks ----

starts <- seq(1, n_prior_samp, by = chunk_size)
chunk_files <- file.path(outdir, sprintf("chunk_%03d.rds", seq_along(starts)))
tStart <- Sys.time()

myChunks <- which((seq_along(starts) - 1) %% n_tasks == task_id - 1)

message(n_prior_samp, " models in ", length(starts), " chunks of up to ", 
        chunk_size, ", on ", n_cores, " cores; task ", task_id, " of ", 
        n_tasks, " runs ", length(myChunks), " chunks")

for(k in myChunks) {
    if(file.exists(chunk_files[k])) {
        message("chunk ", k, " already done, skipping")
        next
    }
    
    ids <- starts[k]:min(starts[k] + chunk_size - 1, n_prior_samp)
    
    # only build the models for this chunk, so memory doesn't scale with
    # the total number of models
    params <- lapply(ids, function(i) {
        as.list(param_samples[i, ]) |> 
            do.call(roleParams, args = _)
    })
    
    t0 <- Sys.time()
    set.seed(1000 + k)
    e <- roleExperiment(params, funs = funList)
    erun <- runRole(e, cores = if(n_cores > 1) n_cores else NULL)
    
    # keep only the summary stats, with model ids that are global rather than
    # numbered within the chunk
    res <- erun@results
    res$model_id <- ids[res$model_id]
    
    # write to a temp file first so a killed job can't leave a partial chunk
    tmpFile <- paste0(chunk_files[k], ".tmp")
    saveRDS(res, tmpFile)
    file.rename(tmpFile, chunk_files[k])
    
    message("chunk ", k, " of ", length(starts), " done in ", 
            format(round(difftime(Sys.time(), t0, units = "mins"), 1)), 
            "; total elapsed ", 
            format(round(difftime(Sys.time(), tStart, units = "mins"), 1)))
    
    rm(params, e, erun, res)
    invisible(gc())
}

# combine ----

if(all(file.exists(chunk_files))) {
    allRes <- do.call(rbind, lapply(chunk_files, readRDS))
    
    # temp file and rename, in case several tasks finish together
    tmpFile <- file.path(outdir, paste0("results_all.rds.tmp", task_id))
    saveRDS(allRes, tmpFile)
    file.rename(tmpFile, file.path(outdir, "results_all.rds"))
    
    nFailed <- length(unique(allRes$model_id[allRes$run_status != "success"]))
    message("all chunks done: ", length(unique(allRes$model_id)), 
            " models, ", nFailed, " failed")
} else {
    message("some chunks are missing, resubmit to finish them")
}
