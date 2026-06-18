library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)

reticulate::source_python("inst/python/msprime_tests.py")

sim_role <- roleR:::get_sim_role()

# pop size ----

ntime <- 6

gens <- seq(0, 1000, length.out = ntime)

tabund <- matrix(c(40, 10), nrow = ntime, ncol = 2, byrow = TRUE) |> 
    as.data.frame() |> 
    setNames(1:2)

tharm <- matrix(c(40, 10), nrow = ntime, ncol = 2, byrow = TRUE) |> 
    as.data.frame() |> 
    setNames(1:2)

tdiv <- matrix(0, nrow = ntime, ncol = 2) |> 
    as.data.frame() |> 
    setNames(1:2)

mabund <- rep(10000, 2)

nwk2 <- "(t1:3, t2:3);"


popExpImm <- replicate(200, {
    sim_role(Jm = sum(mabund), curtime = max(gens), nwk = nwk2, 
                   meta_abund = mabund, local_sad = tabund, 
                   local_hm = tharm, J_harm = 30, 
                   gens = gens, tdiv = tdiv, 
                   alpha = 100, m = 0.01, 
                   sequence_length = 500, mu = 1e-06)[[2]]
}, simplify = FALSE) |> 
    do.call(rbind, args = _) |> 
    mutate(N = ifelse(sp_id == "1", 40, 10) * 100, 
           G_ago = factor(time))

ggplot(popExpImm, aes(pi, fill = factor(N))) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 16) +
    facet_grid(rows = vars(G_ago), scales = "free_y")

ggplot(popExpImm, aes(tajimas_D, fill = factor(N))) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 24) +
    facet_grid(rows = vars(G_ago), scales = "free_y")

summarize(popExpImm, 
          mean_pi = mean(pi, na.rm = TRUE),
          sd_pi = sd(pi, na.rm = TRUE),
          mean_taj_D = mean(tajimas_D, na.rm = TRUE), 
          sd_taj_D = sd(tajimas_D, na.rm = TRUE), 
          .by = c(N, time))


popExpNoImm <- replicate(200, {
    sim_role(Jm = sum(mabund), curtime = max(gens), nwk = nwk2, 
             meta_abund = mabund, local_sad = tabund, 
             local_hm = tharm, J_harm = 30, 
             gens = gens, tdiv = tdiv, 
             alpha = 100, m = 0, 
             sequence_length = 500, mu = 1e-06)[[2]]
}, simplify = FALSE) |> 
    do.call(rbind, args = _) |> 
    mutate(N = ifelse(sp_id == "1", 40, 10) * 100, 
           G_ago = factor(time))

ggplot(popExpNoImm, aes(pi, fill = factor(N))) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 16) +
    facet_grid(rows = vars(G_ago), scales = "free_y")

ggplot(popExpNoImm, aes(tajimas_D, fill = factor(N))) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 24) +
    facet_grid(rows = vars(G_ago), scales = "free_y")

summarize(popExpNoImm, 
          mean_pi = mean(pi, na.rm = TRUE),
          sd_pi = sd(pi, na.rm = TRUE),
          mean_taj_D = mean(tajimas_D, na.rm = TRUE), 
          sd_taj_D = sd(tajimas_D, na.rm = TRUE), 
          .by = c(N, time))


tabund <- matrix(c(25, 25), nrow = ntime, ncol = 2, byrow = TRUE) |> 
    as.data.frame() |> 
    setNames(1:2)

tharm <- matrix(c(25, 25), nrow = ntime, ncol = 2, byrow = TRUE) |> 
    as.data.frame() |> 
    setNames(1:2)

tdiv <- matrix(0, nrow = ntime, ncol = 2) |> 
    as.data.frame() |> 
    setNames(1:2)

mabund <- c(18000, 2000)

popExpMetaSize <- replicate(200, {
    sim_role(Jm = sum(mabund), curtime = max(gens), nwk = nwk2, 
             meta_abund = mabund, local_sad = tabund, 
             local_hm = tharm, J_harm = 30, 
             gens = gens, tdiv = tdiv, 
             alpha = 100, m = 0.01, 
             sequence_length = 500, mu = 1e-06)[[2]]
}, simplify = FALSE) |> 
    do.call(rbind, args = _) |> 
    mutate(N_meta = ifelse(sp_id == "1", 18000, 2000) * 100, 
           G_ago = factor(time))

ggplot(popExpMetaSize, aes(pi, fill = factor(N_meta))) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 16) +
    facet_grid(rows = vars(G_ago), scales = "free_y")

ggplot(popExpMetaSize, aes(tajimas_D, fill = factor(N_meta))) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 24) +
    facet_grid(rows = vars(G_ago), scales = "free_y")

summarize(popExpMetaSize, 
          mean_pi = mean(pi, na.rm = TRUE),
          sd_pi = sd(pi, na.rm = TRUE),
          mean_taj_D = mean(tajimas_D, na.rm = TRUE), 
          sd_taj_D = sd(tajimas_D, na.rm = TRUE), 
          .by = c(N_meta, time))


# re-implementation of island-mainland msprime model
island_pars <- expand.grid(
    N_mainland = mabund[1] * 100,
    N_island = c(tharm[1, 1], tharm[1, 2]) * 100,
    T_found = c(0.000001, gens[-1]),
    n_samples = c(tabund[1, 1], tabund[1, 2]),
    seq_length = 500,
    mu = 1e-06
) |>
    mutate(m_main_to_island = 0.5 * 0.01 * 3000 /
               (2 * N_island)) |> 
    replicate(200, expr = _, simplify = FALSE) |> 
    do.call(rbind, args = _)

sim <- parallel::mclapply(1:nrow(island_pars), 
                          FUN = function(i) {
                              sim_island(
                                  N_mainland = island_pars$N_mainland[i], 
                                  N_island = island_pars$N_island[i], 
                                  m_main_to_island = island_pars$m_main_to_island[i], 
                                  T_found = island_pars$T_found[i], 
                                  n_samples = island_pars$n_samples[i], 
                                  seq_length = island_pars$seq_length[i], 
                                  mu = island_pars$mu[i], 
                                  alpha = 100)
                          },
                          mc.cores = 8)

sim_res <- do.call(rbind, sim) |> 
    mutate(G_ago = factor(T_found))

ggplot(sim_res, aes(pi, fill = factor(N_island))) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 16) +
    facet_grid(rows = vars(G_ago), scale = "free_y") 

ggplot(sim_res, aes(pi, is.nan(tajimas_D))) +
    geom_jitter()

ggplot(popExp, aes(pi, is.nan(tajimas_D))) +
    geom_jitter()

ggplot(sim_res, aes(pi, color = G_ago)) +
    geom_histogram(alpha = 0.25) +
    facet_grid(rows = vars(N_island)) +
    scale_x_continuous(limits = c(-0.05, 0.5))

ggplot(popExp, aes(pi, color = G_ago)) +
    geom_histogram(alpha = 0.25) +
    facet_grid(rows = vars(N)) #+
    # scale_x_continuous(limits = c(-0.05, 0.5))

ggplot(popExp, aes(tajimas_D, color = G_ago)) +
    geom_histogram(alpha = 0.25) +
    facet_grid(rows = vars(N))

# branch length ----

tabund <- (matrix(c(10, 0, 0, 0,
                    10, 10, 2, 10,
                    10, 4, 2, 20), 
                  byrow = TRUE, nrow = 3) * 3) |> 
    as.data.frame() |> 
    setNames(as.character(1:4))

tharm <- apply(1 / tabund, 2, function(x) {
    # y <- numeric(length(x))
    # y[is.finite(x)] <- seq_along(x[is.finite(x)]) / cumsum(x[is.finite(x)])
    # y
    x[!is.finite(x)] <- 2
    seq_along(x) / cumsum(x)
}) |> 
    as.data.frame()


tdiv <- (matrix(c(0, NA, NA, NA, 
                  0, 0.2, 0.9, 0.4, 
                  0, 0.2, 0.9, 0.4), 
                byrow = TRUE, nrow = 3) * 500) |> 
    as.data.frame() |> 
    setNames(as.character(1:4))

mabund <- rep(10000, 4)

nwkShrtRt <- "((t1:3, t2:3):1,(t3:2, t4:2):2);"
nwkLongRt <- "((t1:3, t2:3):9,(t3:2, t4:2):10);"


mspBlExp <- replicate(20, {
    mspShrtRt <- sim_role(Jm = sum(mabund), curtime = 1000, nwk = nwkShrtRt, 
                          meta_abund = mabund, local_sad = tabund, 
                          local_hm = tharm, J_harm = 30, 
                          gens = c(0, 500, 1000), tdiv = tdiv, 
                          alpha = 100, m = 0.05, 
                          sequence_length = 500, mu = 1e-05)
    
    mspLongRt <- sim_role(Jm = sum(mabund), curtime = 1000, nwk = nwkLongRt, 
                          meta_abund = mabund, local_sad = tabund, 
                          local_hm = tharm, J_harm = 30, 
                          gens = c(0, 500, 1000), tdiv = tdiv, 
                          alpha = 100, m = 0.05, 
                          sequence_length = 500, mu = 1e-05)
    
    x <- mutate(tharm * 100, time = c(0, 500, 1000)) |> 
        pivot_longer(!time, names_to = "sp_id", values_to = "hmean") |> 
        right_join(mspShrtRt[[2]][, 1:2])
    
    data.frame(x, 
               shrt = mspShrtRt[[2]]$pi,
               long = mspLongRt[[2]]$pi)
}, simplify = FALSE) |> 
    do.call(rbind, args = _)


ggplot(mspBlExp, aes(shrt, long, color = factor(hmean))) +
    geom_point() +
    scale_color_viridis_d() +
    facet_grid(vars(time))





island_pars <- expand.grid(
    N_mainland = c(100000, 5000),
    N_island = c(20000, 500),
    m_main_to_island = c(0.2, 0),
    T_found = c(10000, 100),
    n_samples = c(1000, 10),
    seq_length = 500,
    mu = c(1e-06, 1e-07)
) |> replicate(10, expr = _, simplify = FALSE) |> 
    do.call(rbind, args = _)

sim <- parallel::mclapply(1:nrow(island_pars), 
                   FUN = function(i) {
                       print(i)
                       sim_island(
                           N_mainland = island_pars$N_mainland[i], 
                           N_island = island_pars$N_island[i], 
                           m_main_to_island = island_pars$m_main_to_island[i], 
                           T_found = island_pars$T_found[i], 
                           n_samples = island_pars$n_samples[i], 
                           seq_length = island_pars$seq_length[i], 
                           mu = island_pars$mu[i])
                   },
                   mc.cores = 8)

sim_res <- do.call(rbind, sim)

ggplot(sim_res, 
       aes(pi, color = factor(N_island))) +
    geom_histogram(alpha = 0.25) +
    facet_grid(rows = vars(mu), cols = vars(N_mainland))

ggplot(sim_res, aes(n_sites, pi, color = mu)) +
    geom_point()

ggplot(sim_res, aes(N_island * mu, pi, color = N_island)) +
    geom_jitter() +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(vars(m_main_to_island)) + 
    geom_abline(intercept = 0, slope = 1) +
    geom_hline(yintercept = 0.75)


sim_one <- parallel::mclapply(1:nrow(island_pars), 
                              FUN = function(i) {
                                  print(i)
                                  sim_single(
                                      N_island = island_pars$N_island[i], 
                                      n_samples = island_pars$n_samples[i], 
                                      seq_length = island_pars$seq_length[i], 
                                      mu = island_pars$mu[i])
                              },
                              mc.cores = 8)

sim_one_res <- do.call(rbind, sim_one)

ggplot(data.frame(sim_res,
                  pi_single = sim_one_res$pi), 
       aes(pi_single, pi, color = N_island)) +
    geom_point() +
    facet_grid(vars(m_main_to_island), vars(N_mainland)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_log10() +
    scale_y_log10()

ggplot(sim_one_res, aes(N_island * mu, pi)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(intercept = 0, slope = 1)




sim_single(N_island = 1000 * 100, T_found = 1000, n_samples = 100, seq_length = 500, mu = 1e-07)


