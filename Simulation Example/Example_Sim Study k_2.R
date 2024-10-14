library(foreach)
library(doSNOW)
library(progress)
library(MatVarMix)
# library(mclust) must be installed

#### Definition of the simulation objects ####

n <- c(250, 500, 1000) # sample sizes.
bal <- c(1, 2) # 1 = equal weights; 2 = unequal weights.
dim <- c(1, 2) # 1 = 2 x 3 matrices; 2 = 7 x 6 matrices.
ov <- c(1, 2) # 1 = close groups, 2 = separated groups.

scenarios <- expand.grid(list(n, bal, dim, ov))
colnames(scenarios) <- c("size", "bal", "dim", "overlap")

### Prepare for parallel computing ###

nThreads <- 24 # number of cores to be set for parallel computing
cluster <- makeCluster(nThreads, type = "SOCK")
registerDoSNOW(cluster)

pb2 <- progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = nrow(scenarios),
  complete = "=", # Completion bar character
  incomplete = "-", # Incomplete bar character
  current = ">", # Current bar character
  width = 100
)

progress <- function(n) {
  pb2$tick()
}
opts <- list(progress = progress)
comb <- function(x, ...) {
  lapply(
    seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
  )
}

#### Simulation and Fitting for K = 2 groups ####

# load here the "All parameters.RData" file containing the parameters before running the simulation #

res_k2 <- foreach(
  l = 1:nrow(scenarios), .combine = "comb", .multicombine = TRUE, .packages = "MatVarMix",
  .init = list(list(), list(), list()), .options.snow = opts
) %dopar% {
  
  nsims <- 500 # number of datasets to simulate
  dat <- iniMV <- fitMV <- vector(mode = "list", length = nsims)

  for (n in 1:nsims) {
    fitMV[[n]] <- iniMV[[n]] <- vector(mode = "list", length = 4) # 2 groups + 2
  } # to be modified if k!=2

  for (j in 1:nsims) {
    check1 <- check2 <- 0

    while (check2 < 1) {
      while (check1 < 1) {
        temp.par <- param[[scenarios[l, 2]]][[scenarios[l, 3]]][[scenarios[l, 4]]]

        dat[[j]] <- r_MVN_MX(num = scenarios[l, 1], pix = temp.par[[1]], M = temp.par[[2]], U = temp.par[[3]], V = temp.par[[4]])

        iniMV.true <- list()
        iniMV.true[["prior"]] <- temp.par[[1]]
        iniMV.true[["M"]] <- temp.par[[2]]
        iniMV.true[["U"]] <- temp.par[[3]]
        iniMV.true[["V"]] <- temp.par[[4]]
        iniMV.true[["class"]] <- dat[[j]]$obs.class

        fitMV.true <- tryCatch(fit_MVN_MX(Y = dat[[j]]$Y, k = 2, init.par = iniMV.true), error = function(e) {
          NA
        }) # to be modified if k!=2
        ari.true <- tryCatch(round(mclust::adjustedRandIndex(dat[[j]]$obs.class, fitMV.true$class), digits = 2), error = function(e) {
          0
        })

        ### Check for ARI Interval ###

        if (l >= 1 & l <= 12) {
          if (ari.true >= 0.80 & ari.true <= 0.85) {
            check1 <- 1
          } 
        } else {
          if (ari.true >= 0.95 & ari.true <= 1) {
            check1 <- 1
          } 
        }
      }

      for (k in 1:4) { # to be modified if k!=2

        iniMV[[j]][[k]] <- tryCatch(init_MVN_MX(Y = dat[[j]]$Y, k = k), error = function(e) {
          NA
        })
        fitMV[[j]][[k]] <- tryCatch(fit_MVN_MX(Y = dat[[j]]$Y, k = k, init.par = iniMV[[j]][[k]]), error = function(e) {
          NA
        })
      }

      ### Check for fitting/computational issues ###

      if (all(!is.na(fitMV[[j]]))) {
        check2 <- 1
      }
    }
  }

  list(dat, fitMV, iniMV)
}
