### For automatically reducing K methods.

pmclust.reduceK <- function(X = NULL, K = 2, MU = NULL,
    algorithm = .PMC.CT$algorithm, RndEM.iter = .PMC.CT$RndEM.iter,
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  ### Run through original pmclust().
  ret <- pmclust(X = X, K = K, MU = MU, algorithm = algorithm,
                 RndEM.iter = RndEM.iter, CONTROL = CONTROL,
                 method.own.X = method.own.X, rank.own.X = rank.own.X,
                 comm = comm)

  ### Repeat if error occurs.
  repeat{
    if(ret$check$convergence == 99 && K > 1){
      ### Drop specific i.k if available or
      ### drop the smallest class or
      ### drop the class with the smallest eta among all small classes or
      ### drop all classes with 0 elements.
      PARAM.new <- ret$param
      if(.pmclustEnv$CONTROL$stop.at.fail && .pmclustEnv$FAIL.i.k > 0){
        i.k <- .pmclustEnv$FAIL.i.k
      } else{
        i.k <- which(ret$n.class == min(ret$n.class))
      }
      if(i.k > 1 && min(ret$n.class) > 0){
        i.k <- i.k[which.min(PARAM.new$ETA[i.k])]
      }
      K <- K - length(i.k)

      ### Initial global storage.
      if(algorithm[1] %in% .PMC.CT$algorithm.gbd){
        PARAM.org <- set.global(K = K)
      } else if(algorithm[1] %in% .PMC.CT$algorithm.dmat){
        PARAM.org <- set.global.dmat(K = K)
      } else{
        comm.stop("The algorithm is not found.")
      }

      ### Replacing PARAM.org by previous PARAM.new.
      PARAM.org$ETA <- PARAM.new$ETA[-i.k] / sum(PARAM.org$ETA[-i.k])
      PARAM.org$log.ETA <- log(PARAM.org$ETA)
      PARAM.org$MU <- matrix(PARAM.new$MU[, -i.k], ncol = K)
      PARAM.org$SIGMA <- PARAM.new$SIGMA[-i.k]

      # Update steps.
      method.step <- switch(algorithm[1],
                            "em" = em.step,
                            "aecm" = aecm.step,
                            "apecm" = apecm.step,
                            "apecma" = apecma.step,
                            "kmeans" = kmeans.step,
                            NULL)
      PARAM.new <- method.step(PARAM.org)

      # Obtain classifications.
      if(algorithm[1] == "kmeans"){
        kmeans.update.class()
      } else{
        em.update.class()
      }

      # Get class numbers.
      N.CLASS <- get.N.CLASS(K)

      # For return.
      ret <- list(algorithm = algorithm[1],
                  param = PARAM.new,
                  class = .pmclustEnv$CLASS.spmd,
                  n.class = N.CLASS,
                  check = .pmclustEnv$CHECK)
    } else{
      break
    }
  }

  ret
} # end of pmclust.reduceK().

