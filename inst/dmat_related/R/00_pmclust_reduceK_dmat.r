### A dmat version for automatically reducing K methods.

pmclust.reduceK.dmat <- function(K = 2, algorithm = .PMC.CT$algorithm){
  # Get an initial start.
  PARAM.org <- set.global.dmat(K = K)
  PARAM.org <- try(initial.em.dmat(PARAM.org), silent = TRUE)

  # Ensure the initial is good. Warning: This may take forever to run!
  repeat{
    if(class(PARAM.org) == "try-error"){
      PARAM.org <- set.global.dmat(K = K)
      PARAM.org <- try(initial.em.dmat(PARAM.org), silent = TRUE)
    } else{
      break
    }
  }

  # Update steps.
  method.step <- switch(algorithm[1],
                        "em.dmat" = em.step.dmat,
                        # "aecm.dmat" = aecm.step.dmat,
                        # "apecm.dmat" = apecm.step.dmat,
                        # "apecma.dmat" = apecma.step.dmat,
                        NULL)
  if(comm.all(is.null(method.step))){
    comm.stop("Algorithm is not found.")
  }
  PARAM.new <- try(method.step(PARAM.org), silent = TRUE)
  em.update.class.dmat()
  N.CLASS <- get.N.CLASS.dmat(K)


  # Reduce K if error occurs.
  repeat{
    if((class(PARAM.new) == "try-error" ||
        .pmclustEnv$CHECK$convergence == 99) &&
       K > 1){
      # Drop specific i.k if available or
      # drop the smallest class or
      # drop the class with the smallest eta among all small classes or
      # drop all classes with 0 elements.
      if(.pmclustEnv$CONTROL$stop.at.fail && .pmclustEnv$FAIL.i.k > 0){
        i.k <- .pmclustEnv$FAIL.i.k
      } else{
        i.k <- which(N.CLASS == min(N.CLASS))
      }
      if(i.k > 1 && min(N.CLASS) > 0){
        i.k <- i.k[which.min(PARAM.new$ETA[i.k])]
      }
      K <- K - length(i.k)
      comm.cat("- Reduce: ", K, "\n")

      # Initial global storage.
      PARAM.org <- set.global.dmat(K = K)

      # Replacing PARAM.org by previous PARAM.new.
      PARAM.org$ETA <- PARAM.new$ETA[-i.k] / sum(PARAM.new$ETA[-i.k])
      PARAM.org$log.ETA <- log(PARAM.org$ETA)
      PARAM.org$MU <- matrix(PARAM.new$MU[, -i.k], ncol = K)
      PARAM.org$SIGMA <- PARAM.new$SIGMA[-i.k]

      # Update steps.
      e.step.dmat(PARAM.org)
      PARAM.new <- try(method.step(PARAM.org), silent = TRUE)
      em.update.class.dmat()
      N.CLASS <- get.N.CLASS.dmat(K)
    } else{
      break
    }
  }

  # For return.
  ret <- list(algorithm = algorithm[1],
              param = PARAM.new,
              class = .pmclustEnv$CLASS.spmd,
              n.class = N.CLASS,
              check = .pmclustEnv$CHECK)

  ret
} # End of pmclust.reduceK.dmat().

