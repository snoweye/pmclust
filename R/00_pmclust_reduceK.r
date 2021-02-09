### For automatically reducing K methods.

### X should be in spmd or gbd and set at .pmclustEnv or so, as used
### in pmclust().
pmclust.reduceK <- function(K = 2, algorithm = .PMC.CT$algorithm){
  if(any(algorithm[1] %in% c("kmeans"))){
    stop("kmeans/pkmeans is not supported in reduceK.")
  }

  if(algorithm[1] %in% .PMC.CT$algorithm.gbd){
    ret <- pmclust.reduceK.spmd(K = K, algorithm = algorithm)
  } else{
    comm.stop("The algorithm is not found.")
  }

  ret
} # End of pmclust.reduceK().


pmclust.reduceK.spmd <- function(K = 2, algorithm = .PMC.CT$algorithm){
  # Get an initial start.
  PARAM.org <- set.global(K = K)
  PARAM.org <- try(initial.em(PARAM.org), silent = TRUE)

  # Ensure the initial is good. Warning: This may take forever to run!
  repeat{
    if(class(PARAM.org) == "try-error"){
      PARAM.org <- set.global(K = K)
      PARAM.org <- try(initial.em(PARAM.org), silent = TRUE)
    } else{
      break
    }
  }

  # Update steps.
  method.step <- switch(algorithm[1],
                        "em" = em.step,
                        "aecm" = aecm.step,
                        "apecm" = apecm.step,
                        "apecma" = apecma.step,
                        NULL)
  if(comm.all(is.null(method.step))){
    comm.stop("Algorithm is not found.")
  }
  PARAM.new <- try(method.step(PARAM.org), silent = TRUE)
  em.update.class()
  N.CLASS <- get.N.CLASS(K)


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
      PARAM.org <- set.global(K = K)

      # Replacing PARAM.org by previous PARAM.new.
      PARAM.org$ETA <- PARAM.new$ETA[-i.k] / sum(PARAM.new$ETA[-i.k])
      PARAM.org$log.ETA <- log(PARAM.org$ETA)
      PARAM.org$MU <- matrix(PARAM.new$MU[, -i.k], ncol = K)
      PARAM.org$SIGMA <- PARAM.new$SIGMA[-i.k]

      # Update steps.
      e.step.spmd(PARAM.org)
      PARAM.new <- try(method.step(PARAM.org), silent = TRUE)
      em.update.class()
      N.CLASS <- get.N.CLASS(K)
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
} # End of pmclust.reduceK.spmd().

