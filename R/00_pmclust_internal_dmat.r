### For general internal methods.

pmclust.internal.dmat <- function(X, K, MU = NULL,
    algorithm = .PMC.CT$algorithm, RndEM.iter = .PMC.CT$RndEM.iter,
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm){
  # Check.
  if(! is.ddmatrix(X)){
    comm.stop("A ddmatrix is required.")
  }
  if(!(algorithm[1] %in% .PMC.CT$algorithm)){
    comm.stop("The algorithm is not found.")
  }

  # Set global variables.
  PARAM.org <- set.global.dmat(K = K, RndEM.iter = RndEM.iter)
  if(!is.null(CONTROL)){
    tmp <- .pmclustEnv$CONTROL[!(names(.pmclustEnv$CONTROL) %in%
                                 names(CONTROL))]
    .pmclustEnv$CONTROL <- c(tmp, CONTROL)
  }

  # Initialization for algorithms.
  if(! is.null(MU)){
    if(algorithm[1] != "kmeans.dmat"){
      PARAM.org <- initial.em.dmat(PARAM.org, MU = MU)
    } else{
      PARAM.org <- initial.center.dmat(PARAM.org, MU = MU)
    }
  } else{
    if(algorithm[1] != "kmeans.dmat"){
      PARAM.org <- initial.RndEM.dmat(PARAM.org)
    } else{
      PARAM.org <- initial.center.dmat(PARAM.org)
    }
  }

  # Update steps.
  method.step <- switch(algorithm[1],
                        "em.dmat" = em.step.dmat,
                        "aecm.dmat" = aecm.step.dmat,
                        "apecm.dmat" = apecm.step.dmat,
                        "apecma.dmat" = apecma.step.dmat,
                        "kmeans.dmat" = kmeans.step.dmat,
                        NULL)
  if(is.null(method.step)){
    comm.stop("Algorithm is not found.")
  }
  PARAM.new <- method.step(PARAM.org)

  # Obtain classifications.
  if(algorithm[1] == "kmeans"){
    kmeans.update.class.dmat()
  } else{
    em.update.class.dmat()
  }

  # Get class numbers.
  N.CLASS <- get.N.CLASS.dmat(K)

  # For return.
  ret <- list(algorithm = algorithm[1],
              param = PARAM.new,
              class = .pmclustEnv$CLASS.spmd,
              n.class = N.CLASS,
              check = .pmclustEnv$CHECK)
  ret
} # End of pmclust.internal().

