### For general methods.

pmclust.internal <- function(X, K,
    algorithm = c("em", "aecm", "apecm", "apecma", "kmeans"),
    MU = NULL, RndEM.iter = 10, CONTROL = NULL,
    method.own.X = c("spmdr", "common", "single"),
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm){
  # Check.
  if(!(algorithm[1] %in% c("em", "aecm", "apecm", "apecma", "kmeans"))){
    comm.stop("The algorithm is not found.")
  }
  if(!(method.own.X[1] %in% c("spmdr", "common", "single"))){
    comm.stop("The method.own.X is not found.")
  }

  # Assign X to .pmclustEnv
  convert.data(X, method.own.X[1], rank.own.X, comm)

  # Set global variables.
  PARAM.org <- set.global(K = K, RndEM.iter = RndEM.iter)
  if(!is.null(CONTROL)){
    tmp <- .pmclustEnv$CONTROL[!(names(.pmclustEnv$CONTROL) %in%
                                 names(CONTROL))]
    .pmclustEnv$CONTROL <- c(tmp, CONTROL)
  }

  # Initialization for algorithms.
  if(! is.null(MU)){
    if(algorithm[1] != "kmeans"){
      PARAM.org <- initial.em(PARAM.org, MU = MU)
    } else{
      PARAM.org <- initial.center(PARAM.org, MU = MU)
    }
  } else{
    if(algorithm[1] != "kmeans"){
      PARAM.org <- initial.RndEM(PARAM.org)
    } else{
      PARAM.org <- initial.center(PARAM.org)
    }
  }

  # Update steps.
  method.step <- switch(algorithm[1],
                        "em" = em.step,
                        "aecm" = aecm.step,
                        "apecm" = apecm.step,
                        "apecma" = apecma.step,
                        "kmeans" = kmeans.step,
                        NULL)
  if(is.null(method.step)){
    comm.stop("Algorithm is not found.")
  }
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
  ret
} # End of pmclust.internal().


pmclust <- function(X, K,
    algorithm = c("em", "aecm", "apecm", "apecma"),
    MU = NULL, RndEM.iter = 10, CONTROL = list(debug = 0),
    method.own.X = c("spmdr", "common", "single"),
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm){
  ret <- pmclust.internal(X, K,
                          algorithm = algorithm[1],
                          MU = MU,
                          RndEM.iter = RndEM.iter,
                          CONTROL = CONTROL,
                          method.own.X = method.own.X[1],
                          rank.own.X = rank.own.X,
                          comm = comm)
  class(ret) <- "pmclust"
  ret
} # end of pmclust().


pkmeans <- function(X, K, MU = NULL, CONTROL = list(debug = 0),
    method.own.X = c("spmdr", "common", "single"),
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm){
  algorithm <- "kmeans"
  ret <- pmclust.internal(X, K,
                          algorithm = algorithm[1],
                          MU = MU,
                          CONTROL = CONTROL,
                          method.own.X = method.own.X[1],
                          rank.own.X = rank.own.X,
                          comm = comm)
  class(ret) <- "pkmeans"
  ret
} # end of pkmeans().

