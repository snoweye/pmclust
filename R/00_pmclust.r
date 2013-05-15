### For general methods.

pmclust <- function(X, K, MU = NULL,
    algorithm = .PMC.CT$algorithm, RndEM.iter = .PMC.CT$RndEM.iter,
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm){

  if(algorithm[1] %in% c("em", "aecm", "apecm", "apecma", "kmeans")){
    ret <- pmclust.internal(X, K,
                            MU = MU,
                            algorithm = algorithm[1],
                            RndEM.iter = RndEM.iter,
                            CONTROL = CONTROL,
                            method.own.X = method.own.X[1],
                            rank.own.X = rank.own.X,
                            comm = comm)
  } else if(algorithm[1] %in% c("em.dmat", "aecm.dmat", "apecm.dmat",
                                "apecma.dmat", "kmeans.dmat")){
    ret <- pmclust.internal.dmat(X, K,
                                 MU = MU,
                                 algorithm = algorithm[1],
                                 RndEM.iter = RndEM.iter,
                                 CONTROL = CONTROL,
                                 method.own.X = method.own.X[1],
                                 rank.own.X = rank.own.X,
                                 comm = comm)
  } else{
    comm.stop("The algorithm is not found.")
  }

  if(algorithm[1] %in% c("kmeans", "kmeans.dmat")){
    class(ret) <- "pkmeans"
  } else{
    class(ret) <- "pmclust"
  }
  ret
} # end of pmclust().


pkmeans <- function(X, K, MU = NULL,
    algorithm = c("kmeans", "kmeans.dmat"),
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm){

  if(algorithm[1] == "kmeans"){
    ret <- pmclust.internal(X, K,
                            MU = MU,
                            algorithm = algorithm[1],
                            CONTROL = CONTROL,
                            method.own.X = method.own.X[1],
                            rank.own.X = rank.own.X,
                            comm = comm)
  } else if(algorithm[1] == "kmeans.dmat"){
    ret <- pmclust.internal.dmat(X, K,
                                 MU = MU,
                                 algorithm = algorithm[1],
                                 CONTROL = CONTROL,
                                 method.own.X = method.own.X[1],
                                 rank.own.X = rank.own.X,
                                 comm = comm)
  } else{
    comm.stop("The algorithm is not found.")
  }
  class(ret) <- "pkmeans"
  ret
} # end of pkmeans().

