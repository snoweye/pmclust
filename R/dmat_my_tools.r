### This function initializes global variables.
set.global.dmat <- function(K = 2, PARAM = NULL,
    method = c("kmeans.dmat"),
    RndEM.iter = 10){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)
  if(! is.ddmatrix(X.dmat)){
    stop("X.dmat is not a ddmatrix.")
  }
  CTXT <- ICTXT(X.dmat)

  ### Get data information.
  N <- nrow(X.dmat)
  p <- ncol(X.dmat)

  ### Set parameters.
  if(is.null(PARAM)){
    PARAM <- list(N = N, p = p, K = K,
                  ETA = NULL, log.ETA = NULL, MU = NULL, SIGMA = NULL,
                  U = rep(list(), K),
                  U.check = rep(TRUE, K),
                  logL = NULL,
                  min.N.CLASS = min(c((p + 1) * p * 0.5 + 1, N / K * 0.2)))
    PARAM$ETA <- rep(1/K, K)
    PARAM$log.ETA <- rep(-log(K), K) 
    PARAM$MU <- matrix(0, p, K)
    PARAM$SIGMA <- rep(list(diag(1.0, p)), K)
  } else{
    PARAM$N <- N
    K <- PARAM$K
  }

  ### Set global storages.
  .pmclustEnv$CONTROL <- list(max.iter = 1000, abs.err = 1e-4, rel.err = 1e-6,
                              debug = 1, RndEM.iter = RndEM.iter,
                              exp.min = log(.Machine$double.xmin),
                              exp.max = log(.Machine$double.xmax),
                              U.max = 1e+4,
                              U.min = 1e-6)

  .pmclustEnv$COMM.SIZE <- spmd.comm.size()
  .pmclustEnv$COMM.RANK <- spmd.comm.rank()

  .pmclustEnv$p.times.logtwopi <- p * log(2 * pi)

  .pmclustEnv$Z.dmat <- ddmatrix(0, N, K)
  .pmclustEnv$Z.colSums <- colSums(.pmclustEnv$Z.dmat)

  .pmclustEnv$W.dmat <- ddmatrix(0, N, K)
  .pmclustEnv$W.dmat.rowSums <- rowSums(.pmclustEnv$W.dmat)

  .pmclustEnv$U.dmat <- ddmatrix(0, N, K)

  .pmclustEnv$CLASS.dmat <- ddmatrix(0, N, 1)

  .pmclustEnv$CHECK <- list(method = method[1], i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)

  ### For semi-supervised clustering.
#  assign.ss.spmd()

  for(i.k in 1:K){
    tmp.U <- decompsigma(PARAM$SIGMA[[i.k]])
    PARAM$U[[i.k]] <- tmp.U$value 
    PARAM$U.check[[i.k]] <- tmp.U$check 
  }

  PARAM
} # End of set.global.dmat().

