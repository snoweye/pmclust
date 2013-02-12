# This file gives a simple initialization.

initial.center.dmat <- function(PARAM, MU = NULL){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)
  if(! is.ddmatrix(X.dmat)){
    stop("X.dmat is not a ddmatrix.")
  }

  if(is.null(MU)){
    N <- nrow(X.dmat)
    id <- spmd.bcast.integer(as.integer(sample(1:N, PARAM$K)))
    PARAM$MU <- t(as.matrix(X.dmat[id, ]))
  } else{
    PARAM$MU <- MU
  }

  for(i.k in 1:PARAM$K){
    # B <- sweep(X.dmat, 2, as.vector(PARAM$MU[, i.k]))			# bug
    # .pmclustEnv$Z.dmat[, i.k] <- -rowSums(B * B)			# bug
    B <- base.pdsweep(dx = X.dmat, vec = PARAM$MU[, i.k],
                      MARGIN = 2L, FUN = "-")
    .pmclustEnv$Z.dmat[, i.k] <- -rowSums(B * B)
  }

  .pmclustEnv$CLASS <- unlist(apply(.pmclustEnv$Z.dmat, 1, which.max))

  PARAM
} # End of initial.center.dmat().

