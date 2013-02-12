### This file contains functions for log density of MVN.
### These will majorly update .pmclustEnv$W.spmd.

logdmvnorm.dmat <- function(PARAM, i.k){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)

#  for(i.k in 1:PARAM$K){
#    U <- chol(PARAM$SIGMA[[i.k]])
    U <- PARAM$U[[i.k]]
    logdet <- sum(log(abs(diag(U)))) * 2
#    B <- t.X.spmd - PARAM$MU[, i.k]
#    A <- backsolve(U, B, upper.tri = TRUE, transpose = TRUE)
#    distval <- colSums(A * A)

    ### SPMD
    # B <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow(X.spmd), ncol(X.spmd))
    # B <- B %*% backsolve(U, diag(1, PARAM$p))
    # distval <- rowSums(B * B)
    # .pmclustEnv$W.spmd[, i.k] <- -(.pmclustEnv$p.times.logtwopi + logdet +
    #                                distval) * 0.5

    ### DMAT
    B <- base.pdsweep(dx = X.dmat, vec = PARAM$MU[, i.k],
                      MARGIN = 2L, FUN = "-")
    C <- backsolve(U, diag(1, PARAM$p))
    B <- B %*% as.ddmatrix(C, bldim = bldim(B), ICTXT = ctxt(B))
    distval <- rowSums(B * B)
    .pmclustEnv$W.dmat[, i.k] <- -(.pmclustEnv$p.times.logtwopi + logdet +
                                   distval) * 0.5
#  }
  invisible()
} # End of logdmvnorm.dmat().

