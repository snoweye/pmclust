### This file gives initializations.

initial.em.dmat <- function(PARAM, MU = NULL){
  if(exists("X.dmat", envir = .pmclustEnv)){
    X.dmat <- get("X.dmat", envir = .pmclustEnv)
  }

  if(! is.ddmatrix(X.dmat)){
    stop("X.dmat is not a ddmatrix.")
  }

  if(is.null(MU)){
    N <- nrow(X.dmat)
    id <- spmd.bcast.integer(as.integer(sample(1:N, PARAM$K)))
    ### WCC: original
    PARAM$MU <- t(as.matrix(X.dmat[id, ]))
    ### WCC: debugging
    # tmp.1 <- X.dmat[id,]
    # tmp.2 <- as.matrix(tmp.1)
    # tmp.3 <- t(tmp.2)
    # PARAM$MU <- tmp.3
  } else{
    PARAM$MU <- MU
  }

  e.step.dmat(PARAM)
  PARAM <- em.onestep.dmat(PARAM)
  PARAM$logL <- logL.step.dmat()
  em.update.class.dmat()

  PARAM
} # End of initial.em.dmat().

initial.RndEM.dmat <- function(PARAM){
  logL.save <- -Inf
  i.iter <- 1

  PARAM.org <- PARAM
  repeat{
    PARAM <- try(initial.em.dmat(PARAM.org))
    if(class(PARAM) == "try-error"){
      comm.cat(PARAM, "\n", quiet = TRUE)
      next
    }

    N.CLASS <- get.N.CLASS.dmat(PARAM$K)
    if(any(N.CLASS < PARAM$min.N.CLASS)){
      if(.pmclustEnv$CONTROL$debug > 0){
        comm.cat("N.CLASS: ", N.CLASS, "\n", quiet = TRUE)
      }
      next
    }

    if(.pmclustEnv$CONTROL$debug > 0){
      comm.cat("Initial: ", format(Sys.time(), "%H:%M:%S"),
               ", iter: ", i.iter, ", logL: ",
                           sprintf("%-20.10f", PARAM$logL), "\n",
               sep = "", quiet = TRUE)
    }

    if(logL.save < PARAM$logL){
      logL.save <- PARAM$logL
      PARAM.save <- PARAM
      PARAM.save$initial.i.iter <- i.iter
    }

    i.iter <- i.iter + 1
    if(i.iter > .pmclustEnv$CONTROL$RndEM.iter){
      break
    }
  }

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat("Using initial iter: ", PARAM.save$initial.i.iter, "\n",
             sep = "", quiet = TRUE)
  }
  PARAM <- initial.em.dmat(PARAM.save, MU = PARAM.save$MU)
  PARAM
} # End of initial.RndEM.dmat().

