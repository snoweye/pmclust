### This file gives initializations.

initial.em.dmat <- function(PARAM, MU = NULL){
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

