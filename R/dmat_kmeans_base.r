### This file provides functions for kmeans.

kmeans.e.step.dmat <- function(PARAM){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)
  for(i.k in 1:PARAM$K){
    B <- sweep(X.dmat, 2, as.vector(PARAM$MU[, i.k]))
    .pmclustEnv$Z.dmat[, i.k] <- sqrt(rowSums(B * B))
  }
  invisible()
} # End of kmeans.e.step.dmat().

kmeans.m.step.dmat <- function(PARAM){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)
  for(i.k in 1:PARAM$K){
    tmp <- colMeans(X.dmat[.pmclustEnv$CLASS.dmat == i.k,])
    PARAM$MU[, i.k] <- as.vector(tmp)
  } 
  PARAM
} # End of kmeans.m.step.dmat().

kmeans.logL.step.dmat <- function(){
  tmp <- apply(.pmclustEnv$Z.dmat, 1, which.min)
  tmp.diff <- sum(.pmclustEnv$CLASS.dmat != tmp)
  .pmclustEnv$CLASS.dmat <- tmp
  as.integer(tmp.diff)
} # End of kmeans.logL.step.dmat().

kmeans.step.dmat <- function(PARAM.org){
  .pmclustEnv$CHECK <- list(method = "kmeans", i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- PARAM.org$N

  ### For debugging.
  if((!is.null(.pmclustEnv$CONTROL$save.log)) && .pmclustEnv$CONTROL$save.log){
    if(! exists("SAVE.iter", envir = .pmclustEnv)){
      .pmclustEnv$SAVE.param <- NULL
      .pmclustEnv$SAVE.iter <- NULL
      .pmclustEnv$CLASS.iter.org <- apply(.pmclustEnv$Z.dmat, 1, which.min)
    }
  }

  repeat{
    ### For debugging.
    if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
        .pmclustEnv$CONTROL$save.log){
      time.start <- proc.time()
    }

    PARAM.new <- kmeans.onestep.dmat(PARAM.org)

    .pmclustEnv$CHECK <- check.kmeans.convergence(PARAM.org, PARAM.new, i.iter)

    if(.pmclustEnv$CHECK$convergence > 0){
      break
    }

    ### For debugging.
    if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
        .pmclustEnv$CONTROL$save.log){
      tmp.time <- proc.time() - time.start

      .pmclustEnv$SAVE.param <- c(.pmclustEnv$SAVE.param, PARAM.new)
      CLASS.iter.new <- apply(.pmclustEnv$Z.dmat, 1, which.min)
      tmp <- sum(CLASS.iter.new != .pmclustEnv$CLASS.iter.org)
      tmp.all <- c(tmp / PARAM.new$N, PARAM.new$logL,
                   PARAM.new$logL - PARAM.org$logL,
                   (PARAM.new$logL - PARAM.org$logL) / PARAM.org$logL)
      .pmclustEnv$SAVE.iter <- rbind(.pmclustEnv$SAVE.iter,
                                     c(tmp, tmp.all, tmp.time))
      .pmclustEnv$CLASS.iter.org <- CLASS.iter.new
    }

    PARAM.org <- PARAM.new
    i.iter <- i.iter + 1
  }

  PARAM.new
} # End of kmeans.step.dmat().

kmeans.onestep.dmat <- function(PARAM){
#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(filename = "kmeans.Rprof", append = TRUE)
#  }

  PARAM <- kmeans.m.step.dmat(PARAM)
  kmeans.e.step.dmat(PARAM)

#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- kmeans.logL.step.dmat()

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat(">>kmeans.onestep: ", format(Sys.time(), "%H:%M:%S"),
             ", iter: ", .pmclustEnv$CHECK$iter, ", logL: ",
                         sprintf("%-30d", PARAM$logL), "\n",
             sep = "", quiet = TRUE)
    if(.pmclustEnv$CONTROL$debug > 10){
      mb.print(PARAM, .pmclustEnv$CHECK)
    }
  }

  PARAM
} # End of kmeans.onestep.dmat().


kmeans.update.class.dmat <- function(){
  .pmclustEnv$CLASS.dmat <- apply(.pmclustEnv$Z.dmat, 1, which.min)
  invisible()
} # End of kmeans.update.class.dmat().

