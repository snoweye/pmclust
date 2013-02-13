### This file contains major functions for EM iterations.

### E-step.
ape.step.dmat <- function(PARAM){
  for(i.k in 1:PARAM$K){
    logdmvnorm.dmat(PARAM, i.k)
  }

  ape.update.expectation.dmat(PARAM)
} # End of ape.step.dmat().

ape.step.dmat.k <- function(PARAM, i.k, update.logL = TRUE){
  logdmvnorm.dmat(PARAM, i.k)
  ape.update.expectation.k.dmat(PARAM, i.k, update.logL)
} # End of ape.step.dmat.k().


### z_nk / sum_k z_n might have numerical problems if z_nk all underflowed.
ape.update.expectation.dmat <- function(PARAM, update.logL = TRUE){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)

  N <- nrow(X.dmat)
  K <- PARAM$K

  .pmclustEnv$W.dmat <- W.plus.y(.pmclustEnv$W.dmat, PARAM$log.ETA, N, K)
  .pmclustEnv$U.dmat <- exp(.pmclustEnv$W.dmat)
  .pmclustEnv$Z.dmat <- .pmclustEnv$U.dmat

  tmp.id <- rowSums(.pmclustEnv$W.dmat < .pmclustEnv$CONTROL$exp.min) == K |
            rowSums(.pmclustEnv$W.dmat > .pmclustEnv$CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.dmat <- .pmclustEnv$W.dmat[tmp.id,]

    if(tmp.flag == 1){
      tmp.scale <- max(tmp.dmat) - .pmclustEnv$CONTROL$exp.max / K
    } else{
      tmp.scale <- unlist(apply(tmp.dmat, 1, max)) -
                   .pmclustEnv$CONTROL$exp.max / K
    }
    .pmclustEnv$Z.dmat[tmp.id,] <- exp(tmp.dmat - tmp.scale)
  }

  .pmclustEnv$W.dmat.rowSums <- rowSums(.pmclustEnv$Z.dmat)
  .pmclustEnv$Z.dmat <- .pmclustEnv$Z.dmat / .pmclustEnv$W.dmat.rowSums

  .pmclustEnv$Z.colSums <- colSums(.pmclustEnv$Z.dmat)
} # End of ape.update.expectation.dmat().

ape.update.expectation.k.dmat <- function(PARAM, i.k, update.logL = TRUE){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)

  N <- nrow(X.dmat)
  K <- PARAM$K

  .pmclustEnv$W.dmat[, i.k] <- W.plus.y.k(.pmclustEnv$W.dmat, PARAM$log.ETA,
                                          N, K, i.k)
  .pmclustEnv$U.dmat[, i.k] <- exp(.pmclustEnv$W.dmat[, i.k])
  .pmclustEnv$Z.dmat <- .pmclustEnv$U.dmat

  tmp.id <- rowSums(.pmclustEnv$W.dmat < .pmclustEnv$CONTROL$exp.min) == K |
            rowSums(.pmclustEnv$W.dmat > .pmclustEnv$CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.dmat <- .pmclustEnv$W.dmat[tmp.id,]

    if(tmp.flag == 1){
      tmp.scale <- max(tmp.dmat) - .pmclustEnv$CONTROL$exp.max / K
    } else{
      tmp.scale <- unlist(apply(tmp.dmat, 1, max)) -
                   .pmclustEnv$CONTROL$exp.max / K
    }
    .pmclustEnv$Z.dmat[tmp.id,] <- exp(tmp.dmat - tmp.scale)
  }

  .pmclustEnv$W.dmat.rowSums <- rowSums(.pmclustEnv$Z.dmat)
  .pmclustEnv$Z.dmat <- .pmclustEnv$Z.dmat / .pmclustEnv$W.dmat.rowSums
  .pmclustEnv$Z.colSums <- colSums(.pmclustEnv$Z.dmat)

  if(update.logL){
    .pmclustEnv$W.dmat.rowSums <- log(.pmclustEnv$W.dmat.rowSums)
    if(tmp.flag){
      .pmclustEnv$W.dmat.rowSums[tmp.id] <- .pmclustEnv$W.dmat.rowSums[tmp.id] +
                                            tmp.scale
    }
  }
} # End of ape.update.expectation.k.dmat().


### APECM-step.
apecm.step.dmat <- function(PARAM.org){
  .pmclustEnv$CHECK <- list(method = "apecm", i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- -.Machine$double.xmax

  ### For debugging.
  if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
      .pmclustEnv$CONTROL$save.log){
    if(! exists("SAVE.iter", envir = .pmclustEnv)){
      .pmclustEnv$SAVE.param <- NULL
      .pmclustEnv$SAVE.iter <- NULL
      .pmclustEnv$CLASS.iter.org <- unlist(apply(.pmclustEnv$Z.dmat, 1,
                                                 which.max))
    }
  }

  repeat{
    ### For debugging.
    if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
        .pmclustEnv$CONTROL$save.log){
      time.start <- proc.time()
    }

    PARAM.new <- try(apecm.onestep.dmat(PARAM.org))
    if(class(PARAM.new) == "try-error"){
      comm.cat("Results of previous iterations are returned.\n", quiet = TRUE)
      .pmclustEnv$CHECK$convergence <- 99
      PARAM.new <- PARAM.org
      break
    }

    .pmclustEnv$CHECK <- check.em.convergence(PARAM.org, PARAM.new, i.iter)
    if(.pmclustEnv$CHECK$convergence > 0){
      break
    }

    ### For debugging.
    if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
        .pmclustEnv$CONTROL$save.log){
      tmp.time <- proc.time() - time.start

      .pmclustEnv$SAVE.param <- c(.pmclustEnv$SAVE.param, PARAM.new)
      CLASS.iter.new <- unlist(apply(.pmclustEnv$Z.dmat, 1, which.max))
      tmp <- as.double(sum(CLASS.iter.new != .pmclustEnv$CLASS.iter.org))
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
} # End of apecm.step.dmat().

apecm.onestep.dmat <- function(PARAM){
#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(filename = "apecm.Rprof", append = TRUE)
#  }

  ### Update ETA
  PARAM <- cm.step.dmat.ETA(PARAM)
  ape.step.dmat(PARAM)

  ### Update MU and SIGMA
  for(i.k in 1:PARAM$K){
    PARAM <- cm.step.dmat.MU.SIGMA.k(PARAM, i.k)
    ape.step.dmat.k(PARAM, i.k,
                    update.logL = ifelse(i.k == PARAM$K, TRUE, FALSE))
  }

#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step.dmat()

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat(">>apecm.onestep: ", format(Sys.time(), "%H:%M:%S"),
             ", iter: ", .pmclustEnv$CHECK$iter, ", logL: ",
                         sprintf("%-30.15f", PARAM$logL), "\n",
             sep = "", quiet = TRUE)
    if(.pmclustEnv$CONTROL$debug > 4){
      logL <- indep.logL.dmat(PARAM)
      comm.cat("  >>indep.logL: ", sprintf("%-30.15f", logL), "\n",
               sep = "", quiet = TRUE)
    }
    if(.pmclustEnv$CONTROL$debug > 20){
      mb.print(PARAM, .pmclustEnv$CHECK)
    }
  }

  PARAM
} # End of apecm.onestep.dmat().

