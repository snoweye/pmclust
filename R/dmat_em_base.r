### This file contains major functions for EM iterations.

### E-step.
e.step.dmat <- function(PARAM, update.logL = TRUE){
  for(i.k in 1:PARAM$K){
    logdmvnorm.dmat(PARAM, i.k)
  }

  update.expectation.dmat(PARAM, update.logL = update.logL)
  invisible()
} # End of e.step.dmat().

### z_nk / sum_k z_n might have numerical problems if z_nk all underflowed.
update.expectation.dmat <- function(PARAM, update.logL = TRUE){
  N <- PARAM$N
  K <- PARAM$K

  .pmclustEnv$U.dmat <- sweep(.pmclustEnv$W.dmat, 2, as.vector(PARAM$log.ETA))
  .pmclustEnv$Z.dmat <- exp(.pmclustEnv$U.dmat)

  tmp.id <- rowSums(.pmclustEnv$U.dmat < .pmclustEnv$CONTROL$exp.min) == K |
            rowSums(.pmclustEnv$U.dmat > .pmclustEnv$CONTROL$exp.max) > 0

comm.print("0")
  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
comm.print("0.1", all.rank = TRUE)
comm.print(str(tmp.id), all.rank = TRUE)
flush(stdout)
Sys.sleep(1)
barrier()
comm.stop()

    tmp.dmat <- .pmclustEnv$U.dmat[tmp.id,]
flush(stdout)
Sys.sleep(1)
barrier()

comm.print("0.2")
    if(tmp.flag == 1){
comm.print("0.3")
      tmp.scale <- max(tmp.dmat) - .pmclustEnv$CONTROL$exp.max / K
    } else{
comm.print("0.4")
      tmp.scale <- apply(tmp.dmat, 1, max) - .pmclustEnv$CONTROL$exp.max / K
    }
comm.print("0.5")
    tmp.scale <- as.vector(tmp.scale)
comm.print("0.6")

    .pmclustEnv$Z.dmat[tmp.id,] <- exp(tmp.dmat - tmp.scale)
  }
comm.print("1")

  .pmclustEnv$W.rowSums <- as.vector(rowSums(.pmclustEnv$Z.dmat))

  .pmclustEnv$Z.dmat <- .pmclustEnv$Z.dmat / .pmclustEnv$W.rowSums
  .pmclustEnv$Z.colSums <- as.vector(colSums(.pmclustEnv$Z.dmat))

  if(update.logL){
    .pmclustEnv$W.rowSums <- log(.pmclustEnv$W.rowSums)
    if(tmp.flag){
      .pmclustEnv$W.rowSums[tmp.id] <-
        .pmclustEnv$W.rowSums[tmp.id] + tmp.scale
    }
  }
  invisible()
} # End of update.expectation.dmat().


### M-step.
m.step.dmat <- function(PARAM){
  if(exists("X.dmat", envir = .pmclustEnv)){
    X.dmat <- get("X.dmat", envir = .pmclustEnv)
  }

  ### MLE For ETA
  PARAM$ETA <- .pmclustEnv$Z.colSums / sum(.pmclustEnv$Z.colSums)
  PARAM$log.ETA <- log(PARAM$ETA)

  p <- PARAM$p
  p.2 <- p * p
  for(i.k in 1:PARAM$K){
    ### MLE for MU
    B <- colSums(X.dmat * as.vector(.pmclustEnv$Z.dmat[, i.k])) /
         .pmclustEnv$Z.colSums[i.k]
    PARAM$MU[, i.k] <- as.vector(B)

    ### MLE for SIGMA
    if(PARAM$U.check[[i.k]]){
      B <- sweep(X.dmat, 2, as.vector(PARAM$MU[, i.k])) *
           (as.vector(sqrt(.pmclustEnv$Z.dmat[, i.k]) /
                      .pmclustEnv$Z.colSums[i.k]))
      tmp.SIGMA <- as.matrix(crossprod(B))
      dim(tmp.SIGMA) <- c(p, p)

      tmp.U <- decompsigma(tmp.SIGMA)
      PARAM$U.check[[i.k]] <- tmp.U$check
      if(tmp.U$check){
        PARAM$U[[i.k]] <- tmp.U$value
        PARAM$SIGMA[[i.k]] <- tmp.SIGMA
      }
    } else{
      if(.pmclustEnv$CONTROL$debug > 2){
        comm.cat("  SIGMA[[", i.k, "]] is fixed.\n", sep = "", quiet = TRUE)
      }
    }
  }

  PARAM
} # End of m.step.dmat().


### log likelihood.
logL.step.dmat <- function(){
  tmp.logL <- sum(.pmclustEnv$W.rowSums)
  tmp.logL
} # End of logL.step.dmat().


### EM-step.
em.step.dmat <- function(PARAM.org){
  .pmclustEnv$CHECK <- list(algorithm = "em", i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- -.Machine$double.xmax

  ### For debugging.
  if((!is.null(.pmclustEnv$CONTROL$save.log)) && .pmclustEnv$CONTROL$save.log){
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

    PARAM.new <- try(em.onestep.dmat(PARAM.org))
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
} # End of em.step.dmat().

em.onestep.dmat <- function(PARAM){
#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(filename = "em.Rprof", append = TRUE)
#  }

  PARAM <- m.step.dmat(PARAM)
  e.step.dmat(PARAM)

#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step.dmat()

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat(">>em.onestep: ", format(Sys.time(), "%H:%M:%S"),
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
} # End of em.onestep.dmat().


### Obtain classifications.
em.update.class.dmat <- function(){
  .pmclustEnv$CLASS.dmat <- apply(.pmclustEnv$Z.dmat, 1, which.max)
  invisible()
} # End of em.update.class.dmat().

