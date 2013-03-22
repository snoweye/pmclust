### This file contains major functions for EM iterations.

### CM-step.
cm.step.dmat.ETA <- function(PARAM){
  ### MLE For ETA
  PARAM$ETA <- .pmclustEnv$Z.colSums / sum(.pmclustEnv$Z.colSums)
  PARAM$log.ETA <- log(PARAM$ETA)
  PARAM
} # End of cm.step.dmat.ETA().

cm.step.dmat.MU <- function(PARAM){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)

  p <- PARAM$p

  for(i.k in 1:PARAM$K){
    ### MLE for MU
    ### bug
    #B <- base.pdsweep(dx = X.dmat, vec = .pmclustEnv$Z.dmat[, i.k],
    #                  MARGIN = 1L, FUN = "*")
    comm.stop("Not implemented yet.")
    B <- NULL
    PARAM$MU[, i.k] <- colSums(B) / .pmclustEnv$Z.colSums[i.k]
  }

  PARAM
} # End of cm.step.dmat.MU().

cm.step.dmat.SIGMA <- function(PARAM){
  X.dmat <- get("X.dmat", envir = .GlobalEnv)

  N <- nrow(X.dmat)

  p <- PARAM$p
  p.2 <- p * p
  for(i.k in 1:PARAM$K){
    if(PARAM$U.check[[i.k]]){
      ### bug
      #B <- base.pdsweep(dx = X.dmat, vec = PARAM$MU[, i.k],
      #                  MARGIN = 2L, FUN = "-") *
      #     sqrt(.pmclustEnv$Z.dmat[, i.k] / .pmclustEnv$Z.colSums[i.k])
      comm.stop("Not implemented yet.")
      B <- NULL
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
        comm.cat("  SIGMA[[", i.k, "]] is fixed.\n", sep = "")
      }
    }
  }

  PARAM
} # End of cm.step.dmat.SIGMA().


### AECM-step.
aecm.step.dmat <- function(PARAM.org){
  .pmclustEnv$CHECK <- list(method = "aecm", i.iter = 0, abs.err = Inf,
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

    PARAM.new <- try(aecm.onestep.dmat(PARAM.org))
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
} # End of aecm.step.dmat().

aecm.onestep.dmat <- function(PARAM){
#  if(COMM.RANK == 0){
#    Rprof(filename = "aecm.Rprof", append = TRUE)
#  }

  PARAM <- cm.step.dmat.ETA(PARAM)
  e.step.dmat(PARAM, update.logL = FALSE)

  PARAM <- cm.step.dmat.MU(PARAM)
  e.step.dmat(PARAM, update.logL = FALSE)

  PARAM <- cm.step.dmat.SIGMA(PARAM)
  e.step.dmat(PARAM, update.logL = TRUE)

#  if(COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step.dmat()

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat(">>aecm.onestep: ", format(Sys.time(), "%H:%M:%S"),
             ", iter: ", .pmclustEnv$CHECK$iter, ", logL: ",
                         sprintf("%-30.15f", PARAM$logL), "\n",
             sep = "", quiet = TRUE)
    if(.pmclustEnv$CONTROL$debug > 4){
      logL <- indep.logL.dmat(PARAM)
      comm.cat("  >>indep.logL: ", sprintf("%-30.15f", logL), "\n", sep = "")
    }
    if(.pmclustEnv$CONTROL$debug > 20){
      mb.print(PARAM, .pmclustEnv$CHECK)
    }
  }

  PARAM
} # End of aecm.onestep.dmat().

