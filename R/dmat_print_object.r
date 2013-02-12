### This files provides functions for output.

mb.print.dmat <- function(PARAM, CHECK){
  if(CHECK$method %in%
     c("em.dmat", "aecm.dmat", "apecm.dmat", "apecma.dmat")){
    ETA <- PARAM$ETA
    SIGMA <- PARAM$SIGMA
  }
  MU <- PARAM$MU

  if(.pmclustEnv$COMM.RANK == 0){
#    cat("=====\n")
    cat("\n")
    cat("Method: ", CHECK$method, "\n", sep = "")
    cat("Convergence: ", CHECK$convergence,
        "  iter: ", CHECK$iter,
        "  abs.err: ", CHECK$abs.err,
        "  rel.err: ", CHECK$rel.err, "\n", sep = "")
    cat("logL: ", PARAM$logL, "\n", sep = "")
    cat("K: ", PARAM$K, "\n", sep = "")
    if(CHECK$method %in%
       c("em.dmat", "aecm.dmat", "apecm.dmat", "apecma.dmat")){
      cat("\nETA:\n")
      print(ETA)
    }
    cat("\nMU:\n")
    print(MU)
    if(CHECK$method %in%
       c("em.dmat", "aecm.dmat", "apecm.dmat", "apecma.dmat")){
      cat("\nSIGMA:\n")
      print(matrix(do.call("c", SIGMA), ncol = PARAM$K))
    }
#    cat("=====\n")
    cat("\n")
  }
  barrier()
} # End of mb.print.dmat().
