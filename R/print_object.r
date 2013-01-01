### This files provides functions for output.

mb.print <- function(PARAM, CHECK){
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
    if(CHECK$method %in% c("em", "aecm", "apecm", "apecma")){
      cat("\nETA:\n")
      print(PARAM$ETA)
    }
    cat("\nMU:\n")
    print(PARAM$MU)
    if(CHECK$method %in% c("em", "aecm", "apecm", "apecma")){
      cat("\nSIGMA:\n")
      print(matrix(do.call("c", PARAM$SIGMA), ncol = PARAM$K))
    }
#    cat("=====\n")
    cat("\n")
  }
  barrier()
} # End of mb.print().
