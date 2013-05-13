### This files provides functions for output.

print.pmclust <- function(x, ...){
  CHECK <- x$check
  PARAM <- x$param

  ETA <- PARAM$ETA
  SIGMA <- matrix(do.call("c", PARAM$SIGMA), ncol = PARAM$K)
  MU <- PARAM$MU

  if(spmd.comm.rank() == 0){
    cat("\n")
    cat("Method: ", CHECK$method, "\n", sep = "")
    cat("Convergence: ", CHECK$convergence,
        "  iter: ", CHECK$iter,
        "  abs.err: ", CHECK$abs.err,
        "  rel.err: ", CHECK$rel.err, "\n", sep = "")
    cat("logL: ", PARAM$logL, "\n", sep = "")
    cat("K: ", PARAM$K, "\n", sep = "")
    cat("\nETA:\n")
    print(ETA)
    cat("\nMU:\n")
    print(MU)
    cat("\nSIGMA:\n")
    print(SIGMA)
    cat("\n")
  }
  barrier()
} # End of print.pmclust().


print.pkmeans <- function(x, ...){
  CHECK <- x$check
  PARAM <- x$param

  MU <- PARAM$MU

  if(spmd.comm.rank() == 0){
    cat("\n")
    cat("Method: ", CHECK$method, "\n", sep = "")
    cat("Convergence: ", CHECK$convergence,
        "  iter: ", CHECK$iter,
        "  abs.err: ", CHECK$abs.err,
        "  rel.err: ", CHECK$rel.err, "\n", sep = "")
    cat("logL: ", PARAM$logL, "\n", sep = "")
    cat("K: ", PARAM$K, "\n", sep = "")
    cat("\nMU:\n")
    print(MU)
    cat("\n")
  }
  barrier()
} # End of print.pkmeans().
