### Convert X.gbd to X.dmat

as.dmat <- function(X.gbd, bldim = .BLDIM, ICTXT = .ICTXT,
    comm = .SPMD.CT$comm){
  X.gbd <- load.balance(X.gbd, comm = comm)

  N.gbd <- nrow(X.gbd)
  p <- ncol(X.gbd)
  N <- spmd.allreduce.integer(N.gbd, integer(1), op = "sum")

  X.dmat <- ddmatrix(0, N, p, bldim = c(N.gbd, p), CTXT = 2)
  X.dmat@Data <- X.gbd
  X.dmat <- redistribute(X.dmat, bldim = bldim, ICTXT = ICTXT)

  X.dmat
} # End of as.dmat().
