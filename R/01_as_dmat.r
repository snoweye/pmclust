### Convert X.spmd to X.dmat

as.dmat <- function(X.spmd, bldim = .BLDIM, ICTXT = .ICTXT,
    comm = .SPMD.CT$comm){
  X.spmd <- load.balance(X.spmd, comm = comm)

  N.spmd <- nrow(X.spmd)
  p <- ncol(X.spmd)
  N <- spmd.allreduce.integer(N.spmd, integer(1), op = "sum")

  X.dmat <- ddmatrix(0, N, p, bldim = c(N.spmd, p), CTXT = 2)
  X.dmat@Data <- X.spmd
  X.dmat <- redistribute(X.dmat, bldim = bldim, ICTXT = ICTXT)

  X.dmat
} # End of as.dmat().
