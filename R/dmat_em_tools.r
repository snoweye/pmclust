### This file contains files for estimating parameters emperically.

### This function collects N.CLASS
get.N.CLASS.dmat <- function(K){
  tabulate(as.vector(.pmclustEnv$CLASS.dmat), nbins = K)
} # End of get.N.CLASS.dmat().


get.CLASS <- function(PARAM){
  A <- exists("CLASS.dmat", envir = .pmclustEnv)
  B <- exists("CLASS.spmd", envir = .pmclustEnv)

  if(A & B){
    comm.stop("CLASS.spmd and CLASS.dmat both exist in .pmclustEnv")
  } else{
    if(A){
      ret <- spmd.allgather.integer(as.integer(.pmclustEnv$CLASS.spmd),
                                    integer(PARAM$N))
      ret <- unlist(ret)
    }
    if(B){
      ret <- as.integer(as.vector(.pmclustEnv$CLASS.dmat))
    }
  }

  ret
} # End of get.CLASS().

