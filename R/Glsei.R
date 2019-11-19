## ==============================================================
## The original version of the code is obtained from 
## R package "limSolve" ver 1.5.5.1 written by 
## Karline Soetaert [aut, cre], Karel Van den Meersche [aut], Dick van Oevelen [aut], LAPACK authors [cph].  
## The code is modified to handle numerical instability.           
## ==============================================================

Glsei <- function (A = NULL, B = NULL, E = NULL, F = NULL, G = NULL, H = NULL, 
                   Wx = NULL, Wa = NULL, type = 2, tol = sqrt(.Machine$double.eps), 
                   tolrank = NULL, fulloutput = FALSE, verbose = TRUE) 
{
  if (is.vector(E) & length(F) == 1) 
    E <- t(E)
  else if (!is.matrix(E) & !is.null(E)) 
    E <- as.matrix(E)
  if (is.vector(A) & length(B) == 1) 
    A <- t(A)
  else if (!is.matrix(A) & !is.null(A)) 
    A <- as.matrix(A)
  if (is.vector(G) & length(H) == 1) 
    G <- t(G)
  else if (!is.matrix(G) & !is.null(G)) 
    G <- as.matrix(G)
  if (!is.matrix(F) & !is.null(F)) 
    F <- as.matrix(F)
  if (!is.matrix(B) & !is.null(B)) 
    B <- as.matrix(B)
  if (!is.matrix(H) & !is.null(H)) 
    H <- as.matrix(H)
  if (is.null(A) && is.null(E)) {
    if (is.null(G)) 
      stop("cannot solve least squares problem - A, E AND G are NULL")
    A <- matrix(data = 0, nrow = 1, ncol = ncol(G))
    B <- 0
  }
  else if (is.null(A)) {
    A <- matrix(data = E[1, ], nrow = 1)
    B <- F[1]
  }
  Neq <- nrow(E)
  Napp <- nrow(A)
  Nx <- ncol(A)
  Nin <- nrow(G)
  if (is.null(Nx)) 
    Nx <- ncol(E)
  if (is.null(Nx)) 
    Nx <- ncol(G)
  if (is.null(Nin)) 
    Nin <- 1
  if (is.null(Neq)) {
    Neq <- 0
    if (verbose & type == 1) 
      warning("No equalities - setting type = 2")
    type = 2
    F <- NULL
  }
  else {
    if (ncol(E) != Nx) 
      stop("cannot solve least squares problem - A and E not compatible")
    if (length(F) != Neq) 
      stop("cannot solve least squares problem - E and F not compatible")
  }
  if (is.null(G)) 
    G <- matrix(data = 0, nrow = 1, ncol = Nx)
  if (is.null(H)) 
    H <- 0
  if (ncol(G) != Nx) 
    stop("cannot solve least squares problem - A and G not compatible")
  if (length(B) != Napp) 
    stop("cannot solve least squares problem - A and B not compatible")
  if (length(H) != Nin) 
    stop("cannot solve least squares problem - G and H not compatible")
  if (!is.null(Wa)) {
    if (length(Wa) != Napp) 
      stop("cannot solve least squares problem - Wa should have length = number of rows of A")
    A <- A * Wa
    B <- B * Wa
  }
  Tol <- tol
  if (is.null(Tol)) 
    Tol <- sqrt(.Machine$double.eps)
  IsError <- FALSE
  if (type == 2) {
    if (!is.null(Wx)) 
      stop("cannot solve least squares problem - weights not implemented for type 2")
    if (!is.null(Wa)) 
      stop("cannot solve least squares problem - weights not implemented for type 2")
    ## library("quadprog")
    dvec <- crossprod(A, B)
    Dmat <- crossprod(A, A)
    Pass <- 0
    tol.d <- tol
    while(Pass==0){
      diag(Dmat) <- diag(Dmat) + tol.d
      Amat <- t(rbind(E, G))
      bvec <- c(F, H)
      sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = Neq)
      sol$IsError <- FALSE
      sol$X <- sol$solution
      if(any(is.na(sol$X))==TRUE){
        Pass <- 0
        tol.d <- 5*tol.d
      }else{
        Pass <- 1
      }
    }
  }
  else stop("cannot solve least squares problem - type unknown")
  X <- sol$X
  X[which(abs(X) < Tol)] <- 0
  if (any(is.infinite(X))) {
    residual <- Inf
    solution <- Inf
  }
  else {
    residual <- 0
    if (Nin > 0) {
      ineq <- G %*% X - H
      residual <- residual - sum(ineq[ineq < 0])
    }
    if (Neq > 0) 
      residual <- residual + sum(abs(E %*% X - F))
    if (residual > Tol) 
      sol$IsError <- TRUE
    solution <- 0
    if (Napp > 0) 
      solution <- sum((A %*% X - B)^2)
  }
  xnames <- colnames(A)
  if (is.null(xnames)) 
    xnames <- colnames(E)
  if (is.null(xnames)) 
    xnames <- colnames(G)
  names(X) <- xnames
  res <- list(X = X, residualNorm = residual, solutionNorm = solution, 
              IsError = sol$IsError, type = "lsei")
  return(res)
}
