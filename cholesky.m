function [] = cholesky()
  A = matriz();
  L=chol(A,"upper");
  U=chol(A,"lower");
A
L
U
