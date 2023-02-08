# Solves the linear knapsack-type problem with a generic
# LP solver for a given budget.
lp_solver = function(reward, cost, budget) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("package `lpSolve` required.")
  }
  K = ncol(reward)
  x.coeffs = c(t(reward))
  A.mat = matrix(0, nrow(reward), length(x.coeffs))
  start = 1
  for (row in 1:nrow(A.mat)) {
    end = start + K -1
    A.mat[row, start:end] = 1
    start = end + 1
  }
  c.coeff = c(t(cost))
  f.con = rbind(A.mat, c.coeff)
  f.dir = rep("<=", nrow(f.con))
  f.rhs = c(rep(1, nrow(A.mat)), budget)
  lp.res = lpSolve::lp("max", x.coeffs, f.con, f.dir, f.rhs)

  list(gain = sum(x.coeffs * lp.res$solution),
       spend = sum(c.coeff * lp.res$solution),
       alloc.mat = matrix(lp.res$solution, nrow = nrow(reward), byrow = TRUE))
}

convex_hull = function(reward, cost, ix.order = NULL) {
  if (is.null(ix.order)) {
    ix.order = matrix(order(col(t(cost)), t(cost)), ncol = ncol(cost), byrow = TRUE) -
      ncol(cost) * (row(cost) - 1)
  }
  ret = convex_hull_rcpp(reward, cost, ix.order - 1)

  lapply(ret[[1]], function(x) x + 1)
}
