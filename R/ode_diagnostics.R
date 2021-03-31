#' Check ODE fit
#'
#' This function give useful diagnostics for how parameter values changed through fits, how objective function values changed, and if fits converged.
#'
#' @param fit1 first fitting using optim
#' @param fit2 second fitting using optim
#' @param Need to fill in params
#'
#' @return Nothing. Will print diagnostics to console.
#'
#' @export


ode_diagnostics = function(fit1,
                           fit2,
                           parms.guess,
                           x0.mat,
                           ode_fun,
                           transf.num=6,
                           dat.real.ls,
                           aug.ls=list()){
  guess.SS = obj_helper(parms=parms.guess,
                        x0.mat=x0.mat,
                        ode_fun = ode_fun,
                        transf.num=6,
                        dat.real.ls = dat.real.ls,
                        aug.ls=aug.ls)
  cat("Parameter changes:\n")
  print(rbind(parms.guess=parms.guess,
              fit1=fit1$par,
              fit2=fit2$par))
  cat("\nSum of squares improvement:\n")
  print(c(guess.SS=guess.SS, fit1.SS = fit1$value, fit2.SS=fit2$value))
  cat("\nconvergence:\n")
  print(c(fit1$convergence, fit2$convergence))
}
