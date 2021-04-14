#' Objective function using sum of squares
#'
#' This function fits the data of a single experiment for given parameters (relying on `exper_pred()`, compares results with actual data, and returns sum of squared CV(!!).
#'
#' Note that at present this function does NOT treat initial conditions as a parameter to fit.
#'
#' @param parms parameters for model implemented between transfers (this is first because `optim()` demands this)
#' @param x0 vector of initial conditions at start of transfer period 1 (e.g. at transfer 0)
#' @param ode_fun The ode function for model implemented between transfers.
#' @param transf.num How many transfers to model (1 is basically no dealing with transfers). With Landis data, this will generally be 6.
#' @param dat.real actual data, in the form of a matrix, with 1st column being transfer number, second being abund species 1, etc. Should not include initial transfer unless we're fitting initial conditions as parameters
#' `dataprep_landis()` preps the Landis data for this, and serves as an example.
#' @param scale.vec atomic or vector with scaling factor for errors. Generally this should be the abundance of each species in the 1-species trials. If atomic, doesn't have meaningful effect.
#' @param transf.dur How many hours are each transfer period? Default is 48 (e.g. Landis data).
#' @param transf.dil At the beginning of each transfer period, what is the fraction of innoculant in the total material. Default is 5/200 (from Landis: 5 microliters innoculant into 195 new flour goo).
#' @param aug.ls Augmentation list. Not currently used, but a way to pass information down the chain of functions.
#'
#' @return Returns single value of the sum of squared error.
#'
#' @export

obj_cv_cont = function(parms, #parameters for model implemented between transfers (first because optim)
                  x0, #initial conditions at start of transfer 1
                  ode_fun, #function for model implemented between transfers
                  transf.num, #how many transfers to model (1 is basically no dealing with transfers)
                  dat.real, #actual data, in the form of a *matrix*,
                  times,
                  #with 1st column being transfer number, second being abund species 1, etc.
                  #should not include initial transfer unless we're fitting initial conditions as parameters
                  scale.vec, # atomic or vector of mean population size within the single-species trials, for use to handle scaling issues.
                  aug.ls=list()
){
  times=dat.real[,1]
  out=exper_pred_cont(parms=parms,
                 x0=x0,
                 ode_fun = ode_fun,
                 aug.ls = aug.ls,
                 times=times,
                 return.all=FALSE)
  ## account for Tilman: if modeling tilman r* model, remove resource column from comparison(!!)
  if(!is.null(aug.ls$Til)){
    out$transf.pred=out$transf.pred[,-ncol(out$transf.pred)]
  }
  SS=0
  for(i.compare in 1:nrow(dat.real)){
    cur.transf=dat.real[i.compare,1]
    err = (dat.real[i.compare,-1]+1)/scale.vec -
      (out$dat.pred[out$dat.pred[,1]==cur.transf,-1]+1)/scale.vec
    SS=SS+sum(err^2)
  }
  return(SS)
}
