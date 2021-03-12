#' Objective function using sum of squares
#'
#' This function fits the data of a single experiment for given parameters (relying on `exper_pred()`, compares results with actual data, and returns sum of squared error.
#'
#' Note that at present this function does NOT treat initial conditions as a parameter to fit.
#'
#' @param parms parameters for model implemented between transfers (this is first because `optim()` demands this)
#' @param x0 vector of initial conditions at start of transfer period 1 (e.g. at transfer 0)
#' @param ode_fun The ode function for model implemented between transfers.
#' @param transf.num How many transfers to model (1 is basically no dealing with transfers). With Landis data, this will generally be 6.
#' @param dat.real actual data, in the form of a matrix, with 1st column being transfer number, second being abund species 1, etc. Should not include initial transfer unless we're fitting initial conditions as parameters
#' `dataprep_landis()` preps the Landis data for this, and serves as an example.
#' @param transf.dur How many hours are each transfer period? Default is 48 (e.g. Landis data).
#' @param transf.dil At the beginning of each transfer period, what is the fraction of innoculant in the total material. Default is 5/200 (from Landis: 5 microliters innoculant into 195 new flour goo).
#' @param aug.ls Augmentation list. Not currently used, but a way to pass information down the chain of functions.
#'
#' @return Returns single value of the sum of squared error.
#'
#' @export

obj_ss = function(parms, #parameters for model implemented between transfers (first because optim)
                  x0, #initial conditions at start of transfer 1
                  ode_fun, #function for model implemented between transfers
                  transf.num, #how many transfers to model (1 is basically no dealing with transfers)
                  dat.real, #actual data, in the form of a *matrix*,
                  #with 1st column being transfer number, second being abund species 1, etc.
                  #should not include initial transfer unless we're fitting initial conditions as parameters
                  transf.dur=48, #how long are transfers? in hours
                  transf.dil = 5/200, #what is the dilution of each transfer. In Landis it's 5 microliter into 195 of new material, so 5/200
                  aug.ls=list() #misc list for anything extra, if needed
){
  reso = 1/transf.dur  #resolution time points to plot per hour.Only need one value now
  return.all = FALSE #  we don't need to get the full time series, just the transfer values.
  out=exper_pred(parms=parms,
                 x0=x0,
                 ode_fun = ode_fun,
                 transf.num=transf.num,
                 transf.dur=transf.dur,
                 transf.dil=transf.dil,
                 aug.ls = aug.ls,
                 reso=reso,
                 return.all=return.all)
  SS=0
  for(i.compare in 1:nrow(dat.real)){
    cur.transf=dat.real[i.compare,1]
    err = (dat.real[i.compare,-1]+1) -
      (out$transf.pred[out$transf.pred[,1]==cur.transf,-1]+1)
    err=err/10
    SS=SS+sum(err^2)
  }
  return(SS)
}
