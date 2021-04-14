#' Simulate an experiment using ODE model
#'
#' This function takes parameters, model, initial conditions, and experiment design (number of transfers, dilution), and simulates trajectory of microbes.
#'
#'
#' @param parms parameters for model implemented between transfers (this is first because `optim()` demands this)
#' @param x0 vector of initial conditions at start of transfer period 1 (e.g. at transfer 0)
#' @param ode_fun The ode function for model implemented between transfers.
#' @param transf.num How many transfers to model (1 is basically no dealing with transfers). With Landis data, this will generally be 6.
#' @param transf.dur How many hours are each transfer period? Default is 48 (e.g. Landis data).
#' @param transf.dil At the beginning of each transfer period, what is the fraction of innoculant in the total material. Default is 5/200 (from Landis: 5 microliters innoculant into 195 new flour goo).
#' @param aug.ls Augmentation list. Currently used to identify modification for Tilman model (see `ode_Til()` and `ode_Til_1spec()`). Argument *$Til* should be a numeric for the number of species, and *$TilR$* is the amount of resouces at the start of each transfer (arbitrary, defaults to 195 for the microliters in Landis 2021).
#' @param reso How many time points per hour should we return for plotting?
#' @param return.all Should we return both "highly" resolved trajectories for plotting (resolution based on `reso`) and predictions at each transfer (`TRUE`), or should we just return the predictions at each transfer (`FALSE`). `FALSE` provides everything needed for model-fitting, and should be faster (although how much is unclear).
#'
#' @return For `return.all=TRUE`, returns a list, with `$transf.pred` being a data frame of predictions at each transfer point, and `$series.tot` being a data frame of predictions every 1/`reso` hours. For `return.all=FALSE`, only returns the data frame of `transf.pred`.
#'
#' @export


exper_pred_cont = function(parms, #parameters for model implemented between transfers (first because optim)
                           x0, #initial conditions at start of transfer 1
                           ode_fun, #function for model implemented between transfers
                           aug.ls=list(), #misc list for anything extra, if needed
                           times, #vector of times to predict at
                           reso = 10,  #resolution time points to plot per hour
                           return.all = TRUE # if TRUE, return list of both the reso-based timeseries and the time based timeseries
){
  times.plot=seq(0, max(times), by= 1/reso)
  if(!is.null(aug.ls$Til)){# carve off final parm for resources
    TilR = parms[length(parms)]
    parms=parms[-length(parms)]
  }
    out.lsoda = ode(x0,times,ode_fun,parms)
    out.plot = ode(x0,times.plot,ode_fun,parms)
    ## slightly hacky way of accomodating the tilman model here
  #make object to return
  res.ls=list(dat.pred=out.lsoda)
  if(return.all==T) {res.ls$series.tot=out.plot}
  return(res.ls)
}
