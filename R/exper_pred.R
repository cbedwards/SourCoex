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
#' @param aug.ls Augmentation list. Not currently used, but a way to pass information down the chain of functions.
#' @param reso How many time points per hour should we return for plotting?
#' @param return.all Should we return both "highly" resolved trajectories for plotting (resolution based on `reso`) and predictions at each transfer (`TRUE`), or should we just return the predictions at each transfer (`FALSE`). `FALSE` provides everything needed for model-fitting, and should be faster (although how much is unclear).
#'
#' @return For `return.all=TRUE`, returns a list, with `$transf.pred` being a data frame of predictions at each transfer point, and `$series.tot` being a data frame of predictions every 1/`reso` hours. For `return.all=FALSE`, only returns the data frame of `transf.pred`.
#'
#' @export


exper_pred = function(parms, #parameters for model implemented between transfers (first because optim)
                      x0, #initial conditions at start of transfer 1
                      ode_fun, #function for model implemented between transfers
                      transf.num, #how many transfers to model (1 is basically no dealing with transfers)
                      transf.dur=48, #how long are transfers? in hours
                      transf.dil = 5/200, #what is the dilution of each transfer. In Landis it's 5 microliter into 195 of new material, so 5/200
                      aug.ls=list(), #misc list for anything extra, if needed
                      reso = 10,  #resolution time points to plot per hour
                      return.all = TRUE # if TRUE, return list of both the timeseries and the transfer numbers
){
  series.tot=NULL #for storing total timeseries
  transf.pred=NULL
  times=seq(0, transf.dur, by= 1/reso)
  x0.cur = x0

  for(i.transf in 1:transf.num){
    out.lsoda=ode(x0.cur,times,ode_fun,parms)
    ## what would we count at the end of this cycle?
    transf.count=tail(out.lsoda[,-1],1)
    ## update times so that they make sense
    out.lsoda[,1]=out.lsoda[,1]+(i.transf-1)*transf.dur
    ## add our total timeseries and cycle counts
    series.tot=rbind(series.tot, out.lsoda)
    transf.pred=rbind(transf.pred, transf.count)
    ## Okay, implement the transfer
    x0.cur=transf.count*transf.dil
  }
  transf.pred = cbind(transf=1:transf.num, transf.pred)
  #make object to return
  res.ls=list(transf.pred=transf.pred)
  if(return.all==T) {res.ls$series.tot=series.tot}
  return(res.ls)
}
