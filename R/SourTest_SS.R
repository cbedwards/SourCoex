#' Test SS objective function
#'
#' This function tests that obj_helper.R is behaving. Takes identical arguments to obj_helper, but parameters should specific values such that microbial abundance is always 0.
#'
#'
#' @param parms0 parameters for model implemented between transfers, WITH VALUES CHOSEN TO ENSURE 0 ABUNDANCE AT ALL MEASUREMENT PERIODS
#' @param x0.mat matrix of initial conditions at start of transfer 1; each row is a different replicate. Note that this MUST be in the same order as dat.real.ls below
#' @param ode_fun The ode function for model implemented between transfers.
#' @param transf.num How many transfers to model (1 is basically no dealing with transfers). With Landis data, this will generally be 6.
#' @param dat.real.ls Empirical data to fit, in the form of a list. Each entry is a *matrix*, with columns of replicate number ('rep'), transfer number ('transf'), and abundance for each species ('abund1', 'abund2', etc)>
#' `dataprep_landis()` preps the Landis data for this, and serves as an example.
#' @param transf.dur How many hours are each transfer period? Default is 48 (e.g. Landis data).
#' @param transf.dil At the beginning of each transfer period, what is the fraction of innoculant in the total material. Default is 5/200 (from Landis: 5 microliters innoculant into 195 new flour goo).
#' @param aug.ls Augmentation list. Used to carry additional information, like identifying the Tilman model.
#'
#' @return Returns vector with the SS calculated by `obj_helper()` (first entry, labeled "objectiveSS"), and calculated by hand assuming real answer is at 0 (second entry, sanitySS). If these match, all is well.
#'
#' @export

SourTest_SS=function(parms.0,#parameter values that should lead to 0 microbes
                     x0.mat, #initial conditions
                     ode_fun, #ode function
                     transf.num, #how many transfers
                     dat.real.ls, #actual data
                     transf.dur=48,
                     transf.dil=5/200,
                     aug.ls=list()
){
  obj.ans = obj_helper(parms=parms.0,
                       x0.mat=x0.mat,
                       ode_fun = ode_fun,
                       transf.num=transf.num,
                       dat.real.ls = dat.real.ls,
                       aug.ls=aug.ls)
  dat.vec = do.call(rbind, dat.real.ls)[,-1]
  sanity.test = sum(dat.vec^2)
  return(c(objectiveSS=obj.ans, sanitySS=sanity.test))
}
