#' Apply objective function to multiple experimental replicates
#'
#' This function applies `SourCoex::obj_ss()` to the supplied list of data sets and initial conditions (e.g. multiple replicates). As such it serves as an intermediate between `optim()` and `obj_ss()`.
#'
#' Note: it may be appropriate to update this with multiple objective function options in the future.
#'
#' Note: this is a function that should be updated when we switch over to using hours.
#'
#' Also note: this function is quite fragile wrt the format of the real data. **CBE should add some safety features and or more flexibility.**
#'
#' @param parms parameters for model implemented between transfers (this is first because `optim()` demands this)
#' @param x0.mat matrix of initial conditions at start of transfer 1; each row is a different replicate. Note that this MUST be in the same order as dat.real.ls below
#' @param ode_fun The ode function for model implemented between transfers.
#' @param transf.num How many transfers to model (1 is basically no dealing with transfers). With Landis data, this will generally be 6.
#' @param dat.real.ls Empirical data to fit, in the form of a list. Each entry is a *matrix*, with columns of replicate number ('rep'), transfer number ('transf'), and abundance for each species ('abund1', 'abund2', etc)>
#' `dataprep_landis()` preps the Landis data for this, and serves as an example.
#' @param scale.vec atomic or vector with scaling factor for errors. Generally this should be the abundance of each species in the 1-species trials. If atomic, doesn't have meaningful effect.
#' @param transf.dur How many hours are each transfer period? Default is 48 (e.g. Landis data).
#' @param transf.dil At the beginning of each transfer period, what is the fraction of innoculant in the total material. Default is 5/200 (from Landis: 5 microliters innoculant into 195 new flour goo).
#' @param aug.ls Augmentation list. Not currently used, but a way to pass information down the chain of functions.
#'
#' @return Returns single value of the sum of the values from objective functions for each replicate.
#'
#' @export


#function to simultaneously fit multiple replicates, potentially with different x0s
obj_helper=function(parms, #parameters for model implemented between transfers (first because optim)
                    x0.mat, #initial conditions at start of transfer 1; each row is a different replicate.
                    ## same order as dat.real.ls below
                    ode_fun, #function for model implemented between transfers
                    transf.num, #how many transfers to model (1 is basically no dealing with transfers)
                    dat.real.ls, #actual data, list of entries, each entry is a *matrix*
                    #with 1st column being transfer number, second being abund species 1, etc.
                    #should not include initial transfer unless we're fitting initial conditions as parameters
                    scale.vec=1, #atomic or vector of scaling factors for error associated with each species.
                    transf.dur=48, #how long are transfers? in hours
                    times=NULL, #vector of times to predict at
                    transf.dil = 5/200, #what is the dilution of each transfer. In Landis it's 5 microliter into 195 of new material, so 5/200
                    aug.ls=list() #misc list for anything extra, if needed
){
  if(nrow(x0.mat) != length(dat.real.ls)){stop("initial densities and real data have mismatched dimensions")}
  ss.tot=0
  for(i.rep in 1:nrow(x0.mat)){
    x0=x0.mat[i.rep,]
    if(!is.null(aug.ls$Til)){
     x0=c(x0,parms[length(parms)])
    }
    dat.real=dat.real.ls[[i.rep]]
    if(sum(is.na(dat.real)) != 0){stop(paste0("NA in `dat.real.ls` entry number ", i.rep,". Check to make sure nothing weird is happing, remove transfer numbers with NAs."))}
    ss.cur=obj_cv(parms=parms,
                  x0=x0,
                  ode_fun = ode_fun,
                  transf.num=transf.num,
                  dat.real=dat.real,
                  scale.vec,
                  transf.dur=transf.dur,
                  times=times,
                  transf.dil=transf.dil,
                  aug.ls=aug.ls)
    # print(ss.cur)
    ss.tot=ss.tot+ss.cur
  }
  return(ss.tot)
}
