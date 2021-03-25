#' Prep Landis et al. 2021 data (ie `sourData`) for model-fitting
#'
#' This function subsets the Landis data for single-species model-fitting, and returns a matrix of appropriate initial conditions, as well as various formattings of the raw data for analysis and plotting. Plays well with `obj_helper()`, `plotter_landis_solo()`. Alternative to `dataprep_landis()` for fitting logistic growth, etc.
#'
#' @param specid Single character of the species to prep.
#' @param data.use Data frame of data to use. Almost assuredly sourData
#' @param aug.ls() Augmentation list. Currently used for setting initial conditions for Tilman experiments
#'
#' @return Returns list with several components, suitable for analysis and plotting
#' \itemize{
#'   \item dat - dataframe of the data for this species, for plotting.
#'   \item dat.real.ls - list form of the data for this species, used to evaluate model fit with `obj_helper()`.
#'   \item x0.mat - matrix of initial conditions, with each row corresponding to one experiment. Same order as dat.real.ls, also used in `obj_helper()`.
#'  }
#' @export

#function to prep fitted landis data for within-species
dataprep_landis_solo = function(specid,#specid should be a single value
                                data.use=sourData,
                                aug.ls=list()){
  if(length(specid)>1){stop("More than one species given. Probably you want a different function?")}
  dat = data.use %>%
    filter(abund=="abs", spec1 == specid, spec2 == specid, transf>0)
  dat=(dat[,c("rep","transf", "abund1")])
  dat$abund1=dat$abund1/1000
  temp=split(dat, f = dat[,"rep"])
  dat.real.ls=lapply(temp,
                     function(x)(return(na.omit(as.matrix(x[,-1]))))
  )

  x0.mat=matrix(c(10),ncol=1, nrow=length(dat.real.ls)) #initial densities
  x0.mat=x0.mat/1000 #scaling to per k.

  return(list(dat=dat,
              dat.real.ls=dat.real.ls,
              x0.mat=x0.mat))
}

