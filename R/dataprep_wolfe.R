#' Prep Wolfe pilot data
#'
#' This function subsets the Landis data, and returns a matrix of appropriate initial conditions, as well as various formattings of the raw data for analysis and plotting. Plays well with `obj_helper()`, `plotter_landis()`.
#'
#' @param specid Vector of the first and second species for the analysis. Order matters.
#' @param data.use data frame of data to use. Almost assuredly sourData
#'
#' @return Returns list with several components, suitable for analysis and plotting
#' \itemize{
#'   \item dat.comp - data of only the 2-species competition experiments. For `plotter_landis()`, joint plot.
#'   \item dat.solo1 - data of only the single-species competition experiments for species 1. For `plotter_landis()`, species 1 plot.
#'   \item dat.solo2 - data of only the single-species competition experiments for species 2. For `plotter_landis()`, species 2 plot.
#'   \item dat.real.ls - combination of dat.comp, dat.solo1, dat.solo2. Used to evaluate model fit with `obj_helper()`.
#'   \item x0.mat - matrix of initial conditions, with each row corresponding to one experiment. Same order as dat.real.ls, also used in `obj_helper()`.
#'   }
#' @export


#function to prep fitted landis data
dataprep_wolfe = function(specid,
                           data.use){
  dat.comp = data.use %>%
    filter(spec == specid)
  dat.comp=(dat.comp[,c("rep", "time","density.cfus")])
  dat.comp = dat.comp %>%
    rename(abund1 = density.cfus)
  dat.comp$abund1=dat.comp$abund1/1000
  # cut out initial conditions
  dat.start = dat.comp[dat.comp$time==0,]
  temp=split(dat.comp, f = dat.comp[,"rep"])
  dat.comp.ls=lapply(temp,
                     function(x)(return(na.omit(as.matrix(x[,-1]))))
  )
  #make initial conditions matrix
  x0.mat=matrix(dat.start$abund1,ncol=1) #initial densities
  x0.mat=x0.mat/1000
  times=sort(unique(dat.comp$time))
  return(list(dat.comp=dat.comp,
              dat.real.ls=dat.comp.ls,
              x0.mat=x0.mat,
              times=times))
}
