#' Prep Landis et al. 2021 data (ie `sour.data`) for model-fitting
#'
#' This function subsets the Landis data, and returns a matrix of appropriate initial conditions, as well as various formattings of the raw data for analysis and plotting. Plays well with `obj_helper()`, `plotter_landis()`.
#'
#' @param specid Vector of the first and second species for the analysis. Order matters.
#' @param data.use data frame of data to use. Almost assuredly sour.data
#'
#' @return Returns list with several components, suitable for analysis and plotting
#' \itemize{
#'   \item dat.comp - data of only the 2-species competition experiments. For `plotter_landis()`, joint plot.
#'   \item dat.solo1 - data of only the single-species competition experiments for species 1. For `plotter_landis()`, species 1 plot.
#'   \item dat.solo2 - data of only the single-species competition experiments for species 2. For `plotter_landis()`, species 2 plot.
#'   \item dat.real.ls - combination of dat.comp, dat.solo1, dat.solo2. Used to evaluate model fit with `obj_helper()`.
#'   \item x0.mat - matrix of initial conditions, with each row corresponding to one experiment. Same order as dat.real.ls, also used in `obj_helper()`.
#' @export


#function to prep fitted landis data
dataprep_landis = function(specid,
                           data.use=sour.data){
  dat.comp = data.use %>%
    filter(abund=="abs", spec1 == specid[1], spec2 == specid[2],transf>0)
  dat.comp=(dat.comp[,c("rep","transf", "abund1", "abund2")])
  temp=split(dat.comp, f = dat.comp[,"rep"])
  dat.comp.ls=lapply(temp,
                     function(x)(return(na.omit(as.matrix(x[,-1]))))
  )

  ## expers of only spec 1
  dat.solo1 = data.use %>%
    filter(abund=="abs", spec1 == specid[1], spec2 == specid[1],transf>0)
  dat.solo1=(dat.solo1[,c("rep","transf", "abund1", "abund2")])
  temp=split(dat.solo1, f = dat.solo1[,"rep"])
  dat.solo1.ls=lapply(temp,
                      function(x)(return(na.omit(as.matrix(x[,-1]))))
  )

  ## expers of only spec 2
  dat.solo2 = data.use %>%
    filter(abund=="abs", spec1 == specid[2], spec2 == specid[2],transf>0)
  #### BIG NOTE: WE HAVE TO SWITCH THE ORDER OF ABUND1 AND ABUND 2
  ## This experiment wasn't designed for a single pairing, so data was formatted such that abund1 is
  ## ## the abundance of the first species *in that competition*. So solos are ALWAYS in the abund1 slot.
  dat.solo2=(dat.solo2[,c("rep","transf", "abund2", "abund1")])
  temp=split(dat.solo2, f = dat.solo2[,"rep"])
  dat.solo2.ls=lapply(temp,
                      function(x)(return(na.omit(as.matrix(x[,-1]))))
  )

  dat.real.ls=c(dat.comp.ls, dat.solo1.ls, dat.solo2.ls)


  #make initial conditions matrix
  x0.1=matrix(c(5,5),nrow=1) #initial densities
  x0.comp=x0.1[rep(1,length(dat.comp.ls)),]
  x0.1=matrix(c(10,0),nrow=1) #initial densities
  x0.solo1=x0.1[rep(1,length(dat.solo1.ls)),]
  x0.1=matrix(c(0,10),nrow=1) #initial densities
  x0.solo2=x0.1[rep(1,length(dat.solo2.ls)),]
  x0.mat=rbind(x0.comp, x0.solo1, x0.solo2)
  return(list(dat.comp=dat.comp,
              dat.solo1=dat.solo1,
              dat.solo2=dat.solo2,
              dat.real.ls=dat.real.ls,
              x0.mat=x0.mat))
}
