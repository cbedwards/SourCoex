#' Plot fitted models and Landis data
#'
#' This function takes (presumably fitted) parameter values, model and experiment information, and raw data, and provides a panel each for experiments of interspecies competition, each of the single-species competition experiments, and parameter values.
#'
#'
#' @param parms vector of parameter values for model implemented between transfers (this is first because `optim()` demands this).
#' @param ode_fun function for the ode model model implemented between transfers.
#' @param parmnames vectofr of characters for the name of the parameters of the ode function. These are generally saved in the ode function files as parnames.XX, e.g. `parnames.LV`, `parnames.log`.
#' @param parmunits Vector of characters for the units of the parameters of the ode function. These are also generally saved in the ode function files, in the form of `units.LV`, `units.log`, etc).
#' @param dat.comp dataframe of the data from the interspecific experiments. Generated in the right format with `dataprep_landis()`.
#' @param dat.solo1 dataframe of the data from the single-species experiment for species 1. Generated in the right format with `dataprep_landis()`.
#' @param dat.solo2 dataframe of the data from the single-species experiment for species 2. Generated in the right format with `dataprep_landis()`.
#' @param specid vector of species names, for generating appropriate titles
#' @param specmap.cur map of species id to species name. For the landis data in `sourData`, the appropriate map is `specMap`.
#' @param reso How many time points per hour should we generate to plot trajectories
#'
#' @return Returns a ggplot2 object: 2x2 grid of plots and parameter table. Actually a 3x2 grid with row for custom title object.
#'
#' @export


plotter_landis=function(parms, # fitted parameters
                         ode_fun, # ode function
                         parmnames, #name of ode function parameters (e.g. parnames.lv)
                         parmunits, #units of ode function parms (e.g. units.LV)
                         dat.comp, #actual data from competition experiments
                         dat.solo1, #actual data from solo competition experiment, spec 1
                         dat.solo2, #actual data from solo competition experiment, spec 2
                         specid, # vector of species names, for mapping
                         aug.ls=list(),
                         specmap.cur=specMap, # map of species id to species name
                         reso=10 #plot points per hour.


){
  name.spec1=specMap[specMap$name.data==specid[1],2][[1]]
  name.spec2=specMap[specMap$name.data==specid[2],2][[1]]
  x0=c(5/1000,5/1000)
  if(!is.null(aug.ls$Til)){x0=c(x0, parms[length(parms)])}
  fit.pred=exper_pred(parms=parms,
                      x0=x0,
                      ode_fun = ode_fun,
                      transf.num=6,
                      return.all=T,
                      aug.ls = aug.ls,
                      reso=reso)
  dat.plot.joint=as.data.frame(fit.pred$series.tot)
  if(!is.null(aug.ls$Til)){
    dat.plot.joint=dat.plot.joint[,-ncol(dat.plot.joint)]
    }
  names(dat.plot.joint)=c("time","spec1","spec2")
  x0=c(10/1000,0)
  if(!is.null(aug.ls$Til)){x0=c(x0, parm[length(parms)])}
  fit.pred=exper_pred(parms=parms,
                      x0=x0,
                      ode_fun = ode_fun,
                      transf.num=6,
                      return.all=T,
                      aug.ls = aug.ls,
                      reso=reso)
  dat.plot.1=as.data.frame(fit.pred$series.tot)
  if(!is.null(aug.ls$Til)){
    dat.plot.1=dat.plot.1[,-ncol(dat.plot.1)]
  }
  names(dat.plot.1)=c("time","spec1","spec2")

  x0=c(0,10/1000)
  if(!is.null(aug.ls$Til)){x0=c(x0, parms[length(parms)])}
  fit.pred=exper_pred(parms=parms,
                      x0=x0,
                      ode_fun = ode_fun,
                      transf.num=6,
                      return.all=T,
                      aug.ls = aug.ls,
                      reso=reso)
  dat.plot.2=as.data.frame(fit.pred$series.tot)
  if(!is.null(aug.ls$Til)){
    dat.plot.2=dat.plot.2[,-ncol(dat.plot.2)]
  }
  names(dat.plot.2)=c("time","spec1","spec2")




  dat.comp.mean = dat.comp %>%
    group_by(transf) %>%
    summarise(abund1 = mean(abund1, na.rm=T),
              abund2 = mean(abund2, na.rm=T))
  dat.solo1.mean = dat.solo1 %>%
    group_by(transf) %>%
    summarise(abund1 = mean(abund1, na.rm=T),
              abund2 = mean(abund2, na.rm=T))
  dat.solo2.mean = dat.solo2 %>%
    group_by(transf) %>%
    summarise(abund1 = mean(abund1, na.rm=T),
              abund2 = mean(abund2, na.rm=T))

  gg.joint = ggplot(data=dat.plot.joint) +
    # fitted curves
    geom_path(aes(x=time, y=spec1), linetype="dashed")+
    geom_path(aes(x=time, y=spec2), col="blue", linetype="dashed")+
    # actual data
    geom_point(data=dat.comp, aes(x=transf*48, y = abund1))+
    geom_point(data=dat.comp, aes(x=transf*48, y = abund2), col='blue')+
    # data means
    geom_point(data=dat.comp.mean, aes(x=transf*48, y = abund1),
               col="black", alpha=.5, size=4, shape=17)+
    geom_point(data=dat.comp.mean, aes(x=transf*48, y = abund2),
               col='blue', alpha=.5, size=4, shape=17)+
    # plotting details
    xlab("time (hrs)")+
    ylab("CFUs")+
    ggtitle("Both species")
  gg.spec1 = ggplot(data=dat.plot.1) +
    #fitted curve
    geom_path(aes(x=time, y=spec1), linetype="dashed")+
    geom_path(aes(x=time, y=spec2), col="blue", linetype="dashed")+
    #actual data
    geom_point(data=dat.solo1, aes(x=transf*48, y = abund1))+
    geom_point(data=dat.solo1.mean, aes(x=transf*48, y = abund1),
               col="black", alpha=.5, size=4, shape=17)+
    xlab("time (hrs)")+
    ylab("CFUs")+
    ggtitle(paste0(name.spec1," only"))
  gg.spec2 = ggplot(data=dat.plot.2) +
    # fitted curve
    geom_path(aes(x=time, y=spec1), linetype="dashed")+
    geom_path(aes(x=time, y=spec2), col="blue", linetype="dashed")+
    # actual data
    geom_point(data=dat.solo2, aes(x=transf*48, y = abund1), col='blue')+
    geom_point(data=dat.solo2.mean, aes(x=transf*48, y = abund1),
               col='blue', alpha=.5, size=4, shape=17)+
    xlab("time (hrs)")+
    ylab("CFUs")+
    ggtitle(paste0(name.spec2," only"))



  title.gg <- ggplot() +
    labs(title = paste0(name.spec1, " vs ",
                        name.spec2,
                        ""),
         subtitle = "Based on Landis et al. 2021")+
    theme_void()+
    theme(plot.title=element_text(size=22))

  par.df=data.frame(parameter=parmnames,
                    fitted.value=format(round(parms,4),
                                        big.mark=",",
                                        decimal.mark=".",
                                        scientific = FALSE),
                    units=parmunits)
  rownames(par.df)=NULL

  tt.lv =  ttheme_minimal(
    col.just="right",
    core=list(bg_params = list(fill = blues9[c(1,1,2,2,3,3)], col=NA),
              fg_params=list(fontface=3, cex=2, parse=TRUE)),
    colhead=list(fg_params=list(col="navyblue", fontface=4L, cex=2)),
    rowhead=list(fg_params=list(col="navyblue", fontface=3L, cex=2)))

  gg.par = tableGrob(par.df, theme=tt.lv, rows=NULL)

  return(plot_grid(title.gg, NULL,
                   gg.spec1, gg.spec2,
                   gg.joint, gg.par, nrow=3, rel_heights = c(.15, 1, 1))
  )
}
