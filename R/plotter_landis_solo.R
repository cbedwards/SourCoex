#' Plot fitted models and Landis data (single-species)
#'
#' This function takes (presumably fitted) parameter values, model and experiment information, and raw data, and provides a panel for data + trajectory, and parameter values.
#'
#'
#' @param parms vector of parameter values for model implemented between transfers (this is first because `optim()` demands this).
#' @param ode_fun function for the ode model model implemented between transfers.
#' @param parmnames vectofr of characters for the name of the parameters of the ode function. These are generally saved in the ode function files as parnames.XX, e.g. `parnames.LV`, `parnames.log`.
#' @param parmunits Vector of characters for the units of the parameters of the ode function. These are also generally saved in the ode function files, in the form of `units.LV`, `units.log`, etc).
#' @param dat dataframe of the data from the experiments. Generated in the right format with `dataprep_landis.solo()`.
#' @param specid vector of species names, for generating appropriate titles
#' @param specmap.cur map of species id to species name. For the landis data in `sourData`, the appropriate map is `specMap`.
#' @param reso How many time points per hour should we generate to plot trajectories
#'
#' @return Returns a ggplot2 object: 2x2 grid of plots and parameter table. Actually a 3x2 grid with row for custom title object.
#'
#' @export

#function to plot fitted landis data
plotter_landis_solo=function(parms, # fitted parameters
                             ode_fun, # ode function
                             parmnames, #name of ode function parameters (e.g. parnames.lv)
                             parmunits, #units of ode function parms (e.g. units.LV
                             dat,
                             specid, # vector of species names, for mapping
                             aug.ls=list(), # for tilman etc
                             specmap.cur=specMap, # map of species id to species name
                             reso=10, #plot points per hour.
                             noparms=FALSE #if true, don't plot the parms panel (useful for feeding into shiny)
){
  specname=specMap[specMap$name.data==specid,2][[1]]
  x0=10/1000
  if(!is.null(aug.ls$Til)){ #doing tilman model
    x0=c(x0,parms[length(parms)])
  }
  fit.pred=exper_pred(parms=parms,
                      x0=x0,
                      ode_fun = ode_fun,
                      transf.num=6,
                      return.all=T,
                      reso=reso,
                      aug.ls=aug.ls)
  dat.plot=as.data.frame(fit.pred$series.tot)
  if(!is.null(aug.ls$Til)){ #doing tilman model
    names(dat.plot)=c("time","abund","R")
  }else{
    names(dat.plot)=c("time","abund")
  }

  dat.mean = dat %>%
    group_by(transf) %>%
    summarise(abund1 = mean(abund1, na.rm=T))

  gg.plot = ggplot(data=dat.plot) +
    # fitted curves
    geom_path(aes(x=time, y=abund), linetype="dashed")+
    # actual data
    geom_point(data=dat, aes(x=transf*48, y = abund1))+
    # data means
    geom_point(data=dat.mean, aes(x=transf*48, y = abund1),
               col="black", alpha=.5, size=4, shape=17)+
    # plotting details
    xlab("time (hrs)")+
    ylab("Thousands of CFUs")+
    ggtitle("")
  require(cowplot)
  require(gridExtra)
  title.gg <- ggplot() +
    labs(title = paste0(specname,
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
  if(noparms==FALSE){
    out=plot_grid(title.gg, NULL,
                  gg.plot, gg.par, nrow=2, rel_heights = c(.15, 2))
  }else{
    out=plot_grid(title.gg,
              gg.plot, nrow=2, rel_heights = c(.15, 2))
  }
    return(out)
}
