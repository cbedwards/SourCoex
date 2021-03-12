# Step 1: LV model as example ode_func

library(deSolve)
library(tidyverse)
library(SourCoex)
library(here)
library(cowplot)
library(gridExtra)
library(beepr)


ode_LV=function(t,y,parms) {
  # This function is based on the Lotka-Volterra competition model
  #state variables:
  x1=max(y[1],0); x2=max(y[2],0);
  # Parameters:
  r1=parms[1];
  r2=parms[2];
  a12 = parms[3]
  a21 = parms[4]
  k1=parms[5]
  k2=parms[6]
  # Model:
  dx1= r1*x1*(1 - (x1 + a12*x2)/k1)
  dx2= r2*x2*(1 - (x2 + a21*x1)/k2)
  dY=c(dx1,dx2);
  return(list(dY));
}
# for an ode function: make vector of parameter names in order
# and make vector of units in order, using plotmath notation.
parmnames.LV=c("r[1]", "r[2]", "alpha[12]","alpha[21]","k[1]","k[2]")
units.LV=c("growth per capita per hour",
           "growth per capita per hour",
           "unitless (multiplier)",
           "unitless (multiplier)",
           "individuals",
           "individuals"
)

ode_log=function(t,y,parms) {
  # Yes, this can be solved analytically, but I want to make sure the method is working.
  #state variables:
  x=max(y[1],0);
  # Parameters:
  r=parms[1];
  k=parms[2]
  # Model:
  dx= r*x*(1 - x/k)
  dY=c(dx);
  return(list(dY));
}

# for an ode function: make vector of parameter names in order
# and make vector of units in order, using plotmath notation.
parmnames.log=c("r","k")
units.log=c("growth per capita per hour",
            "individuals"
)


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

#function to prep fitted landis data for within-species
dataprep_landis.solo = function(specid,#specid should be a single value
                                data.use=sour.data){
  if(length(specid)>1){stop("More than one species given. Probably you want a different function?")}
  dat = data.use %>%
    filter(abund=="abs", spec1 == specid, spec2 == specid, transf>0)
  dat=(dat[,c("rep","transf", "abund1")])
  temp=split(dat, f = dat[,"rep"])
  dat.real.ls=lapply(temp,
                     function(x)(return(na.omit(as.matrix(x[,-1]))))
  )

  x0.mat=matrix(c(10),ncol=1, nrow=length(dat.real.ls)) #initial densities
  return(list(dat=dat,
              dat.real.ls=dat.real.ls,
              x0.mat=x0.mat))
}


#function to plot fitted landis data
plotter_landis=function(parms, # fitted parameters
                        ode_fun, # ode function
                        parmnames, #name of ode function parameters (e.g. parnames.lv)
                        parmunits, #units of ode function parms (e.g. units.LV
                        dat.comp, #actual data from competition experiments
                        dat.solo1, #actual data from solo competition experiment, spec 1
                        dat.solo2, #actual data from solo competition experiment, spec 2
                        specid, # vector of species names, for mapping
                        specmap.cur=spec.map, # map of species id to species name
                        reso=10 #plot points per hour.
){
  name.spec1=spec.map[spec.map$name.data==specid[1],2][[1]]
  name.spec2=spec.map[spec.map$name.data==specid[2],2][[1]]

  fit.pred=exper_pred(parms=parms,
                      x0=c(5,5),
                      ode_fun = ode_fun,
                      transf.num=6,
                      return.all=T,
                      reso=reso)
  dat.plot.joint=as.data.frame(fit.pred$series.tot)
  names(dat.plot.joint)=c("time","spec1","spec2")

  fit.pred=exper_pred(parms=parms,
                      x0=c(10,0),
                      ode_fun = ode_fun,
                      transf.num=6,
                      return.all=T,
                      reso=reso)
  dat.plot.1=as.data.frame(fit.pred$series.tot)
  names(dat.plot.1)=c("time","spec1","spec2")

  fit.pred=exper_pred(parms=parms,
                      x0=c(0,10),
                      ode_fun = ode_fun,
                      transf.num=6,
                      return.all=T,
                      reso=reso)
  dat.plot.2=as.data.frame(fit.pred$series.tot)
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
    #actual data
    geom_point(data=dat.solo1, aes(x=transf*48, y = abund1))+
    geom_point(data=dat.solo1.mean, aes(x=transf*48, y = abund1),
               col="black", alpha=.5, size=4, shape=17)+
    xlab("time (hrs)")+
    ylab("CFUs")+
    ggtitle(paste0(name.spec1," only"))
  gg.spec2 = ggplot(data=dat.plot.2) +
    # fitted curve
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
                        ", Lotka-Volterra model"),
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

#function to plot fitted landis data
plotter_landis.solo=function(parms, # fitted parameters
                             ode_fun, # ode function
                             parmnames, #name of ode function parameters (e.g. parnames.lv)
                             parmunits, #units of ode function parms (e.g. units.LV
                             dat,
                             specid, # vector of species names, for mapping
                             specmap.cur=spec.map, # map of species id to species name
                             reso=10 #plot points per hour.
){
  specname=spec.map[spec.map$name.data==specid,2][[1]]

  fit.pred=exper_pred(parms=parms,
                      x0=c(10),
                      ode_fun = ode_fun,
                      transf.num=6,
                      return.all=T,
                      reso=reso)
  dat.plot=as.data.frame(fit.pred$series.tot)
  names(dat.plot)=c("time","abund")

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
    ylab("CFUs")+
    ggtitle("")
  require(cowplot)
  require(gridExtra)
  title.gg <- ggplot() +
    labs(title = paste0(specname,
                        ", logistic model"),
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
                   gg.plot, gg.par, nrow=2, rel_heights = c(.15, 2))
  )
}



##########################
# For flexibility: write around ode_fun, define from library of models
ode_fun=ode_LV

times=seq(0,48,by=.2)
x0=c(x1=5,x2=5)
parms=c(r1=.2,r2=.2,a12=.5,a21=.5,k1=66000, k2=66000)
out.lsoda=lsoda(x0,times,ode_fun,parms)
matplot(out.lsoda[,1],out.lsoda[,2:3],type="l",xlab="time t",ylab="x1,x2",main="lsoda solutions");
plot(out.lsoda[,2:3],type='l',xlab='u',ylab='v',main="")


# Step 2: exper_pred to predict multiple transfers

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
    out.lsoda=lsoda(x0.cur,times,ode_fun,parms)
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

out=exper_pred(parms=parms,
               x0=x0,
               ode_fun = ode_fun,
               transf.num=6,
               return.all=T,
               reso=10)
matplot(out$series.tot[,1],out$series.tot[,2:3],type="l",xlab="time t",ylab="x1,x2",main="LV solutions");

## test logistic

out=exper_pred(parms=c(.2, 1000),
               x0=10,
               ode_fun = ode_log,
               transf.num=6,
               return.all=T,
               reso=10)
plot(out$series.tot)


# Step 3: obj_fun to implement above, find difference between data and predictions

# objective function for a single replicate
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

#function to simultaneously fit multiple replicates, potentially with different x0s
obj_helper=function(parms, #parameters for model implemented between transfers (first because optim)
                    x0.mat, #initial conditions at start of transfer 1; each row is a different replicate.
                    ## same order as dat.real.ls below
                    ode_fun, #function for model implemented between transfers
                    transf.num, #how many transfers to model (1 is basically no dealing with transfers)
                    dat.real.ls, #actual data, list of entries, each entry is a *matrix*
                    #with 1st column being transfer number, second being abund species 1, etc.
                    #should not include initial transfer unless we're fitting initial conditions as parameters
                    transf.dur=48, #how long are transfers? in hours
                    transf.dil = 5/200, #what is the dilution of each transfer. In Landis it's 5 microliter into 195 of new material, so 5/200
                    aug.ls=list() #misc list for anything extra, if needed
){
  if(nrow(x0.mat) != length(dat.real.ls)){stop("initial densities and real data have mismatched dimensions")}
  ss.tot=0
  for(i.rep in 1:nrow(x0.mat)){
    x0=x0.mat[i.rep,]
    dat.real=dat.real.ls[[i.rep]]
    if(sum(is.na(dat.real)) != 0){stop(paste0("NA in `dat.real.ls` entry number ", i.rep,". Check to make sure nothing weird is happing, remove transfer numbers with NAs."))}
    ss.cur=obj_ss(parms=parms,
                  x0=x0,
                  ode_fun = ode_fun,
                  transf.num=transf.num,
                  dat.real=dat.real,
                  transf.dur=transf.dur,
                  transf.dil=transf.dil,
                  aug.ls=aug.ls)
    # print(ss.cur)
    ss.tot=ss.tot+ss.cur
  }
  return(ss.tot)
}

prepped.solo=dataprep_landis.solo("1")

## test on logistic
obj_helper(parms=c(.2, 1000),
           x0=prepped.solo$x0.mat,
           ode_fun = ode_log,
           transf.num=6,
           dat.real.ls = prepped.solo$dat.real.ls)




parm_check = function(parms.cur, #current estimate
                      parm.lower, #lower bounds allowed
                      parm.upper, #upper bounds allowed
                      parmnames,
                      tol=0.001){ #name of parameters
  ## function to check if any parameters are at boundaries
  if(any(abs((parms.cur - parm.lower)/((parm.lower+parm.upper)/2))<tol)){
    warning(paste0("One or more parameters very near lower bounds:",
                   "\n   ",paste0((parmnames[abs((parms.cur - parm.lower)/((parm.lower+parm.upper)/2))<tol]), collapse="; ")
    ))
  }
  if(any(abs((parms.cur - parm.upper)/((parm.lower+parm.upper)/2))<tol)){
    warning(paste0("One or more parameters very near upper bounds:",
                   "\n   ",paste0((parmnames[abs((parms.cur - parm.upper)/((parm.lower+parm.upper)/2))<tol]), collapse="; ")
    ))
  }
}



########################################
# step 4: use optim to find best parms (logistic)
# #########################################
#
# Species 1
parms.guess=c(r=.1, k=60000)
specid="1"
dat.prepped=dataprep_landis.solo(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
modfit=plotter_landis.solo(parms=fit2$par,
                    ode_fun=ode_log,
                    parmnames = parmnames.log,
                    parmunits = units.log,
                    dat=dat.prepped$dat,
                    specid=specid)
modfit
ggsave(here("sandbox/figs",paste0("logistic-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)


### Species 2
parms.guess=c(r=.1, k=10000)
specid="2"
dat.prepped=dataprep_landis.solo(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
modfit=plotter_landis.solo(parms=fit2$par,
                           ode_fun=ode_log,
                           parmnames = parmnames.log,
                           parmunits = units.log,
                           dat=dat.prepped$dat,
                           specid=specid)
modfit
ggsave(here("sandbox/figs",paste0("logistic-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)



### Species 3
parms.guess=c(r=.1, k=2000)
specid="3"
dat.prepped=dataprep_landis.solo(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
modfit=plotter_landis.solo(parms=fit2$par,
                           ode_fun=ode_log,
                           parmnames = parmnames.log,
                           parmunits = units.log,
                           dat=dat.prepped$dat,
                           specid=specid)
modfit
ggsave(here("sandbox/figs",paste0("logistic-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species 4
parms.guess=c(r=.1, k=3000)
specid="4"
dat.prepped=dataprep_landis.solo(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
modfit=plotter_landis.solo(parms=fit2$par,
                           ode_fun=ode_log,
                           parmnames = parmnames.log,
                           parmunits = units.log,
                           dat=dat.prepped$dat,
                           specid=specid)
modfit
ggsave(here("sandbox/figs",paste0("logistic-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species A
parms.guess=c(r=.1, k=3000)
specid="A"
dat.prepped=dataprep_landis.solo(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
modfit=plotter_landis.solo(parms=fit2$par,
                           ode_fun=ode_log,
                           parmnames = parmnames.log,
                           parmunits = units.log,
                           dat=dat.prepped$dat,
                           specid=specid)
modfit
ggsave(here("sandbox/figs",paste0("logistic-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)


### Species B
parms.guess=c(r=.1, k=30000)
specid="B"
dat.prepped=dataprep_landis.solo(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
modfit=plotter_landis.solo(parms=fit2$par,
                           ode_fun=ode_log,
                           parmnames = parmnames.log,
                           parmunits = units.log,
                           dat=dat.prepped$dat,
                           specid=specid)
modfit
ggsave(here("sandbox/figs",paste0("logistic-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)


### Species C
parms.guess=c(r=.1, k=30000)
specid="C"
dat.prepped=dataprep_landis.solo(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
modfit=plotter_landis.solo(parms=fit2$par,
                           ode_fun=ode_log,
                           parmnames = parmnames.log,
                           parmunits = units.log,
                           dat=dat.prepped$dat,
                           specid=specid)
modfit
ggsave(here("sandbox/figs",paste0("logistic-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species D
parms.guess=c(r=.1, k=100000)
specid="D"
dat.prepped=dataprep_landis.solo(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_log,
           dat.real.ls=dat.prepped$dat.real.ls
           ##Methods stuff
)
modfit=plotter_landis.solo(parms=fit2$par,
                           ode_fun=ode_log,
                           parmnames = parmnames.log,
                           parmunits = units.log,
                           dat=dat.prepped$dat,
                           specid=specid)
modfit
ggsave(here("sandbox/figs",paste0("logistic-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

########################################
# step 4: use optim to find best parms (LV)
# #########################################

convergence.lv.ls=list()

## Note: define boundaries for plausible.
## Reminders:
##   r=0 means they cannot increase at low density with no competition. Reasonable lower bound, might put it slightly higher
##   r=1 means they double every hour
##   alpha = 0 means no effect of other species. Necessary lower bound (maybe better if slightly above this)
##   alpha = 10 means each interspec counts as 10 intraspec. This bound could be moved.
##   k = 5 means that the carrying capacity is the starting density. I've decided to push that higher.
##   I think a nice general rule is the upper bounds on k should be 1 order of magnitude beyond highest observed density
##   (this should probably be joint density)
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=1, r2=1, a12=5, a21=5, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("1","2")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
plotter_landis(parms=fit2$par,
               ode_fun=ode_fun,
               parmnames = parmnames.LV,
               parmunits=units.LV,
               dat.comp=dat.prepped$dat.comp,
               dat.solo1 = dat.prepped$dat.solo1,
               dat.solo2 = dat.prepped$dat.solo2,
               specid = specid,
               specmap.cur = spec.map)


mod.fit = plotter_landis(parms=fit2$par,
               ode_fun=ode_fun,
               parmnames = parmnames.LV,
               parmunits=units.LV,
               dat.comp=dat.prepped$dat.comp,
               dat.solo1 = dat.prepped$dat.solo1,
               dat.solo2 = dat.prepped$dat.solo2,
               specid = specid,
               specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)

convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)

########################
# Species 1 and 3
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=1, r2=1, a12=5, a21=5, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("1","3")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)

convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)

########################
# Species 1 and 4
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=1, r2=1, a12=5, a21=5, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("1","4")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)

########################
# Species 2 and 3
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=1, r2=1, a12=5, a21=5, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("2","3")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)

########################
# Species 2 and 4
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=1, r2=1, a12=5, a21=5, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("2","4")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)
###############################
# Species 3 and 4
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=1, r2=1, a12=5, a21=5, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("3","4")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)

###############################
# Species A and B
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=5, a21=5, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("A","B")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)

###############################
# Species A and C
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("A","C")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
# beep(4)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)


###############################
# Species A and D
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("A","D")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
# beep(4)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)

###############################
# Species B and C
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("B","C")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
# beep(4)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)
###############################
# Species B and D
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("B","D")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
# beep(4)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)
###############################
# Species C and D
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("C","D")
dat.prepped=dataprep_landis(specid)
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
parm_check(parms.cur=fit1$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
# beep(4)
fit2
parm_check(parms.cur=fit2$par,
           parm.lower=parm.lower,
           parm.upper=parm.upper,
           parmnames=parmnames.LV
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.LV,
                         parmunits=units.LV,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map)
mod.fit
ggsave(here("sandbox/figs",paste0("LV-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)
convergence.lv.ls[[length(convergence.lv.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)

saveRDS(convergence.lv.ls, here("sandbox","lv-convergence-list.RDS"))

# Step 5: report and plot results







#
# Systems test: does giving LV with known result give right answer?
#   Approach: generate right answer first.

# parm.sim=c(r1=.1,r2=.1,a12=1.1,a21=.7,k1=66000, k2=66000)
# out.sim=exper_pred(parm.sim,
#                    x0=c(5,5),
#                    ode_fun=ode_fun,
#                    reso=10,
#                    transf.num=6)
# matplot(out.sim$series.tot[,1],out.sim$series.tot[,2:3],type="l",xlab="time t",ylab="x1,x2",main="LV solutions");
#
# dat.real.sim=out.sim$transf.pred[c(1,3,6),]
#
# ## simple case: two experiments with identical results
# dat.ls.sim=list(dat.real.sim, dat.real.sim)
# x0.mat.sim=x0.mat[1:2,]
#
# parm.lower=c(r1=0, r2=0, a12 = 0, a21=0, k1 = 5, k2=5)
# parm.upper=c(r1=1, r2=1, a12=10, a21=10, k1=10^8, k2=10^8)
#
# fit1=optim(par=parm.sim*1.1,
#            fn=obj_helper,
#            ## ad'l args
#            x0.mat=x0.mat.sim,
#            transf.num=6,
#            ode_fun=ode_fun,
#            dat.real.ls=dat.ls.sim,
#            ## to keep values in sane ranges, so deSolve can handle them
#            method="L-BFGS-B",
#            lower=parm.lower,
#            upper=parm.upper)
#
# format(fit1$par, scientific = F)
#
# (parm.sim-fit1$par)/parm.sim
