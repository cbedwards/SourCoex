library(SourCoex)
library(here)







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




prepped.solo=dataprep_landis.solo("1")

## test on logistic
obj_helper(parms=c(.2, 1000),
           x0=prepped.solo$x0.mat,
           ode_fun = ode_log,
           transf.num=6,
           dat.real.ls = prepped.solo$dat.real.ls)







########################################
# step 4: use optim to find best parms (logistic)
# #########################################
#
# Species 1 #############
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


### Species 2 #############
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



### Species 3 #############
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

### Species 4 #############
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

### Species A #############
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


### Species B #############
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


### Species C #############
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

### Species D #############
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

## multi-species ---------------

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
##
# 1 and 2 ########################
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


# 1 and 3 ########################
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


# 1 and 4 ########################

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


# 2 and 3 ########################
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


# 2 and 4 ########################
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

# 3 and 4 ###############################

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


# A and B ###############################
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


# A and C ###############################
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



# A and D ###############################
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


# B and C ###############################
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

# B and D ###############################

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

#  C and D ###############################

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


# Cross-taxa competition ----------------

# A and 1 #############################

parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("A","1")
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


# A and 2 ###############################

parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("A","2")
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

# A and 3 ###############################
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("A","3")
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



# A and 4 ###############################
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("A","4")
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



# B and 1 ###############################
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("B","1")
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



# B and 2 ###############################
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("B","2")
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



# B and 3 ###############################
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("B","3")
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




# B and 4 ###############################
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("B","4")
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



# C and 1 ###############################
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("C","1")
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


# C and 2 ###############################
#
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("C","2")
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


# C and 3 ###############################
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("C","3")
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





# C and 4 ###############################

parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("C","4")
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




# D and 1 ###############################

parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("D","1")
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




# D and 2 ###############################

parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("D","2")
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




# D and 3 ###############################

parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("D","3")
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



# D and 4 ###############################
parm.lower=c(r1=0.01, r2=0.01, a12 = 0, a21=0, k1 = 10^3, k2=10^3)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.2,r2=.2,a12=.5,a21=.9,k1=2*10^5, k2=6000)
ode_fun=ode_LV
specid=c("D","4")
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
