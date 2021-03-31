library(SourCoex)
library(here)







##########################
# For flexibility: write around ode_fun, define from library of models
ode_fun=ode_Til

times=seq(0,1,by=.01)
x0=c(x1=5/1000,x2=5/1000,R=195)
parms=c(r1=.02,r2=.02,d1=.01,d2=.01, R=195)
out.lsoda=lsoda(x0,times,ode_fun,parms)
matplot(out.lsoda[,1],out.lsoda[,2:3],type="l",xlab="time t",ylab="x1,x2",main="lsoda solutions");
plot(out.lsoda[,2:3],type='l',xlab='u',ylab='v',main="")


# Step 2: exper_pred to predict multiple transfers


out=exper_pred(parms=parms,
               x0=x0,
               ode_fun = ode_fun,
               transf.num=6,
               return.all=T,
               reso=10,
               aug.ls=list(Til=2))
matplot(out$series.tot[,1],out$series.tot[,2:3],type="l",xlab="time t",ylab="x1,x2",main="LV solutions");

## test 1-spec

out=exper_pred(parms=c(.002, .005, 195),
               x0=c(10/1000, 195),
               ode_fun = ode_Til_1spec,
               transf.num=6,
               return.all=T,
               reso=10,
               aug.ls=list(Til=1)
)
plot(out$series.tot)


# Step 3: obj_fun to implement above, find difference between data and predictions

# objective function for a single replicate




prepped.solo=dataprep_landis_solo("1",aug.ls=list(Til=1))

## test on logistic
SourTest_SS(parms=c(.2, 1000, 195),
           x0.mat=prepped.solo$x0.mat,
           ode_fun = ode_Til_1spec,
           transf.num=6,
           dat.real.ls = prepped.solo$dat.real.ls,
           aug.ls=list(Til=1))







########################################
# step 4: use optim to find best parms (single-species)
# #########################################
#NOTE: we are, I think, overfitting with 1 species (3 time points, 3 parameters). So don't be surprised if things get weird.
#
# Species 1 #############
parms.guess=c(r=.05, d=.01, R=50)
parm.lower=c(r1=0, d1=0.00001, R = .0001)
parm.upper=c(r1=100, d1=100, R = 10^10)
specid="1"
dat.prepped=dataprep_landis_solo(specid, aug.ls=list(Til=1))
##quickcheck: is our guess reasonable?
modfit=plotter_landis_solo(parms=parms.guess,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
modfit=plotter_landis_solo(parms=fit2$par,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_Til_1spec,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=1)
)

ggsave(here("sandbox/figs",paste0("Til-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)


### Species 2 #############
parms.guess=c(r=.005, d=.05, R=30)
parm.lower=c(r1=0, d1=0.00001, R = .0001)
parm.upper=c(r1=100, d1=100, R = 10^10)
specid="2"
dat.prepped=dataprep_landis_solo(specid, aug.ls=list(Til=1))
##quickcheck: is our guess reasonable?
modfit=plotter_landis_solo(parms=parms.guess,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
modfit=plotter_landis_solo(parms=fit2$par,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit
ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_Til_1spec,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=1)
)

ggsave(here("sandbox/figs",paste0("Til-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species 3 #############
parms.guess=c(r=.01, d=.05, R=15)
parm.lower=c(r1=0, d1=0.00001, R = .0001)
parm.upper=c(r1=100, d1=100, R = 10^10)
specid="3"
dat.prepped=dataprep_landis_solo(specid, aug.ls=list(Til=1))
##quickcheck: is our guess reasonable?
modfit=plotter_landis_solo(parms=parms.guess,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
modfit=plotter_landis_solo(parms=fit2$par,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit
ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_Til_1spec,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=1)
)
ggsave(here("sandbox/figs",paste0("Til-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species 4 #############
parms.guess=c(r=.005, d=.05, R=28)
parm.lower=c(r1=0, d1=0.00001, R = .0001)
parm.upper=c(r1=100, d1=100, R = 10^10)
specid="4"
dat.prepped=dataprep_landis_solo(specid, aug.ls=list(Til=1))
##quickcheck: is our guess reasonable?
modfit=plotter_landis_solo(parms=parms.guess,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
modfit=plotter_landis_solo(parms=fit2$par,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit
ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_Til_1spec,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=1)
)
ggsave(here("sandbox/figs",paste0("Til-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species A #############
parms.guess=c(r=.02, d=.05, R=30)
parm.lower=c(r1=0, d1=0.00001, R = .0001)
parm.upper=c(r1=100, d1=100, R = 10^10)
specid="A"
dat.prepped=dataprep_landis_solo(specid, aug.ls=list(Til=1))
##quickcheck: is our guess reasonable?
modfit=plotter_landis_solo(parms=parms.guess,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
modfit=plotter_landis_solo(parms=fit2$par,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit
ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_Til_1spec,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=1)
)
ggsave(here("sandbox/figs",paste0("Til-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species B #############
parms.guess=c(r=.0008, d=.02, R=600)
parm.lower=c(r1=0, d1=0.00001, R = .0001)
parm.upper=c(r1=100, d1=100, R = 10^10)
specid="B"
dat.prepped=dataprep_landis_solo(specid, aug.ls=list(Til=1))
##quickcheck: is our guess reasonable?
modfit=plotter_landis_solo(parms=parms.guess,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
modfit=plotter_landis_solo(parms=fit2$par,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit
ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_Til_1spec,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=1)
)
ggsave(here("sandbox/figs",paste0("Til-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species C - fitting issues #############
parms.guess=c(r=.00045, d=.05, R=600)
parm.lower=c(r1=0, d1=0.00001, R = 1)
parm.upper=c(r1=100, d1=100, R = 10^10)
specid="C"
dat.prepped=dataprep_landis_solo(specid, aug.ls=list(Til=1))
##quickcheck: is our guess reasonable?
modfit=plotter_landis_solo(parms=parms.guess,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           # method="L-BFGS-B",
           # upper=parm.upper,
           # lower=parm.lower
)
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           # method="L-BFGS-B",
           # upper=parm.upper,
           # lower=parm.lower
)
modfit=plotter_landis_solo(parms=fit2$par,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit
ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_Til_1spec,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=1)
)
ggsave(here("sandbox/figs",paste0("Til-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

### Species D #############
parms.guess=c(r=.008, d=.0002, R=100)
parm.lower=c(r1=0, d1=0.00001, R = .0001)
parm.upper=c(r1=100, d1=100, R = 10^10)
specid="D"
dat.prepped=dataprep_landis_solo(specid, aug.ls=list(Til=1))
##quickcheck: is our guess reasonable?
modfit=plotter_landis_solo(parms=parms.guess,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_Til_1spec,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=1),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2
modfit=plotter_landis_solo(parms=fit2$par,
                           ode_fun=ode_Til_1spec,
                           parmnames = parmnames.Til.1spec,
                           parmunits = units.Til.1spec,
                           dat=dat.prepped$dat,
                           specid=specid,
                           aug.ls=list(Til=1))
modfit
ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_Til_1spec,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=1)
)
ggsave(here("sandbox/figs",paste0("Til-spec-",specid,".jpg")),
       modfit,
       width=18, height=9)

## multi-species ---------------
##
## Unlike before, we now have 9 time points (3 per time series) and 5 parameters. So this is plausible to fit.

convergence.til.ls=list()

## Note: define boundaries for plausible.

# 1 and 2 ########################
#
# r=.05, d=.01, R=50
parm.lower=c(r1=0, r2=0, d1 = 0, d2=0, R=0.01)
parm.upper=c(r1=5, r2=5, d1 = 5, d2=5, R=10^10)
## Hard to guess what is a good starting point
parms.guess=c(r1=.0025,r2=.002,d1=.0001, d2=.01, R = 50)
ode_fun=ode_Til
specid=c("1","2")
dat.prepped=dataprep_landis(specid)

plotter_landis(parms=parms.guess,
               ode_fun=ode_fun,
               parmnames = parmnames.Til,
               parmunits=units.Til,
               dat.comp=dat.prepped$dat.comp,
               dat.solo1 = dat.prepped$dat.solo1,
               dat.solo2 = dat.prepped$dat.solo2,
               specid = specid,
               specmap.cur = spec.map,
               aug.ls=list(Til=2))

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=2),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=2),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2

## Quick diagnostics
ode_diagnostics(fit1,fit2,
                parms.guess=parms.guess,
                x0.mat=dat.prepped$x0.mat,
                ode_fun=ode_fun,
                transf.num=6,
                dat.real.ls=dat.prepped$dat.real.ls,
                aug.ls=list(Til=2)
)

mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.Til,
                         parmunits=units.Til,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map,
                         aug.ls=list(Til=2))
mod.fit
ggsave(here("sandbox/figs",paste0("Til-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)

convergence.til.ls[[length(convergence.Til.ls)+1]] = list(specid = specid,
                                                        convergence = fit2$convergence,
                                                        message = fit2$message)
# 1 and 3 ########################
#
# r=.05, d=.01, R=50
parm.lower=c(r1=0, r2=0, d1 = 0, d2=0, R=0.01)
parm.upper=c(r1=5, r2=5, d1 = 5, d2=5, R=10^10)
## Hard to guess what is a good starting point
parms.guess=c(r1=.0025,r2=.002,d1=.0001, d2=.01, R = 50)
ode_fun=ode_Til
specid=c("1","3")
dat.prepped=dataprep_landis(specid)

plotter_landis(parms=parms.guess,
               ode_fun=ode_fun,
               parmnames = parmnames.Til,
               parmunits=units.Til,
               dat.comp=dat.prepped$dat.comp,
               dat.solo1 = dat.prepped$dat.solo1,
               dat.solo2 = dat.prepped$dat.solo2,
               specid = specid,
               specmap.cur = spec.map,
               aug.ls=list(Til=2))
obj_helper(parms=parms.guess,
           ode_fun=ode_fun,
           dat.solo1 = dat.prepped$dat.solo1,
           dat.solo2 = dat.prepped$dat.solo2,
           specid = specid,
           specmap.cur = spec.map,
           aug.ls=list(Til=2))
fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=2),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=2),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2





mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.Til,
                         parmunits=units.Til,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map,
                         aug.ls=list(Til=2))
mod.fit
ggsave(here("sandbox/figs",paste0("Til-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)

convergence.til.ls[[length(convergence.til.ls)+1]] = list(specid = specid,
                                                          convergence = fit2$convergence,
                                                          message = fit2$message)




# 1 and 4 ########################
#
# r=.05, d=.01, R=50
parm.lower=c(r1=0, r2=0, d1 = 0, d2=0, R=0.01)
parm.upper=c(r1=5, r2=5, d1 = 5, d2=5, R=10^10)
## Hard to guess what is a good starting point
parms.guess=c(r1=.0025,r2=.002,d1=.0001, d2=.01, R = 50)
ode_fun=ode_Til
specid=c("1","4")
dat.prepped=dataprep_landis(specid)

plotter_landis(parms=parms.guess,
               ode_fun=ode_fun,
               parmnames = parmnames.Til,
               parmunits=units.Til,
               dat.comp=dat.prepped$dat.comp,
               dat.solo1 = dat.prepped$dat.solo1,
               dat.solo2 = dat.prepped$dat.solo2,
               specid = specid,
               specmap.cur = spec.map,
               aug.ls=list(Til=2))

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=2),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=2),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2





mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.Til,
                         parmunits=units.Til,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map,
                         aug.ls=list(Til=2))
mod.fit
ggsave(here("sandbox/figs",paste0("Til-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)

convergence.til.ls[[length(convergence.til.ls)+1]] = list(specid = specid,
                                                          convergence = fit2$convergence,
                                                          message = fit2$message)


# 2 and 3 ########################
#
# r=.05, d=.01, R=50
parm.lower=c(r1=0, r2=0, d1 = 0, d2=0, R=0.01)
parm.upper=c(r1=5, r2=5, d1 = 5, d2=5, R=10^10)
## Hard to guess what is a good starting point
parms.guess=c(r1=.0025,r2=.002,d1=.015, d2=.001, R = 45)
ode_fun=ode_Til
specid=c("2","3")
dat.prepped=dataprep_landis(specid)

plotter_landis(parms=parms.guess,
               ode_fun=ode_fun,
               parmnames = parmnames.Til,
               parmunits=units.Til,
               dat.comp=dat.prepped$dat.comp,
               dat.solo1 = dat.prepped$dat.solo1,
               dat.solo2 = dat.prepped$dat.solo2,
               specid = specid,
               specmap.cur = spec.map,
               aug.ls=list(Til=2))

fit1=optim(par=parms.guess,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=2),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit1
fit2=optim(par=fit1$par,
           fn=obj_helper,
           ## ad'l args
           x0.mat=dat.prepped$x0.mat,
           transf.num=6,
           ode_fun=ode_fun,
           dat.real.ls=dat.prepped$dat.real.ls,
           aug.ls=list(Til=2),
           ##Methods stuff
           method="L-BFGS-B",
           upper=parm.upper,
           lower=parm.lower
)
fit2





mod.fit = plotter_landis(parms=fit2$par,
                         ode_fun=ode_fun,
                         parmnames = parmnames.Til,
                         parmunits=units.Til,
                         dat.comp=dat.prepped$dat.comp,
                         dat.solo1 = dat.prepped$dat.solo1,
                         dat.solo2 = dat.prepped$dat.solo2,
                         specid = specid,
                         specmap.cur = spec.map,
                         aug.ls=list(Til=2))
mod.fit
ggsave(here("sandbox/figs",paste0("Til-specs-",paste0(specid, collapse="-"),".jpg")),
       mod.fit,
       width=18, height=13)

convergence.til.ls[[length(convergence.til.ls)+1]] = list(specid = specid,
                                                          convergence = fit2$convergence,
                                                          message = fit2$message)




##############LV below this --------


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
parm.lower=c(r1=1/10^10, r2=1/10^10, a12 = 0, a21=0, k1 = 1/1000, k2=1/1000)
parm.upper=c(r1=5, r2=5, a12=100, a21=100, k1=10^6, k2=10^6)
## Hard to guess what is a good starting point
parms.guess=c(r1=.001,r2=.001,a12=.5,a21=.9,k1=50, k2=50)
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
