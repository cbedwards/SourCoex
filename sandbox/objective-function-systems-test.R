

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


SourTest_SS(parms=c(.2, 1000, 195),
           x0.mat=prepped.solo$x0.mat,
           ode_fun = ode_Til_1spec,
           transf.num=6,
           dat.real.ls = prepped.solo$dat.real.ls,
           aug.ls=list(Til=1))
