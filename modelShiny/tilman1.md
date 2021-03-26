## Tilman single-species model
 
This is a modified version of the Tilman R* model. This model explicitly models resources available - in the Tilman model, this is assumed to be a resource pool with some constant inflow rate (e.g. nitrogen int eh soil, with some constant amount of nitrogen being added per year through stone weathering). We have tweaked this - we still explicitly model a resource pool, but we have zero resource inflow. Instead, we have some amount of resource that microbes have access to at the start of each "transfer period", which they then consume by growing. THis gives us the following set of equations:

\\(\frac{dN}{dt} = N(rR - d)\\)
\\(\frac{dR}{dt} = - NrR)\\)

where N is the population size of the microbe, and R is the amount of resources available. Layered on top of this is an additional parameter, $R_0$, which is the amount of resource available at the beginning of the experiment and at the start of each transfer period.


### Parameters
- `r` is the per capita conversion rate of resources into more microbes. It's easy to think of this as an intrinsic growth rate, but it's actually not - since it's being multiplied by the amount of resouces, with near-infinite resources, the population grows near-infinitely fast.
- `d` is the per capita mortality rate.
- `$R_0$` is the amount of resources available per transfer. I'm fairly sure this is a distinct parameter from `r`, but it's possible that the equations can be rescaled to avoid using this. Regardless, large values mean the population can grow quite large before running out of resources -- but also means the population will grow faster! 

(Note that I should swap back into using `a` instead of `r` for the growth/conversion rate. This both clarifies that it's not a growth rate, and avoids the "rR" terms, which look weird).