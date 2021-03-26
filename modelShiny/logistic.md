## Logistic Growth
 
The dashed line currently shows the trajectory of a population following logistic growth. The logistic growth equation is defined by 

\\(\frac{dN}{dt} = rN(1-\frac{N}{K})\\)

where N is the population size.


### Parameters
 - `r` is the "intrinsic growth rate" - this captures how fast the population would grow if it weren't limited by competition (e.g. exponential growth).
 - `k` is the "carrying capacity" - as the population size approaches the carrying capacity, competition reduces the effective growth rate, until when the popution size *is* `k`, the population is steady.

###  Intuitions
One of the key assumptions of the logistic growth equation, and its multi-species cousin, the Lotka-Volterra equation, is that there is an assumed constant supply of resources, such that there is some specific populations size that can persist indefinitely.
