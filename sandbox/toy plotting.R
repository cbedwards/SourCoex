## do some initial plotting
##

library(tidyverse)
specmat=unique(sour.data[,c("spec1","spec2")])
specmat=specmat[specmat$spec1 != specmat$spec2,]
gglist=NULL
for(i in 1:nrow(specmat)){
  specid=specmat[i,]
  dat.abs = sour.data %>%
    filter(abund=="rel", spec1 == specid[[1]], spec2 == specid[[2]])
  # View(dat.abs)

  gglist[[i]]=ggplot(dat.abs, aes(transf, abund1, color = rep)) +
    geom_line()+
    ggtitle(paste(spec.map$name.full[spec.map$name.data==specid[[1]]],
                  " vs \n",
                  spec.map$name.full[spec.map$name.data==specid[[2]]])
    )+
    ylim(0,1)+
    xlab("Transfer number (each is 48 hrs)")+
    ylab("Abund spec 1")+
    theme(legend.position="none")

  # print(gglist[[i]])
}

library(cowplot)
library(gridExtra)
library(sjPlot)
length(gglist)
grid.arrange(gglist[[1]],
             gglist[[2]],
             gglist[[3]],
             gglist[[4]],
             gglist[[5]],
             gglist[[6]],
             gglist[[7]],
             gglist[[8]],
             gglist[[9]],
             gglist[[10]],
             gglist[[11]],
             gglist[[12]],
             gglist[[13]],
             gglist[[14]],
             gglist[[15]],
             gglist[[16]],
             gglist[[17]],
             gglist[[18]],
             gglist[[19]],
             gglist[[20]],
             gglist[[21]],
             gglist[[22]],
             gglist[[23]],
             gglist[[24]],
             gglist[[25]],
             gglist[[26]],
             gglist[[27]],
             gglist[[28]],
             nrow=6)

##############
## Absolute densities
##############
specid="4"
dat.abs = sour.data %>%
  filter(abund=="abs", spec1 == specid, spec2 == specid)
# View(dat.abs)

ggplot(dat.abs, aes(transf, log(abund1), color = rep)) +
  geom_line()+
  geom_point()+
  ggtitle(paste(spec.map$name.full[spec.map$name.data==specid],
                " vs \n",
                spec.map$name.full[spec.map$name.data==specid])
  )+
  # ylim(0,1)+
  xlab("Transfer number (each is 48 hrs)")+
  ylab("Abund spec 1")+
  theme(legend.position="none")


####
#### Abs densities, sum total
####


specmat=unique(sour.data[,c("spec1","spec2")])
specmat=specmat[specmat$spec1 != specmat$spec2,]
gglist=NULL
for(i in 1:nrow(specmat)){
  specid=specmat[i,]
  dat.abs = sour.data %>%
    filter(abund=="rel", spec1 == specid[[1]], spec2 == specid[[2]])
  # View(dat.abs)

  gglist[[i]]=ggplot(dat.abs, aes(transf, abund1, color = rep)) +
    geom_line()+
    ggtitle(paste(spec.map$name.full[spec.map$name.data==specid[[1]]],
                  " vs \n",
                  spec.map$name.full[spec.map$name.data==specid[[2]]])
    )+
    ylim(0,1)+
    xlab("Transfer number (each is 48 hrs)")+
    ylab("Abund spec 1")+
    theme(legend.position="none")

  # print(gglist[[i]])
}

library(cowplot)
library(gridExtra)
library(sjPlot)
length(gglist)
grid.arrange(gglist[[1]],
             gglist[[2]],
             gglist[[3]],
             gglist[[4]],
             gglist[[5]],
             gglist[[6]],
             gglist[[7]],
             gglist[[8]],
             gglist[[9]],
             gglist[[10]],
             gglist[[11]],
             gglist[[12]],
             gglist[[13]],
             gglist[[14]],
             gglist[[15]],
             gglist[[16]],
             gglist[[17]],
             gglist[[18]],
             gglist[[19]],
             gglist[[20]],
             gglist[[21]],
             gglist[[22]],
             gglist[[23]],
             gglist[[24]],
             gglist[[25]],
             gglist[[26]],
             gglist[[27]],
             gglist[[28]],
             nrow=6)
