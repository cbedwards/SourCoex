## do some initial plotting
##

library(tidyverse)
library(here)

## First, simple plot of data for species 2 and species 4 (their abs is pretty similar)
##
dat.cur=sour.data %>%
  filter(abund=="abs", spec1 == "2", spec2=="4")

gg=ggplot(dat.cur, aes(x=transf))+
  geom_point(aes(y=abund1, x=transf+.02), color='indianred', alpha=.5)+
  geom_point(aes(y=abund2,x=transf-.02), color='cornflowerblue', alpha=.5)+
  ggtitle(paste(spec.map$name.full[spec.map$name.data=="2"],
                " vs \n",
                spec.map$name.full[spec.map$name.data=="4"]))+
  xlab("transfer number")+
  ylab("abundance (CFUs)")

ggsave(here("sandbox/figs","example-data-specs2-4.jpg"),
       gg,
       width=8, height=5,
       units="in")




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
gg.fin=grid.arrange(gglist[[1]],
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
ggsave(here("sandbox/figs","rel-abund.jpg"),
       gg.fin,
       width=20, height=15,
       units="in")

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
gglist=list()
for(i in 1:nrow(specmat)){
  specid=specmat[i,]
  dat.abs = sour.data %>%
    filter(abund=="abs", spec1 == specid[[1]], spec2 == specid[[2]]) %>%
    mutate(combined=abund1+abund2)
  # View(dat.abs)

  gglist[[i]]=ggplot(dat.abs, aes(transf, combined, color = rep)) +
    geom_line()+
    ggtitle(paste(spec.map$name.full[spec.map$name.data==specid[[1]]],
                  " and \n",
                  spec.map$name.full[spec.map$name.data==specid[[2]]])
    )+
    # ylim(0,1)+
    xlab("Transfer number (each is 48 hrs)")+
    ylab("Total abundance")+
    theme(legend.position="none")

  # print(gglist[[i]])
}

library(cowplot)
library(gridExtra)
library(sjPlot)
length(gglist)
gg.fin=grid.arrange(gglist[[1]],
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
ggsave(here("sandbox/figs","total-abund.jpg"),
       gg.fin,
       width=20, height=15,
       units="in")


####
#### Abs densities, log sum total
####


specmat=unique(sour.data[,c("spec1","spec2")])
specmat=specmat[specmat$spec1 != specmat$spec2,]
gglist=list()
for(i in 1:nrow(specmat)){
  specid=specmat[i,]
  dat.abs = sour.data %>%
    filter(abund=="abs", spec1 == specid[[1]], spec2 == specid[[2]]) %>%
    mutate(lcombined=log(abund1+abund2))
  # View(dat.abs)

  gglist[[i]]=ggplot(dat.abs, aes(transf, lcombined, color = rep)) +
    geom_line()+
    ggtitle(paste(spec.map$name.full[spec.map$name.data==specid[[1]]],
                  " and \n",
                  spec.map$name.full[spec.map$name.data==specid[[2]]])
    )+
    # ylim(0,1)+
    xlab("Transfer number (each is 48 hrs)")+
    ylab("Log total abundance")+
    theme(legend.position="none")

  # print(gglist[[i]])
}

library(cowplot)
library(gridExtra)
library(sjPlot)
length(gglist)
gg.fin=grid.arrange(gglist[[1]],
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
ggsave(here("sandbox/figs","total-log-abund.jpg"),
       gg.fin,
       width=20, height=15,
       units="in")
