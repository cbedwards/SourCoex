library(readr)
library(ggplot2)
library(viridis)
library(here)
library(SourCoex)

theme.mine=theme(plot.title = element_text(face="bold", size=24),
                 text=element_text(size=rel(5.5)),
                 legend.text=element_text(size=rel(3)),
                 strip.text = element_text(size = rel(3))
)



raw = read_excel("inst/extdata/Wolfe-sour-pilotgrowth.xlsx")
names(raw)
dat=raw
names(dat)=c("spec","strain","time","rep","dilut","num.cfus","density.cfus")
dat$time=gsub(" .*","",dat$time)
dat$time=as.numeric(dat$time)

unique(dat$spec)

# raw data
gp = ggplot(data = dat)+
  geom_point(aes(x = time, y = (density.cfus), color=rep))+
  geom_path(aes(x = time, y = (density.cfus), color=rep, group=rep))+
  facet_wrap(. ~ spec, scales="free")+
  scale_color_viridis()+
  theme.mine+
  xlab("hours")+
  ylab("density (cfus)")
ggsave(here("sandbox/figs/wolfe-pilot","rawdata.jpg"),
       gp,
       device="jpeg",
       width=14, height=12,
       units="in"
       )

#log scale plot of raw data
gp2 = gp + scale_y_log10()
ggsave(here("sandbox/figs/wolfe-pilot","rawdata-log.jpg"),
       gp2,
       device="jpeg",
       width=14, height=12,
       units="in"
)


#### Quick and dirty model fitting

exper_pred_cont(parms=c(r=1.5, k=600000),
                x0=180000,
                ode_fun=ode_log,
                times=unique(dat$time),
                )

