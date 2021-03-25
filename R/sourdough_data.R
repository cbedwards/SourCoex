## but I'm not sure how it knows about the package otherwise.
## leaving it in for now.

## code to prepare `wine_data` dataset goes here
# retrieve paths to datafiles
filename <- here("inst/extdata","Figure_3_-_Source_Data_1.xlsx")

# read the two .csv files
spec.map=read_excel(#this tab was added by me, modified from the Readme tab for easier formatting.
  filename,
  sheet="speciesMap"
)

dat.xfer1RA= read_excel(path=filename,
                        sheet="xfer1RA"
)
dat.xfer3RA= read_excel(path=filename,
                        sheet="xfer3RA"
)
dat.xfer6RA= read_excel(path=filename,
                        sheet="xfer6RA"
)
dat.xfer1AA= read_excel(path=filename,
                        sheet="xfer1AA"
)
dat.xfer3AA= read_excel(path=filename, ##Note: Dryad includes a random "2" about 40 cols in
                        ## I have removed this in the excel file.
                        sheet="xfer3AA"
)
dat.xfer6AA= read_excel(path=filename,
                        ##Note: Dryad includes a random "2" about 40 cols in
                        ## I have removed this in the excel file.
                        sheet="xfer6AA"
)

#helper function: restructuring from the specific structure given in Dryad.
restruc_ = function(data, #data is the original data from Landis et al sheet.
                    abund, #abund is the type of abundance ("rel" or "abs" for relative or absolut)
                    transf #which transfer was it? Numeric.
){
  temp=rbind(
    select(data, Sp1, Sp2, starts_with("R1"))%>%
      rename(spec1=1, spec2=2, abund1=3, abund2=4) %>%
      mutate(rep = 1),
    select(data, Sp1, Sp2, starts_with("R2"))%>%
      rename(spec1=1, spec2=2,abund1=3, abund2=4) %>%
      mutate(rep = 2),
    select(data, Sp1, Sp2, starts_with("R3"))%>%
      rename(spec1=1, spec2=2, abund1=3, abund2=4) %>%
      mutate(rep = 3),
    select(data, Sp1, Sp2, starts_with("R4"))%>%
      rename(spec1=1, spec2=2, abund1=3, abund2=4) %>%
      mutate(rep = 4),
    select(data, Sp1, Sp2, starts_with("R5"))%>%
      rename(spec1=1, spec2=2, abund1=3, abund2=4) %>%
      mutate(rep = 5)
  )
  temp = cbind(temp,
               abund=abund,
               transf=transf)
  temp = temp %>%
    relocate(abund, spec1, spec2, rep, transf, abund1, abund2) %>%
    arrange(abund, spec1, spec2, rep)
  return(temp)
}

sour.data=rbind(restruc_(dat.xfer1RA, abund="rel", transf=1),
                restruc_(dat.xfer3RA, abund="rel", transf=3),
                restruc_(dat.xfer6RA, abund="rel", transf=6),
                restruc_(dat.xfer1AA, abund="abs", transf=1),
                restruc_(dat.xfer3AA, abund="abs", transf=3),
                restruc_(dat.xfer6AA, abund="abs", transf=6)
)

## turn relative measures from percents to proportions (math is easier)

sour.data[sour.data$abund=="rel",c("abund1","abund2")]=
  sour.data[sour.data$abund=="rel",c("abund1","abund2")]/100

## transf1, relative
## First, be dumb, make it all interspec ratios

transf0.rel=data.frame(abund="rel",
                       transf=0,
                       sour.data[sour.data$transf==1 & sour.data$abund =="rel",
                                 c("spec1","spec2", "rep")],
                       abund1=.5,
                       abund2=.5)
transf0.rel=transf0.rel[,names(sour.data)]
## where spec1 = spec2, set abund1 = 1 and abund2 = 0
transf0.rel[transf0.rel$spec1==transf0.rel$spec2, c("abund1")] = 1
transf0.rel[transf0.rel$spec1==transf0.rel$spec2, c("abund2")] = 0

## to make absolute, just change "abund" label and multiply relative abundance by initial densities
## Note: our measures elsewhere are in CFU per microliter, so we are expressing it that way here.
## 2000 CFUs innoculant into 200 microliters, so 10 CFUs/microliter
transf0.abs=transf0.rel
transf0.abs$abund="abs"
transf0.abs[,c("abund1", "abund2")]=transf0.abs[,c("abund1", "abund2")]*10

sour.data=rbind(sour.data, transf0.rel, transf0.abs)

sour.data$rep = as.factor(sour.data$rep)

# save the wine_data dataframe as an .rda file in WineReviews/data/
usethis::use_data(sour.data, spec.map, overwrite = TRUE)
