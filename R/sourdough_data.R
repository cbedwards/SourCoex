library(tidyverse)
library(readxl) ## apparently this is bad practice?
## but I'm not sure how it knows about the package otherwise.
## leaving it in for now.

## code to prepare `wine_data` dataset goes here
# retrieve paths to datafiles
filename <- system.file(
  "extdata",
  "Figure_3_-_Source_Data_1.xlsx",
  package = "SourCoex"
)

# read the two .csv files
spec.map= read_excel(#this tab was added by me, modified from the Readme tab for easier formatting.
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
  temp=rbind(select(data, Sp1, Sp2, starts_with("R1"))%>%
          rename(spec1=1, spec2=2, abund1=3, abund2=4),
        select(data, Sp1, Sp2, starts_with("R1"))%>%
          rename(spec1=1, spec2=2, abund1=3, abund2=4),
        select(data, Sp1, Sp2, starts_with("R3"))%>%
          rename(spec1=1, spec2=2, abund1=3, abund2=4),
        select(data, Sp1, Sp2, starts_with("R4"))%>%
          rename(spec1=1, spec2=2, abund1=3, abund2=4),
        select(data, Sp1, Sp2, starts_with("R5"))%>%
          rename(spec1=1, spec2=2, abund1=3, abund2=4))
  temp = cbind(temp,
               abund="rel",
               transf=1)
  return(temp)
}

sour.data=rbind(restruc_(dat.xfer1RA, abund="rel", transf=1),
               restruc_(dat.xfer3RA, abund="rel", transf=3),
               restruc_(dat.xfer6RA, abund="rel", transf=6),
               restruc_(dat.xfer1AA, abund="abs", transf=1),
               restruc_(dat.xfer3AA, abund="abs", transf=3),
               restruc_(dat.xfer6AA, abund="abs", transf=6)
)
sour.data=sour.data %>%
  relocate(abund, transf, spec1, spec2, abund1, abund2)
# save the wine_data dataframe as an .rda file in WineReviews/data/
usethis::use_data(sour.data, spec.map, overwrite = TRUE)
