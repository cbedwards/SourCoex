library(readxl)

## code to prepare `wine_data` dataset goes here
# retrieve paths to datafiles
first.file <- system.file(
  "extdata",
  "Figure_3_-_Source_Data_1.xlsx",
  package = "SourCoex"
)

# read the two .csv files
dat.readme= read_excel(
  first.file
)

dat.xfer1RA= read_excel(path=first.file,
                         sheet="xfer1RA"
)
dat.xfer1RA.ls=list(time="first transfer",
                    abund="relative",
                    data=dat.xfer1RA)
dat.xfer3RA= read_excel(path=first.file,
                        sheet="xfer3RA"
)
dat.xfer3RA.ls=list(time="third transfer",
                    abund="relative",
                    data=dat.xfer3RA)

dat.xfer6RA= read_excel(path=first.file,
                        sheet="xfer6RA"
)
dat.xfer6RA.ls=list(time="sixth transfer",
                    abund="relative",
                    data=dat.xfer6RA)

dat.xfer1AA= read_excel(path=first.file,
                        sheet="xfer1AA"
)
dat.xfer1AA.ls=list(time="first transfer",
                    abund="absolute",
                    data=dat.xfer1AA)

dat.xfer3AA= read_excel(path=first.file, ##Note: Dryad includes a random "2" about 40 cols in
                        ## I have removed this in the excel file.
                        sheet="xfer3AA"
)
dat.xfer3AA.ls=list(time="third transfer",
                    abund="absolute",
                    data=dat.xfer3AA)

dat.xfer6AA= read_excel(path=first.file,
                        ##Note: Dryad includes a random "2" about 40 cols in
                        ## I have removed this in the excel file.
                        sheet="xfer6AA"
)
dat.xfer6AA.ls=list(time="sixth transfer",
                    abund="absolute",
                    data=dat.xfer6AA)

sour_data = list(RA.xfer.1=dat.xfer1RA.ls,
                 RA.xfer.3=dat.xfer3RA.ls,
                 RA.xfer.6=dat.xfer6RA.ls,
                 AA.xfer.1=dat.xfer1AA.ls,
                 AA.xfer.3=dat.xfer3AA.ls,
                 AA.xfer.6=dat.xfer6AA.ls,
                 readme=dat.readme)

# save the wine_data dataframe as an .rda file in WineReviews/data/
usethis::use_data(sour_data, overwrite = TRUE)
