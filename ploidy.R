# Ploidy 

setwd("/Users/eilishmcmaster/Documents/ZierObco")

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/ea46cc026bb56cafd339f5af383c94f46e0de2dd/read_dart_counts_csv_faster_new.r?raw=TRUE")
# devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")
counts2 <- read_dart_counts_csv_faster('ZierObco/dart_raw/Report_DZ22-7321_SNPcount_3.csv', # import readcount data 
                                       minAlleleCount=1, 
                                       minGenotypeCount=0)

# plot histograms of the readcound data for all samples 
# run by species dms
# minor allele frequency 0.05
par(mfrow = c(2, 2),mai=c(0.5,0.5,0.2,0.2)) 
doitall(dms_zi, counts2, 0.05, "Z. ingramii")
doitall(dms_zod, counts2, 0.05, "Z. odorifera")
doitall(dms_zo1, counts2, 0.05, "Z. obcordata (1)")
doitall(dms_zo2, counts2, 0.05, "Z. obcordata (2)")