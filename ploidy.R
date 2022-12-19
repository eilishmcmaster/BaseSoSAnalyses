# Ploidy 

setwd("/Users/eilishmcmaster/Documents/ZierObco")

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/ea46cc026bb56cafd339f5af383c94f46e0de2dd/read_dart_counts_csv_faster_new.r?raw=TRUE")
# devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")
counts2 <- read_dart_counts_csv_faster('ZierObco/dart_raw/Report_DZ22-7321_SNPcount_3.csv', # import readcount data 
                                       minAlleleCount=1, 
                                       minGenotypeCount=0)


# plot a histogram for each sample in dms_zod

test2 <- count_subsetter(dms_zod, counts2, 0.05)

tr <-  t(test2$c1)
o <- lapply(split(tr,rownames(tr)), as.list)

tr2 <-  t(test2$c2)
o2 <- lapply(split(tr2,rownames(tr2)), as.list)

nn <- mergeLists_internal(o, o2)

minor <- lapply(nn, sapply, function(x) min(x)/sum(x))
a <- do.call(rbind, minor) #make matrix
major <- lapply(nn, sapply, function(x) max(x)/sum(x))
b <- do.call(rbind, major) #make matrix
c <- cbind(a,b)

# only for six samples at a time
par(mfrow = c(4, 5), mai=c(0.5,0.2,0.2,0.2))  # Set up a 2 x 2 plotting space

# Create the loop.vector (all the columns)
loop.vector <- 1:nrow(c)

for (i in loop.vector) { # Loop over loop.vector
  
  # store data in column.i as x
  x <- c[i,]
  
  # Plot histogram of x
  hist(x,breaks=20,
       main = paste(rownames(c)[i]),
       xlab = "",#"MAF reads/ total reads",
       ylab="",
       xlim = c(0, 1))
}

hist(c)


# plot histograms of the readcound data for all samples at the same time
# run by species dms
# minor allele frequency 0.05
par(mfrow = c(2, 2),mai=c(0.5,0.5,0.2,0.2)) 
doitall(dms_zi, counts2, 0.05, "Z. ingramii")
doitall(dms_zod, counts2, 0.05, "Z. odorifera")
doitall(dms_zo1, counts2, 0.05, "Z. obcordata (1)")
doitall(dms_zo2, counts2, 0.05, "Z. obcordata (2)")