
library(openxlsx) #reading and writing xlsx
library(stringr) #replacing strings while data wrangling 
library(dplyr) #data wrangling 
library(data.table) #data wrangling 
library(ggfortify) # calculating glm and pca
library(ggpubr) # ggarrange
library(pracma) # maths
library(RRtools) #Jason's package for dart data
library(ggthemes) #themes for ggplots
library(RColorBrewer) #used for making colour scemes for plots
library(ozmaps) #draws australia coastlines and state boundaries
library(adegenet) #essential for processing dart data
library(ggrepel) #used for plotting labels on ggplote

topskip   <- 6
nmetavar  <- 18
RandRbase <- "" #main directory 
species <- "LantCama" #species name
dataset <- "DLan22-7500" #dart order
missingness <- 0.3

# source my custom functions
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")

setwd("/Users/eilishmcmaster/Documents/LantCama")

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/ea46cc026bb56cafd339f5af383c94f46e0de2dd/read_dart_counts_csv_faster_new.r?raw=TRUE")

counts2 <- read_dart_counts_csv_faster('LantCama/dart_raw/Report_DLan22-7500_3_moreOrders_SNPcount_3.csv', # import readcount data 
                                       minAlleleCount=1, 
                                       minGenotypeCount=0)

library(readxl)
m2 <- custom.read(species, dataset) #read custom metadata csv

#can deal with single sample groups and NA in the species_col function 
read_histogram_function <- function(meta, counts, filter_reads, species_col) {
  species <- unique(meta[[species_col]][!is.na(meta[[species_col]])])
  
  # filter the reads
  combined_reads <- counts$c1 + counts$c2
  counts$c1[combined_reads < filter_reads] <- NA
  counts$c2[combined_reads < filter_reads] <- NA
  
  # get the proportions for all (rows are samples)
  c3_min <- pmin(t(counts$c1), t(counts$c2), na.rm = TRUE) / t(combined_reads)
  c3_max <- pmax(t(counts$c1), t(counts$c2), na.rm = TRUE) / t(combined_reads)
  c3 <- cbind(c3_min, c3_max) 
  
  for (i in seq_along(species)) {
    print(paste("Running", species[i], "now"))
    samples <- meta$sample[meta[[species_col]] == species[i]]
    c3_species <- c3[row.names(c3) %in% samples, ] 
    
    par(mfrow = c(4, 5), mai = c(0.5, 0.2, 0.2, 0.2))  # Set up a 2 x 2 plotting space
    
    hist(c3_species, main = species[i], xlab = "", ylab = "", breaks = 50, col = "red", xaxt = 'n')
    axis(side = 1, at = c(0, 0.25,  0.5,  0.75, 1), labels = c(0, 0.25,  0.5,  0.75, 1))
    
    if (class(c3_species) %in% "array" || class(c3_species) %in% "matrix"){
      loop.vector <- 1:nrow(c3_species)
      for (i in loop.vector) { # Loop over loop.vector
        
        # store data in row.i as x
        x <- c3_species[i,]
        if (sum(x, na.rm=TRUE) > 0) { # skip empties
          # Plot histogram of x
          hist(x, breaks = 50,
               main = paste(rownames(c3_species)[i]),
               xlab = "",#"MAF reads/ total reads",
               ylab = "",
               xlim = c(0, 1),
               xaxt = 'n')
          axis(side = 1, at = c(0, 0.25,  0.5,  0.75, 1), labels = c(0, 0.25,  0.5,  0.75, 1))
        }
      }
    }
  }
}

        


start <- Sys.time()
pdf(file="lantana_all_10read_nomaf_pminmethod.pdf")
read_histogram_function(m3, counts2, 10, species_col = "sp") #needs meta, analysis column, counts data, and minimum number of reads per cell
dev.off()
fin <- Sys.time()

