
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
species <- "PolyMoor" #species name
dataset <- "DPolys23-8127" #dart order
missingness <- 0.3

# source my custom functions
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")

setwd('/Users/eilishmcmaster/Documents/PolyMoor/')

devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/ea46cc026bb56cafd339f5af383c94f46e0de2dd/read_dart_counts_csv_faster_new.r?raw=TRUE")

counts2 <- read_dart_counts_csv_faster('PolyMoor/dart_raw/Report_DPolys23-8127_SNPcount_3.csv', # import readcount data 
                                       minAlleleCount=0, 
                                       minGenotypeCount=0)

library(readxl)
m2 <- custom.read(species, dataset) #read custom metadata csv

#can deal with single sample groups and NA in the species_col function 
read_histogram_function <- function(meta, counts, filter_reads, species_col) {
  species <- unique(meta[[species_col]][!is.na(meta[[species_col]])])
  
  # anything that is NA must be 0 for the dividing etc
  counts$c1[is.na(counts$c1)] <- 0
  counts$c2[is.na(counts$c2)] <- 0
  
  # filter the reads
  combined_reads <- counts$c1 + counts$c2
  counts$c1[combined_reads < filter_reads] <- 0
  counts$c2[combined_reads < filter_reads] <- 0
  combined_reads <- counts$c1 + counts$c2 #update readcount
  
  # get the proportions for all (rows are samples)
  c3_min <- pmin(t(counts$c1), t(counts$c2), na.rm = TRUE) / t(combined_reads)
  c3_min[is.infinite(c3_min)] <- 1 # 1/0 is Inf, so making results where there are reads for one but not the other =1 AF
  c3_min[is.nan(c3_min)] <- NA # 0/0 is NaN, so removing results where there were no reads
  
  c3_max <- pmax(t(counts$c1), t(counts$c2), na.rm = TRUE) / t(combined_reads)
  c3_max[is.infinite(c3_max)] <- 1 # 1/0 is Inf, so making results where there are reads for one but not the other =1 AF
  c3_max[is.nan(c3_max)] <- NA # 0/0 is NaN, so removing results where there were no reads
  
  # removes full homozygotes otherwise theres waaaay too many to plot
  c3_min[c3_min==0] <- NA
  c3_min[c3_min==1] <- NA
  c3_max[c3_max==0] <- NA
  c3_max[c3_max==1] <- NA
  
  c3 <- cbind(c3_min, c3_max) 
  
  for (i in seq_along(species)) {
    print(paste("Running", species[i], "now"))
    samples <- meta$sample[meta[[species_col]] == species[i]]
    c3_species <- c3[row.names(c3) %in% samples, ] 
    
    if(isTRUE(nrow(c3_species)==0)){
      next
    }
    
    par(mfrow = c(4, 5), mai = c(0.5, 0.2, 0.2, 0.2))  # Set up a 2 x 2 plotting space
    par(cex.axis = 0.5)     # Adjust the value (0.8 in this example) to change the label size
    par(cex.main = 0.7)     # Adjust the value (0.9 in this example) to change the title size
    hist(c3_species, main = species[i], xlab = "", ylab = "", breaks = 50, col = "red", xaxt = 'n')
    axis(side = 1, at = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1), labels = c("0", "1/4", "1/3", "1/2", "2/3", "3/4", "1"), las=2)
    
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
          axis(side = 1, at = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1), labels = c("0", "1/4", "1/3", "1/2", "2/3", "3/4", "1"), las=2)
        }
      }
    }
  }
}




start <- Sys.time()
pdf(file="polystichum_all_10read_nomaf_pminmethod.pdf")
read_histogram_function(m2, counts2, 10, species_col = "sp") #needs meta, analysis column, counts data, and minimum number of reads per cell
dev.off()
fin <- Sys.time()


start <- Sys.time()
pdf(file="polystichum_all_0read_nomaf_pminmethod.pdf")
read_histogram_function(m2, counts2, 0, species_col = "sp") #needs meta, analysis column, counts data, and minimum number of reads per cell
dev.off()
fin <- Sys.time()

start <- Sys.time()
pdf(file="polystichum_all_20read_nomaf_pminmethod.pdf")
read_histogram_function(m2, counts2, 20, species_col = "sp") #needs meta, analysis column, counts data, and minimum number of reads per cell
dev.off()
fin <- Sys.time()

########################



read_histogram_function2 <- function(meta, counts, filter_reads, species_col, dms=NULL, remove_by_dms=NULL) {
  
  species <- unique(meta[[species_col]][!is.na(meta[[species_col]])])
  
  # anything that is NA must be 0 for the dividing etc
  counts$c1[is.na(counts$c1)] <- 0
  counts$c2[is.na(counts$c2)] <- 0
  
  # filter the reads
  combined_reads <- counts$c1 + counts$c2
  counts$c1[combined_reads < filter_reads] <- 0
  counts$c2[combined_reads < filter_reads] <- 0
  
  # get the proportions for all (rows are samples)
  c3_min <- pmin(t(counts$c1), t(counts$c2), na.rm = TRUE) / t(combined_reads)
  c3_min[is.infinite(c3_min)] <- 1 # 1/0 is Inf, so making results where there are reads for one but not the other =1 AF
  c3_min[is.nan(c3_min)] <- NA # 0/0 is NaN, so removing results where there were no reads
  
  c3_max <- pmax(t(counts$c1), t(counts$c2), na.rm = TRUE) / t(combined_reads)
  c3_max[is.infinite(c3_max)] <- 1 # 1/0 is Inf, so making results where there are reads for one but not the other =1 AF
  c3_max[is.nan(c3_max)] <- NA # 0/0 is NaN, so removing results where there were no reads
  
  # removes full homozygotes otherwise theres waaaay too many to plot
  c3_min[c3_min==0] <- NA
  c3_min[c3_min==1] <- NA
  c3_max[c3_max==0] <- NA
  c3_max[c3_max==1] <- NA
  
  c3 <- cbind(c3_min, c3_max) 
  
  out_data  <- list() # place to put the data 
  
  for (i in seq_along(species)) {
    print(paste("Running", species[i], "now"))
    samples <- meta$sample[meta[[species_col]] == species[i]] # get the NSW ID for that species samples
    if(isTRUE(remove_by_dms)){
      samples <- samples[samples %in% dms$sample_names]
    }
    
    c3_species <- c3[row.names(c3) %in% samples, ] # get the readcount df with those samples
    
    if (class(c3_species) %in% "array" || class(c3_species) %in% "matrix") {
      c3_species <- c3_species[, colSums(!is.na(c3_species) & c3_species != "") > 0] 
      out_data[[paste0(species[i])]] <- data.frame(c3_species)
    } else {
      c3_species <- c3_species[!is.na(c3_species)] # remove empty columns
      out_data[[paste0(species[i])]] <- as.vector(c3_species)
    }
  }
  
  return(out_data)
}


test <- read_histogram_function2(m2, counts2, 10, species_col="sp", dms=dms, remove_by_dms = TRUE) #needs meta, analysis column, counts data, and minimum number of reads per cell


####

whole_sp_plots <- function(data, species, max){
  plots <- list()
  for(i in seq_along(species)){
    x <- na.omit(unlist(list(data[[species[i]]])))
    
    p <- ggplot(data.frame(x = x), aes(x = x)) +
      geom_histogram(bins = 30, fill = "lightpink", color="black", size=0.2) + 
      scale_x_continuous(limits = c(0, 1),
                         breaks = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1),
                         labels = c("0", "1/4", "1/3", "1/2", "2/3", "3/4", "1")) +
      labs(x = element_blank(), y = element_blank(), title=paste0(species[i]))+
      theme_few()+
      # scale_y_continuous(expand=c(0,0), limits = c(0, max[i]))+
      geom_vline(xintercept = c(0.25, 0.333, 0.5, 0.666, 0.75), 
                 color = c("red","blue","black","blue","red"),
                 alpha = 0.3) +
      theme(axis.text.x = element_text(size=8, colour = c("black","red","blue","black","blue","red","black"), angle=90),
            axis.text.y = element_text(size=8),
            title = element_text(size=8, face="italic"),
            plot.margin = margin(0, 0, 0, 0))
    
    plots[[i]] <- p
  }
  return(plots)
}

# counts3 <- count_subsetter(dms, counts2, 0)


z <- whole_sp_plots(test, c("Polystichum moorei", "Polystichum whiteleggei", "Polystichum proliferum"), c(125,225,150))
sp_hist_plots <- ggarrange(z[[1]],z[[2]],z[[3]], align="hv", ncol=3,
                           labels=c("A","B","C"), font.label = list(size = 10, color = "black", face = "bold", family = NULL)) %>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count")

z[[1]]+ geom_vline(xintercept = 0.14)
# sp_hist_plots
ggsave("PolyMoor/outputs/plots/species_ploidy_hist2.png", plot = sp_hist_plots, width = 150, height = 60, dpi = 300, units = "mm")


###
specific_sample_plots <- function(data, samples){
  plots <- list()
  for(i in samples){
    print(i)
    x <- data[i,] %>% as.numeric()
    z <- data.frame(x=x)
    
    p <- ggplot(z, aes(x)) +
      geom_histogram(bins = 30, fill = "lightblue", color="black", size=0.2) + #usually use 50 bins
      scale_x_continuous(limits = c(0, 1),
                         breaks = c(0,0.25, 0.333, 0.5, 0.666, 0.75, 1),
                         labels = c("0", "1/4 ", "1/3", "1/2", "2/3", " 3/4", "1")) +
      labs(x = element_blank(), y = element_blank(), title=paste0(i))+
      theme_few()+
      # scale_y_continuous(expand=c(0,0), limits = c(0, max[i]))+
      geom_vline(xintercept = c(0.25, 0.333, 0.5, 0.666, 0.75), 
                 color = c("red","blue","black","blue","red"),
                 alpha = 0.3) +
      theme(axis.text.x = element_text(size=8, colour = c("black","red","blue","black","blue","red","black"), angle=90),
            axis.text.y = element_text(size=8),
            title = element_text(size=8),
            plot.margin = margin(0,0,0,0))
    
    
    plots[[paste0(i)]] <- p
    rm(x)
  }
  return(plots)
}

moor_samples <- specific_sample_plots(test[["Polystichum moorei"]],
                                      c("NSW1092305","NSW1159561","NSW1096976")) #"NSW1159620"

n1_samples <- specific_sample_plots(test[["Polystichum moorei"]],
                                    c("NSW1096976","NSW1096959","NSW1092300"))

white_samples <- specific_sample_plots(test[["Polystichum whiteleggei"]],
                                       c("NSW1159613","NSW1092324","NSW1159569"))

prof_samples <- specific_sample_plots(test[["Polystichum proliferum"]],
                                      c("NSW1159486","NSW1159487","NSW1159480"))

example_hist <- ggarrange(moor_samples[[1]],moor_samples[[2]],moor_samples[[3]],
                          white_samples[[1]],white_samples[[2]],white_samples[[3]],
                          prof_samples[[1]],prof_samples[[2]],prof_samples[[3]],
                          align="hv", ncol=3, nrow=3,
                          labels=c("A","","","","","","B","","","C","",""),
                          font.label = list(size = 10, color = "black", face = "bold", family = NULL))%>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count"
  )

all_hist <- ggarrange(z[[1]], moor_samples[[1]],moor_samples[[2]],moor_samples[[3]],
                      z[[2]], white_samples[[1]],white_samples[[2]],white_samples[[3]],
                      z[[3]], prof_samples[[1]],prof_samples[[2]],prof_samples[[3]],
                      align="hv", ncol=4, nrow=3,
                      labels=c("A","","","","B","","","","C","",""),
                      font.label = list(size = 10, color = "black", face = "bold", family = NULL))%>%
  annotate_figure(.,
                  bottom = "Allele frequency",
                  left="Count"
  )

example_hist
ggsave("PolyMoor/outputs/plots/example_ploidy_hist.png", plot = example_hist, width = 140, height = 160, dpi = 300, units = "mm")
ggsave("PolyMoor/outputs/plots/all_ploidy_hist.png", plot = all_hist, width = 190, height = 170, dpi = 300, units = "mm")

# ggsave("LantCama/outputs/specific_examples.png", plot = last_plot(), width = 190, height = 160, dpi = 300, units = "mm")


#####################

moor_samples <- m2[(m2$sp %in% c("Polystichum moorei")), ] %>% .$sample
dms_pm <- remove.by.list(dms, moor_samples) %>% remove.poor.quality.snps(., min_repro=0.96, max_missing=0.2)

count_subsetter2 <- function(dms, count){
  # ds <- dms$gt
  # keepers <- get_minor_allele_frequencies(ds)
  # ds <- ds[,which(keepers>=min)]
  # cat("Are there any NAs in the altcount data? ", any(is.na(ds)),"\n")
  # cat("Loci with NAs:")
  # print(table(apply(ds, 2, function(x) any(is.na(x)))))
  # 
  samples_tk <- dms$sample_names
  
  s_tk_location <- which(count$sample_names %in% samples_tk)
  
  count$c1 <- count$c1[,s_tk_location]
  count$c2 <- count$c2[,s_tk_location]
  count$sample_names <- colnames(count$c1)
  
  rownames(count$c1) <- count$locus_labels
  rownames(count$c2) <- count$locus_labels
  
  
  count$meta <- count$meta[,s_tk_location]
  count$sample_qc <- count$sample_qc[,s_tk_location]
  
  return(count)
}


count_pm <- count_subsetter(dms_pm, counts2, 0)

start <- Sys.time()
pdf(file="polystichum_moorei_0read_nomaf_pminmethod.pdf")
read_histogram_function(m2, count_pm, 0, species_col = "sp") #needs meta, analysis column, counts data, and minimum number of reads per cell
dev.off()
fin <- Sys.time()
