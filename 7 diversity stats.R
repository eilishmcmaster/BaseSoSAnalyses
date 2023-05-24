# Diversity statistics 
# We use diveRsity package to calculate these stats. It is CRUCIAL that only one genetic group is 
# considered at a time, because this method is based upon shared common alleles. These groups should be species 
# or genetic groups identified by PCA and UMAP. They are first made into group specific dms data, then stats are 
# calculated. Within a genetic group there can be 1 or more sites (collection sites).
# This method will calculate the diversity statistics by site.
# 
# This method requires a genpop file to be written, however the base method does not allow single site groups. I have 
# made the function `single_site_genepop_basicstats` to overcome this. The other function `multi_site_genepop_basicstats`
# does the exact same thing that Jason's original functions did, but processes the dataframe within the function so the
# final output is more manageable. 


# New methods ####################################################################
# specify loci missingness threshold and minor allele frequency 
max_missing <- 0.3
maf_stats <- 0.05

# Single site example
# this is an outgroup species with only one collection site of 6 samples 
# make a dms data for only the species Z. odorifera and filter missingness and reproducability 
dms_zod <- remove.by.list(dms, m2[(m2$sp %in% "Z. odorifera"),] %>%.$sample) %>%
  remove.poor.quality.snps(., min_repro=0.96,max_missing=max_missing)

# calculate stats for that dms
zo <- single_site_genepop_basicstats(dms_zod, maf_stats, "Zo")


# Multi site example 
dms_zo1<- remove.by.list(dms, m2[(m2$sp %in% "Z. obcordata (1)"),] %>%.$sample) %>%
  remove.poor.quality.snps(., min_repro=0.96,max_missing=max_missing)

zo1<- multi_site_genepop_basicstats(dms_zo1, maf_stats,"obcordata1", dms_zo1$meta$site)

# combine stats 
stats <- rbind(zo, zo1)

## multiple species at the same time (species are calculated separately)

# dms has all species under `sp`, the maf filter is 0.05, and each species is filtered for loci with 0.3 missingness separately 

species_stats <- multispecies_stats(dms, 0.05) 

# Multiple species with multiple sites simultaneously 

site_stats <- species_site_stats(dms, 0.05, "genetic_group2", "site2")
# dms= dms with samples of interest 
# 0.05 = MAF threshold
# "genetic_group2" = the genetic grouping to use, make sure each group is a single structure group or results will be biased
# site2 = the column name of the sites

# Original method ####################################################################
# multi site groups only 
# same as new multi site example 

dms_subset <- dms_zo1 # species specific dms as calculated above 
sppop_freq <- as.data.frame(table(dms_subset$meta$site)) # get the numbers of samples per site 

gp   <- dart2genepop(dms_subset, RandRbase, species, dataset, dms_subset$meta$site, maf_val=0.05)
#note that in dart2genepop, you can set your own min. allele frequency (maf_val). now it is set as maf_val=0.05
#the dart2genepop() also removes any loci that are missing for a population

#make sure that populations consisting of 1 individual have been removed as they can stuff up the scripts
library(diveRsity)
bs <- basicStats(infile = gp, outfile = NULL, #this step takes ages to run
                 fis_ci = FALSE, ar_ci = TRUE, 
                 ar_boots = 999, 
                 rarefaction = FALSE, ar_alpha = 0.05)

population_names  <- as.data.frame(names(bs$main_tab))#ls() rearranges the names 
colnames(population_names)[1] <- "sample"
help <- merge(population_names, m2[,c("site", "sample")], by="sample", all.x=TRUE, all.y=FALSE)

npop <- length(population_names$sample)
result <- as.data.frame(mat.or.vec(npop,11))
measurement_names <- rownames(bs$main_tab[[1]])
population_names  <- names(bs$main_tab) #ls() rearranges the names 
rownames(result) <- population_names
colnames(result) <- measurement_names

for (r in 1:npop) {
  popstats <- bs$main_tab[[r]][,"overall"] ##extract from a list
  result[r,] <- popstats}

result <- as.data.frame(result)
result$sample <- rownames(result)

div_stats <-result[,c("size","obs_het","exp_het","fis","sample")]
# div_stats$sppops <- rownames(div_stats)
div_stats <- merge(help, div_stats, by="sample")
div_stats$sample <- NULL
# div_stats <- merge(div_stats, sppop_freq, by.x="site", by.y="Var1")
# colnames(div_stats) <- c("Subpopulation","Ho","He","Fis","n")

div_stats


# Venn diagram of alleles ########################################################
# This function makes a venn diagrams of the alleles genetic groups have. Not super polished, takes a while.
library(ggVennDiagram)
# variables are dms, grouping variable, and MAF
venn <- venner(dms, dms$meta$analyses[,"sp"], 0.05)

big_venn <- ggVennDiagram(venn, label_alpha = 0, edge_size = 0, label_size=2)+ #label="percent"
  scale_fill_gradient(low="lightblue",high = "white", trans="log")+theme(legend.position = "none")

# venn3 <- venner(dms, dms$meta$analyses[,"venn"], 0.05)
# 
# little_venn <- ggVennDiagram(venn3, label_alpha = 0, edge_size = 0, label_size=4)+ #label="percent"
#   scale_fill_gradient(low="lightblue",high = "white", trans="log")+theme(legend.position = "none")

big_venn

# save(dms, venner, matcher2, filter, file = "venner_example.RData")

