# This method is used to visualise suspected hybrids and their parent species
# This makes:
#   (1) A plot of PC1 vs sample heterozygosity. You should see parent 1 and 2 at opposite ends of the PC1 axis, with hybrids in the middle, and parent 1 and 2 at the same level oh Ho (y axis) and hybrids with clearly higher heterozygosity.
#   (2) A plot of genotypes at sites that are segregating in the two parent species. You should see fixed loci in both parents (fixed genotypes are 0 or 2) and heterozygous loci in the hybrids -- crosses would inherit one allele from each parent. 
# 
# Here I have already made the dms in my environment. In m2$sp3 I have my putative species and hybrid descriptions. 

#### Filter dms so that only putative parent species and hybrids exist ####

x <- m2[!(m2$sp3 %in% "Z. obcordata (2)"),] %>% .$sample
dms <- remove.by.list(dms2, x) %>% 
  remove_missing_loci_by_pop(., "sp3", 0.1) %>%
  remove.sample.by.pop.missingness(., "sp3", 0.3, 0.3)


#### Part 1: PCA PC1 vs sample heterozygosity ####
# get individual heterozygosity
Ho_sample <- rowMeans(dms$gt == 1, na.rm = TRUE)
names(Ho_sample) <- dms$sample_names

# do PCA
gen_d5 <- new("genlight", dms[["gt"]]) # convert df to genlight object for glPca function
gen_pca <- glPca(gen_d5, parallel=TRUE, nf=5) # do pca -- this method allows the input to have NAs 
g_pca_df <- gen_pca[["scores"]] # extract PCs 
g_pca_df2 <- merge(g_pca_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=FALSE) # add metadata 

# Match vector names to Row.names column in df
matching_rows <- match(names(Ho_sample), g_pca_df2$Row.names)

# Add vector as a new column in the df
g_pca_df2$sample_Ho <- NA
g_pca_df2$sample_Ho[matching_rows] <- Ho_sample

pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:5]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)") # create names for axes

# PCA plots 
pca_plot <- ggplot(g_pca_df2, aes(x=PC1, y=PC2, colour=sp2, label=site))+ 
  geom_point()+theme_bw()+xlab(pcnames[1])+ylab(pcnames[2])+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right")+
  guides(colour = guide_legend(title.position = "top", direction = "vertical"))

ggplot(g_pca_df2, aes(x=PC1, y=sample_Ho, colour=sp2, label=site))+ 
  geom_point()+theme_bw()+xlab(pcnames[1])+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right")+
  guides(colour = guide_legend(title.position = "top", direction = "vertical"))


#### Part 2: Tile plot of fixed loci ####
# Identify putative parental pop and hybrid samples based on PCA
species1 <- g_pca_df2$Row.names[which(g_pca_df2$PC1 > 10)]
species2 <- g_pca_df2$Row.names[which(g_pca_df2$PC1 < -0)]
hybrids <- g_pca_df2$Row.names[which(g_pca_df2$PC1 < 10 & g_pca_df2$PC1 > 0)]

# Get genotypes
ind_species1 <- match(species1, as.matrix(rownames(dms$gt)))
ind_species2 <- match(species2, as.matrix(rownames(dms$gt)))
ind_hybrids <- match(hybrids, as.matrix(rownames(dms$gt)))

# Filter the genotype table to remove loci missing in >6 samples (~20% of the dataset)
miss_thresh <- 6
imiss <- which(colSums(is.na(dms$gt)) >= miss_thresh)
dms$gt <- dms$gt[, -imiss]

# Identify reciprocally fixed loci in the parental pop samples
f_species1 <- colSums(dms$gt[ind_species1,], na.rm=TRUE) / (2 * colSums(!is.na(dms$gt[ind_species1,])))
f_species2 <- colSums(dms$gt[ind_species2,], na.rm=TRUE) / (2 * colSums(!is.na(dms$gt[ind_species2,])))

# Keep only loci where the absolute difference in allele frequencies between f_species1 and f_species2 is greater than 0.7
diagloc <- which(abs(f_species1 - f_species2) > 0.7)

# Filter the genotype table to only include the diagnostic loci
diaggt <- dms$gt[, diagloc]


# Create a vector of allele frequencies in f_species1 for the loci in diagloc
allele_freq <- f_species1[diagloc]

# Order the diaggt dataframe by allele frequency in f_species1
diaggt <- diaggt[, order(allele_freq)]

# Display the ordered dataframe
ordered_loci <- colnames(diaggt)

# Filter the genotype table to remove loci missing in >6 samples (~20% of the dataset)
miss_thresh <- 0.05
imiss <- which((colSums(is.na(diaggt))/(2*nrow(diaggt))) >= miss_thresh)
if(length(imiss)>=1){
  diaggt <- diaggt[, -imiss]
}

# Choose 100 of those loci at random
diagsamp <- sample(diaggt, 100, replace = FALSE)

# Plot an alignment (table) of those loci
# Order the table by locus and cluster
dialn.t <- data.frame(dms$meta$site, dms$meta$analyses[,"sp3"], diaggt)
colnames(dialn.t)[1:2] <- c("site", "sp3")
dialn.t <- arrange(dialn.t, desc(factor(sp3)))
dialn.t <- tibble::rownames_to_column(dialn.t, "irn")

#make long 
dialn.m <- reshape2::melt(dialn.t, id.vars=c("irn", "site", "sp3"))

ggplot(dialn.m, aes(fill=as.numeric(value), y=factor(irn, levels = irn_order), x=variable, label=sp3))+
  geom_tile() + 
  theme_bw() +
  ylab(element_blank()) +
  xlab(element_blank()) +
  scale_fill_gradient2(low="#66ccee", mid="#009988", high="#CCBB44", na.value = "white", midpoint=1, limits=c(0,2)) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")+
  facet_grid(rows=factor(dialn.m$sp3), scales = "free_y", space = "free_y")

