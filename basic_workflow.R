library(adegenet)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggtree)
library(igraph)
library(openxlsx)
library(RColorBrewer)
library(RRtools)
library(RSplitsTree)
library(scales)
library(stringr)
library(tanggle)
library(tidyr)
library(vegan)

source('https://github.com/eilishmcmaster/SoS_functions/blob/33bd7065baa91b7e2a1800b6481d63573fb38d28/dart2svdquartets.r?raw=TRUE')
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")


topskip   <- 6
nmetavar  <- 18
RandRbase <- "" #main directory
species <- "DodoVisc" #species name
dataset <- "DDo23-8381" #dart order
basedir <- ""

d1        <- new.read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=FALSE)
meta      <- read.meta.data.full.analyses.df(d1, basedir, species, dataset)

d2        <- dart.meta.data.merge(d1, meta) 

# d2 <- remove.by.list(d2, d2$sample_names[d2$meta$analyses[,'situ']!='herbarium'])

d3        <- remove.poor.quality.snps(d2, min_repro=0.96, max_missing=0.2)%>% remove.fixed.snps()
d4        <- sample.one.snp.per.locus.random(d3, seed=12345) 

d5 <- remove.by.missingness(d4, 0.2)

dms <- d5

length(dms$locus_names)
length(dms$sample_names)

dms_maf5 <- remove.by.maf(dms, 0.05)

m2 <- dms$meta$analyses %>% as.data.frame()
m2$lat <- as.numeric(m2$lat)
m2$long <- as.numeric(m2$long)

# 
# # Define colors and shapes for species (instead of sites)
# species_levels <- factor(m2$pop[which(!is.na(m2$toothbrush))])  # Use species column "sp"
# species_colours <- named_list_maker(levels(species_levels), palette = 'Set3',4)
# species_colours['G. DodoVisc (Enmore)'] <- '#e0e06c'
# pop_colours <- species_colours
# names(pop_colours) <- str_extract(names(species_colours), "\\(([^)]+)\\)") %>%
#   str_remove_all("[()]")
# 
# #create annotations
# enmore_cols <- named_list_maker(m2[which(m2$pop=="Enmore"), 'site'], palette = "YlGn", num = 3)
# torrington_cols <- named_list_maker(m2[which(m2$pop=="Torrington"), 'site'], palette = "RdPu", num = 3)
# chambigne_cols <- named_list_maker(m2[which(m2$pop=="Chambigne"), 'site'], palette = "Blues", num = 3)
# gf_cols <- named_list_maker(m2[which(m2$pop=="Guy Fawkes"), 'site'], palette = "Purples", num = 3)
# site_cols <- c(enmore_cols, torrington_cols, chambigne_cols, gf_cols)
# 
# 
# species_shapes <- setNames(1:length(levels(species_levels)), levels(species_levels))  # Assign shapes


#### RAW GENOTYPE HEATMAP ####
submat <- dms$gt[,sample(1:ncol(dms$gt), 1000, replace = FALSE)] %>% as.matrix
submat[is.na(submat)] <- -9

row_ann <- rowAnnotation(Group = m2$pop[match(rownames(submat), m2$sample)])

# Define a set of discrete color values
discrete_colors <- c("-9" = "white", "0" = "lightblue", "1" = "coral4", "2" = "orange")

# Apply the discrete color scheme to the heatmap
gt <- Heatmap(submat, 
              col = discrete_colors,
              use_raster = FALSE,
              right_annotation = row_ann,
              cluster_columns = TRUE, cluster_rows = TRUE,
              column_names_gp = gpar(fontsize = 0),
              row_names_gp = gpar(fontsize = 6))

# Set up the PNG device
pdf("DodoVisc/outputs/2025/DodoVisc_gt.pdf", width = 25, height = 20)

# Draw the plot
draw(gt)

# Turn off the PNG device
dev.off()


#### KINSHIP ####
kin <- individual_kinship_by_pop(dms, RandRbase, species, dataset, dms$meta$analyses[,"pop"], maf=0.05, mis=0.2, as_bigmat=TRUE)

kin2 <- as.data.frame(kin) %>%mutate_all(~replace(.,.<0.35355339, 0)) #removes all of the pairwise connections that are k<0.45
kin3 <- kin2

network <- graph_from_adjacency_matrix(as.matrix(kin3), mode="undirected", diag=F,weighted=T) #makes the network based on k>0.45
plot(network)

ceb <- cluster_fast_greedy(network) # make the bubbles for the network plot
ceb_net_plot <- plot(ceb, network, vertex.label.color="transparent", vertex.size=2, edge.width=0.4) #make the network plot
ceb_net_plot

##### GET CLONES ####

clones <- as.data.frame(cbind(genet=ceb$membership, sample=ceb$names)) #get the clones from the network as a df
clones_out <- merge(clones, m2[,c("sample","lat","long","site","pop")], by="sample", all.x = TRUE) #add some metadata

clones_out <- clones_out %>%
  arrange(pop, genet)
clones_out$genet_original <- clones_out$genet
old_genet_unique <- unique(clones_out$genet)
new_genet_unique <- 1:length(old_genet_unique)
clones_out$genet <- new_genet_unique[match(clones_out$genet,old_genet_unique)]

write.xlsx(clones_out, "DodoVisc/outputs/clones_out.xlsx", rowNames=FALSE, colNames=TRUE)

# Remove duplicate clones from working data
#close clones with highest quality data (least missingess)
missingness_gt <- as.data.frame(rowSums(is.na(dms[["gt"]]))/ncol(dms[["gt"]])) #get missingness of all samples
colnames(missingness_gt) <- "missingness"
clone_missing <- merge(clones_out, missingness_gt, by.x="sample", by.y=0) # merge the clone data
clone_missing <- clone_missing[order(clone_missing$missingness),] #order by missingness low to high

unique_genets <- distinct(clone_missing, genet, .keep_all=TRUE) #keeps the top result for each genet (lowest missingness)

#make a list removing duplicate clones
non_clones <- dms$sample_names[which(!dms$sample_names %in% clones_out$sample)]
approved_clones <- unique_genets$sample

clones_dereplicated <- c(non_clones, approved_clones) # list of sample to keep

# remove clones from dms
dms_no_clones <- remove.by.list(dms, clones_dereplicated)

##### GET FAMILIES ####

fam_network <- graph_from_adjacency_matrix(as.matrix(kin), mode="undirected", diag=F,weighted=T) #makes the network based on k>0.45
plot(fam_network)

fam_ceb <- cluster_fast_greedy(fam_network) # make the bubbles for the network plot
fam_ceb_net_plot <- plot(fam_ceb, fam_network, vertex.label.color="transparent", vertex.size=2, edge.width=0.4) #make the network plot
fam_ceb_net_plot

families <-as.data.frame(cbind(family=fam_ceb$membership, sample=fam_ceb$names)) #get the clones from the network as a df
families <- merge(families, m2[,c("sample","lat","long","site","pop",'situ')], by="sample") #add some metadata
families <- families[order(as.numeric(families$family)),] #order the table by genet

write.xlsx(families, "DodoVisc/outputs/families.xlsx", rowNames=FALSE, colNames=TRUE)

##### SAMPLES SUMMARY ####

unfiltered_samples_summary <- d2$meta$analyses %>% 
  as.data.frame %>%
  group_by(pop, site) %>%
  summarise(
    n = n()
  )

filtered_samples_summary <- dms$meta$analyses %>% 
  as.data.frame %>%
  group_by(pop, site) %>%
  summarise(
    n = n()
  )

noclones_samples_summary <- dms_no_clones$meta$analyses %>% 
  as.data.frame %>%
  group_by(pop, site) %>%
  summarise(
    n = n()
  )

combined_sample_summary <- merge(unfiltered_samples_summary, filtered_samples_summary, by=c('pop', 'site'), all=TRUE)
combined_sample_summary2 <- merge(combined_sample_summary, noclones_samples_summary, by=c('pop', 'site'), all=TRUE)

colnames(combined_sample_summary2) <- c('Population','Site','Unfiltered QC (n)','Filtered QC (n)', 'Clones removed (n)')
View(combined_sample_summary2)

write.xlsx(combined_sample_summary2, paste0(species,'/outputs/DodoVisc_sample_summary.xlsx'))

##### KINSHIP HEATMAP ####

col_fun2 = colorRamp2(c(0,0.2,0.34,0.35355339), c("grey95", "red",'blue','black'))

kinship_df <- kin %>% as.data.frame()
kinship_df$sample <- rownames(kinship_df)

hm_sites2 <- merge(kinship_df, m2[,c("sample","site","pop", "lat", "long")],
                   by="sample", all.x=TRUE, all.y=FALSE)
hm_sites2 <- hm_sites2[match(rownames(kinship_df),hm_sites2$sample),]
rownames(hm_sites2) <- hm_sites2[,"sample"]

hm_sites <- hm_sites2
hm_sites2[,"sample"] <- NULL

site_ann <- HeatmapAnnotation(Site = hm_sites2$site,
                              col=list(Site=site_cols),
                              annotation_name_gp = gpar(fontsize = 10),
                              annotation_name_side="right")

pop_ann <- HeatmapAnnotation(Population = hm_sites2$pop, col = list(Population = pop_colours),
                             annotation_name_gp = gpar(fontsize = 10),
                             annotation_name_side="right")

site_ann2 <- rowAnnotation(Site = hm_sites2$site,
                           col=list(Site=site_cols),
                           annotation_name_gp = gpar(fontsize = 0))


pop_ann2 <- rowAnnotation(Population = hm_sites2$pop, col = list(Population = pop_colours),
                          annotation_name_gp = gpar(fontsize = 0))


hma <- Heatmap( as.matrix(hm_sites2[ , c(1:(nrow(hm_sites2)))]), 
                col=col_fun2, 
                bottom_annotation=c(site_ann, pop_ann),
                right_annotation=c(site_ann2, pop_ann2),
                name = "Kin", #title of legend
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_names_max_width = unit(15, "cm"),
                border_gp = gpar(col = "black", lty = 1),
                column_order=order(hm_sites2$pop),
                row_order=order(hm_sites2$pop)
)

draw(hma, merge_legend = TRUE)

# Set up the PNG device
png('DodoVisc/outputs/DodoVisc_kin.png', width = 18, height = 14, units = 'cm', res = 300)

# Draw the plot
draw(hma, merge_legend = TRUE)

# Turn off the PNG device
dev.off()


#### DIVERSITY STATS #####

stats_0.3 <- fastDiversity::faststats(dms_no_clones$gt, genetic_group_variable = dms_no_clones$meta$analyses[,'pop'],
                                      site_variable = dms_no_clones$meta$analyses[,'pop'], max_missingness = 0.3, maf=0.05)

write.xlsx(stats_0.3, paste0(species,'/outputs/DodoVisc_stats.xlsx'))

m2$ind_ho <- fastDiversity::individual_Ho(dms$gt, genetic_group_variable = dms$meta$analyses[,'sp'])

#### PCA ####
##### PCA OF ALL SAMPLES ####
# Filter by MAF and convert to genlight
gen_d5 <- new("genlight", dms_maf1$gt)
gen_pca <- glPca(gen_d5, parallel = TRUE, nf = 7)

# Extract PCA scores
g_pca_df <- gen_pca[["scores"]]

# Merge with metadata (fix m2 to d2$meta$analyses)
g_pca_df2 <- merge(g_pca_df, d2$meta$analyses, by.x = "row.names", by.y = "sample", all.y = FALSE, all.x = FALSE)

# Create axis labels with % variance
pcnames <- paste0(colnames(g_pca_df), " (", 
                  round(gen_pca[["eig"]][1:6] / sum(gen_pca[["eig"]]) * 100, 2), 
                  "%)")

pca1_hull <- g_pca_df2 %>% 
  filter(species %like% 'Grevillea DodoVisc')%>%
  group_by(pop) %>% 
  slice(chull(PC1, PC2)) 

# PCA Plot 1: PC1 vs PC2
pca_plot1 <- ggplot(g_pca_df2, aes(x = PC1, y = PC2, colour = pop, shape = pop)) +
  geom_vline(xintercept = 0, alpha = 0.2) + geom_hline(yintercept = 0, alpha = 0.2) +
  xlab(pcnames[1]) + ylab(pcnames[2]) +
  geom_point(size = 1.5) +
  theme_few() +
  labs(colour = "Species", shape = "Species") +
  theme(legend.key.size = unit(0, 'lines'), 
        legend.position = "right",
        legend.text = element_text(face = "italic"),
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8)) +
  guides(         color = guide_legend(title.position = "left", ncol=3),
                  shape = guide_legend(title.position = "left", ncol=3)) +
  geom_shape(data = pca1_hull, alpha = 00, expand = 0.03, radius = 0.03,size=0.3,color="black",
             aes(x=PC1, y=PC2, group=pop), linetype='dotted')+
  scale_colour_manual(values = species_colours) +
  scale_shape_manual(values = species_shapes)

print(pca_plot1)

# PCA Plot 2: PC3 vs PC4
pca_plot2 <- ggplot(g_pca_df2, aes(x = PC3, y = PC4, colour = pop, shape = pop)) +
  geom_vline(xintercept = 0, alpha = 0.2) + geom_hline(yintercept = 0, alpha = 0.2) +
  xlab(pcnames[3]) + ylab(pcnames[4]) +
  geom_point(size = 1.5) +
  theme_few() +
  labs(colour = "Species", shape = "Species") +
  theme(legend.key.size = unit(0, 'lines'), 
        legend.position = "right",
        legend.text = element_text(face = "italic"),
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8)) +
  guides(colour = guide_legend(title.position = "top")) +
  scale_colour_manual(values = species_colours) +
  scale_shape_manual(values = species_shapes)

print(pca_plot2)

# PCA Plot 3: PC5 vs PC6
pca_plot3 <- ggplot(g_pca_df2, aes(x = PC5, y = PC6, colour = pop, shape = pop)) +
  geom_vline(xintercept = 0, alpha = 0.2) + geom_hline(yintercept = 0, alpha = 0.2) +
  xlab(pcnames[5]) + ylab(pcnames[6]) +
  geom_point(size = 1.5) +
  theme_few() +
  labs(colour = "Species", shape = "Species") +
  theme(legend.key.size = unit(0, 'lines'), 
        legend.position = "right",
        legend.text = element_text(face = "italic"),
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8)) +
  guides(colour = guide_legend(title.position = "top")) +
  scale_colour_manual(values = species_colours) +
  scale_shape_manual(values = species_shapes)

print(pca_plot3)

all3_pca_plots <- ggarrange(pca_plot1, pca_plot2, pca_plot3, labels="AUTO",
                            font.label=list(face="plain"),
                            common.legend = TRUE, ncol=3, legend = "bottom")
all3_pca_plots

ggsave("DodoVisc/outputs/DodoVisc_pca.png",
       all3_pca_plots, width = 22, height = 10, units = "cm", dpi=300, bg='white')


##### PCA PER POPULATION ####
run_pca_plot <- function(pop_name, dms) {
  # Filter and clean
  dms_filtered <- dms %>%
    remove.by.list(dms$sample_names[dms$meta$analyses[,'pop'] == pop_name]) %>%
    remove.poor.quality.snps(min_repro = 0.96, max_missing = 0.3) %>%
    remove.fixed.snps() %>%
    remove.by.maf(0.05)%>% 
    remove.by.missingness(0.3)
  
  # Create genlight and run PCA
  gen_d5 <- new("genlight", dms_filtered$gt)
  gen_pca <- glPca(gen_d5, parallel = TRUE, nf = 7)
  
  # Merge PCA scores with metadata
  g_pca_df <- gen_pca[["scores"]]
  g_pca_df2 <- merge(g_pca_df, dms$meta$analyses, by.x = "row.names", by.y = "sample")
  
  # Axis labels
  pcnames <- paste0(colnames(g_pca_df), " (", 
                    round(gen_pca[["eig"]][1:6] / sum(gen_pca[["eig"]]) * 100, 2), 
                    "%)")
  
  # Plot
  ggplot(g_pca_df2, aes(x = PC1, y = PC2, colour = site, shape = site)) +
    geom_vline(xintercept = 0, alpha = 0.2) + 
    geom_hline(yintercept = 0, alpha = 0.2) +
    geom_point(size = 1.5) +
    xlab(pcnames[1]) + ylab(pcnames[2]) +
    theme_few() +
    labs(colour = "Site", shape = "Site", title=pop_name) +
    theme(legend.key.size = unit(0, 'lines'), 
          legend.position = "bottom",
          legend.direction = 'vertical',
          legend.text = element_text(face = "plain"),
          axis.title = element_text(size = 10), 
          axis.text = element_text(size = 8), 
          plot.title = element_text(size = 11))
}

# Example usage:
pca_enmore <- run_pca_plot("Enmore", dms_no_clones)
pca_torrington <- run_pca_plot("Torrington", dms_no_clones)
pca_guyfawkes <- run_pca_plot("Guy Fawkes", dms_no_clones)

# Then just call `pca_enmore` etc. to display each
ggarrange(pca_enmore, pca_torrington, pca_guyfawkes, ncol=3, labels='AUTO', font.label=list(face="plain"))


#### SPLITSTREE ####

# Prepare distance matrix
# splitstree(dist(dms_maf1$gt), paste0(species,'/outputs/DodoVisc.nex'))

# Read the network
Nnet <- phangorn::read.nexus.networx(paste0(species,'/outputs/DodoVisc.nex'))
# Prepare plotting data
x <- data.frame(x = Nnet$.plot$vertices[,1], 
                y = Nnet$.plot$vertices[,2], 
                sample = rep(NA, nrow(Nnet$.plot$vertices)))

x[Nnet$translate$node, "sample"] <- Nnet$translate$label

# Merge with metadata
x <- merge(x, m2, by.x = "sample", by.y = "sample", all.x = TRUE, all.y = FALSE)
x_tips <- x[!is.na(x$sample) & !is.na(x$site), ]

# Adjust axis limits
net_x_axis <- max(x_tips$x) - min(x_tips$x)
net_y_axis <- max(x_tips$y) - min(x_tips$y)

# Update tip labels to show site only
Nnet$translate$label <- gsub("\\\\n", "\n", x_tips[match(Nnet$tip.label, x_tips$sample), "pop_label"])

hull <- x %>% group_by(pop) %>% 
  slice(chull(x, y)) %>%
  filter(!is.na(pop))

# Create the plot with species-based colors and shapes
splitstree_plot <- ggplot(Nnet, aes(x = x, y = y)) +
  geom_shape(data = hull, alpha = 0.2, expand = 0.01, radius = 0.01,
             aes(fill = pop), color = 'transparent') +
  geom_splitnet(layout = "slanted", size = 0.2) +
  geom_point(data = x_tips, aes(x = x, y = y, colour = pop, shape = situ), size = 1) +
  scale_fill_manual(values = species_colours, na.translate = FALSE, guide = "none") +    
  scale_colour_manual(values = species_colours, na.translate = FALSE, guide = "none") +    
  scale_shape_manual(values = c(16, 3, 17), guide = guide_legend("Location")) +            
  geom_tiplab2(size = 3, hjust = -0.15, angle = 0, fontface = "italic", lineheight = 0.8) +
  theme_void() +
  expand_limits(x = c(min(x_tips$x) - 0.2 * net_x_axis, max(x_tips$x) + 0.3 * net_x_axis),
                y = c(min(x_tips$y) - 0.1 * net_y_axis, max(x_tips$y) + 0.1 * net_y_axis)) +
  theme(
    legend.text = element_text(face = "plain"),
    legend.key.size = unit(0, 'lines')
  ) +
  coord_fixed()

splitstree_plot

# Save the plot
ggsave(paste0(species,'/outputs/DodoVisc_splitstree.pdf'), 
       splitstree_plot, width=20, height=20, units='cm', dpi=300)

# Display the plot
pca_splits_plot <- ggarrange(all3_pca_plots, splitstree_plot, nrow=2, labels=c("","D"), font.label=list(face="plain"))

# Save the plot
ggsave(paste0(species,'/outputs/DodoVisc_pca_splitstree.png'), 
       pca_splits_plot, width=22, height=18, units='cm', dpi=300, bg='white')

#### FST ANALYSIS ####
# Remove samples with NA in site column (not species, since all are Grevillea DodoVisc)
# Filter sites with n < 2
sitepop_freq <- as.data.frame(table(dms_no_clones$meta$analyses[,'site']))
not_n1_sites <- as.vector(sitepop_freq[sitepop_freq$Freq < 3, 1])
keep_samples <- dms_no_clones$sample_names[which(dms_no_clones$meta$analyses[,'site'] %in% sitepop_freq[sitepop_freq$Freq >= 3, 1])]
fst_dms <- remove.by.list(dms_no_clones, keep_samples)

# Check filtered data
print("After filtering sites with n < 2:")
print(length(fst_dms$sample_names))
print(table(fst_dms$meta$analyses[,'site']))

# Calculate FST and geographic distance
gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pFst <- population.pw.Fst(fst_dms, fst_dms$meta$analyses[,'site'], RandRbase, species, dataset, maf_val=0.05, miss_val=0.3)
pS <- population.pw.spatial.dist(fst_dms, fst_dms$meta$analyses[,'site'])

# Make self comparisons NA
diag(pFst$Fst) <- NA
diag(pS$S) <- NA

##### FST MANTEL ####
# Mantel test
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 10000, na.rm = TRUE)
print(man)

# Mantel plot preparation
Fst_sig <- as.data.frame(as.table(pFst$Fst))
colnames(Fst_sig) <- c("Site1", "Site2", "Fst")
dist_long <- as.data.frame(as.table(pS$S))
Fst_sig$Geo_dist <- dist_long$Freq
Fst_sig$Geo_dist2 <- Fst_sig$Geo_dist / 1000  # Convert to km

# Aggregate metadata
meta_agg_df <- data.frame(lat = fst_dms$meta$lat, 
                          long = fst_dms$meta$long,
                          site = fst_dms$meta$analyses[,"site"],
                          pop = fst_dms$meta$analyses[,"pop"])
meta_agg <- meta_agg_df %>%
  group_by(pop, site) %>%
  summarise(lat = mean(lat, na.rm = TRUE),
            long = mean(long, na.rm = TRUE),
            .groups = 'drop')

# Add metadata for sites
Fst_sig2 <- merge(Fst_sig, distinct(meta_agg[, c("pop", "site")]), by.x = "Site1", by.y = "site", all.y = FALSE)
Fst_sig2 <- merge(Fst_sig2, distinct(meta_agg[, c("pop", "site")]), by.x = "Site2", by.y = "site", all.y = FALSE)
Fst_sig2$same_sp <- ifelse(Fst_sig2$pop.x == Fst_sig2$pop.y, "Within group", "Between group")

# Mantel plot
fstp1 <- ggplot(Fst_sig2, aes(x = Geo_dist2, y = Fst, color = same_sp)) +
  geom_point(size = 1, alpha = 0.3) +
  labs(x = "Distance (km)", y = "FST", colour = "Comparison") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom") +
  ggtitle("Mantel Test: Fst vs Geographic Distance",
          subtitle = paste("r =", round(man$statistic, 3), "p =", sprintf("%.3f", man$signif)))

# Display and save plot
print(fstp1)
ggsave("DodoVisc/outputs/DodoVisc_fst_mantel.png", 
       fstp1, width = 15, height = 10, units = "cm", dpi = 600)

##### FST HEATMAP ####
# Print Mantel result
print(paste("Mantel statistic r is", round(man$statistic, 3), ", P =", man$signif))

# Make heatmaps
# Geo distance
geo_d <- pS$S
mat_dist <- geo_d / 1000  # Convert to km

# FST
mat_fst <- pFst$Fst

# Order heatmaps consistently
order_hm <- Heatmap(mat_fst, cluster_rows = TRUE, cluster_columns = TRUE)
od <- colnames(mat_fst)[column_order(order_hm)]
mat_fst <- mat_fst[od, od]
mat_dist <- mat_dist[od, od]

mat_dist <- merge(mat_dist, meta_agg[,c('pop','site')], by.x=0, by.y="site", all.y=FALSE) #add aggregated df to mat2 (fst)
rownames(mat_dist) <- mat_dist$Row.names
mat_dist$Row.names <- NULL
mat_dist <- mat_dist[match(colnames(mat_dist)[1:nrow(mat_dist)],rownames(mat_dist)),]

# Specify heatmap colors
library(circlize)
gene_col <- colorRamp2(c(0, 0.5, 1), c("#8DD3C7", "white", "#FB8072"))
palette <- colorRamp2(c(0, max(geo_d/1000, na.rm = TRUE)), c("white", "#80B1D3"))


# Define species colors (replace named_list_maker if needed)
row_ann <- rowAnnotation(Group = mat_dist$pop,
                         col=list(Group=species_colours),  show_legend = FALSE,
                         na_col="white",
                         annotation_legend_param = list(labels_gp=gpar(fontface="italic",fontsize=8),
                                                        title_gp=gpar(fontsize=10)),
                         annotation_name_gp = gpar(fontsize = 0),
                         annotation_name_side="top")

col_ann <- HeatmapAnnotation(Group = mat_dist$pop, 
                             col = list(Group = species_colours),  show_legend = FALSE,
                             annotation_name_gp = gpar(fontsize = 0),
                             annotation_legend_param = list(labels_gp=gpar(fontface="italic", fontsize=8),
                                                            title_gp=gpar(fontsize=10)),
                             annotation_name_side="left",
                             na_col = "white")

# Geo heatmap
geo <- Heatmap(as.matrix(mat_dist[,1:nrow(mat_dist)]), 
               rect_gp = gpar(type = "none"),
               bottom_annotation = col_ann,
               width = nrow(mat_dist) * unit(6, "mm"),
               height = nrow(mat_dist) * unit(6, "mm"),
               col = palette,
               na_col = "white",
               column_names_gp = gpar(fontsize = 8),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               name = "Distance (km)",
               heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                           labels_gp = gpar(fontsize = 6)),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if (i > j) {
                   grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                   grid.text(sprintf("%.0f", mat_dist[i, j]), x, y, gp = gpar(fontsize = 6))
                 }
               })

geo

# FST heatmap
gene <- Heatmap(as.matrix(mat_fst), 
                right_annotation = c(row_ann),
                rect_gp = gpar(type = "none"),
                width = nrow(mat_fst) * unit(6, "mm"),
                height = nrow(mat_fst) * unit(6, "mm"),
                col = gene_col,
                na_col = "grey",
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 0),
                border_gp = gpar(col = "black", lty = 1),
                name = "FST",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                            labels_gp = gpar(fontsize = 5)),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if (i <= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    grid.text(sprintf("%.2f", mat_fst[i, j]), x, y, gp = gpar(fontsize = 6))
                  }
                })

gene
# Draw combined heatmap
gene_width <- nrow(mat_fst) * unit(6, "mm")
draw(geo + gene, ht_gap = -gene_width, merge_legend=TRUE)

# Set the file name and parameters
filename <- "DodoVisc/outputs/DodoVisc_fst.png"
width <- 8
height <- 6
dpi <- 600
units <- "cm"

# Set up the PNG device
png(filename, width = width, height = height, units = units, res = dpi)

# Draw the plot
draw(geo + gene, ht_gap = -gene_width)

# Turn off the PNG device
dev.off()


#### LEA ####
library(LEA)

# # Run once, takes time
# nd_lea <- dart2lea(dms_maf5, RandRbase, species, dataset)
# kvalrange <- 1:20
# snmf1 <- snmf(nd_lea, project='new', K=kvalrange, entropy = TRUE, repetitions = 5,  CPU=8)
# lea_samples<- dms_maf5$sample_names
# save(snmf1,lea_samples, file='DodoVisc/popgen/DodoVisc_snmf.RData')

load(file='DodoVisc/popgen/DodoVisc_snmf.RData')

K_chosen <- 10
best = which.min(cross.entropy(snmf1, K = K_chosen))
plot(snmf1, col = "blue", pch = 19, cex = 1.2)
entropy <- t(summary(snmf1)$crossEntropy) %>% as.data.frame()
entropy$k <- 1:20

entropy_plot <- ggplot(entropy, aes(x=k, y=mean))+
  geom_point(color='blue')+
  labs(y="Cross entropy")+
  theme_bw()
entropy_plot

ggsave('DodoVisc/outputs/paper/entropy_plot.png', entropy_plot, width = 12, height = 8,dpi=300, units = "cm")

qmatrix_df <- as_tibble(Q(snmf1, K = K_chosen, run=which.min(cross.entropy(snmf1, K = K_chosen)))) %>%
  mutate(sample = lea_samples) %>%
  pivot_longer(-sample, names_to = "lea_cluster", values_to = "proportion")

qmatrix_df2 <- merge(qmatrix_df, m2, by='sample')

qmatrix_df2 <- qmatrix_df2 %>%
  mutate(lea_cluster = gsub("V", "", lea_cluster))

# Order samples by latitude
qmatrix_df2 <- qmatrix_df2 %>%
  arrange(lat) %>%  # Arrange by latitude
  mutate(sample = factor(sample, levels = unique(sample)))  # Preserve order in factor


qmatrix_df2$pop <- factor(qmatrix_df2$pop)

lea_plot <- ggplot(qmatrix_df2, aes(x = sample, y = proportion, fill = as.factor(lea_cluster))) +
  geom_bar(stat = "identity", width = 1, alpha=0.8) +
  facet_grid(pop ~ ., scales = "free", space = "free") +  # Swap facet direction
  theme_few() +
  paletteer::scale_fill_paletteer_d("ggthemes::Classic_20") +  
  theme(
    legend.position = 'none',
    axis.title.x = element_text(size = 10, angle = 0),
    axis.text = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    plot.margin = margin(2, 2, 0, 2, unit = "pt"),
    strip.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0, face='italic', size=8) # Adjust facet labels
  ) +
  scale_y_continuous(limits = c(0, 1.001), expand = c(0, 0)) +
  labs(y = "Ancestry Proportion", fill = "") +
  coord_flip()  # Flips x and y axes

lea_plot

gt::info_paletteer(color_pkgs = "dichromat")
gt::info_paletteer(color_pkgs = "ggsci")
gt::info_paletteer(color_pkgs = "rcartocolor")
gt::info_paletteer(color_pkgs = "ggthemes")

lea_plot2 <- ggplot(qmatrix_df2, aes(x = sample, y = proportion, fill = as.factor(lea_cluster))) +
  geom_bar(stat = "identity", width = 1, alpha=0.8) +
  facet_grid(. ~pop, scales = "free", space = "free") +   # Swap facet direction
  theme_few() +
  paletteer::scale_fill_paletteer_d("ggthemes::Classic_20") +
  theme(
    legend.position = 'none',
    axis.title.y = element_text(size = 10, angle = 90),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    plot.margin = margin(2, 5, 2, 2, unit = "pt"),
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, face='italic', size=8) # Adjust facet labels
  ) +
  scale_y_continuous(limits = c(0, 1.001), expand = c(0, 0)) +
  labs(y = "Ancestry Proportion", fill = "")

lea_plot2

# #### COMBINED PLOT ####
# library(multipanelfigure)
# 
# combined_plots2 <- multi_panel_figure(
#   width = c(7,3.2,0.01,3.2,7),   # Adjust these dimensions as needed
#   height = c(9,10),
#   unit = "cm",
#   panel_label_type = "upper-roman"
# )
# 
# # Fill the panels with the respective plots
# combined_plots2 %<>%
#   fill_panel(splitstree_plot+theme(legend.position = 'top'), column = 1:2, row = 1, label = "A") %<>%
#   fill_panel(draw(geo + gene, ht_gap = -gene_width, merge_legend=TRUE), column = 4:5, row = 1, label = "B")%<>%
#   fill_panel(pca_enmore , column = 1, row = 2, label = "C") %<>%
#   fill_panel(pca_torrington, column = 2:4, row = 2, label = "D") %<>%
#   fill_panel(pca_guyfawkes, column = 5, row = 2, label = "E") 
# 
# combined_plots2
# 
# ggsave("DodoVisc/outputs/DodoVisc_combined.png",
#        combined_plots2, width = 23, height = 20, units = "cm", dpi = 600)

#### MAKE ADMIXTURE INPUT FILES ####
igds_file <- dart2gds(dms, paste0(out_dir_species,'/ADMIXTURE_'), species, dataset)

igds <- snpgdsOpen(igds_file)

snpgdsGDS2BED(igds, bed.fn=paste0(out_dir_species,'/ADMIXTURE_',species, "_", dataset,"_bed"))


fam <- data.frame(fam = 1:length(dms$sample_names), ind = dms$sample_names, 
                  fat = rep(0, length(dms$sample_names)), mot = rep(0, length(dms$sample_names)), 
                  sex = rep(0, length(dms$sample_names)),
                  pheno = rep(-9, length(dms$sample_names)), stringsAsFactors = FALSE)

write.table(fam, paste0(out_dir_species,'/ADMIXTURE_',species, "_", dataset,"_bed.fam"), sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

snpgdsClose(igds)

cd /home/eilish/projects/admixture/DodoVisc_May25

# #### LINUX COMMAND ADMIXTURE ####
# for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; \
# do admixture/dist/admixture_linux-1.3.0/admixture --cv=20 ADMIXTURE_*.bed $K -j20 | tee ADMIXTURE_log${K}.out; done
# grep -h CV ADMIXTURE_log*.out > ADMIXTURE_CV_results.tsv

#### STRUCTURE INPUT FILES ####
dms_structure <- remove.loci.randomly(dms, 1000)
dart2struct(dms_structure, RandRbase, species, dataset)
write.table(data.frame(index=1:length(dms_structure$sample_names), sample_names=dms_structure$sample_names),
            file=paste(RandRbase, species, "/popgen/", dms_structure$treatment,
                       "/struct/", species, "_", dataset, "_ids.txt",
                       sep = ""))

write_structure_mainparams(length(dms_structure$sample_names), length(dms_structure$locus_names), 'DodoVisc/popgen/raw_SNPFilt_1SNPperClone_SNPFilt/struct/mainparams')

#### RAXML INPUT FILE ####
df <- dms$gt

# Function to summarize combinations of values in columns
summarize_column_combinations <- function(df) {
  # Get the unique sets of values (including NA) for each column
  value_combinations <- apply(df, 2, function(col) {
    # Sort unique values (including NA) as a character string
    paste(sort(unique(col)), collapse = ",")
  })
  
  # Count occurrences of each unique combination
  table(value_combinations)
}

# Apply the function
comb_summary <- summarize_column_combinations(df)
View(comb_summary)


### loci that are 0,1 or 1,2 -- considered invariant by raxml

# Function to identify combinations and columns
identify_combinations <- function(df) {
  # Get the unique sets of values for each column
  combinations <- sapply(df, function(col) {
    paste(sort(unique(col)), collapse = ",")
  })
  
  # Group columns by their unique combinations
  split(names(combinations), combinations)
}

combination_groups <- apply(df, 2, function(x){sort(unique(x))})
# Apply the function
combination_groups <- identify_combinations(df)

out <- c()
for(i in 1:length(combination_groups)){
  v1 <- identical(combination_groups[[i]], c(0,1))
  v2 <- identical(combination_groups[[i]], c(1,2))
  
  if(isTRUE(v1)|isTRUE(v2))
  {out <- c(out, i)}
}

dms_raxml <- remove.snps.from.dart.data(dms, snps_to_remove = out, input_as_names = FALSE)

dart2svdquartets(dms_raxml, RandRbase, species, paste0(dataset,'_raxml'), add_pop=TRUE, pop=dms_raxml$sample_names)

# read what we just wrote
genotype_matrix <- read.nexus.data('DodoVisc/popgen/raw_SNPFilt_1SNPperClone/svdq/DodoVisc_raxml.nex')

genotype_matrix_ape <- ape::nexus2DNAbin(genotype_matrix)

ape::write.FASTA(genotype_matrix_ape, file='DodoVisc/popgen/raw_SNPFilt_1SNPperClone/svdq/DodoVisc_raxml.fasta')
