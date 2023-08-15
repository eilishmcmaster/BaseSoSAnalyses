# FST and geographic distance 
# FST is a measure of the similarity of allele frequencies between populations. If allele frequencies are identical (0)
# then populations are expected to have complete gene flow with ongoing migrations between populations. If allele
# frequencies are completely different (1) the populations have no migration/gene flow. 
# Often gene flow correlates with geographic distance, as it becomes harder for pollen or seeds to disperse between groups. 

# FST and distance calculations #####################################################################
# remove sites where n=1
sppop_freq <- as.data.frame(table(dms$meta$site))
not_n1_sites <- as.vector(sppop_freq[sppop_freq$Freq<=1,1]) #remove groups where n<=1
not_n1_samples <- dms$sample_names[which(!(dms$meta$site %in% not_n1_sites))]
fst_dms <- remove.by.list(dms, not_n1_samples)

library(SNPRelate)
library(geosphere)

# calculate FST and geodist
gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pFst      <- population.pw.Fst(fst_dms, fst_dms$meta$site, RandRbase,species,dataset, maf_val=0.05, miss_val=0.3) #calculates genetic distance 
pS        <- population.pw.spatial.dist(fst_dms, fst_dms$meta$site) #calculates geographic distance between populations

####plot IBD plot

library(reshape2) #for melting data
library(vegan) #for mantel test

# Make self comparisons NA
diag(pFst$Fst) <- NA
diag(pS$S) <- NA

#Mantel test 
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 10000, na.rm = TRUE) #mantel test, finds if matrices are signficantly similar
man

# mantel plot
Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$Geo_dist2 <-Fst_sig$Geo_dist/1000 

# adding metadata for sites
Fst_sig2 <- merge(Fst_sig, distinct(m2[,c("site","sp")]), by.x="Var1", by.y="site", all.y=FALSE)
Fst_sig2 <- merge(Fst_sig2, distinct(m2[,c("site","sp")]), by.x="Var2", by.y="site", all.y=FALSE)
Fst_sig2$same_sp <- ifelse(Fst_sig2$sp.x == Fst_sig2$sp.y, "Intraspecific", "Interspecific")

library(ggforce)
fstp1 <- ggplot(Fst_sig2, aes(x= Geo_dist2, y=Fst, color=same_sp))+geom_point(size=1, alpha=0.3)+
  labs(x="Distance (km)", y="FST", colour="Comparison")+
  facet_zoom(x=Geo_dist2<25, zoom.size=1)+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")
fstp1

ggsave("BossFrag/outputs/paper/BossFrag_manning_fst.tiff",
       fstp1, width = 15, height = 15, units = "cm", dpi=600)

paste("Mantel statistic r is", round(man$statistic, 3), ", P =", man$signif)

# Make heatmaps
# geo dist
geo_d <-pS$S #this is a square matrix
mat <- geo_d/1000 # convert to km 

#FST
mat2 <-pFst$Fst
agg <- unique(m2[, c("site", "sp")]) # create aggregated df of species and site 
mat2 <- merge(mat2, agg, by.x=0, by.y="site", all.y=FALSE) #add aggregated df to mat2 (fst)
rownames(mat2) <- mat2$Row.names

mat2$Row.names <- NULL
mat2 <- mat2[match(colnames(mat2)[1:nrow(mat2)],rownames(mat2)),]

row_group_ann <- rowAnnotation(Group = mat2$sp,
                               col=list(Group=species_colours),
                               na_col="white",
                               annotation_legend_param = list(labels_gp=gpar(fontface="italic",fontsize=8),
                                                              title_gp=gpar(fontsize=10)),
                               annotation_name_gp = gpar(fontsize = 8),
                               annotation_name_side="top")

bottom_group_ann <- HeatmapAnnotation(Group = mat2$sp, col = list(Group = species_colours),
                                      annotation_name_gp = gpar(fontsize = 0),
                                      annotation_legend_param = list(labels_gp=gpar(fontface="italic", fontsize=8),
                                                                     title_gp=gpar(fontsize=10)),
                                      annotation_name_side="left",
                                      na_col = "white")

# specify fst heatmap colours 
gene_col <-  colorRamp2(c(0,0.5,1), c("#8DD3C7", "white", "#FB8072"))


#specify geo heatmap colours
palette <-  colorRamp2(c(0, max(mat)), c("white", "#80B1D3"))
# make geo heatmap 
geo <- Heatmap(mat, col=palette,na_col="white",
               bottom_annotation = bottom_group_ann,
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               row_names_max_width = unit(15, "cm"),
               border_gp = gpar(col = "black", lty = 1),
               name="Distance (km)",
               heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                           labels_gp = gpar(fontsize = 8)),
               cluster_rows = TRUE, 
               cluster_columns = TRUE
               # column_order=column_order(gene)
)

# make fst heatmap
gene <- Heatmap(mat2[,1:nrow(mat2)], right_annotation = row_group_ann,
                bottom_annotation = bottom_group_ann,
                col=gene_col,na_col="#8DD3C7",
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                row_names_max_width = unit(15, "cm"),
                border_gp = gpar(col = "black", lty = 1),
                # column_order=rownames(arrange(mat2, sp)),
                # row_order=rownames(arrange(mat2,sp)),
                column_order=column_order(geo),
                row_order=row_order(geo),
                name="FST",
                heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                            labels_gp = gpar(fontsize = 8)),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", mat2[,1:nrow(mat2)][i, j]), x, y, gp = gpar(fontsize = 4))})

# combine plots
fst_plots <- geo+gene
draw(fst_plots)


# Set the file name and parameters
filename <- "BossFrag/outputs/paper/fst_plot.tiff"
width <- 21
height <- 12
dpi <- 300
units <- "cm"

# Set up the PNG device
tiff(filename, width = width, height = height, units = units, res = dpi)

# Draw the plot
draw(fst_plots, row_order = row_order(gene), merge_legend = TRUE)

# Turn off the PNG device
dev.off()
