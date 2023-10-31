# FST and geographic distance 
# FST is a measure of the similarity of allele frequencies between populations. If allele frequencies are identical (0)
# then populations are expected to have complete gene flow with ongoing migrations between populations. If allele
# frequencies are completely different (1) the populations have no migration/gene flow. 
# Often gene flow correlates with geographic distance, as it becomes harder for pollen or seeds to disperse between groups. 

# FST and distance calculations #####################################################################
# remove sites where n=1
insitu_samples <- m2[(m2$sp %in% c("Pherosphaera fitzgeraldii"))& 
                       !(m2$pop_large %in% c("Mt Tomah Ex-situ","ANBG Ex-Situ")), ] %>% .$sample
dms_pf_insitu <- remove.by.list(dms_pf, insitu_samples)

sppop_freq <- as.data.frame(table(dms_pf_insitu$meta$site))
not_n1_sites <- as.vector(sppop_freq[sppop_freq$Freq<=1,1]) #remove groups where n<=1
not_n1_samples <- dms_pf_insitu$sample_names[which(!(dms_pf_insitu$meta$site %in% not_n1_sites))]
fst_dms <- remove.by.list(dms_pf_insitu, not_n1_samples)


# calculate FST and geodist
gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pFst      <- population.pw.Fst(fst_dms, fst_dms$meta$site, RandRbase,species,dataset, maf_val=0.05, miss_val=0.2) #calculates genetic distance 
pS        <- population.pw.spatial.dist(fst_dms, fst_dms$meta$site) #calculates geographic distance between populations


####plot IBD plot

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

meta_agg <- m2 %>%
  group_by(genetic_group, site) %>%
  summarize(lat = mean(lat, na.rm=TRUE),
            long = mean(long,na.rm=TRUE),
            .groups = 'drop')

# adding metadata for sites
Fst_sig2 <- merge(Fst_sig, distinct(meta_agg[,c("site","genetic_group")]), by.x="Var1", by.y="site", all.y=FALSE)
Fst_sig2 <- merge(Fst_sig2, distinct(meta_agg[,c("site","genetic_group")]), by.x="Var2", by.y="site", all.y=FALSE)
Fst_sig2$same_sp <- ifelse(Fst_sig2$genetic_group.x == Fst_sig2$genetic_group.y, "Within group", "Between group")

fstp1<- ggplot(Fst_sig2, aes(x= Geo_dist2, y=Fst, color=same_sp))+geom_point(size=1, alpha=0.3)+
  labs(x="Distance (km)", y="FST", colour="Comparison")+
  facet_zoom(x=Geo_dist2<2, zoom.size=1)+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom")
fstp1

ggsave("PherFitz/outputs/plots/PherFitz_manning_fst.png",
       fstp1, width = 15, height = 15, units = "cm", dpi=600)

paste("Mantel statistic r is", round(man$statistic, 3), ", P =", man$signif)


# Make heatmaps
# geo dist
geo_d <-pS$S #this is a square matrix
mat <- geo_d/1000 # convert to km 

#FST
mat2 <-pFst$Fst
agg <- unique(m2[, c("site", "pop_large")]) # create aggregated df of pop_largeecies and site
mat2 <- merge(mat2, agg, by.x=0, by.y="site", all.y=FALSE) #add aggregated df to mat2 (fst)
rownames(mat2) <- mat2$Row.names

mat2$Row.names <- NULL
mat2 <- mat2[match(colnames(mat2)[1:nrow(mat2)],rownames(mat2)),]

order_hm <- Heatmap(mat,
                    cluster_rows = TRUE,
                    cluster_columns = TRUE)
od <- colnames(mat)[column_order(order_hm)]

mat = mat[od, od]
mat2 = mat2[od, c(od,"pop_large")]

# specify fst heatmap colours 
gene_col <-  colorRamp2(c(0,0.5,1), c("#8DD3C7", "white", "#FB8072"))


#specify geo heatmap colours
palette <-  colorRamp2(c(0, max(mat, na.rm=TRUE)), c("white", "#80B1D3"))


row_Subpopulation_ann <- rowAnnotation(Subpopulation = mat2$pop_large,
                                       col=list(Subpopulation=pop_colours),
                                       na_col="white",
                                       annotation_legend_param = list(labels_gp=gpar(fontsize=6),#fontface="italic",
                                                                      title_gp=gpar(fontsize=8)),
                                       annotation_name_gp = gpar(fontsize = 0),
                                       annotation_name_side="top")

bottom_Subpopulation_ann <- HeatmapAnnotation(Subpopulation = mat2$pop_large, col = list(Subpopulation = pop_colours),
                                              annotation_name_gp = gpar(fontsize = 0),
                                              show_legend = FALSE,
                                              annotation_name_side="right",
                                              na_col = "white")

geo <- Heatmap(mat,rect_gp = gpar(type = "none"),
               width = nrow(mat)*unit(4, "mm"),
               height = nrow(mat)*unit(4, "mm"),
               col=palette,na_col="white",
               bottom_annotation = bottom_Subpopulation_ann,
               column_names_gp = gpar(fontsize = 6),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               name="Distance (km)",
               heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                           labels_gp = gpar(fontsize = 6)),
               # cluster_rows = TRUE, 
               # cluster_columns = TRUE,
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(i >= j) {
                   grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                 }
               }
)

# make fst heatmap
gene <- Heatmap(as.matrix(mat2[,1:nrow(mat2)]), rect_gp = gpar(type = "none"),
                width = nrow(mat2)*unit(4, "mm"),
                height = nrow(mat2)*unit(4, "mm"),
                right_annotation = row_Subpopulation_ann,
                col=gene_col,na_col="grey",
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 0),
                border_gp = gpar(col = "black", lty = 1),
                name="FST",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                            labels_gp = gpar(fontsize = 6)),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(i <= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    grid.text(sprintf("%.2f", mat2[,1:nrow(mat2)][i, j]), x, y, gp = gpar(fontsize = 4))
                  }
                })

gene_width <- nrow(mat2)*unit(4, "mm")

draw(geo + gene, ht_gap = -gene_width)

# Set the file name and parameters
filename <- "/Users/eilishmcmaster/Documents/PherFitz/PherFitz/outputs/plots/fst_plot_overlay.png"
width <- 16
height <- 12
dpi <- 600
units <- "cm"

# Set up the PNG device
png(filename, width = width, height = height, units = units, res = dpi)

# Draw the plot
draw(geo + gene, ht_gap = -gene_width)

# Turn off the PNG device
dev.off()
