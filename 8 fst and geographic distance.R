# FST and geographic distance 
# FST is a measure of the similarity of allele frequencies between populations. If allele frequencies are identical (0)
# then populations are expected to have complete gene flow with ongoing migrations between populations. If allele
# frequencies are completely different (1) the populations have no migration/gene flow. 
# Often gene flow correlates with geographic distance, as it becomes harder for pollen or seeds to disperse between groups. 

# FST and distance calculations #####################################################################
library(SNPRelate)

gds_file <- dart2gds(dms, RandRbase, species, dataset)

pFst      <- population.pw.Fst(dms, dms$meta$site, RandRbase,species,dataset, maf_val=0.05, miss_val=1) #calculates genetic distance 

library(geosphere)
pS        <- population.pw.spatial.dist(dms, dms$meta$site) #calculates geographic distance between populations

####plot IBD plot

library(reshape2) #for melting data
library(vegan) #for mantel test

#tiff("E:/test/test fst plot.tiff", units="in", width=10, height=5, res=300)
Fst_sig <- cbind(melt(pS$S), unlist(as.list(pFst$Fst)))
colnames(Fst_sig)[3] <- "Geo_dist"
colnames(Fst_sig)[4] <- "Fst"
Fst_sig$Geo_dist2 <-Fst_sig$Geo_dist/1000 


# Mantel test  #####################################################################
# tests the correlation of the dist and fst matrices 

man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 999, na.rm = TRUE) #mantel test for IBD

ggplot(Fst_sig, aes(x= Geo_dist2, y=Fst))+geom_point()+labs(x="Distance (km)", y="Fst", title="Pairwise Fst plots")+
  annotation_custom(textGrob(paste("Mantel statistic r is", round(man$statistic, 3), ", P =", man$signif),
                             x=0.7,  y=0.1, gp=gpar(col="red", fontsize=8, fontface="italic")))+
  theme_bw()


# Basic FST and Geodist plots   #####################################################################

geo_d <-pS$S #this is a square matrix
geo_d <- geo_d/1000 # convert to km 
# geo_d[upper.tri(geo_d)] <- NA #makes the upper triangular part of the matrix into nothing
rownames(geo_d) <- colnames(pS$S) #make sure rownames are the same as colnames

sppop_freq <- as.data.frame(table(dms$meta$site))
n1_sp <- as.vector(sppop_freq[sppop_freq$Freq<=1,1]) #remove groups where n<=1

# geodist 
dimnames <- list (var1 = colnames(pS$S), var2 = colnames(pS$S)) 
mat <- matrix(geo_d, ncol=length(colnames(geo_d)), nrow=length(colnames(geo_d)), dimnames = dimnames)
mat <- mat[!(rownames(mat) %in% n1_sp),!(colnames(mat) %in% n1_sp)]

#FST
genetic_d <-pFst$Fst
rownames(genetic_d) <- colnames(pFst$Fst)

dimnames2 <- list (var1 = colnames(pFst$Fst), var2 = colnames(pFst$Fst))
mat2 <- matrix(genetic_d, ncol=length(colnames(geo_d)), nrow=length(colnames(geo_d)), dimnames = dimnames)
mat2 <- mat2[!(rownames(mat2) %in% n1_sp),!(colnames(mat2) %in% n1_sp)] # make fst matrix

agg <- unique(m2[, c("site","sp")]) # create aggregated df of species and site 

mat2 <- merge(mat2, agg, by.x=0, by.y="site", all.y=FALSE) #add aggregated df to mat2 (fst)
rownames(mat2) <- mat2$Row.names

mat2$Row.names <- NULL
mat2 <- mat2[match(colnames(mat2)[1:nrow(mat2)],rownames(mat2)),]

# make plots with custom functions
# geo plot variables are (matrix, axis text size)
geo <- geo_heat_function2(mat,10)

# fst plot variables are (matrix, geo plot, annotation1, annotation2, internal text size, axis text size)
gene <- gene_heat_function2(mat2[,1:nrow(mat2)], geo, NULL, NULL, 8,10) 

#combine
fst_plots2 <- geo+gene

#plot
draw(fst_plots2, row_order = row_order(geo))

# Advanced plots   #####################################################################
# this section is pretty finicky -- trying to get the annotations to be in the right order can be hard. 
# also ensure that the colnames of the distance plot match the fst plot. 

#make annotations
#can be anything, doesnt have to be species 
clust_ann2 <- HeatmapAnnotation(Species = mat2$sp,col=list(Species=sp_colours))
row_clust_ann2 <- rowAnnotation(Species = mat2$sp,col=list(Species=sp_colours),
                                show_legend=FALSE, annotation_label=" ")

# specify fst heatmap colours 
gene_col <-  colorRamp2(c(0,0.5,1), c("#8DD3C7", "white", "#FB8072"))
# make fst heatmap
gene <- Heatmap(mat2[,1:nrow(mat2)], bottom_annotation = clust_ann2, right_annotation = row_clust_ann2,
                col=gene_col,
                row_names_gp = gpar(fontsize = 11),
                column_names_gp = gpar(fontsize = 11),
                row_names_max_width = unit(15, "cm"),
                border_gp = gpar(col = "black", lty = 1),
                column_order=rownames(arrange(mat2, sp)),
                row_order=rownames(arrange(mat2,sp)),
                name="Pairwise Fst",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", mat2[,1:nrow(mat2)][i, j]), x, y, gp = gpar(fontsize = 9))})

#specify geo heatmap colours
palette <-  colorRamp2(c(0, max(mat)), c("white", "#80B1D3"))
# make geo heatmap 
geo <- Heatmap(mat, col=palette,
               row_names_gp = gpar(fontsize = 11),
               column_names_gp = gpar(fontsize = 11),
               row_names_max_width = unit(15, "cm"),
               border_gp = gpar(col = "black", lty = 1),
               name="Distance (km)",
               column_order=column_order(gene))

# combine plots
fst_plots <- geo+gene

#plot
draw(fst_plots , row_order = row_order(gene),merge_legend=TRUE) 
