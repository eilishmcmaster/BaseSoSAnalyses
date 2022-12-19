# make a kinship heatmap 

# Original method with SNPrelate package###################################################################

# claculate kinship -- also done in `3 find clones.R`, does not need to be repeated
# # VERY important that the population groups are true genetic groups and not conglomerates of multiple genetic groups
# # this can be species, subpops, or sites 
# kin <- individual_kinship_by_pop(dms2, RandRbase, species, dataset, dms2$meta$analyses[,"sp"], maf=0.1, mis=0.2, as_bigmat=TRUE)

# NEW method with popkin package ###################################################################




# Make figure ###################################################################

library(heatmaply)
library(circlize)
library(ComplexHeatmap)
library(tidyr)

col_fun2 = colorRamp2(c(0,0.25,0.45), c("white", "red","black")) # make colour palette to fill heatmap

clonekin <- as.data.frame(kin) # make it a separate df to mess around with 
# clonekin <- clonekin[rownames(clonekin) %in% clones_dereplicated, colnames(clonekin) %in% clones_dereplicated]
clonekin$sample <- rownames(clonekin)

# add metadata, can be different variables to what I used 
hm_sites2 <- merge(clonekin,clones_out[,c("sample","genet")], by="sample", all.x=TRUE, all.y=FALSE)
hm_sites2 <- merge(hm_sites2, m2[,c("sample","site", "sp")],
                   by="sample", all.x=TRUE, all.y=FALSE)

# reorder rows to match columns IMPORTANT 
hm_sites2 <- hm_sites2[match(rownames(clonekin),hm_sites2$sample),]

# move sample id back to rownames
rownames(hm_sites2) <- hm_sites2[,"sample"]
hm_sites2[,"sample"] <- NULL

# replace blanks with NA string
hm_sites2$genet[hm_sites2$genet==""]<-"NA"

# add 0 to the start of single digit genet ID numbers (5 -> 05)
hm_sites2$genet <- sprintf("%02d", as.numeric(hm_sites2$genet))

# make colour range for genet IDs 
clone_colours <- named_list_maker(hm_sites2[order(as.numeric(hm_sites2$genet)),"genet"], "YlGnBu", 8)
clone_colours["NA"] <- "white"

#create annotations -- go on the bottom of the heatmap
site_ann <- HeatmapAnnotation(Site = hm_sites2$site,
                             col=list(Site=site_colours))

sp_ann <- HeatmapAnnotation(Species = hm_sites2$sp,
                            col=list(Species=sp_colours))

# clust_ann <- HeatmapAnnotation(Cluster = hm_sites2$cluster,col=list(Cluster=cluster_colours))
clone_ann <- HeatmapAnnotation(Genet = hm_sites2$genet,col=list(Genet=clone_colours))

# make heatmap 
hma <- Heatmap( as.matrix(hm_sites2[ , c(1:(nrow(hm_sites2)))]), # specify kinship part of dataframe
                col=col_fun2, # colours to fill heatmap
                bottom_annotation=c(clone_ann), # which annotations to plot 
                name = "Kinship", #title of legend
                row_names_gp = gpar(fontsize = 4), #sample ID fontsize
                column_names_gp = gpar(fontsize = 4),
                row_names_max_width = unit(15, "cm"),
                border_gp = gpar(col = "black", lty = 1),
                column_order=rownames(arrange(hm_sites2, site)), #order of samples, can be automatically grouped if this is removed
                row_order=rownames(arrange(hm_sites2,site))
)

draw(hma, merge_legend = TRUE) #draw plot 