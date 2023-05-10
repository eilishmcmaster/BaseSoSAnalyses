# make a kinship heatmap 

# Original method with SNPrelate package###################################################################

# claculate kinship -- also done in `3 find clones.R`, does not need to be repeated
# # VERY important that the population groups are true genetic groups and not conglomerates of multiple genetic groups
# # this can be species, subpops, or sites 
# kin <- individual_kinship_by_pop(dms2, RandRbase, species, dataset, dms2$meta$analyses[,"sp"], maf=0.1, mis=0.2, as_bigmat=TRUE)

# NEW method with dist  ###################################################################
# i prefer dist over popkin for initial anaalyses as it is fast and simple

#### Kinship by distance ####

library(heatmaply)
library(circlize)
library(ComplexHeatmap)
library("tidyr")

col_fun2 = colorRamp2(c(0,0.5,1), c("white", "red","black"))


kinship_df <- dist_kinship_matrix(dms$gt) %>% as.data.frame(.)
# kin <- as.matrix(dist(dms$gt, diag=TRUE))
# kinship_df <- 1- (kin/max(kin)) %>% as.data.frame(.)

kinship_df$sample <- rownames(kinship_df)

hm_sites2 <- merge(kinship_df, m2[,c("sample","site","sp", "lat", "long","site_number", "genetic_group")],
                   by="sample", all.x=TRUE, all.y=FALSE)
hm_sites2 <- hm_sites2[match(rownames(kinship_df),hm_sites2$sample),]
rownames(hm_sites2) <- hm_sites2[,"sample"]


hm_sites2[,"sample"] <- NULL


#create annotations

site_ann <- HeatmapAnnotation(Location = hm_sites2$site,
                              col=list(Location=site_colours))

site_numeric_ann <- HeatmapAnnotation(Location = hm_sites2$site_number,
                                      col=list(Location=site_numeric_colours),
                                      annotation_name_gp = gpar(fontsize = 10),
                                      annotation_name_side="left")

sp_ann <- HeatmapAnnotation(Species = hm_sites2$sp,
                            col=list(Species=species_colours),
                            na_col="white",
                            annotation_legend_param = list(labels_gp=gpar(fontface="italic")),
                            annotation_name_gp = gpar(fontsize = 10),
                            annotation_name_side="left")


group_ann2 <- HeatmapAnnotation(Group = hm_sites2$genetic_group, col = list(Group = group_colours),
                                annotation_name_gp = gpar(fontsize = 10),
                                annotation_name_side="left",
                                na_col = "white")

lat_ann <- rowAnnotation(Latitude = hm_sites2$lat, col = list(Latitude = colorRamp2(c(min(hm_sites2$lat), mean(hm_sites2$lat),
                                                                                      max(hm_sites2$lat)), c("blue", "white", "red"))))

long_ann <- rowAnnotation(Longitude = hm_sites2$long, 
                          col = list(Longitude = colorRamp2(c(min(hm_sites2$long),
                                                              mean(hm_sites2$long),
                                                              max(hm_sites2$long)), c("blue", "white", "red"))),
                          annotation_name_gp = gpar(fontsize = 10),
                          annotation_name_side="top")


hma <- Heatmap( as.matrix(hm_sites2[ , c(1:(nrow(hm_sites2)))]), 
                col=col_fun2, 
                bottom_annotation=c(sp_ann, group_ann2),
                right_annotation = c(long_ann),
                name = "Genetic similarity", #title of legend
                row_names_gp = gpar(fontsize = 0),
                column_names_gp = gpar(fontsize = 0),
                row_names_max_width = unit(15, "cm"),
                border_gp = gpar(col = "black", lty = 1),
                # column_order=order(hm_sites2$lat),
                # row_order=order(hm_sites2$lat)
)


draw(hma, merge_legend = TRUE)
# 
# # Popkin method ###########################################################################
# library(popkin)
# library(BEDMatrix)
# library(openxlsx)
# library(lfa)
# 
# dms_pv1 <- remove.by.list(dms, m2[(m2$sp2 %in% "P. linifolia"),] %>%.$sample) 
# 
# 
# X <- t(dms_pv1$gt)
# subpops <- dms_pv1$sample_names
# # subpops <- dms_pv1$meta$analyses[,"sp2"]
# kinship <- popkin(X, subpops) # kinship matrix, can be used same as kin in heatmap builder
# plot_popkin(
#   inbr_diag(kinship),
#   labs = NULL,
#   # shared bottom and left margin value, to make space for labels
#   mar = 1
# )
# 
# # compute pairwise FST matrix from kinship matrix
# pairwise_fst <- pwfst(kinship)
# # fancy legend label
# leg_title <- expression(paste('Pairwise ', F[ST]))
# # NOTE no need for inbr_diag() here!
# plot_popkin(
#   pairwise_fst,
#   labs = subpops
# )
# 
