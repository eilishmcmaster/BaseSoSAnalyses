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

### Pairwise kinship boxplots

hm_x <- hm_sites[ , c(rownames(hm_sites), "sample")]
hm_x[lower.tri(hm_x, diag=TRUE)] <- NA
hm_x2 <- reshape2::melt(hm_x, na.rm=TRUE)
hm_x3 <- merge(hm_x2, m2[,c("sample","pop_large")], by.x="sample", by.y="sample", all.y=FALSE)
hm_x3 <- merge(hm_x3, m2[,c("sample","pop_large","genetic_group2")], by.x="variable", by.y="sample", all.y=FALSE)
hm_x3_filtered <- subset(hm_x3, (pop_large.x == pop_large.y))

kin_boxplot <- ggplot(hm_x3_filtered, aes(x=value, y=pop_large.x))+geom_boxplot(outlier.shape = NA)+
  theme_few()+geom_point(colour="blue", size=0.5)+
  scale_x_continuous(limits = c(0,0.5), expand=c(0,0))+
  geom_vline(xintercept = 0.45, colour="red", linetype="dashed")+
  geom_vline(xintercept = 0.25, colour="orange", linetype="dashed")+
  geom_vline(xintercept = 1/8, colour="pink", linetype="dashed")+
  geom_vline(xintercept = 1/16, colour="lightblue", linetype="dashed")+
  ylab(element_blank())+xlab("Pairwise kinship (k)")+
  facet_grid(genetic_group2~., scales = "free", drop=TRUE, space = "free")+
  theme(strip.text.y = element_text(angle = 0, face="italic"))+
  annotate(geom = "rect", xmin = 0.45, xmax = 0.5, ymin = -Inf, ymax = +Inf,alpha = 0.2, fill='red') # clones

kin_boxplot

# ggsave("PherFitz/outputs/plots/kin_boxplot.png", plot = kin_boxplot, width = 150, height = 150, dpi = 300, units = "mm")

# Table 1 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025716/


# Calculate mean, standard deviation, minimum, maximum, and count using aggregate()
summary_stats <- aggregate(value ~ pop_large.x+ genetic_group2, data = hm_x3_filtered, simplify=TRUE,
                           FUN = function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x)))%>%
  cbind(.[[ncol(.)]])

summary_stats$value <- NULL

# Display the summary statistics
print(summary_stats)
# write.xlsx(summary_stats, file="PherFitz/outputs/average_ibd_kin_per_site.xlsx", rowNames=FALSE)
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
