# dimensionality reduction methods 
# any groups identified here should inform the groupings used for kinship analysis, so this process is not linear.

# PCA reduces dimensions by consolidating variables (driving variable focus)
# PCA calculations
gen_d5 <- new("genlight", dms[["gt"]]) #convert df to genlight object for glPca function
gen_pca <- glPca(gen_d5, parallel=TRUE, nf=5) #do pca -- this method allows the input to have NAs 
g_pca_df <- gen_pca[["scores"]] #extract PCs 
g_pca_df2 <- merge(g_pca_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=FALSE) # add metadata 

pcnames <- paste0(colnames(g_pca_df)," (",
                  paste(round(gen_pca[["eig"]][1:5]/sum(gen_pca[["eig"]]) *100, 2)),
                  "%)") #create names for axes

# PCA plots 
sp_colours <- named_list_maker(dms$meta$analyses[,"sp"], "Spectral", 9)

pca_plot <- ggplot(g_pca_df2, aes(x=PC1, y=PC2, colour=sp))+ 
  geom_point()+theme_bw()+
  labs(colour="Species")+xlab(pcnames[1])+ylab(pcnames[2])+
  theme(legend.key.size = unit(0, 'lines'), legend.position = "right")+
  guides(colour = guide_legend(title.position = "top", direction = "vertical"))+
  scale_colour_manual(values=sp_colours)
pca_plot


# custom functions to plot with lat long gradients
lat_p <- pca_grad_funct(g_pca_df2, PC1, PC2,1,2, lat,"Latitude", "yellow", "red")
long_p <- pca_grad_funct(g_pca_df2, PC1, PC2,1,2, long,"Longitude", "lightblue", "darkblue")
lat_long <- ggarrange(lat_p, long_p)
lat_long


# UMAP reduces dimensions by retaining relationships between points (sample focussed)
# UMAP calculations
library(umap)
custom.config = umap.defaults
custom.config$random_state = 666666

#run umap
umer <- umap(d5_no0, config=custom.config) # run umap

umap_df <- umer$layout %>% as.data.frame() #extract output vectors
umap_df2 <- merge(umap_df, m2, by.x=0, by.y="sample", all.y=FALSE, all.x=TRUE) #add metadata

# UMAP plots 
umap_plot <- ggplot(umap_df2, aes(x=V1, y=V2, colour=factor(sp)))+
  geom_point()+theme_bw()+labs(colour="Species")+
  scale_colour_manual(values=sp_colours)

umap_plot

# make groups based on UMAP -- works for very clear groupings
# this script can also be used on PCA output dataframe if there are distinct groups
library("dbscan")

#clustering with dbscan
cl <- hdbscan(as.matrix(umer$layout), minPts = 5) # makes HDB object (list of lissts )
clustered <- cbind(umer$layout, cluster=(cl$cluster)) # extracts df of umap V1, V2, and cluster -- sample rownames
clustered2 <- merge(umap_df2,clustered, by=c("V1","V2"), all=TRUE) # matches cluster to metadata by UMAP position V1 and V2
if(min(clustered2$cluster)==0){
  clustered2$cluster <- clustered2$cluster+1
}
#export 
#would recommend adding to meta file so you don't have to run every time 
dm <- merge(dm, clustered2[,c('Row.names','cluster')], by.x='sample', by.y='Row.names', all.x=TRUE) # add cluster to