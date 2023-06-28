# find and remove clones using SNPrelate kinship 


#calculate kinship by population 
# VERY important that the population groups are true genetic groups and not conglomerates of multiple genetic groups
# this can be species, subpops, or sites 
kin <- individual_kinship_by_pop(dms2, RandRbase, species, dataset, dms2$meta$analyses[,"sp"], maf=0.05, mis=0.2, as_bigmat=TRUE)


# Finding the clones
#https://kateto.net/netscix2016.html
kin2 <- as.data.frame(kin) %>%mutate_all(~replace(.,.<0.45, 0)) #VERY IMPORTANT, removes all of the pairwise connections that are k<0.45 
diag(kin2) <- 0
kin3 <- kin2[rowSums(kin2)>0, rowSums(kin2)>0] #removes all of the individuals with no connections
diag(kin2) <- 0.5

kin3<- as.data.frame(kin3) %>%mutate_all(~replace(.,.>0, 1)) # replaces all remaining pairwise connections with 1 (all are equal above k>0.45)

# cluster the clones 
library(igraph)
network <- graph_from_adjacency_matrix(as.matrix(kin3), mode="undirected", diag=F,weighted=T) #makes the network based on k>0.45 
plot(network)

ceb <- cluster_fast_greedy(network) # makes cluster groupings

#visualise groups
plot(ceb, network, vertex.label.color="transparent", vertex.size=2, edge.width=0.4) #make the network plot

# get clones 
clones <-as.data.frame(cbind(genet=ceb$membership, sample=ceb$names)) #get the clones from the network as a df
clones_out <- merge(clones, m2[,c("sample","lat","long")], by="sample") #add some metadata
clones_out <- clones_out[order(as.numeric(clones_out$genet)),] #order the table by genet 


# Export data if you want 
# write.table(clones_out, paste0(species,"/outputs/clones_out.tsv"), sep="\t", row.names=FALSE, col.names=TRUE)
clones_out
# write.xlsx(clones_out, paste0(species,"/outputs/20220929_clones_out.xlsx"),
# asTable = FALSE, overwrite = TRUE)


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
dms <- remove.by.list(dms2, clones_dereplicated)

#this dataframe is used for some pca but takes a while to make
#better to only make if youre going to use it
# d5_no0 <- as.data.frame(dms[["gt"]]) %>% mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) 

