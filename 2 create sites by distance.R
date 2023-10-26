# 2 classify sites based on geographic distance 

# This method calculates pairwise distances between samples and then groups them by a specified distance. 
# This can then be added to metadata file. Useful for species where metadata is messy. 

library(geosphere)
# calculate pairwise distances in metres
pS <- population.pw.spatial.dist(dms2, dms2$sample_names)

# convert distance to km
distance_mat <-pS$S/1000

#group everything within X metres
distance_mat3 <- as.data.frame(distance_mat) %>%mutate_all(~replace(.,.>0.4, NA)) # example 400 metres (0.4 km)

# remove all connections that are more than 400 m
distance_mat3 <- as.data.frame(distance_mat3) %>%mutate_all(~replace(.,.>0, 1)) 

diag(distance_mat3) <- 1 #restore diagonal, not necessary


library(igraph)
# create network based on matrix 
dist_network <- graph_from_adjacency_matrix(as.matrix(distance_mat3), mode="undirected", diag=FALSE)

#create clusters from network
ceb <- cluster_fast_greedy(dist_network)

#view clusters
plot(ceb, dist_network, vertex.label.color="transparent", vertex.size=2, edge.width=0.4)

# get sites for all samples 
sites <-as.data.frame(cbind(sites=as.numeric(ceb$membership), sample=ceb$names))


# Export the data 
# sites_out <- merge(sites, m2, by="sample") #merge new sites with existing m2 metadata 
clipr::write_clip(sites) # write dataframe to clipboard
# could also write an xlsx 


## for when you dont have dms, just a dataframe:
require(geosphere)

S <- mat.or.vec(nrow(m2),nrow(m2) )
for (i in 1:nrow(m2)) {
  for (j in 1:nrow(m2)) {
    if (i > j) {
      LLi <- m2[i, c("long","lat")]
      LLj <- m2[j, c("long","lat")]
      Dij <- distCosine(as.numeric(LLi), as.numeric(LLj))
      S[i, j] <- Dij
      S[j, i] <- Dij
    }
  }
}
colnames(S) <- m2$sample
rownames(S) <- m2$sample

distance_mat <-S/1000

distance_mat2 <- ifelse(distance_mat <= 1, 1, 0)


dist_network <- igraph::graph_from_adjacency_matrix(as.matrix(distance_mat2), mode="undirected", diag=F)

ceb <- cluster_fast_greedy(dist_network)

ceb_net_plot <- plot(ceb, dist_network, vertex.label.color="transparent",vertex.size=2, edge.width=0.4)
ceb_net_plot
  
new_sites <-as.data.frame(cbind(sites=as.numeric(ceb$membership), sample=ceb$names))

