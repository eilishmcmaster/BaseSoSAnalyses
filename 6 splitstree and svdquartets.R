# trees of various kinds

# Splitstree ###################################################################
# exporting data to be viewed in splitstree desktop gui 
library(RSplitsTree)

d5_splitstree <- as.matrix(dms$gt) # dist object
rownames(d5_splitstree) <- paste(rownames(dms$gt), dms$meta$analyses[,"sp"], dms$meta$analyses[,"genetic_group"], sep=" ") # make leaf names include species
splitstree(dist(d5_splitstree), paste0(species,'/outputs/all_pultenaea_nexus_file.nex'))


# Visualise splitstre in R #######################################################
library(tanggle)
library(RSplitsTree)

splitstree(dist(dms$gt), 'PherFitz/outputs/all_pher_nexus_file_for_R.nex')

#need to open and save the file in Splitstree app for it to open here, IDK why
Nnet <- phangorn::read.nexus.networx('PherFitz/outputs/all_pher_nexus_file_for_R.nex')
# mat2[match(colnames(mat2)[1:nrow(mat2)],rownames(mat2)),]

x <- data.frame(x=Nnet$.plot$vertices[,1], y=Nnet$.plot$vertices[,2], 
                sample=rep(NA, nrow(Nnet$.plot$vertices)))


x[Nnet$translate$node,"sample"] <- Nnet$translate$label
x <- merge(x, m2, by="sample", all.x=TRUE, all.y=FALSE)

net_x_axis <- max(x$x)-min(x$x)
net_y_axis <- max(x$y)-min(x$y)

Nnet$translate$label <-  x[match(Nnet$tip.label, x$sample), "site"] %>% .[1:length(Nnet$tip.label)]


splitstree_plot <- ggplot(Nnet, mapping = aes_(~x, ~y), layout = "slanted", mrsd = NULL, 
                          as.Date = FALSE, yscale = "none", yscale_mapping = NULL, 
                          ladderize = FALSE, right = FALSE, branch.length = "branch.length", 
                          ndigits = NULL)+
  geom_splitnet(layout = "slanted", size=0.2)+
  geom_point(data=x, aes(x, y, colour=genetic_group2))+
  scale_colour_manual(values=genetic_group2_colours, na.translate=FALSE,
                      guide = guide_legend("Genetic group"))+
  # geom_tiplab2(size=2, hjust=-0.2)+
  theme_void()+
  expand_limits(x=c(min(x$x)-0.05*net_x_axis, max(x$x)+0.05*net_x_axis),
                y=c(min(x$y)-0.05*net_y_axis, max(x$y)+0.05*net_y_axis))+
  theme(legend.text = element_text(face="italic"), legend.position = "top")+coord_fixed()

splitstree_plot

# SVDquartets ###################################################################
# export data to calculate ML trees in SVDquartet
# output must be uploaded to server to run SVDquartets in PAUP
library(dartR)
# make genlight file -- uses the original data before clones removed BEWARE
bab <- gl.read.dart("/Users/eilishmcmaster/Documents/ZierObco/ZierObco/dart_raw/Report_DZ22-7321_SNP_mapping_2.csv",
                    ind.metafile="/Users/eilishmcmaster/Documents/ZierObco/ZierObco/meta/ZierObco_DZ22-7321_meta.csv",
                    mono.rm = TRUE)
# make nexus file for SVDquartets
gl2svdquartets(bab, outfile="20221213_zieria.nex",
               outpath="/Users/eilishmcmaster/Documents/ZierObco/ZierObco/outputs/")

# to run on server 
# eilish@rcn02:~/projects/zieria_phylo_Dec22$ paup all_zieria_nexus_file.nex -L log.log
