# trees of various kinds

# Splitstree ###################################################################
# exporting data to be viewed in splitstree desktop gui 
library(RSplitsTree)
# d5_no0 <- as.data.frame(dms[["gt"]]) %>% mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) # calculated previously
d5_splitstree <-d5_no0
rownames(d5_splitstree) <- paste(rownames(d5_no0), dms$meta$analyses[,"sp"], sep=" ") # make leaf names include species
splitstree(dist(d5_splitstree), paste0(species,'/outputs/all_zieria_nexus_file.nex'))

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