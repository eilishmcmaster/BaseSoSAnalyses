# migration on map 

################# Migration and FST ###################

make_genepop_file <- function(dms, maf,missing, group, grouping){

  dmsx <- dms %>%
    remove.poor.quality.snps(., min_repro=0.96,max_missing=missing) %>%
    remove.by.maf(., maf)

  ds <- dmsx$gt

  print(paste("Loci:", ncol(ds)))
  print(paste("Samples:",nrow(ds)))


  # if(ncol(ds) >=50){# make into genepop format
  old <- c("0","1","2", NA)
  new <- c("0101","0102","0202","0000")
  ds[ds %in% old] <- new[match(ds, old, nomatch = 0000)]

  # populations
  pops <- unique(grouping) # get population names

  # write genepop file
  gf <- paste0(species, "/popgen/genepop_",group,".gen") # make genepop file path

  cat(paste0("genepop file: ",species, " with MAF ", paste0(maf)), # first line of genepop file
      file=gf,sep="\n")
  cat(colnames(ds),file=gf,sep="\n", append=TRUE) # one loci name per line

  remove <- c() # vector for populations excluded from the analysis

  for (i in 1:length(pops)){ #loop for making the population groups
    if (length(which(grouping %in% pops[i]))<=1){ # find if the population is n=1
      cat("Removing population ", pops[i], " due to n=1")
      remove <- c(remove, pops[i]) # add the pop name to remove vector
    }else{
      cat("pop",file=gf,sep="\n", append=TRUE) # add the data to the genepop file
      df <- ds[which(grouping %in% pops[i]),]
      for (j in 1:nrow(df)){
        cat(c(paste0(pops[i],","),df[j,], "\n"),file=gf,sep="\t", append=TRUE)
      }
    }

  } #end of pops loop
  return(gf)
}

# here i am grouping by subpopulation (dms$meta$analyses[,"pop_large_short"]) but you can use dms$meta$site or sometyhing else
c <- make_genepop_file(dms, maf=0.05,missing=0.2, "pfitz", dms$meta$analyses[,"pop_large_short"])


#Migration was estimated using the divMigrate method (19) with Jost’s D metric of differentiation (47), as implemented in the R package diveRsity (46). This approach uses allele frequency differences between population pairs to estimate rates of migration in each direction; note that these rates are relative to other population pairs in the same data set and cannot be compared across data sets.

v <- diveRsity::divMigrate(infile=c, outfile=NULL, stat="gst",plot_network=TRUE, filter_threshold = 0.2, boots=1000, para=TRUE)

# Save a single object to a file because it takes a while to make
# saveRDS(v, "/Users/eilishmcmaster/Documents/PherFitz/PherFitz/outputs/plots/gst_m_10000_reps.rds")
# Restore it under a different name
# v <- readRDS("/Users/eilishmcmaster/Documents/PherFitz/PherFitz/outputs/plots/gst_m_10000_reps.rds")

d_mig <- v$gRelMig #v$dRelMig

# from is rows, to is columns
colnames(d_mig) <- unique(dms$meta$analyses[,"pop_large_short"]) 
rownames(d_mig) <- unique(dms$meta$analyses[,"pop_large_short"]) 

# Filter by significance -- significance is if there is a significant difference in the directions
# Test for overlap of the estimated 95% confidence intervals. Where there is no overlap, the directional gene flow components are said to be significantly different (asymmetric).
mig_sig <- v$gRelMigSig #v$dRelMigSig
colnames(mig_sig) <- unique(dms$meta$analyses[,"pop_large_short"])
rownames(mig_sig) <- unique(dms$meta$analyses[,"pop_large_short"])

d_mig[mig_sig>0.01] <- NA

# qgraph::qgraph(d_mig,legend = TRUE, edge.labels = TRUE, 
# curve = 2.5, mar = c(2, 2, 5, 5))

long_mig <- melt(d_mig)

meta_agg <- m2 %>%
  group_by(pop_large_short,pop_large, genetic_group) %>%
  summarize(lat = mean(lat, na.rm=TRUE),
            long = mean(long,na.rm=TRUE),
            .groups = 'drop')%>%
  subset(.,pop_large_short!="Ex_situ_PF")

long_mig <- merge(long_mig, distinct(meta_agg[, c("pop_large_short","pop_large", "genetic_group","lat","long")]), by.x = "Var1", by.y = "pop_large_short", all.y = FALSE)
long_mig <- merge(long_mig, distinct(meta_agg[, c("pop_large_short","pop_large", "genetic_group","lat","long")]), by.x = "Var2", by.y = "pop_large_short", all.y = FALSE)
long_mig<- distinct(long_mig)  

colnames(long_mig)[1:3] <- c("from", "to", "m")

long_mig <- long_mig[long_mig$from!="Ex_situ_PF"&long_mig$to!="Ex_situ_PF"&!is.na(long_mig$m),]
pop_Fst_sig2 <- pop_Fst_sig2[order(pop_Fst_sig2$Fst),]

# Common theme settings
common_theme <- theme(
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),  # Added legend text size
  strip.text.x = element_text(size = 8)
)


# Plot 2: long_mig
long_mig <- long_mig[order(long_mig$m), ]

gst_no_map <- ggplot() + coord_cartesian() + coord_fixed() +
  theme_few() +
  geom_curve(data = long_mig[long_mig$m > 0.2, ], #
             aes(x = long.x, y = lat.x,
                 xend = long.y, yend = lat.y, colour = m), 
             size = 0.5, na.rm = TRUE, curvature = 0.3, 
             arrow = arrow(angle = 20, ends = "first", type = "open", length = unit(2, "mm"))) +
  scale_color_gradient(low = "white", high = "red") + # midpoint = 0.5, mid = "blue",
  geom_point(data = meta_agg, mapping = aes(x = long, y = lat), colour = "black") +
  xlim(lims[1], lims[2]) +
  ylim(lims[3], lims[4]) +
  labs(x = "Longitude", y = "Latitude", colour = "Migration (m)") + guides(size = "none", alpha = "none") +
  ggrepel::geom_label_repel(data = meta_all_fitz, aes(x = long, y = lat, label = pop_large_short),
                            min.segment.length = 0.25, color = "black", fill = "white", size = 3,
                            segment.colour = "white", alpha = 0.9, label.size = 0, nudge_y = 0.003) +
  theme(legend.position = "bottom") 


gst_no_map