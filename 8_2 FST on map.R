# calculate fst 
# calculate FST and geodist
pop_gds_file <- dart2gds(fst_dms, RandRbase, species, dataset)
pop_pFst <- population.pw.Fst(fst_dms, fst_dms$meta$analyses[,'pop_large'], RandRbase, species, dataset, maf_val = 0.05, miss_val = 0.2)
pop_pS <- population.pw.spatial.dist(fst_dms, fst_dms$meta$analyses[,'pop_large'])

meta_agg <- m2 %>%
  group_by(pop_large, genetic_group) %>%
  summarize(lat = mean(lat, na.rm=TRUE),
            long = mean(long,na.rm=TRUE),
            .groups = 'drop')


####plot IBD plot

# Make self comparisons NA
diag(pop_pFst$Fst) <- NA
diag(pop_pS$S) <- NA

# Mantel test
pop_man <- mantel(xdis = pop_pS$S, ydis = pop_pFst$Fst, permutations = 10000, na.rm = TRUE)
pop_man

# mantel plot
pop_Fst_sig <- cbind(melt(pop_pS$S), unlist(as.list(pop_pFst$Fst)))
colnames(pop_Fst_sig)[3] <- "Geo_dist"
colnames(pop_Fst_sig)[4] <- "Fst"
pop_Fst_sig$Geo_dist2 <- pop_Fst_sig$Geo_dist / 1000
pop_Fst_sig <- pop_Fst_sig[!is.na(pop_Fst_sig$Geo_dist),]

# adding metadata for pop_larges
pop_Fst_sig2 <- merge(pop_Fst_sig, distinct(meta_agg[, c("pop_large", "genetic_group","lat","long")]), by.x = "Var1", by.y = "pop_large", all.y = FALSE)
pop_Fst_sig2 <- merge(pop_Fst_sig2, distinct(meta_agg[, c("pop_large", "genetic_group","lat","long")]), by.x = "Var2", by.y = "pop_large", all.y = FALSE)
pop_Fst_sig2$same_sp <- ifelse(pop_Fst_sig2$genetic_group.x == pop_Fst_sig2$genetic_group.y, "Within group", "Between group")
pop_Fst_sig2<- distinct(pop_Fst_sig2)  




# removed samples i dont want
meta_all_fitz1 <- subset(m2, sp == "Pherosphaera fitzgeraldii" & !(site %like% "Ex_situ"))  # all collected fitz

# get average lat long per population
meta_all_fitz <- meta_all_fitz1 %>%
  group_by(pop_large, pop_large_short) %>%
  summarize(lat = mean(lat, na.rm=TRUE),
            long = mean(long,na.rm=TRUE),
            .groups = 'drop')


# get map bounds
bound <- c(
  left = min(meta_all_fitz$long, na.rm = T) - 0.05, bottom = min(meta_all_fitz$lat, na.rm = T) - 0.05,
  right = max(meta_all_fitz$long, na.rm = T) + 0.05, top = max(meta_all_fitz$lat, na.rm = T) + 0.05)

# get stadia map
ggmap::register_stadiamaps("your-api")

map <- ggmap::get_stadiamap(bbox = bound, zoom=14, scale=4, 
                            maptype="stamen_terrain",
                            color="bw") %>%
  ggmap::ggmap()


lims <- c(150.27103, 150.378419148936, -33.7563, -33.6996072727273)

fst_map <- map+theme_bw()+
  geom_segment(data=pop_Fst_sig2[pop_Fst_sig2$Fst<=0.5,],
               aes(x=long.x, y=lat.x,
                   xend = long.y, yend = lat.y,
                   size = 1, colour=Fst, alpha=0.5-Fst), lineend = "round")+
  scale_color_gradient(low = "#8DD3C7", high = "white")+
  geom_text(data=pop_Fst_sig2[pop_Fst_sig2$Fst<=0.5,],
            aes(label = round(Fst, 2), x = (long.x + long.y) / 2, y = (lat.y+lat.x) / 2),
            size = 2,
            color = "black")+
  geom_point(data=meta_all_fitz, mapping=aes(x=long, y=lat))+
  xlim(lims[1], lims[2])+
  ylim(lims[3], lims[4])+   labs(x = "Longitude", y = "Latitude", colour="FST") +guides(size = "none", alpha="none")+
  ggrepel::geom_label_repel(data = meta_all_fitz,aes(x = long, y = lat, label=pop_large_short),
                            min.segment.length=0.25, color="black",fill="white",size=3, segment.colour="white", alpha=0.9, label.size=0, nudge_y=0.004)