## advanced admix method

## uses a lot of vectorisation
## plots multiple at the same time

# Only run the first time #########################################################
library(LEA)
library(adegenet)

# initially we use a range of Ks to identify the ideal K 
kvalrange <- 2:10 # range of Ks to calculate with 

# make lea files
nd_lea <- dart2lea(dms, RandRbase, species, dataset)

# run admixture analysis 
snmf1 <- snmf(nd_lea, K=kvalrange, entropy = TRUE, repetitions = 3, project = "new", CPU=8) 

# save all of the outputs to csvs
for (k in kvalrange){
  ce           <- cross.entropy(snmf1, K=k) #snmf_project
  Rbest        <- which.min(ce) # run with the lowest cross entropy criterion
  qdf <- Q(snmf1, K=k, run=Rbest) # get admixture coefficients (Q) for each sample
  qdf <- cbind(as.data.frame(dms$sample_names),qdf) # add the metadata
  colnames(qdf)[1] <- "sample_names"
  # write the entropy and admixture results to a CSV so that these dont have to be run again 
  write.table(qdf, paste0(species, paste0("/outputs/20230509_k",k,"_admix.csv")), sep=",", row.names=FALSE, col.names=TRUE)
}

# save the entropy
entropy <- t(summary(snmf1)$crossEntropy)
write.table(t(summary(snmf1)$crossEntropy), paste0(species, "/outputs/cross_entropy_admix.csv"),
            sep=",", row.names=FALSE, col.names=TRUE)



## make the plots from the saved CSV files ######################################
#### read in the csv files with the k of interest (in this case 5:8) ####
kvalrange <- 5:8
qdfs_wide <- list()
for (k in kvalrange) {
  qdfs_wide[[k]] <- read.csv(paste0(species, "/outputs/20230509_k", k, "_admix.csv"))
  colnames(qdfs_wide[[k]])[2:(k+1)] <- 1:k # rename columns
}

# If you are doing a facet grid my site, create a vector with the site order you need
# latitude order but parviflora is grouped
site_order <- c("Slippy Downs", "Evans Head", "Snapper Head", "Woombah", "Pillar Valley",
                "Pebbly Beach", "Barcoongere", "Woolgoolga headland", "Bare Bluff",
                "Green Bluff", "Macauley", "Coffs Creek", "Crescent Head", "Limeburners Creek",
                "Herbarium", "Boat Harbour", "Hickson Street walk", "Scenic drive",
                "Pinny Beach Headland", "Terrigal Haven", "Huskisson", "Hammondville", "Pitt Town", 
                "Wianamatta NR", "Shanes Park", "Kemps Creek",
                "Menai")


#### order samples by similarity -- to be used for ordering WITHIN facets ####
d <-  dist(dms$gt, diag=TRUE)#dist(qdf7[, -1])

# Perform hierarchical clustering
hc <- hclust(d)

# Get the order of the samples based on the clustering
order <- hc$order %>% dms$sample_names[.]

# Reorder the data frame by the sample order
library(tidyr)
for (k in kvalrange) {
  qdfs_wide[[k]] <- qdfs_wide[[k]][match(order,qdfs_wide[[k]][,"sample_names"]),]
  qdfs_wide[[k]][,"order"] <- 1:nrow(qdfs_wide[[k]])
  qdfs_wide[[k]] <- merge(qdfs_wide[[k]], m2, by.x="sample_names", by.y="sample", all.x=TRUE, all.y=FALSE)
}

# make the qdfs long instead of wide
qdfs <- qdfs_wide
for (k in kvalrange) {
  qdfs[[k]] <- pivot_longer(qdfs[[k]], cols=2:(k+1), names_to="population", values_to="Q") 
}

# get the colours for the plots using the largest K
getPalette2 <- colorRampPalette(brewer.pal(n=9, paste("Spectral"))) #YlOrRd
colz <- getPalette2(length(unique(qdfs[[max(kvalrange)]]$population))) #get palette for bargraphs

#### functions for making plots ####
# for the first plot (has facet labels) -- best for LARGE groups eg species level, genetic group level
admix_plot_list <- function(df, group){
  ggplot(df, 
         aes(x=reorder_within(sample_names, order, {{group}}), y=Q, fill=population))+
    geom_bar(position="stack", stat="identity")+
    theme_few()+
    labs(y="Admixture\ncoefficient (Q)", x=element_blank(), fill="Source\npopulation")+
    scale_fill_manual(values=colz)+
    facet_grid(cols=vars({{group}}), scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4),
          strip.text.x = element_text(angle = 90, size=8, face="italic"),
          strip.background = element_blank(),
          panel.spacing = unit(0.07, "lines"),
          legend.position = "bottom",
          legend.text=element_text(size=8),
          legend.key.width=unit(0.5, "cm"))+
    scale_y_continuous(limits = c(0,1.001), expand=c(0,0))+
    scale_x_discrete(labels=NULL)+
    guides(fill = guide_legend(nrow = 1))
  
}

# for the subsequent plots (no facet labels)
admix_plot_list2 <- function(df, group){
  ggplot(df, 
         aes(x=reorder_within(sample_names, order, {{group}}), y=Q, fill=population))+
    geom_bar(position="stack", stat="identity")+
    theme_few()+
    labs(y="Admixture\ncoefficient (Q)", x=element_blank(), fill="Source\npopulation")+
    scale_fill_manual(values=colz)+
    facet_grid(cols=vars({{group}}), scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4),
          strip.text.x = element_blank(),
          strip.background = element_blank(),
          panel.spacing = unit(0.07, "lines"),
          legend.position = "bottom",
          legend.text=element_text(size=8),
          legend.key.width=unit(0.5, "cm"))+
    scale_y_continuous(limits = c(0,1.001), expand=c(0,0))+
    scale_x_discrete(labels=NULL)+
    guides(fill = guide_legend(nrow = 1))
  
}

# for the first sit plot (has facet labels) -- best for SMALL groups eg sites
admix_plot_siteorder <- function(df, site_order){
  ggplot(df, 
         aes(x=reorder_within(sample_names, order, site), #, labels=NULL),
             y=Q, fill=population))+
    geom_bar(position="stack", stat="identity")+
    theme_few()+
    labs(y="Admixture\ncoefficient (Q)", x=element_blank(), fill="Source\npopulation")+
    scale_fill_manual(values=colz)+
    facet_grid(~factor(site, levels={{site_order}}), scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4),
          strip.text.x = element_text(size = 6), strip.background = element_blank())+
    scale_y_continuous(limits = c(0,1.001), expand=c(0,0))+
    scale_x_discrete(labels=NULL)+
    theme(strip.text.x = element_text(angle = 90, size=8), panel.spacing = unit(0.07, "lines"))
}

# for subsequent site plots
admix_plot_siteorder2 <- function(df, site_order){
  ggplot(df, 
         aes(x=reorder_within(sample_names, order, site), #, labels=NULL),
             y=Q, fill=population))+
    geom_bar(position="stack", stat="identity")+
    theme_few()+
    labs(y="Admixture\ncoefficient (Q)", x=element_blank(), fill="Source\npopulation")+
    scale_fill_manual(values=colz)+
    facet_grid(~factor(site, levels={{site_order}}), scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4),
          strip.text.x = element_blank(),
          panel.spacing = unit(0.07, "lines"),
          strip.background = element_blank())+
    scale_y_continuous(limits = c(0,1.001), expand=c(0,0))+
    scale_x_discrete(labels=NULL)+
    guides(fill = guide_legend(nrow = 1))
}


#### Making admix plots with LARGE facet grid groups eg species level, genetic group level ####
# make the plots
plots <- list()
for (k in kvalrange) {
  if (k== min(kvalrange)){
    p <- admix_plot_list(qdfs[[k]], genetic_group2)
    plots[[k]] <- p
  }
  else{
    p <- admix_plot_list2(qdfs[[k]], genetic_group2)
    plots[[k]] <- p
  }
}
# plot the plots
ggarrange(plots[[5]],plots[[6]],plots[[7]], plots[[8]], ncol=1, heights=c(1.6,1,1,1),
          common.legend = TRUE, legend.grob = get_legend(plots[[8]]), legend="bottom")


#### Making admix plots with SMALL facet grid groups eg site level ####
plots2 <- list()
for (k in kvalrange) {
  if (k== min(kvalrange)){
    p <- admix_plot_siteorder(qdfs[[k]], site_order)
    plots2[[k]] <- p
  }
  else{
    p <- admix_plot_siteorder2(qdfs[[k]], site_order)
    plots2[[k]] <- p
  }
}

ggarrange(plots2[[5]],plots2[[6]],plots2[[7]], plots2[[8]], ncol=1, heights=c(1.6,1,1,1),
          common.legend = TRUE, legend.grob = get_legend(plots2[[8]]), legend="bottom")


### scatterpie ####

library(scatterpie)
library(ozmaps)
library(ggmap)
library(ggsn)

# This plots pies for each site with average admixtre proportions. 
ag_gps <- aggregate(cbind(lat, long,`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`)~site, data = qdfs_wide[[8]], mean)

divxlims <- c(min(ag_gps$long, na.rm=TRUE)-0.2,max(ag_gps$long, na.rm=TRUE)+0.2) #find the min / max longitude
divylims <- c(min(ag_gps$lat, na.rm=TRUE)-0.2,max(ag_gps$lat, na.rm=TRUE)+0.2) #find the min / max latitude


ggplot(ozmaps::abs_ste) + geom_sf(fill="#f9f9f9", colour="grey") +
  coord_sf(xlim = divxlims, ylim = divylims) + labs(y=element_blank(), x=element_blank(), fill="Source\npopulation")+
  geom_scatterpie(aes(x=long, y=lat, group ="site", r = 0.15),
                  data =ag_gps,
                  cols=colnames(ag_gps)[4:(8+3)],  alpha=1, size=0.01, colour=NA)+#
  theme_few()+theme(axis.text.x = element_text(angle=90) )+
  scale_fill_manual(values=colz)+
  ggsn::scalebar(x.min=divxlims[1]+0.01, x.max=divxlims[2]-0.01, # make scalebar, might need adjustments 
                 y.min=divylims[1]+0.01, y.max=divylims[2]-0.01,
                 transform=TRUE,
                 dist = 100,dist_unit = "km",  model = 'WGS84',
                 location="bottomright", border.size=0.2, st.size = 2, st.bottom = FALSE)


