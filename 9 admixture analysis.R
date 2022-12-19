# Admixture analysis 
# Admixture analysis allows the prediction of ancestral genetic groups. I find this often aligns with the UMAP results.
# It takes a long time to do these calculations, so after the initial run I save the results to a csv that I use to make
# plots.
# You will do an initial run of snmf with the kvalrange to identify the optimum number of ancestral populations. 
# Then you can just run one K at a time (ie K=3 instead of K=2:16) in future.
# The ideal number of populations (K) is based upon cross entropy. The ideal K minimises cross entropy and |K|.
# https://rdrr.io/bioc/LEA/man/crossEntropy.html


# Only run the first time #########################################################
library(LEA)
library(adegenet)

# initially we use a range of Ks to identify the ideal K 
kvalrange <- 2:16 # range of Ks to calculate with 

# make lea files
nd_lea <- dart2lea(dms, RandRbase, species, dataset)

# run admixture analysis 
snmf1 <- snmf(nd_lea, K=kvalrange, entropy = TRUE, repetitions = 10, project = "new", CPU=4) 


ce           <- cross.entropy(snmf1, K=2) #snmf_project
Rbest        <- which.min(ce) # run with the lowest cross entropy criterion


# ONLY REQUIRED FOR INITIAL RUN 
# uses a function to find the first K that minimises entropy 
library(nlraa) # get the plateau point
nl <- nls(formula = mean ~ SSquadp3xs(kvalrange, a, b, xs),
          data = cbind(as.data.frame(t(summary(snmf1)$crossEntropy)),as.data.frame(kvalrange)))
plateau <- round(as.data.frame(summary(nl)$coefficients)[3,1]) # plateau is optimium K 

qdf <- Q(snmf1, K=plateau, run=Rbest) # get admixture coefficients (Q) for each sample
qdf <- cbind(as.data.frame(dms$sample_names),qdf) # add the metadata
colnames(qdf)[1] <- "sample_names"

# write the entropy and admixture results to a CSV so that these dont have to be run again 
write.table(qdf, paste0(species, "/outputs/20221118_qdf_admix.csv"), sep=",", row.names=FALSE, col.names=TRUE)
entropy <- t(summary(snmf1)$crossEntropy)
write.table(t(summary(snmf1)$crossEntropy), paste0(species, "/outputs/cross_entropy_admix.csv"),
            sep=",", row.names=FALSE, col.names=TRUE)


# Subsequent analyses / plots  #########################################################
# uses the CSV written from the initial run 

## Make entropy plot   #########################################################
entropy <- fread(paste0(species, "/outputs/cross_entropy_admix.csv"))
kvalrange <- 2:16
# entropy plot 
entropy_plot <- ggplot(entropy, aes(x=kvalrange, y=mean ))+geom_point()+theme_few()+
  labs(x="K", y="Mean\ncross entropy")+geom_vline(xintercept=plateau, col='blue', linetype="dotted")#+geom_line(y=fitted(nl),col='red')
entropy_plot

## Make admix plot  #########################################################
# import qdf data 
qdf <- read.csv(paste0(species, "/outputs/20221118_qdf_admix.csv"))

#get optimum k (plateau)
plateau <- ncol(qdf)-1
colnames(qdf)[2:(plateau+1)] <- 1:plateau # rename columns

# add metadata 
qdf <- merge(qdf, m2, by.x="sample_names", by.y="sample", all.x=TRUE, all.y=FALSE)

# make data wide --> long 
qdf2 <- melt(qdf, measure.vars=colnames(qdf)[2:(plateau+1)], id.vars=colnames(qdf)[-c(2:4)], 
             variable.name="population", value.name="Q") # make data long 


getPalette2 <- colorRampPalette(brewer.pal(n=9, paste("YlOrRd"))) 
colz <- getPalette2(length(unique(qdf2$population))) #get palette for bargraphs

# plot using custom functions
# function variables are (data, grouping variable, x axis text size, x axis text rotation)
site_admix <- admix_plotter(qdf2, site, 10, 90)
site_admix
species_admix <- admix_plotter(qdf2, sp, 10, 90)
species_admix



# Map with admix pie plots #######################################################

library(scatterpie)
library(ozmaps)
library(ggmap)

# This plots pies for each site with average admixtre proportions. 
qdf_keyspecies  <- merge(qdf, ag_gps, by="site")

divxlims <- c(min(qdf_keyspecies$long.y, na.rm=TRUE)-0.2,max(qdf_keyspecies$long.y, na.rm=TRUE)+0.2) #find the min / max longitude
divylims <- c(min(qdf_keyspecies$lat.y, na.rm=TRUE)-0.2,max(qdf_keyspecies$lat.y, na.rm=TRUE)+0.2) #find the min / max latitude

scatterpie_plot_function(qdf_keyspecies, "site", colz, divxlims, divylims, 0.02)+
  ggsn::scalebar(x.min=divxlims[1]+0.01, x.max=divxlims[2]-0.01, # make scalebar, might need adjustments 
                 y.min=divylims[1]+0.01, y.max=divylims[2]-0.01,
                 transform=TRUE,
                 dist = 10,dist_unit = "km",  model = 'WGS84',
                 location="bottomright", border.size=0.2, st.size = 2, st.bottom = FALSE)
