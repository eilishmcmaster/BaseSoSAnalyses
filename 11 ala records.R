# plotting ALA records

# install.packages("galah")
library(galah)
library(ggrepel)
library(ggmap)
library(openxlsx)

# meta <- read.xlsx("/Users/eilishmcmaster/Documents/ZierObco/ZierObco/meta/ZierObco_DZ22-7321_meta.xlsx") # get meta for lat long limits

galah_config(email = "##") # your email here ##

result <- galah_call() |>
  galah_identify("Zieria obcordata", "Zieria ingramii", "Zieria odorifera") |> # list of species to get records 
  atlas_occurrences()

result2 <- result[result$dataResourceName!="iNaturalist Australia",] # remove inaturalist (contributed by the public) 

divxlims <- c(min(meta$long)-0.3,
              max(meta$long)+0.3) #find the min / max longitude

divylims <- c(min(meta$lat)-0.2,
              max(meta$lat)+0.2) #find the min / max latitude

library("ozmaps")
ggplot(ozmaps::abs_ste) + geom_sf(fill="#f9f9f9", fill="grey") +
  coord_sf(xlim = divxlims, ylim = divylims) + 
  labs(y=element_blank(), x=element_blank(), colour="ALA records")+
  geom_point(result2, mapping=aes(x=decimalLongitude, y=decimalLatitude, colour=scientificName), alpha=0.5)+
  theme_few()+theme(axis.text.x = element_text(angle=90))


