# Maps 


# Basic grey map #################################################################
library("ggsn")
# ag_gps <- aggregate(cbind(lat, long)~site, data = m2, mean) # made previously


ggplot(ozmaps::abs_ste) + geom_sf(fill="#f9f9f9", fill="grey") +
  coord_sf(xlim = divxlims, ylim = divylims) + 
  labs(y=element_blank(), x=element_blank(), fill="sp")+
  geom_point(data = ag_gps, 
             mapping = aes(x = long, y = lat, fill=sp),colour="black",
             size=3, pch=21)+
  theme_few()+
  theme(legend.key.size = unit(0, 'lines'), axis.text.x = element_text(angle=90),
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,-5,-5,-5))+
  guides(fill = guide_legend(title.position = "top", ncol=1))+
  scale_colour_manual(values=sp_colours, aesthetics="fill")+
  ggsn::scalebar(x.min=divxlims[1]+0.1, x.max=divxlims[2]-0.1, 
                 y.min=divylims[1]+0.01, y.max=divylims[2]-0.01,
                 transform=TRUE,
                 dist = 10,dist_unit = "km",  model = 'WGS84',
                 location="bottomright", border.size=0.2, st.size = 3, st.bottom = FALSE)



# Satellite map ##################################################################
# Satellite maps require some setting up. You need to register for a google key to use google satellite images.
# It is worth it tho, once its set up it makes nice maps. Info below: 
#   https://cran.r-project.org/web/packages/ggmap/readme/README.html
#   https://developers.google.com/maps/documentation/geocoding/get-api-key

library(ggmap)
library(ggrepel)

# make an aggregated dataframe with site and species 
ag_gps <- aggregate(cbind(lat, long)~site+sp, data = g_pca_df2, mean)

# sbbox <- c(150.474,-31.488)

# x and y coordinate limits for the plot 
divxlims <- c(min(ag_gps$long)-0.3,
              max(ag_gps$long)+0.3) #find the min / max longitude

divylims <- c(min(ag_gps$lat)-0.2,
              max(ag_gps$lat)+0.2) #find the min / max latitude

# input your google API here ###
register_google(key="###")

# get the map from google
# map <- get_googlemap(sbbox, zoom = 7, maptype = "satellite")
map2 <- get_googlemap(sbbox, zoom = 7, maptype = "hybrid")


# plot the map with aggregated points on top
sat_map <- ggmap(map2, alpha=0.5)+geom_point(ag_gps,
                                             mapping =aes(x = long, y = lat,  fill=sp),
                                             pch=21,color="white", size=3)+labs(x=NULL, y=NULL, fill="Region")+
  scale_color_manual(values=sp_colours, aesthetics="fill")+
  labs(shape="Site")+
  ggsn::scalebar(x.min=divxlims[1]+0.01, x.max=divxlims[2]-0.01,
                 y.min=-32.5, y.max=-29.6,
                 transform=TRUE,
                 dist = 50,dist_unit = "km",  model = 'WGS84',
                 location="topleft", border.size=0.2, st.size = 2.5, st.bottom = FALSE,
                 box.color="white", st.color="white")+theme_bw()+
  theme(legend.position="bottom",
        legend.key.size = unit(0.1, "cm"))+
  guides(fill = guide_legend(direction="vertical", nrow=2))+
  ylim(divylims)+xlim(divxlims)

ggarrange(sat_map, pca_plot, align="hv", labels=c("A","B"), widths=c(1,1.5), 
          legend.grob=get_legend(pca_plot),legend = "right", common.legend = TRUE)


# scatterpie on the satellite map
 ggmap(map2, alpha=0.5) +
  geom_point(ag_gps,mapping =aes(x = long, y = lat),pch=21,color="white", size=3.4)+
  geom_scatterpie(aes(x=long.y, y=lat.y, group ="site", r = 0.15),data =qdf_keyspecies,
                  cols=colnames(qdf_keyspecies)[3:(plateau+2)],  alpha=1, size=0.01, colour=NA)+#
  theme_few()+theme(axis.text.x = element_text(angle=90), legend.position="none")+
  scale_fill_manual(values=colz)+
  labs(x=NULL, y=NULL)+
  ggsn::scalebar(x.min=divxlims[1]+0.01, x.max=divxlims[2]-0.01,
                 y.min=-32.5, y.max=-29.6,
                 transform=TRUE,
                 dist = 50,dist_unit = "km",  model = 'WGS84',
                 location="topleft", border.size=0.2, st.size = 2.5, st.bottom = FALSE,
                 box.color="white", st.color="white")+
  ylim(divylims)+xlim(divxlims)