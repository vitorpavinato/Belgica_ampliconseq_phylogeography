g###############################################
# Making Google maps for publication (almost) #
###############################################

# Vitor Pavinato
# vitor.pavi@gmail.com

# Code based on "http://yarkerconsulting.com/index.php/blog/15-google-maps-and-r"

## ---- rcode_0_A
sessionInfo()$running
sessionInfo()$platform
R.version.string
.Platform$GUI
## ---- rcode_0_B
install.packages(c("RgoogleMaps","ggplot2","ggmap"))
library("RgoogleMaps")
library("ggplot2")
library("ggmap")
library("grid")

## ---- rcode_0end

#############################
## Belgica WGS collections ##
#############################

## Prepare the points to be ploted
## ---- rcode_1
populations <- c("D1", "HP")
longitude <- c(-64.223267,-64.085483)
latitude <- c(-64.725550,-64.764817)

stdata<-data.frame(populations=populations, lat=latitude, lon=longitude)
print(stdata)
## ---- rcode_1end

## Making the first map
## ---- rcode_2
MyMap <- MapBackground(lat=latitude,lon=longitude,zoom=2)
PlotOnStaticMap(MyMap,latitude,longitude,pch=16,col=c("#ff7f00","#377eb8"))
## ---- rcode_2end

## Making a better plot
## ---- rcode_3
# Base map
basemap <- get_map(location=c(lon=mean(longitude),lat=mean(latitude)), zoom = 8, maptype="terrain", source="google",crop=TRUE)
map0 <- ggmap(basemap)
map0
## ---- rcode_3end

## Base map with the corrected axis names
## ---- rcode_4
map1 <- ggmap(basemap, extent="panel", base_layer=ggplot(data=stdata, aes(x=longitude, y=latitude)))
map1
## ---- rcode_4end

## Final map
## ---- rcode_5
pdf("figures/wgs_marker_disc_collections_sites.pdf")
mapf<- map1 + 
        geom_point(aes(colour=factor(stdata$populations)),colour=c("#ff7f00","#377eb8"),size=2,alpha=3/5) +
        labs(x="Longitude", y="Latitude") +
        annotate("segment", x = -64.223267, xend = -65, y = -64.725550, yend = -64.8, colour = "grey40", size=.3) +
        annotate("text", x = -65, y = -64.82, label = "Dream Island site 1 - D1", size=3, colour="gray40") + 
        annotate("segment", x = -64.085483, xend = -64.4, y = -64.764817, yend = -65, colour = "grey40", size=.3) +
        annotate("text", x = -64.4, y = -65.02, label = "Humble Island Penguin - HP", size=3, colour="gray40") 
print(mapf)  
dev.off()
## ---- rcode_5end
