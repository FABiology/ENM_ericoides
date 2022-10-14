{memory.limit(size = 2*memory.limit())  #Increased limit of available memory
Packages <- c("ecospat", "landscapetools", "raster", "red", "rgdal", "sf")
lapply(Packages, library, character.only = T)}

#Calculating EOO and AOO from species distribution modeling (SDM) using thresholds ####
#Threshold funcion (MTP - Minimum Training Presence, P10 - 10%-training-presence and user)

raster_threshold <- function(input_raster, points = NULL, type = NULL, threshold = NULL, binary = F) {
  if (!is.null(points)) {
    pointVals <- raster::extract(input_raster, points)
    if (type == "mtp") {
      threshold <- min(na.omit(pointVals))
    } else if (type == "p10") {
      if (length(pointVals) < 10) {
        p10 <- floor(length(pointVals) * 0.9)
      } else {
        p10 <- ceiling(length(pointVals) * 0.9)
      }
      threshold <- rev(sort(pointVals))[p10]
    }
  }
  raster_thresh <- input_raster
  raster_thresh[raster_thresh < threshold] <- NA
  if (binary) {
    raster_thresh[raster_thresh >= threshold] <- 1
  }
  return(raster_thresh)
}

{Records = read.table("./sample.csv", h = T, sep = ";", dec = ".")[-1]
sdm <- raster("./Present.asc")}

#Binarizing through ecospat (MPA - Minimal Predicted Area) ####
mpa <- ecospat.mpa(sdm, Records, 0.9)

mpa_bin <- ecospat.binary.model(sdm, mpa)
#plot(mpa_bin)
#points(Records, pch = 19, cex = 0.5)

user_threshold <- raster_threshold(sdm, threshold = 0.194, binary = T)
#plot(user_threshold)
#points(Records, pch = 19, cex = 0.5)

#Extract occupied patches of a species in geographic space
occupied <- ecospat.occupied.patch(mpa_bin, Records, 500000) #For present SDM
occupied <- ecospat.occupied.patch(user_threshold, Records, 500000) #For projections SDM
#par(mfrow=c(1,2))
plot(occupied)
#points(Records, col="red", cex=0.5, pch=19)

occupied_bin <- raster_threshold(occupied, threshold = 0, binary = T)
#plot(occupied_bin)

#Save binarized shapefile
{conv_poly <- rasterToPolygons(occupied_bin)
st_poly <- st_as_sf(conv_poly)
dis <- st_union(st_poly)
st_write(dis, "Present.shp", driver="ESRI Shapefile", dsn = "D:/MaxEnt_Bin")}

#writeOGR(conv_poly, dsn = "./Bin", layer = "capitatus_occp", driver = "ESRI Shapefile")
