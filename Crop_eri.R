{require(rgdal)
require(raster)}

#Importing environmental variables ####
#Present
Lista <- list.files("./ericoides/w/p", full.names = T, pattern = ".tif$")
Env <- stack(Lista)
plot(Env)

Base <- readOGR("./ericoides/Area/area.shp")
Base

Knife <- mask(crop(Env, Base), Base)
plot(Knife)

for (i in c(1:19)) {
  writeRaster(Knife[[i]], paste0("./ericoides/CROP/Present/", "bio_", sprintf("%02d", 1:19)[i]), format = "GTiff")
}


#Future rcp85
ListaF85 = list.files("./ericoides/w/f", full.names = T, pattern = ".tif$")
EnvF85 = stack(ListaF85)
plot(EnvF85)

KnifeF85 <- mask(crop(EnvF85, Base), Base)
plot(KnifeF85)

for (i in c(1:19)) {
  writeRaster(KnifeF85[[i]], paste0("./ericoides/CROP/F85/", "bio_", sprintf("%02d", 1:19)[i]), format = "GTiff")
}


#Past Mid-Holocene
Lista_mid = list.files("./ericoides/w/mid", full.names = T, pattern = ".tif$")
Env_mid = stack(Lista_mid)
plot(Env_mid)

Knife_mid <- mask(crop(Env_mid, Base), Base)
plot(Knife_mid)

for (i in c(1:19)) {
  writeRaster(Knife_mid[[i]], paste0("./ericoides/CROP/MID/", "bio_", sprintf("%02d", 1:19)[i]), format = "GTiff")
}

#Past Last Glacial Maximum
Lista_lgm = list.files("./ericoides/w/lgm", full.names = T, pattern = ".tif$")
Env_lgm = stack(Lista_lgm)
plot(Env_lgm)

Knife_lgm <- mask(crop(Env_lgm, Base), Base)
plot(Knife_lgm)

for (i in c(1:19)) {
  writeRaster(Knife_lgm[[i]], paste0("./ericoides/CROP/LGM/", "bio_", sprintf("%02d", 1:19)[i]), format = "GTiff")
}

#Past Last Inter-Glacial
Lista_lig = list.files("./ericoides/w/lig", full.names = T, pattern = ".tif$")
Env_lig = stack(Lista_lig)
plot(Env_lig)

Knife_lig <- mask(crop(Env_lig, Base), Base)
plot(Knife_lig)

for (i in c(1:19)) {
  writeRaster(Knife_lig[[i]], paste0("./ericoides/CROP/LIG/", "bio_", sprintf("%02d", 1:19)[i]), format = "GTiff")
}


#Resample 
x = raster("./bio3.asc")
#res(x)
y = raster("./altitude.asc")
#res(y)
z = resample(x, y, "ngb", "./bio3.asc", overwrite = T)
#z = raster("./ericoides/CROP/bio6.asc")
#res(z)
#plot(z)

#Soil
Lista_soil <- list.files("./SoilGrids/Brazil", full.names = T, pattern = ".tif$")
Env_soil <- stack(Lista_soil)
plot(Env_soil)

Base <- readOGR("./ericoides/Area/Brazil.shp")
Base

Knife_soil <- mask(crop(Env_soil, Base), Base)
plot(Knife_soil)

for (i in c(1:19)) {
  writeRaster(Knife_soil[[i]], paste0("./ericoides/CROP/Soil/", "bio_", sprintf("%02d", 1:19)[i]), format = "GTiff")
}


Crop = list.files("./ericoides/CROP/PaleoClim/TIF/Present", full.names = T, pattern = ".tif$")
EnvC = stack(Crop)
plot(EnvC)

#writeRaster(env2, "./WC_SA/", format = "GTiff", bylayer = T, suffix = "names")