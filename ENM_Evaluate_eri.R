#Don't forget to set your working directory!

#only need to run this code the first time you install these packages on your machine/user account.

#install.packages("rJava", dependencies = TRUE)
#install.packages("ENMeval", dependencies = TRUE)
#install.packages("raster", dependencies = TRUE)

##installing rJava can be tricky. If you're experiencing problems installing either rJava or
##ENMeval, it's the installation of rJava that is very likely your problem. Keeping in mind
##that these instructions are for a PC, visit this website for troubleshooting:
##https://cimentadaj.github.io/blog/2018-05-25-installing-rjava-on-windows-10/installing-rjava-on-windows-10/

{memory.limit(size = 2*memory.limit())  #Increased limit of available memory
Packages <- c("ENMeval", "MASS", "raster", "rJava")
lapply(Packages, library, character.only = T)}

#Find where the java directory is for the dismo package using the command below. Then, 
#you need to put the file maxent.jar into the indicated directory. You'll need to do
#this second step by hand (not using r)

#system.file("java", package="dismo")

#Occurrences ####
occ = read.table("./sample.csv", h = T, sep = ",", dec = ".")[,-1]

#Bioclimatic variables ####
{lista = list.files("./Bio", full.names = T, pattern = ".asc$")
env = stack(lista)
env}

#make a bias file ####
{occur.ras <- rasterize(occ, env, 1)
#plot(occur.ras)

presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras)), 
              lims = c(extent(env)[1], extent(env)[2], extent(env)[3], extent(env)[4]))
dens.ras <- raster(dens, env)
dens.ras2 <- resample(dens.ras, env)
plot(dens.ras2)

writeRaster(dens.ras2, "biasfile.asc", overwrite = T)}

#Check how many potential background points you have available in the dataset. How many places/grid cells there are on the 
#landscape where pseudo absences/background points could be selected from.

length(which(!is.na(values(subset(env, 1)))))

#If this number is far in excess of 10,000, then use 10,000 background points. If this number is comparable to, or 
#smaller than 10,000, then use 5,000, 1,000, 500, or even 100 background points. The number of available non-NA spaces 
#should be well in excess of the number of background points used.

#For the evaluation below, we need to convert the bias object into another format.
#The code is set up to sample 5,000 background points. It would be better if we
#could sample 10,000 background points, but there are not enough places available.
#If we could change it to 10,000 background points we would change the ", 5000," to ", 10000,"

bg <- xyFromCell(dens.ras2, sample(which(!is.na(values(subset(env, 1)))), 10000, 
                                   prob=values(dens.ras2)[!is.na(values(subset(env, 1)))]))

#Evaluation ####
#If you have fewer than 50 occurrence points, you will need to use the "jackknife" method of model validation.
enmeval_results <- ENMevaluate(occ, env, bg, partitions = "jackknife", algorithm = 'maxent.jar', quiet = T,
                               parallel = T, parallelType = "doParallel", numCores = 12,
                               tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 0.1:5))

write.csv(enmeval_results@results, "enmeval_results.csv")
