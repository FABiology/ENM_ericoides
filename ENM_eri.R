#Instalation ####
remotes::install_github("Model-R/modleR", build = T, build_opts = c("--no-resave-data", "--no-manual"))

#Packages ####
{memory.limit(size = 2*memory.limit())  #Increased limit of available memory
Packages <- c("modleR", "raster", "rgdal", "rJava", "doParallel")
lapply(Packages, library, character.only = T)
numCores <- detectCores() #Detect the number of CPU cores on the current host. It is important for parallel computing
registerDoParallel(numCores)} #Register the parallel backend with the foreach package

#Occurrences ####
{records = read.table("sample.csv", h = T, sep = ";", dec = ".")
str(records)
head(records)
tail(records)
records$sp = as.character(records$sp)
species = unique(records$sp)
species}

#Bioclimatic variables ####
{lista = list.files("./Present", full.names = T, pattern = ".tif$")
env = stack(lista)
env}
#plot(env)

#Results directory ####
result_folder <- "./Results"

#Model parameters ####
(ini <- Sys.time())
system.time(foreach(i = 1:length(species), .packages = "modleR", .combine = c) %dopar% {
  occs <- records[records$sp %in% species[i], -1]
  
  sdmdata <- setup_sdmdata(
    species_name = species[i],
    models_dir = result_folder,
    occurrences = occs,
    predictors = env,
    buffer_type = "maximum",
    #dist_buf = 20,
    dist_min = 0.05,
    write_buffer = T,
    clean_dupl = T,
    clean_nas = T,
    clean_uni = T,
    #select_variables = T,
    #percent = 0.5,
    #cutoff = 0.7,
    plot_sdmdata = T,
    n_back = 1000,
    partition_type = "crossvalidation",
    cv_partitions = 10,
    cv_n = 10,
    seed = 50
  )
})

#Modelling and projection ####
pfiles = "./Projection/"

system.time(foreach(i = 1:length(species), .packages = "modleR", .combine = c) %dopar% {
  sp <- species[i]
  occs <- records[records$sp %in% species[i], -1]
  
  do_many(species_name = sp,
          predictors = env,
          models_dir = result_folder,
          project_model = T,
          proj_data_folder = pfiles,
          mask = NULL,
          conf_mat = T,
          write_png = T,
          bioclim = T,
          maxent = T,
          rf = F,
          svmk = F,
          svme = F,
          brt = F,
          glm = F,
          domain = F,
          mahal = F,
          equalize = T,
          write_bin_cut = F)


  final_model(species_name = sp,
              #algorithms = "bioclim",
              select_partitions = T,
              select_par = "AUC",
              select_par_val = 0.8,
              consensus_level = 0.5,
              models_dir = result_folder,
              #proj_dir = "F26_CC",
              which_models = c("raw_mean", "bin_consensus"),
              uncertainty = T,
              scale_models = T,
              write_final = T,
              overwrite = T)

  final_model(species_name = sp,
              #algorithms = "bioclim",
              select_partitions = T,
              select_par = "AUC",
              select_par_val = 0.8,
              consensus_level = 0.5,
              models_dir = result_folder,
              proj_dir = "F85",
              which_models = c("raw_mean", "bin_consensus"),
              uncertainty = T,
              scale_models = T,
              write_final = T,
              overwrite = T)

  final_model(species_name = sp,
              #algorithms = "bioclim",
              select_partitions = T,
              select_par = "AUC",
              select_par_val = 0.8,
              consensus_level = 0.5,
              models_dir = result_folder,
              proj_dir = "MID",
              which_models = c("raw_mean", "bin_consensus"),
              uncertainty = T,
              scale_models = T,
              write_final = T,
              overwrite = T)

  final_model(species_name = sp,
              #algorithms = "bioclim",
              select_partitions = T,
              select_par = "AUC",
              select_par_val = 0.8,
              consensus_level = 0.5,
              models_dir = result_folder,
              proj_dir = "LGM",
              which_models = c("raw_mean", "bin_consensus"),
              uncertainty = T,
              scale_models = T,
              write_final = T,
              overwrite = T)
  
  final_model(species_name = sp,
              #algorithms = "bioclim",
              select_partitions = T,
              select_par = "AUC",
              select_par_val = 0.8,
              consensus_level = 0.5,
              models_dir = result_folder,
              proj_dir = "LIG",
              which_models = c("raw_mean", "bin_consensus"),
              uncertainty = T,
              scale_models = T,
              write_final = T,
              overwrite = T)
  
  ensemble_model(species_name = sp,
                occurrences = occs,
                which_final = "raw_mean",
                models_dir = result_folder,
                #proj_dir = "F26_CC",
                write_ensemble = T,
                overwrite = T)

  ensemble_model(species_name = sp,
                 occurrences = occs,
                 which_final = "raw_mean",
                 models_dir = result_folder,
                 proj_dir = "F85",
                 write_ensemble = T,
                 overwrite = T)

  ensemble_model(species_name = sp,
                 occurrences = occs,
                 which_final = "raw_mean",
                 models_dir = result_folder,
                 proj_dir = "MID",
                 write_ensemble = T,
                 overwrite = T)
  
  ensemble_model(species_name = sp,
                 occurrences = occs,
                 which_final = "raw_mean",
                 models_dir = result_folder,
                 proj_dir = "LGM",
                 write_ensemble = T,
                 overwrite = T)
  
  ensemble_model(species_name = sp,
                 occurrences = occs,
                 which_final = "raw_mean",
                 models_dir = result_folder,
                 proj_dir = "LIG",
                 write_ensemble = T,
                 overwrite = T)
})

stopImplicitCluster()


----------------------------------------------------------------------------------------------------------------------------


library(ENMeval)

#Niche Overlap (similarity)
y = list.files("./Piptolepis_ericoides", full.names = T, pattern = ".tif$")
y = list.files("./Models", full.names = T, pattern = ".asc$")
x = stack(y)

D = calc.niche.overlap(x, overlapStat = "D")
I = calc.niche.overlap(x, overlapStat = "I")
D
I

#t.test(D,I)

library(ENMTools)

#Environmental Layers Correlation
z = list.files("./Present", full.names = T, pattern = ".tif$")
s = stack(z)

cor_plot = raster.cor.plot(s)
#cor_matrix = raster.cor.matrix(s)

#write.csv(cor_matrix, "./cor_matrix.csv")
