#Installing Rquired packages:
installed.packages("rgbif")
install.packages('dplyr')
install.packages("magrittr")
install.packages('raster')
install.packages('ggplot2')
install.packages('sp')
install.packages('ggplot')
install.packages('CoordinateCleaner')
install.packages('rgdal')
install.packages('countrycode')
install.packages("rnaturalearthdata")
install.packages('usdm')
install.packages('dismo')
install.packages('sdm')
install.packages('tidyverse')
#Loading required packages:
library(rgbif)
library(dplyr)
library(tidyverse)
library(magrittr)
library(raster)
library(ggplot2)
library(ggplot)
library(sp)
library(CoordinateCleaner)
library(rgdal)
library(countrycode)
library(rnaturalearthdata)
library(usdm)
library(dismo)
library(sdm)
installAll()
#Getting occurrence for our target species
#First, we are going to check on the name. Usually, one species has more then one name, so we need
#to certificate that we are searching for the right name
name = name_suggest(q = 'Panthera onca')
#It appears that there is indeed Panthera onca but also names for the subspecies. 
#We need to make sure that we want to take the unique Key for  our target species
name$data$key[1]#getting the key for our target species


#OK, now we can search for the species that match this key:
my_data=occ_search(taxonKey = name$data$key[1],hasCoordinate = TRUE,limit = 10000)
my_data=my_data$data
my_data
#select columns of interest
#First, let us access all column names
names(my_data)


my_data <- my_data %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
                basisOfRecord, institutionCode, datasetName)
#remove records without coordinates
my_data <- my_data%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

#convert country code from ISO2c to ISO3c
my_data$countryCode <-  countrycode(my_data$countryCode, origin =  'iso2c', destination = 'iso3c')

#Now, using the clean coordinates package, we are going to flag some problems in our database

my_data <- data.frame(my_data)
View(my_data)
flags <- clean_coordinates(x = my_data,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros", "countries",'seas')) # most test are on by default
flags
View(flags)
#Exclude problematic records
data_clean <- my_data[flags$.summary,]
#Checking temporal issues with data:
temp_flags <- cf_age(x = data_clean,
                lon = "decimalLongitude",
                lat = "decimalLatitude",
                taxon = "species",
                min_age = "year",
                max_age = "year",
                value = "flagged")

temp_flags

data_clean[!temp_flags, "year"]
data_clean <- data_clean[temp_flags, ]
#OK, now let us visualize the distribution of coordinate Uncertainty
hist(data_clean$coordinateUncertaintyInMeters, breaks = 20)#As you can see by the plot,
#the majority of observations are falling between 0 and 1000km.

data_clean <- data_clean %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 1000 | is.na(coordinateUncertaintyInMeters))

#Plotting again:
hist(data_clean$coordinateUncertaintyInMeters, breaks = 20)#Now, the majority of observations are falling
#between 0 and 200km
#Remove unsuitable data sources
table(data_clean$basisOfRecord)

data_clean <- dplyr::filter(data_clean, basisOfRecord == "HUMAN_OBSERVATION" |
                   basisOfRecord == "OBSERVATION" |
                   basisOfRecord == "PRESERVED_SPECIMEN" |
                     basisOfRecord=='MACHINE_OBSERVATION' |
                     basisOfRecord =='MATERIAL_SAMPLE')
table(data_clean$basisOfRecord)
#In the next step we will remove records with suspicious individual counts. 
#GBIF includes few records of absence (individual count = 0) and suspiciously high occurrence counts, which might indicate inappropriate data or data entry problems.
table(data_clean$individualCount)
data_clean <- data_clean%>%
  filter(individualCount > 0 | is.na(individualCount))%>%
  filter(individualCount < 99 | is.na(individualCount)) # high counts are not a problem
#Checking again:
table(data_clean$individualCount)
#Excluding very old records might be a good idea, specially if these records are before second world war
#Due to low quality technology.
#Age of records
table(data_clean$year)
#Excluding records before WWII
data_clean <- data_clean%>%
  filter(year> 1945)
table(data_clean$year)
#Besides this problems treated here, we also have to check for taxonomic problems.
#First, let us check if every observation has Felidae as family
table(data_clean$family)
#Great! Every record has only one value (Felidae)
#Now, let's check on taxonRank level:
table(data_clean$taxonRank)
#I am good with that

#Now, we are going to check if our data are falling inside the IUCN Polygon
#Setting the working directory for our polygon
setwd('C:/Users/User/Documents/GitHub/Species_Distribution_Model')
polygon = readOGR('data_0.shp')
plot(polygon)

polygon$species <- 'Panthera onca'
polygon@data <- polygon@data[,'species',drop=F]




range_flags <- cc_iucn(x = data_clean,
                       range = polygon,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       value = "flagged")
range_flags
data_final <- data_clean[range_flags, ]


#Identify data set with ddmm(degree minute) to dd.dd conversion error
out.ddmm <- cd_ddmm(data_final, lon = "decimalLongitude", lat = "decimalLatitude",
                    ds = "species", diagnostic = T, diff = 1,
                    value = "dataset")


#Test for rasterized sampling
#We are going to identify datasets with a significant proportion of coordinates that have been collected in large scale lattice designs. 
#These records might have a low precision and might therefore be problematic for some analyses.
par(mfrow = c(2,2), mar = rep(2, 4))
out.round <- cd_round(data_final, lon = "decimalLongitude",
                      lat = "decimalLatitude",
                      ds = "species",
                      value = "dataset",
                      T1 = 11,
                      graphs = T)

#Our data has some problems with autocorrelation between 15 - 20 decimal latitude. We will deal 
#with that now.


#we are going to deal with spatial autocorrelation.
# It is important to deal with this problem because it can cause inflation to the niche model
#First, we are going to transform our dataframe into spatial dataframe
#Saving our dataframe first:
write.csv(data_final,'pantera_data.csv')

#Transforming into spatialpoints dataframe
coordinates(data_final) <- ~decimalLongitude+decimalLatitude
plot(polygon)
points(data_final)
#Applying species thinning
install.packages('spThin') 
library(spThin)


df =data.frame(data_final@coords)
df$species='Panthera onca'
thined=thin(
  loc.data=df ,
  lat.col = "decimalLatitude",
  long.col = "decimalLongitude",
  spec.col = "species",
  thin.par=100,
  reps=5,
  locs.thinned.list.return = TRUE,
  write.files = FALSE, 
  write.log.file = FALSE
)
par(mfrow = c(2,2), mar = rep(2, 4))
plotThin( thined )
#We are going to use the df number 5
panthera=thined[[5]]
coordinates(panthera) <- ~Longitude+Latitude
#Plotting on the polygon

plot(polygon)
points(panthera)
#Great!! Now we have a spatial dataframe with unique coordinates and also corrected for
#spatial auto correlation

#Applying  a CRS:

proj4string(panthera) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

#Now, lets save this as a shape file
shapefile(x = panthera, file = "C:/Users/User/Documents/GitHub/Species_Distribution_Model/panthera_onca.shp")

#OK! We have our spatial rarefied shape file, now we need to gather our independent variables(environmental layers)
#We are going to use the function getdata from Raster Package
bioclim = getData(name = 'worldclim',var = 'bio',res = 5)
bioclim@layers
bio_names=c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10',
            'bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')

setwd('C:/Users/User/Documents/GitHub/Species_Distribution_Model/BioClim')
i=1
for( i in 1:length(bioclim@layers)){
  writeRaster(bioclim@layers[[i]],filename =paste(bio_names[[i]],'.img',sep = ""))
}  

#Now, before starting the modelling procedure, we need couple more step to be done:
#Crop the variables and Retain environmental variables that are not correlated
setwd('C:/Users/User/Documents/GitHub/Species_Distribution_Model/BioClim')
a = list.files(pattern = '.img$')
var = raster::stack(a)
polygon
#cropping...
b = crop(x=var, y =polygon)
var_panthera= mask(x=b, mask=polygon)
par(mfrow = c(1,2))
plot(var_panthera[[1]])
plot(polygon)
#OK, all of our variables are now masked by our polygon, lets select only those with lower VIF


my_VIF =vifstep(x = var_panthera, th=11)

my_VIF@corMatrix#Correlation matrix between variables. It seems that var 14 and 15 have high correlation
#We are going to retain bio15 for ecological meaning 
x = raster::unstack(var_panthera)
var_panthera=stack(c(x[[5]],x[[6]],x[[7]],x[[10]],x[[11]],x[[12]],x[[13]],x[[14]],x[[18]]))
#Great! Let's save this raster stack
writeRaster(var_panthera,'panthera_var_stack.img')

#OK, let's begin the modeling procedure 
#Preparing species occurrence data
setwd('C:/Users/User/Documents/GitHub/Species_Distribution_Model')
a = list.files(pattern='.shp$')


panthera_onca$species <- 1###########Create column where all values are number 1(1 means presence of species)
panthera_onca@data <- panthera_onca@data[,'species',drop=F]####Drop all columns and keep only the species column

#For our models to produce satisfactory outputs, first we need to select the best parameters.
#We are going to produce several models with different numbers of pseudo absence and do a sensibility test
#Our goal is to choose the model with high sensibility

setwd('C:/Users/User/Documents/GitHub/Species_Distribution_Model/Settings Selection')
#Creating SDM data
sdmData_1=sdmData(species~.,train=panthera_onca,
          predictors=var_panthera,bg=list(n=128,method='vRandom',remove=T))#Creating a model with proportion of 1 to 1 presence pseudo absence
write.sdm(x = sdmData_1,filename = 'data_1.sdd')

sdmData_2=sdmData(species~.,train=panthera_onca,
                  predictors=var_panthera,bg=list(n=256,method='vRandom',remove=T))#Creating a model with proportion of 1 to 2 presence pseudo absence
write.sdm(x = sdmData_2,filename = 'data_2.sdd')

sdmData_10=sdmData(species~.,train=panthera_onca,
                  predictors=var_panthera,bg=list(n=1280,method='vRandom',remove=T))#Creating a model with proportion of 1 to 10 presence pseudo absence
write.sdm(x = sdmData_10,filename = 'data_10.sdd')

#If you want to load our sdm data, please follow these steps:
setwd('C:/Users/User/Documents/GitHub/Species_Distribution_Model/Settings Selection')
a = list.files(pattern = '.sdd$')
sdmData_1 = read.sdm(a[[1]])
sdmData_10= read.sdm(a[[2]])
sdmData_2= read.sdm(a[[3]])

#Now that we have all sdm data, lets test Which model has a better performance
#Creating our SDM model:
#First, for cross validation method
m1_CV<-sdm(species~.,data = sdmData_1, methods=c('glm','brt','maxent','svm'),
        replication='cross-validation',cv.folds=5,n = 10, parallelSettings = list(ncore=4))

m2_CV<-sdm(species~.,data = sdmData_2, methods=c('glm','brt','maxent','svm'),
        replication='cross-validation',cv.folds=5,n = 10, parallelSettings = list(ncore=4))

m10_CV<-sdm(species~.,data = sdmData_10, methods=c('glm','brt','maxent','svm'),
        replication='cross-validation',cv.folds=5,n = 10, parallelSettings = list(ncore=4))


#Saving those models:
write.sdm(x = m1_CV,filename = 'm1_CV.sdm')
write.sdm(x = m2_CV,filename = 'm2_CV.sdm')
write.sdm(x = m10_CV,filename = 'm10_CV.sdm')

#Now, we are going to create models using bootstreping partitioning

m1_boot<-sdm(species~.,data = sdmData_1, methods=c('glm','brt','maxent','svm'),
           replication='boot',n = 10, parallelSettings = list(ncore=4))

m2_boot<-sdm(species~.,data = sdmData_2, methods=c('glm','brt','maxent','svm'),
           replication='boot',n = 10, parallelSettings = list(ncore=4))

m10_boot<-sdm(species~.,data = sdmData_10, methods=c('glm','brt','maxent','svm'),
            replication='boot',n = 10, parallelSettings = list(ncore=4))

#Saving those models:
write.sdm(x = m1_boot,filename = 'm1_boot.sdm')
write.sdm(x = m2_boot,filename = 'm2_boot.sdm')
write.sdm(x = m10_boot,filename = 'm10_boot.sdm')

#OK, analyzing the results, we can see that none of the models did well.
#One possible reason for that is that our calibration area are to small. Let's take a look
#First, let's create a shape file for our background data

m1_boot@data@info@coords#This is all data points(presence and background)
m1_boot@data@species$species@background#This is the index of only the background data
background = as.data.frame(m1_boot@data@info@coords)#Creating a dataframe with all datapoints
background_sliced=background[c(m1_boot@data@species$species@background),]#Slicing our dataframe to only return
#background data
names(background_sliced)
coordinates(background_sliced) <- ~coords.x1+coords.x2
proj4string(background_sliced) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

#Let's plot:
plot(polygon)
points(panthera)
points(background_sliced,pch = 0)
#As we can see, the number of background its all over the place. This is a problem because the model
#Will be truncated
#Lets create models but with half of the number of background!

sdmData0.5=sdmData(species~.,train=panthera,
                  predictors=var_panthera,bg=list(n=64,method='vRandom',remove=T))#Creating a model with proportion of 1 to 0.5 presence pseudo absence
write.sdm(x = sdmData0.5,filename = 'data_05.sdd')
#OK, creating our models for CV method
m0.5_CV<-sdm(species~.,data = sdmData0.5, methods=c('glm','brt','maxent','svm'),
           replication='cross-validation',cv.folds=5,n = 10, parallelSettings = list(ncore=4))
#Model with bootstrapping method
m0.5_boot<-sdm(species~.,data = sdmData0.5, methods=c('glm','brt','maxent','svm'),
             replication='boot',n = 10, parallelSettings = list(ncore=4))


m0.5_CV
m0.5_boot

m1_boot
m2_boot
m10_boot

m1_CV
m2_CV
m10_CV
#Now, we will get the mean sensibility for all models
my_models=c(m0.5_CV,m0.5_boot,m1_boot,m2_boot,m10_boot,m1_CV,m2_CV,m10_CV)
i=1
for (i in 1:length(my_models)){
  print(mean(getEvaluation(my_models[[i]],wtest='test.dep',stat='sensitivity',opt=2)$sensitivity))
}
#According to our print, the model m1_boot had the best sensibility. So, this is the model that
#we are going to use.

#Let's create some model statistics for our selected model
m1_boot
my_stat = getEvaluation(m1_boot,wtest='test.dep',stat=c('sensitivity','specificity','TSS','AUC'),opt=2)
mean(my_stat$AUC)
mean(my_stat$sensitivity)
mean(my_stat$specificity)
mean(my_stat$TSS)
#Our model got reasonable results, not great, but reasonable
#Let's create our ensemble model
#First, we are going to select only those models ID that had a TSS higher then the mean TSS value
ids = my_stat$modelID[which(my_stat$TSS>mean(my_stat$TSS))]

en = ensemble(m1_boot, var_panthera, filename="panthera_ensemble.img",setting=list(method='mean-weighted',stat='TSS'))
plot(en)
points(panthera)

#OK, now we are going to create our presence-absence model:
list_th=list()
i = 1
for(i in 1:40){
  list_th[[i]]=getEvaluation(m1_boot,w=i,wtest='test.dep',stat=2)[[2]][[2]]
}
mean(unlist(list_th))
binary<-en
binary[]<-ifelse(en[]>=0.4993123,1,0)
plot(binary)
points(panthera)
writeRaster(x = binary,'panthera_binary.img')

#Generating ROC Curves:
roc(m1_boot,smooth = TRUE)
#Generating response curves
rcurve(m1_boot,id=1:40,smooth = T)
#Generating variable importance
my_var_imp=list()
i=1
for(i in 1:40){
  my_var_imp[[i]]=getVarImp(m1_boot,id=i,wtest='test.dep')@varImportance$AUCtest
}

my_var_imp

