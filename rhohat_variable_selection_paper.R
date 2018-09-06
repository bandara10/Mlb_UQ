#### Load libraries ####

# Do list. Rasters in 300 dpi. and rhohat in 300 dpi. 
library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatstat.utils)

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data")

#Step 1: load koala data from BoalaBASE#####

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata=mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]
myalldata=mydata
myalldata[!duplicated(myalldata$X &myalldata$Y), ]

all.locations= subset(myalldata,yearnew >=2005, select=X:Y)#>2005 :mydata not myalldata.
pp=SpatialPoints(all.locations)
selected.locations=all.locations
### Step 01: Bring in distance rasters #####

myfullstack.b <- list.files(path="C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/warton_data_allclimate/distance_vars",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = stack(myfullstack.b)
plot(myfullstack.b,2)
names(myfullstack.b)
myextent=c(387900, 553100, 6862400, 7113600) 
studyarea.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data\\warton_data_allclimate\\ss_area.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 

myextent=studyarea.shp
habitat.rr=crop(myfullstack.b, myextent, snap="near")
names(habitat.rr)
habitat.r<- subset(habitat.rr, c(1:2)) # select my variables

set.seed(123)
selected.locations=gridSample(all.locations, habitat.r, n = 10)



####### Select records based on distance######
# 
# source("Lib_DistEstimatesToItself.r")## set a minimum distance between koalas
# all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
# select.locations = subset(all.locations, all.locations$disttoitself > 200)
# selected.locations = select.locations[,1:2]
# selected.locations=SpatialPoints(selected.locations)
# proj4string(selected.locations) <- CRS("+init=epsg:28356") 
######### Bring in MGA56 square study area boundary map:######


# detail study area

dstudyarea.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 

plot(dstudyarea.shp, axes = TRUE)
plot(studyarea.shp, add = TRUE)
points(x = mydata$X, y = mydata$Y) #  points(x = selected.locations$X, y = selected.locations$Y)


# Make a ppp object for spatstat:=================================================================================================================================
# Create an observation window which is an intersection of the square study area boundaries and the detailed study area:

source("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\owin2sp_source.r")

studyarea.shp= readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data\\warton_data_allclimate\\ss_area.shp") # warton data folder.
dstudyarea.shp <- dstudyarea.shp[, -(1:3)] 

dat.w <- as(as(dstudyarea.shp, "SpatialPolygons"), "owin") # warton data climate
dat.spw <- owin2SP(dat.w)   # ignore error message:undefined columns selected

# Set projection of the owin as GDA94 / SA Lambert:
proj4string(dat.spw) <- CRS("+init=epsg:28356") 
dat.spw <- raster::intersect(x = dstudyarea.shp, y = dat.spw)

###### Convert the sp object into a window:#####
dat.w <- as.owin(dat.spw)
plot(dat.ppp)

# #new
rdat.ppp <- rescale(X = dat.ppp, s = 1000, unitname = "kilometres")
plot(rdat.w, axes = TRUE)

# Select only those koalas within the owin:
id <- inside.owin(x = selected.locations$X, y = selected.locations$Y, w = dat.w)
selected.locations <- selected.locations[id,]

# acsel=selected.locations
plot(dstudyarea.shp, axes = TRUE)
points(x = selected.locations[,1], y = selected.locations[,2])

# Make a ppp object:
dat.ppp <- ppp(x = selected.locations[,1], y = selected.locations[,2], window = dat.w)

plot(dat.ppp, axes = TRUE)



# Rhohat distance_vars Convert rasters into a spatstat image object=================================================================================================================================

#step 1distance_motorwayandlink==================================================================================================================================
setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/Paper one abundance/Figures_paper on abundace")

tiff("roads_rhohat.tiff",width=20,height=10,units="in",res=300) # 8:4 for three mpas. one map 8:10 # 
par(mfrow=c(3,3),oma=c(1,1,1,1),mai = c(.5, .5, .5, .5),mar=c(4,5,4,4))


#1.distance_primaryandlink==================================================================================================================================

habitat.rr$distance_primaryandlink <- habitat.rr$distance_primaryandlink / 1000
windows(); plot(habitat.rr$distance_primaryandlink)

dprl.im <- as.im(habitat.rr$distance_primaryandlink)
plot(dprl.im )

rdprl.im <- rescale.im(X = dprl.im, s = 1000, unitname = "kilometres")
windows(); plot(rdprl.im)

# What is the nature of the association between koala sight locations and total phosphorous?
rdprl.rho <- rhohat(object = rdat.ppp, covariate = rdprl.im)

windows(); plot(rdprl.rho, cex.axis=2,cex.lab=2,cex.main=2, cex.sub=2, ylab = "Number of koalas per square km", xlab = "Distance (km)", main = "distance_primary", xlim = c(0, 5))
# p=predict(dprl.rho) ; plot(p)

#2.distance_secondaryandlink==================================================================================================================================
dsecLink.im <- as.im(habitat.rr$distance_secondaryandlink)
plot(dsecLink.im)

# What is the nature of the association between koala sight locations and total phosphorous?
dsecLink.rho <- rhohat(object = dat.ppp, covariate = dsecLink.im)
plot(dsecLink.rho ,ylab = "intensity",xlab = "distance (m)",ylim= c(0,8e-07), cex.axis=2,cex.lab=2,cex.main=2, cex.sub=2,main = "distance_secondary")

#3.distance_tertiaryandlink==================================================================================================================================
dTerLink.im <- as.im(habitat.rr$distance_tertiaryandlink)
plot(dTerLink.im)

# What is the nature of the association between koala sight locations and total phosphorous?
dTerLink.rho <- rhohat(object = dat.ppp, covariate = dTerLink.im)
plot(dTerLink.rho,ylab = "intensity",xlab = "distance (m)", cex.axis=2,cex.lab=2,cex.main=2, main = " distance_tertiary")


#4.distance_residentil==================================================================================================================================
dresid.im <- as.im(habitat.rr$distance_residentil)
plot(dresid.im)
windows(); plot(rdprl.im)
# What is the nature of the association between koala sight locations and total phosphorous?
dresid.rho <- rhohat(object = dat.ppp, covariate = dresid.im)
plot(dresid.rho, ylab = "intensity",xlab = "distance (m)", cex.axis=2,cex.lab=2,cex.main=2, main = "distance_residential")

#5.distance_unclassified==================================================================================================================================
dUnclasi.im <- as.im(habitat.rr$distance_unclassified)
plot(dUnclasi.im)

# What is the nature of the association between koala sight locations and total phosphorous?
dUnclasi.rho <- rhohat(object = dat.ppp, covariate = dUnclasi.im)
plot(dUnclasi.rho, ylab = "intensity", xlab = "distance (m)", cex.axis=2,cex.lab=2,cex.main=2,main = "distance_unclassified")

### motorway 

dmwl.im <- as.im(habitat.rr$distance_motorwayandlink)
plot(dmwl.im )

# What is the nature of the association between koala sight locations and total phosphorous?
dmwl.rho <- rhohat(object = dat.ppp, covariate = dmwl.im)
plot(dmwl.rho, ylab = "intensity",xlab = "distance (m)", cex.axis=2,cex.lab=2,cex.main=2, cex.sub=2, main = "distance_motorway")
graphics.off()


#Dis_habitat_suitable_1==================================================================================================================================
jpeg("habitat_per1.jpeg",width=25,height=6,units="in",res=800) # 8:4 for three mpas. one map 8:10 # graphics.off()
par(mfrow=c(1,3),oma=c(1,1,1,1),mai = c(.5, .5, .5, .5),mar=c(4,5,4,4))


#####
dHS1.im <- as.im(habitat.rr$Dis_habitat_suitable_1)
plot(dHS1.im)

# What is the nature of the association between koala sight locations and total phosphorous?
dHS1.rho <- rhohat(object = dat.ppp, covariate = dHS1.im)
plot(dHS1.rho, ylab= "intensity", xlab = "distance(m)",cex.axis=2,cex.lab=2,cex.main=2, cex.sub=2, main = "Distance_highly_suitable_habitat")

#Dis_habitat_suitable_2==================================================================================================================================
dHS2.im <- as.im(habitat.rr$Dis_habitat_suitable_2)
plot(dHS2.im)

# What is the nature of the association between koala sight locations and total phosphorous?
dHS2.rho <- rhohat(object = dat.ppp, covariate = dHS2.im)
plot(dHS2.rho,  ylab= "intensity", xlab = "distance(m)", cex.axis=2,cex.lab=2,cex.main=2, cex.sub=2,main = "Distance_suitable_habitat")

#Dis_habitat_suitable_3==================================================================================================================================
dHS3.im <- as.im(habitat.rr$Dis_habitat_suitable_3)
plot(dHS3.im)

# What is the nature of the association between koala sight locations and total phosphorous?
dHS3.rho <- rhohat(object = dat.ppp, covariate = dHS3.im)
plot(dHS3.rho, ylab= "intensity", xlab = "distance(m)", cex.axis=2,cex.lab=2,cex.main=2, cex.sub=2,main = "Distance_low_suitable_suitable")

graphics.off()
###
#s1_secondry_dist==================================================================================================================================
jpeg("habitat_perS_road.jpeg",width=12,height=8,units="in",res=300) # 8:4 for three mpas. one map 8:10 # graphics.off()
par(mfrow=c(2,3),oma=c(.5,.5,.5,.5))


ds1_secondaryR.im <- as.im(habitat.rr$s1_secondry_dist)
plot(ds1_secondaryR.im)

# What is the nature of the association between koala sight locations and total phosphorous?
ds1_secondaryR.rho <- rhohat(object = dat.ppp, covariate = ds1_secondaryR.im)
plot(ds1_secondaryR.rho, xlab = "s1_secondry_dist(units)", main = "")

#s1_unclassified_dist==================================================================================================================================
ds1_unclassi.im <- as.im(habitat.rr$s1_unclassified_dist)
plot(ds1_secondaryR.im)

# What is the nature of the association between koala sight locations and total phosphorous?
ds1_unclassi.rho <- rhohat(object = dat.ppp, covariate = ds1_unclassi.im)
plot(ds1_unclassi.rho, xlab = "s1_unclassified_dist(units)", main = "")

#s2_residential_dist==================================================================================================================================
ds2_resid.im <- as.im(habitat.rr$s2_residential_dist)
plot(ds1_secondaryR.im)

# What is the nature of the association between koala sight locations and total phosphorous?
ds2_resid.rho <- rhohat(object = dat.ppp, covariate = ds2_resid.im)
plot(ds2_resid.rho, xlab = "s2_residential_dist(units)", main = "")


#s3_secondry_dist==================================================================================================================================
ds3_secondary.im <- as.im(habitat.rr$s3_secondry_dist)
plot(ds3_secondary.im)

# What is the nature of the association between koala sight locations and total phosphorous?
ds3_secondary.rho <- rhohat(object = dat.ppp, covariate = ds3_secondary.im)
plot(ds3_secondary.rho, xlab = "s3_secondry_dist(units)", main = "")

#s3_unclassified_dist==================================================================================================================================
ds3_unclass.im <- as.im(habitat.rr$s3_unclassified_dist)
plot(ds3_unclass.im)

# What is the nature of the association between koala sight locations and total phosphorous?
ds3_unclass.rho <- rhohat(object = dat.ppp, covariate = ds3_unclass.im)
plot(ds3_unclass.rho, xlab = "s3_unclassified_dist(units)", main = "")

graphics.off()


###Step 02:  Bring in climate rasters #####

myfullstack.c <- list.files(path="C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/warton_data_allclimate/climate_vars",pattern="\\.tif$",full.names=TRUE)
myfullstack.c = stack(myfullstack.c)
plot(myfullstack.c,2)
names(myfullstack.c)
# myextent=c(387900, 553100, 6862400, 7113600) 
myextent=studyarea.shp
habitat.rr=crop(myfullstack.c, myextent, snap="near")

habitat.r<- subset(habitat.rr, c(1:2)) # select my variables
plot(habitat.r)

#Annual_Mean_Temperature==================================================================================================================================
amt.im <- as.im(habitat.rr$Annual_Mean_Temperature)
plot(amt.im )

# What is the nature of the association between koala sight locations and total phosphorous?
amt.rho <- rhohat(object = dat.ppp, covariate = amt.im)
plot(amt.rho, xlab = "Annual_Mean_Temperature (units)", main = "")

plot(amt.rho, xlim = c(100, 120), xlab = "Annual_Mean_Temperature (metres)", main = "")


#Annual_Precipitation==================================================================================================================================
apt.im <- as.im(habitat.rr$Annual_Precipitation)
plot(apt.im )

# What is the nature of the association between koala sight locations and total phosphorous?
apt.rho <- rhohat(object = dat.ppp, covariate = apt.im)
plot(apt.rho, xlab = "Annual_Precipitation (units)", main = "")

#Max_Temperature_of_Warmest_Month"    (Temperature annual range)==================================================================================================================================
mtwm.im <- as.im(habitat.rr$Max_Temperature_of_Warmest_Month)
plot(mtwm.im )

# What is the nature of the association between koala sight locations and total phosphorous?
mtwm.rho <- rhohat(object = dat.ppp, covariate = mtwm.im)
plot(mtwm.rho, xlim = c(140, 200),xlab = "Max_Temperature_of_Warmest_Month (units)", main = "")

#Mean_Diurnal_Range_    (Temperature annual range)==================================================================================================================================
meanDiuR.im <- as.im(habitat.rr$Mean_Diurnal_Range_)
plot(meanDiuR.im )

# What is the nature of the association between koala sight locations and total phosphorous?
meanDiuR.rho <- rhohat(object = dat.ppp, covariate = meanDiuR.im)
plot(meanDiuR.rho, xlab = "Mean_Diurnal_Range_(units)", main = "")

#Mean_Temperature_of_Coldest_Quarter   ==================================================================================================================================
meanTcoldestM.im <- as.im(habitat.rr$Mean_Temperature_of_Coldest_Quarter)
plot(meanDiuR.im )

# What is the nature of the association between koala sight locations and total phosphorous?
meanTcoldestM.rho <- rhohat(object = dat.ppp, covariate = meanTcoldestM.im)
plot(meanTcoldestM.rho, xlab = "Mean_Temperature_of_Coldest_Quarter(units)", main = "")


#Mean_Temperature_of_Driest_Quarter   ==================================================================================================================================
meanTDquarter.im <- as.im(habitat.rr$Mean_Temperature_of_Driest_Quarter)
plot(meanTDquarter.im )

# What is the nature of the association between koala sight locations and total phosphorous?
meanTDquarter.rho <- rhohat(object = dat.ppp, covariate = meanTDquarter.im)
plot(meanTDquarter.rho, xlab = "Mean_Temperature_of_Driest_Quarter(units)", main = "")

#Mean_Temperature_of_Warmest_Quarter ==================================================================================================================================
meanTWquarter.im <- as.im(habitat.rr$Mean_Temperature_of_Warmest_Quarter)
plot(meanTWquarter.im )

# What is the nature of the association between koala sight locations and total phosphorous?
meanTWquarter.rho <- rhohat(object = dat.ppp, covariate = meanTWquarter.im)
plot(meanTWquarter.rho, xlab = "Mean_Temperature_of_Warmest_Quarter(units)", main = "")

#Mean_Temperature_of_Wettest_Quarter ==================================================================================================================================
meanTWetestQ.im <- as.im(habitat.rr$Mean_Temperature_of_Wettest_Quarter)
plot(meanTWetestQ.im )

# What is the nature of the association between koala sight locations and total phosphorous?
meanTWetestQ.rho <- rhohat(object = dat.ppp, covariate = meanTWetestQ.im)
plot(meanTWetestQ.rho, xlab = "Mean_Temperature_of_Wettest_Quarter(units)", main = "")

#Min_Temperature_of_Coldest_Month ==================================================================================================================================
minTcoldestM.im <- as.im(habitat.rr$Min_Temperature_of_Coldest_Month)
plot(minTcoldestM.im )

# What is the nature of the association between koala sight locations and total phosphorous?
minTcoldestM.rho <- rhohat(object = dat.ppp, covariate = minTcoldestM.im)
plot(minTcoldestM.rho, xlab = "Min_Temperature_of_Coldest_Month(units)", main = "")

#Precipitation_of_Coldest_Quarter ==================================================================================================================================
PreColdest_Q.im <- as.im(habitat.rr$Precipitation_of_Coldest_Quarter)
plot(PreColdest_Q.im )

# What is the nature of the association between koala sight locations and total phosphorous?
PreColdest_Q.rho <- rhohat(object = dat.ppp, covariate = PreColdest_Q.im)
plot(PreColdest_Q.rho, xlab = "Precipitation_of_Coldest_Quarter(units)", main = "")


#Precipitation_of_Driest_Month ==================================================================================================================================
PreDreistM.im <- as.im(habitat.rr$Precipitation_of_Driest_Month)
plot(PreDreistM.im )

# What is the nature of the association between koala sight locations and total phosphorous?
PreDreistM.rho <- rhohat(object = dat.ppp, covariate = PreDreistM.im)
plot(PreDreistM.rho, xlab = "Precipitation_of_Driest_Month(units)", main = "")

#Precipitation_of_Driest_Quarter ==================================================================================================================================
PreDreistQ.im <- as.im(habitat.rr$Precipitation_of_Driest_Quarter)
plot(PreDreistQ.im )

# What is the nature of the association between koala sight locations and total phosphorous?
PreDreistQ.rho <- rhohat(object = dat.ppp, covariate = PreDreistQ.im)
plot(PreDreistQ.rho, xlab = "Precipitation_of_Driest_Quarter(units)", main = "")

#Precipitation_of_Warmest_Quarter ==================================================================================================================================
PreWarmestQ.im <- as.im(habitat.rr$Precipitation_of_Warmest_Quarter)
plot(PreWarmestQ.im )

# What is the nature of the association between koala sight locations and total phosphorous?
PreWarmestQ.rho <- rhohat(object = dat.ppp, covariate = PreWarmestQ.im)
plot(PreWarmestQ.rho, xlab = "Precipitation_of_Warmest_Quarter(units)", main = "")

#Precipitation_of_Wettest_Month ==================================================================================================================================
PreWettestM.im <- as.im(habitat.rr$Precipitation_of_Wettest_Month)
plot(PreWettestM.im )

# What is the nature of the association between koala sight locations and total phosphorous?
PreWettestM.rho <- rhohat(object = dat.ppp, covariate = PreWettestM.im)
plot(PreWettestM.rho, xlab = "Precipitation_of_Wettest_Month(units)", main = "")

#Precipitation_of_Wettest_Quarter ==================================================================================================================================
PreWettestQ.im <- as.im(habitat.rr$Precipitation_of_Wettest_Quarter)
plot(PreWettestQ.im )

# What is the nature of the association between koala sight locations and total phosphorous?
PreWettestQ.rho <- rhohat(object = dat.ppp, covariate = PreWettestQ.im)
plot(PreWettestQ.rho, xlab = "Precipitation_of_Wettest_Quarter(units)", main = "")

#Precipitation_Seasonality ==================================================================================================================================
PreSeason.im <- as.im(habitat.rr$Precipitation_Seasonality)
plot(PreSeason.im )

# What is the nature of the association between koala sight locations and total phosphorous?
PreSeason.rho <- rhohat(object = dat.ppp, covariate = PreSeason.im)
plot(PreSeason.rho, xlab = "Precipitation_Seasonality(units)", main = "")

#Temperature_Annual_Range ==================================================================================================================================
TempAnnuRange.im <- as.im(habitat.rr$Temperature_Annual_Range)
plot(TempAnnuRange.im )

# What is the nature of the association between koala sight locations and total phosphorous?
TempAnnuRange.rho <- rhohat(object = dat.ppp, covariate = TempAnnuRange.im)
plot(TempAnnuRange.rho, xlab = "Temperature_Annual_Range(units)", main = "")

#Temperature_Seasonality ==================================================================================================================================
Tempseasonality.im <- as.im(habitat.rr$Temperature_Seasonality)
plot(Tempseasonality.im )

# What is the nature of the association between koala sight locations and total phosphorous?
Tempseasonality.rho <- rhohat(object = dat.ppp, covariate = Tempseasonality.im)
plot(Tempseasonality.rho, xlab = "Temperature_Seasonality(units)", main = "")



### Step 03: Soil rasters #####

myfullstack.d <- list.files(path="C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/warton_data_allclimate/soil_habitat_vars",pattern="\\.tif$",full.names=TRUE)
myfullstack.d = stack(myfullstack.d)
plot(myfullstack.d,2)
names(myfullstack.d)
# myextent=c(387900, 553100, 6862400, 7113600) 
myextent=studyarea.shp
habitat.rr=crop(myfullstack.d, myextent, snap="near")

plot(habitat.rr)

#aspect91 ==================================================================================================================================
aspect91.im <- as.im(habitat.rr$aspect91)
plot(aspect91.im )

# What is the nature of the association between koala sight locations and total phosphorous?
aspect91.rho <- rhohat(object = dat.ppp, covariate = aspect91.im)
plot(aspect91.rho, xlab = "aspect91(units)", main = "")

#awc ==================================================================================================================================
awc.im <- as.im(habitat.rr$awc)
plot(awc.im )

# What is the nature of the association between koala sight locations and total phosphorous?
awc.rho <- rhohat(object = dat.ppp, covariate = awc.im)
plot(awc.rho, xlab = "awc(units)", main = "")

#clay ==================================================================================================================================
clay.im <- as.im(habitat.rr$clay)
plot(clay.im )

# What is the nature of the association between koala sight locations and total phosphorous?
clay.rho <- rhohat(object = dat.ppp, covariate = clay.im)
plot(clay.rho, xlab = "clay(units)", main = "")

#elev ==================================================================================================================================
elev.im <- as.im(habitat.rr$elev)
plot(elev.im )

# What is the nature of the association between koala sight locations and total phosphorous?
elev.rho <- rhohat(object = dat.ppp, covariate = elev.im)
plot(elev.rho, xlab = "elev(units)", main = "")

#fpcnew ==================================================================================================================================
fpcnew.im <- as.im(habitat.rr$fpcnew)
plot(fpcnew.im )

# What is the nature of the association between koala sight locations and total phosphorous?
fpcnew.rho <- rhohat(object = dat.ppp, covariate = fpcnew.im)
plot(fpcnew.rho,xlim = c(10, 90), xlab = "fpcnew(units)", main = "")

##habit1decimal ==================================================================================================================================
habit1decimal.im <- as.im(habitat.rr$habit1decimal)
plot(habit1decimal.im )

# What is the nature of the association between koala sight locations and total phosphorous?
habit1decimal.rho <- rhohat(object = dat.ppp, covariate = habit1decimal.im)
plot(habit1decimal.rho, xlab = "habit1decimal(units)", main = "") 

##habit2decimal ==================================================================================================================================
habit2decimal.im <- as.im(habitat.rr$habit2decimal)
plot(habit2decimal.im )

# What is the nature of the association between koala sight locations and total phosphorous?
habit2decimal.rho <- rhohat(object = dat.ppp, covariate = habit2decimal.im)
plot(habit2decimal.rho, xlab = "habit2decimal(units)", main = "") 

##habit3decimal ==================================================================================================================================
habit3decimal.im <- as.im(habitat.rr$habit3decimal)
plot(habit3decimal.im )

# What is the nature of the association between koala sight locations and total phosphorous?
habit3decimal.rho <- rhohat(object = dat.ppp, covariate = habit3decimal.im)
plot(habit3decimal.rho, xlab = "habit3decimal(units)", main = "") 

##hpop ==================================================================================================================================
hpop.im <- as.im(habitat.rr$hpop)
plot(hpop.im )

# What is the nature of the association between koala sight locations and total phosphorous?
hpop.rho <- rhohat(object = dat.ppp, covariate = hpop.im)
plot(hpop.rho,xlim=c(3000,7000), xlab = "hpop(units)", main = "") 

##lot_density ==================================================================================================================================
lot_density.im <- as.im(habitat.rr$lot_density)
plot(lot_density.im )

# What is the nature of the association between koala sight locations and total phosphorous?
lot_density.rho <- rhohat(object = dat.ppp, covariate = lot_density.im)
plot(lot_density.rho,xlim = c(0, 1000),ylim = c(0, 0.000001), xlab = "lot_density(units)", main = "") 

##nitro ==================================================================================================================================
nitro.im <- as.im(habitat.rr$nitro)
plot(nitro.im )

# What is the nature of the association between koala sight locations and total phosphorous?
nitro.rho <- rhohat(object = dat.ppp, covariate = nitro.im)
plot(nitro.rho, xlab = "nitro(units)", main = "") 

##sbd ==================================================================================================================================
sbd.im <- as.im(habitat.rr$sbd)
plot(sbd.im )

# What is the nature of the association between koala sight locations and total phosphorous?
sbd.rho <- rhohat(object = dat.ppp, covariate = sbd.im)
plot(sbd.rho, xlab = "sbd(units)", main = "") 

##tpo ==================================================================================================================================
tpo.im <- as.im(habitat.rr$tpo)
plot(tpo.im )

# What is the nature of the association between koala sight locations and total phosphorous?
tpo.rho <- rhohat(object = dat.ppp, covariate = tpo.im)
plot(tpo.rho, xlab = "tpo(units)", main = "") 

##twi ==================================================================================================================================
twi.im <- as.im(habitat.rr$twi)
plot(twi.im )

# What is the nature of the association between koala sight locations and total phosphorous?
twi.rho <- rhohat(object = dat.ppp, covariate = twi.im)
plot(twi.rho, xlab = "twi(units)", main = "") 

####### density variables has no relationship.

myfullstack.d <- list.files(path="C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/warton_data_allclimate/density_varProxcy",pattern="\\.tif$",full.names=TRUE)
myfullstack.d = stack(myfullstack.d)
plot(myfullstack.d,1)
names(myfullstack.d)
# myextent=c(387900, 553100, 6862400, 7113600) 
myextent=studyarea.shp
habitat.rr=crop(myfullstack.d, myextent, snap="near")
##hpop ==================================================================================================================================
hpop.im <- as.im(habitat.rr$hpop)
plot(hpop.im )

# What is the nature of the association between koala sight locations and total phosphorous?
hpop.rho <- rhohat(object = dat.ppp, covariate = hpop.im)
plot(hpop.rho,xlim=c(3000,7000), xlab = "hpop(units)", main = "") 

##lot ==================================================================================================================================
lot.im <- as.im(habitat.rr$lot_density)
plot(lot.im )

# What is the nature of the association between koala sight locations and total phosphorous?
lot.rho <- rhohat(object = dat.ppp, covariate = lot.im)
plot(lot.rho, xlab = "lot(units)", main = "") 


##roadsM ==================================================================================================================================
roadsM.im <- as.im(habitat.rr$roads_motor)
plot(roadsM.im )

# What is the nature of the association between koala sight locations and total phosphorous?
roadsM.rho <- rhohat(object = dat.ppp, covariate = roadsM.im)
plot(roadsM.rho, xlab = "roadsM(units)", main = "") 

##roadsM ==================================================================================================================================
roads_other.im <- as.im(habitat.rr$roads_other)
plot(roads_other.im )

# What is the nature of the association between koala sight locations and total phosphorous?
roads_other.rho <- rhohat(object = dat.ppp, covariate = roads_other.im)
plot(roads_other.rho, xlab = "twi(units)", main = "") 


#List enveromental variables which has a ralationship.
# Annual_Mean_Temperature
# Annual_Precipitation
# Max_Temperature_of_Warmest_Month
# Mean_Diurnal_Range_
# Mean_Temperature_of_Coldest_Quarter
# Mean_Temperature_of_Driest_Quarter
# Mean_Temperature_of_Warmest_Quarter
# Mean_Temperature_of_Wettest_Quarter
# Min_Temperature_of_Coldest_Month
# Precipitation_of_Coldest_Quarter
# Precipitation_of_Driest_Month
# Precipitation_of_Driest_Quarter
# Precipitation_of_Warmest_Quarter
# Precipitation_of_Wettest_Month
# Precipitation_of_Wettest_Quarter
# Precipitation_Seasonality
# Temperature_Annual_Range
# Temperature_Seasonality


#list of soil variables which has relatiohsip
#clay
#aspect91
#elev
# habit1decimal
# habit2decimal
# habit3decimal
# lot_density
# nitro
# tpo


#check coliniarity# put all variables together to a folder

setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/vif_check_env_soil")

myfullstack.v <- list.files(path="C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/vif_check_env_soil",pattern="\\.tif$",full.names=TRUE)
myfullstack.v = stack(myfullstack.v)
plot(myfullstack.v,1)
names(myfullstack.v)
# myextent=c(387900, 553100, 6862400, 7113600) 
myextent=studyarea.shp
habitat.rr=crop(myfullstack.v, myextent, snap="near")
plot(habitat.rr$Annual_Precipitation)
vifstep(habitat.rr, th=10)




# Variables      VIF
# 1                                clay 2.100519
# 2                                elev 7.426360
# 3                       habit1decimal 1.108261
# 4                       habit2decimal 1.128784
# 5                       habit3decimal 1.046693
# 6                         lot_density 1.104974
# 7  Mean_Temperature_of_Driest_Quarter 9.072536
# 8                               nitro 2.998084
# 9       Precipitation_of_Driest_Month 5.241322
# 10          Precipitation_Seasonality 3.219922
# 11            Temperature_Seasonality 4.268737
# 12                                tpo 2.525569



##plot distance variables 
tiff("distance-vars2.tiff",width=13,height=12,units="in",res=300) # 8:4 for three mpas. one map 8:10 # 
par(mfrow=c(3,3),oma=c(.1, .1, .1, .1),mai = c(.5, .5, .5, .5),mar=c(4,5,4,7))
plot(myfullstack.b$distance_primaryandlink,cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2, main = "distance_primary",legend.width=2)
plot(myfullstack.b$distance_secondaryandlink,cex.lab=1.5, cex.axis=1.5,cex.lab=2,cex.main=2,main = "distance_secondary",legend.width=2 )
plot(myfullstack.b$distance_tertiaryandlink, cex.lab=1.5, cex.axis=1.5,cex.lab=2,cex.main=2,main = "distance_tertiary",legend.width=2)
plot(myfullstack.b$distance_residentil, cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2,main = "distance_residential",legend.width=2)
plot(myfullstack.b$distance_unclassified, cex.lab=1.5, cex.axis=1.5,cex.lab=2,cex.main=2,main = "distance_unclassified",legend.width=2)
plot(myfullstack.b$distance_motorwayandlink, cex.lab=1.5, cex.axis=1.5,cex.lab=2,cex.main=2, main = "distance_motorway",legend.width=2)
# graphics.off()
# dev.off()
# 
# jpeg("habitat_vars.jpeg",width=15,height=5.5,units="in",res=600) # 8:4 for three mpas. one map 8:10 # 
# par(mfrow=c(1,3),oma=c(.1, .1, .1, .1),mai = c(.5, .5, .5, .5),mar=c(4,5,4,7))

plot(myfullstack.b$Dis_habitat_suitable_1,cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2, main = "distance_high_suitable",legend.width=2)
plot(myfullstack.b$Dis_habitat_suitable_2,cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2, main = "distance_suitable",legend.width=2)
plot(myfullstack.b$Dis_habitat_suitable_3,cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2, main = "distance_low suitable",legend.width=2)
graphics.off()
dev.off()

levelplot(d, par.settings = BuRdTheme)



setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/Data_raw_spatial/Habitat percentages rasters")

habitat1 <- raster("habit11.tif")
habitat2 <- raster("habit21.tif")
habitat3 <- raster("habit31.tif")
habitat0 <- raster("habit01.tif")

#Rastervis package.
habit <- stack(habitat1,habitat2,habitat3)# BTCTheme viridisTheme plasmaTheme
levelplot(habit, par.settings = BuRdTheme)#RdBuTheme YlOrRdTheme BuRdTheme GrTheme
setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/Paper one abundance/Figures_paper on abundace")



jpeg("habitat_vars.jpeg",width=18,height=5.5,units="in",res=600) # 8:4 for three mpas. one map 8:10 # 
par(mfrow=c(1,4),oma=c(.1, .1, .1, .1),mai = c(.5, .5, .5, .5),mar=c(4,5,4,6))


plot(habitat1, cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2, main = "highly_suitable",legend.width=4)
plot(habitat2, cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2, main = "suitable",legend.width=4)
plot(habitat3, cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2, main = "low_suitable",legend.width=4)
plot(habitat0, cex.lab=1.5, cex.axis=1.5,cex.lab=2, cex.main=2, main = "unsuitable",legend.width=4)

graphics.off()
dev.off()
