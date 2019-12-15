##############################################################################
#
# Suspended Sediment Yield Top-Kriging
#
# Anatoly Tsyplenkov
# atsyplenkov@gmail.com
#
##############################################################################

library(rtop)
library(rgdal)
library(tidyverse)
set.seed(1)

source("/WORK/99_SOFT/rtop/R_rtop_Tutorial/Efficiencies.R")

# 0) Read data sets
observations = readOGR("data/spatial/utm", "terek_ws-utm")
predictionLocations = readOGR("data/spatial/utm", "terek_predictions-utm")

# 1) Standardize values
# (to remove their dependency on catchment size)
plot(observations$AREA, observations$R,
     xlab = "Area (km^2)",
     ylab = "Suspended Sediment Discharge, kg/s", 
     log = "xy")
# Linear regression of log-transformed variables
lm.area <- lm(log(observations$R)~log(observations$AREA))
summary(lm.area)
# Variable to predict
observations$obs <- observations$R/(observations$AREA^lm.area$coefficients[2])
# the variable "obs" will be used by TK routines (see below)
# Dependence on catchment area has been removed:
plot(observations$AREA, observations$obs,
     xlab = "Area (km^2)",
     ylab = "Standardized SS discharge (kg/s/km^2)",
     log = "xy")


# 2) Set parameters and create rtop object
vic <- 6 # neighbourhood (no. of considered catchments)
# Create an rtop object
rtopObj.R <- createRtopObject(observations = observations,
                              predictionLocations = observations,
                              formulaString = obs~1,    # variable to krige
                              params=list(gDist=TRUE, # how to compute distances between areas
                                          rresol=500, # catchment area discretization
                                          nmax=vic,   # no. of neighbours
                                          cv=TRUE,    # cross-validation mode: ON
                                          partialOverlap=T)) # for inaccurate catchment boundaries
#
# 3) Build the empirical variogram
rtopObj.R <- rtopVariogram(rtopObj.R)    # create variogram for the rtop object
rtopObj.R <- rtopFitVariogram(rtopObj.R) # fit theoretical variogram 
rtopObj.R <- checkVario(rtopObj.R)       # check variogram fitting
#
# 4) Perform the prediction in cross-validation
rtopObj.R <- rtopKrige(rtopObj.R)
#
# 5) Analyze predictions
# Prediction of standardized quantiles:
Predictions <- rtopObj.R$predictions$var1.pred # in cross-validation
# TK predictions in cross-validation in cubic meters per second
Predicted.R <- Predictions*(observations$AREA^lm.area$coefficients[2])
# Nash and Sutcliffe Efficiency
NSE(observations$R, Predicted.R) 
# Scatterplot
plot(observations$R, Predicted.R,
     xlab="At-site SS Discharge (kg/s)",
     ylab="Predictions (kg/s)",log = "xy")
abline(a = 0,b = 1)
#
# 6) Predict to ungauged basins



##############################################################################

# Create an rtop-object and fit a variogram model
rtopObj = createRtopObject(observations, predictionLocations,
                           formulaString = obs~1,
                           params = list(gDist = TRUE, rresol = 25))

rtopObj = rtopFitVariogram(rtopObj)

# Exploratory data analyses and visualization of model fit
rtopObj = checkVario(rtopObj, acor = 0.000001, acomp =
                       data.frame(acl1 = c(2,2,2,2,3,3,3,4,4), acl2 = c(2,3,4,5,3,4,5,4,5)))

# Cross-validation
rtopObj = rtopKrige(rtopObj,cv=TRUE)
predictions = rtopObj$predictions
summary(predictions)
sstot = sum((predictions$observed-mean(predictions$observed))^2)
rtopsserr = sum((predictions$observed- predictions$var1.pred)^2)
rtoprsq = 1-rtopsserr/sstot
rtoprsq

# Predictions at new locations - visualize the result on a stream network
rtopObj = rtopKrige(rtopObj)
rnet = readOGR("data/spatial/utm", "terek_river-network-utm")
pred = rtopObj$predictions
# rnet$pred = pred$var1.pred[match(rnet$EZGA, pred$EZGID)]
spplot(pred, "var1.pred", col.regions = bpy.colors())
at = seq(0,max(rnet$pred,na.rm = TRUE),0.01)
cols = bpy.colors(length(at))
cobs = observations@data[,c("XSTATION", "YSTATION", "obs")]
names(cobs) = c("x","y","obs")
coordinates(cobs) = ~x+y
cobs$class = findInterval(cobs$obs, at)

spplot(rnet,"pred",col.regions = bpy.colors(), at = at, panel = function(x,y, ...){
  panel.polygonsplot(x,y, ...)
  sp.points(cobs[,"obs"], cex=1, pch = 16, col = cols[cobs$class])
})
