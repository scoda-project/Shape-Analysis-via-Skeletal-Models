#####################################################################################################
#####################################################################################################
# Libraries

library(shapes)
library(rgl)
library(matlib)
library(RiemBase)
library(plotrix)
library(Rvcg)
library(Matrix)
library(Morpho)
library(fields)
library(ggplot2)
library(rotations)

#####################################################################################################
#####################################################################################################

#clear the environment
remove(list=ls())

#####################################################################################################
#####################################################################################################
#set working directory to file location

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#####################################################################################################
#####################################################################################################
# load functions
if(TRUE){
  source("subFunctions/MathFunctions.R")
  source("subFunctions/euclideanization.R")
  source("subFunctions/ray_triangle_intersection.R")
  source("subFunctions/cutAndStrechSpokes.R")
  source("subFunctions/frameGenerator.R")
  source("subFunctions/normalsOfSkeletalSheetByTriangles.R")
  source("subFunctions/normalsOfSkeletalSheetBySpline.R")
  source("subFunctions/rotateFrameForwardAndBackward.R")
  source("subFunctions/meanFrames.R")
  source("subFunctions/readSrepsData.R")
  source("subFunctions/kurtosisTestFunction.R")
}


#####################################################################################################
#####################################################################################################
# choose the type of study

#choose the type of study "sizeAndShapeAnalysis" or "shapeAnalysis"
typeOfStudy<-"shapeAnalysis"              #removing scale
# typeOfStudy<-"sizeAndShapeAnalysis"     #preserving scale

# choose euclideanization method typeOfStudy4directions as "PNS" or "tangent space"
# PNS takes more time!
# typeOfStudy4directions<-"tangent space"
typeOfStudy4directions<-"PNS"

#choose type of mean direction
typeOfMeanDirection<-"Frechet"
# typeOfMeanDirection<-"PNS"

# choose typeOfTest as "Parametric" or "Permutation"
typeOfTest<-"Parametric"    # Fast Parametric test is based on normality assumption
# typeOfTest<-"Permutation" # nPerm default is 10000 permutations

##############################################################################
##############################################################################
# load SPHARM-PDMs

G1_path <- "../Files/SPHARM_Obj_G1"
G1_list <- list.files(path = G1_path, pattern = "\\.obj$")
nSamplesG1 <- length(G1_list)
nSamplesG1

G2_path <- "../Files/SPHARM_Obj_G2"
G2_list <- list.files(path = G2_path, pattern = "\\.obj$")
nSamplesG2 <- length(G2_list)
nSamplesG2

# Extracting boundary points of the SPHARM-PDM meshes
# Each mesh consists of 1002 boundary points

#G1 case
nSamplesG1 <- length(G1_list)
data_list <- vector(mode="list",length=nSamplesG1)
data_obj1 <- vcgImport(paste0(G1_path,"/",G1_list[1]))
boundaryPDM_G1 <- array(NA, dim=c(nrow(vert2points(data_obj1)),3,nSamplesG1))
for(i in 1:nSamplesG1){
  data_obj <- vcgImport(paste0(G1_path,"//",G1_list[i]))
  boundaryPDM_G1[,,i] <- vert2points(data_obj)
}

#G2 case
nSamplesG2 <- length(G2_list)
cont_data_list <- vector(mode="list",length=nSamplesG2)
boundaryPDM_G2 <- array(NA, dim=c(nrow(vert2points(data_obj1)),3,nSamplesG2))
for(i in 1:nSamplesG2){
  data_obj <- vcgImport(paste0(G2_path,"/",G2_list[i] ))
  boundaryPDM_G2[,,i] <- vert2points(data_obj)
}

# Plot a sample from G1
if(TRUE){
  sampleNumber<-1
  open3d()
  verts <- rbind(t(as.matrix(boundaryPDM_G1[,,sampleNumber])),1)
  trgls <- as.matrix(data_obj1$it)
  tmesh <- tmesh3d(verts, trgls)
  shade3d(tmesh, col="white",alpha=1)  #surface mesh
  wire3d(tmesh, col="lightgrey")  #surface mesh
  # plot only 300 points
  plot3d(boundaryPDM_G1[1:300,,sampleNumber],type="s",radius = 0.4 ,col = "black",expand = 10,box=FALSE,add = TRUE)
  rglwidget()
}

# plot all SPHARM-PDMs
if(TRUE){
  open3d()
  for (i in 1:nSamplesG1) {
    plot3d(boundaryPDM_G1[,,i],type="p" ,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in 1:nSamplesG1) {
    plot3d(boundaryPDM_G2[,,i],type="p" ,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  rglwidget()
}


# Align all boundary points by Generalized Procrustes Analysis (GPA)
# ⚠️ Note that we DO NOT need aligned meshes for DSRep model fitting or analysis
# In this R script the act of alignment is only to provide better visualization
# Thus, the whole procedure can be done without alignments
# Besides, the LPDSRep analysis supports both shapeAnalysis and 
# sizeAndShapeAnalysis regardless of the act of alignment

if(TRUE){
  all_boundaryPoints<-abind(boundaryPDM_G1,boundaryPDM_G2)
  if(typeOfStudy=="shapeAnalysis"){
    proc<-procGPA(all_boundaryPoints, scale = TRUE) 
  }else if(typeOfStudy=="sizeAndShapeAnalysis") {
    proc<-procGPA(all_boundaryPoints, scale = FALSE) 
  }else{
    stop("Please specify the type of study!")
  }
  aligned_SpharmPDM_G1<-proc$rotated[,,1:nSamplesG1]
  aligned_SpharmPDM_G2<-proc$rotated[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)]
}


# plot aligned PDMs
if(TRUE){
  open3d()
  for (i in 1:nSamplesG1) {
    plot3d(aligned_SpharmPDM_G1[,,i],type="p" ,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in 1:nSamplesG2) {
    plot3d(aligned_SpharmPDM_G2[,,i],type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(main = "Proc aligned SPHARM-PDMs of G1 and G2 in blue and red")
  rglwidget()
}

#####################################################################################################
#####################################################################################################
# SPHARM-PDM and DSRep of an eccentric ellipsoid

meshPoints <- read.csv(file = paste("../Files/ellipsoid.csv",sep = ""),check.names = FALSE, header=TRUE, sep=",")
# Creating a matrix of the mesh points with 3 columns 
ellipsoid_SpharmPDM<-matrix(meshPoints[[1]], ncol = 3, byrow = TRUE)

open3d()
verts <- rbind(t(as.matrix(ellipsoid_SpharmPDM)),1)
trgls <- as.matrix(data_obj1$it)
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="white",alpha=0.2)  #surface mesh
wire3d(tmesh, col="lightgrey")  #surface mesh
rglwidget()


# Ellipsoid DSRep
srepsDataEllipsoid<- read.csv(file=paste("../Files/ellipsoid_Skeleton.csv",sep = ""), header=TRUE, sep=",")

upSpoeksNumber<-max(srepsDataEllipsoid$SpokesNumber[which(srepsDataEllipsoid$srepNumber==1 & srepsDataEllipsoid$Spoke=='up')])
downSpoeksNumber<-max(srepsDataEllipsoid$SpokesNumber[which(srepsDataEllipsoid$srepNumber==1 & srepsDataEllipsoid$Spoke=='down')])
crestSpoksNumber<-max(srepsDataEllipsoid$SpokesNumber[which(srepsDataEllipsoid$srepNumber==1 & srepsDataEllipsoid$Spoke=='crest')])
nTotalRadii <- upSpoeksNumber + downSpoeksNumber + crestSpoksNumber
skelPointNo <- nTotalRadii-downSpoeksNumber
skelRange<-c(1:downSpoeksNumber,(2*downSpoeksNumber+1):nTotalRadii)

# BoundaryPDM and SkeletalPDM
tempEllipsoid<-readSrepsData(srepsData = srepsDataEllipsoid)
SkeletalPDMEllipsoid<-tempEllipsoid$SkeletalPDM[,,1]
BoundaryPDMEllipsoid<-tempEllipsoid$BoundaryPDM[,,1]
boundaryPlusSkeletal_Ellipsoid<-tempEllipsoid$boundaryPlusSkeletal


# plot ellipsoid DSRep
if(TRUE){
  srep1<-rbind(SkeletalPDMEllipsoid,BoundaryPDMEllipsoid)
  plot3d(SkeletalPDMEllipsoid[skelRange,],type="s", size=0.5,col = "red",expand = 10,box=FALSE,add = TRUE)
  plot3d(BoundaryPDMEllipsoid[103:122,],type="s", radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  rglwidget()
}


############################################################################################
################ Fit DSRep via Boundary Deformation Method #################################
############################################################################################
# Fit DSRep via thin plate spline from Morpho::tps3d

skelPlusBoundary_ellipsoid<-rbind(SkeletalPDMEllipsoid,BoundaryPDMEllipsoid)


skelPlusBoundary_G1<-array(NA,dim = c(dim(skelPlusBoundary_ellipsoid),nSamplesG1))
for (i in 1:nSamplesG1) {                       # Matrix info of vertices
  skelPlusBoundary_G1[,,i]<-Morpho::tps3d(x = skelPlusBoundary_ellipsoid,
                                          refmat = ellipsoid_SpharmPDM,   # Reference matrix
                                          tarmat = aligned_SpharmPDM_G1[,,i]) # Target matrix
}
skelPlusBoundary_G2<-array(NA,dim = c(dim(skelPlusBoundary_ellipsoid),nSamplesG2))
for (i in 1:nSamplesG2) {
  skelPlusBoundary_G2[,,i]<-Morpho::tps3d(x = skelPlusBoundary_ellipsoid,
                                          refmat = ellipsoid_SpharmPDM,
                                          tarmat = aligned_SpharmPDM_G2[,,i])
}



#####################################################################################################
#####################################################################################################
# Extract DSRep info including the spokes lengths and global directions

# Separating Skeletal and Boundary data from skelPlusBoundary_G1
SkeletalPDMG1<-skelPlusBoundary_G1[1:nTotalRadii,,]
BoundaryPDMG1<-skelPlusBoundary_G1[(nTotalRadii+1):dim(skelPlusBoundary_G1)[1],,]
SkeletalPDMG2<-skelPlusBoundary_G2[1:nTotalRadii,,]
BoundaryPDMG2<-skelPlusBoundary_G2[(nTotalRadii+1):dim(skelPlusBoundary_G1)[1],,]

#####################################################################################################
#####################################################################################################
# Plot a fitted DSRep

#choose sampleNo between 1 to nSamplesG1 to see other ds-reps
sampleNo<-1
#plot
if(TRUE){
  open3d()
  srep1<-rbind(SkeletalPDMG1[,,sampleNo],BoundaryPDMG1[,,sampleNo])
  plot3d(SkeletalPDMG1[skelRange,,sampleNo],type="s", size=0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  # plot mesh and normal vectors of a sample
  verts <- rbind(t(as.matrix(aligned_SpharmPDM_G1[,,sampleNo])),1)
  trgls <- as.matrix(data_obj1$it)
  tmesh <- tmesh3d(verts, trgls)
  shade3d(tmesh, col="white",alpha=0.2)  #surface mesh
  wire3d(tmesh, col="lightgrey")  #surface mesh
  # decorate3d(xlab = "x", ylab = "y", zlab = "z",
  #            box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
}

#####################################################################################################
#####################################################################################################
# Check correspondence

# Plot to see spoke correspondence
# We see spokes' tip and tails in separated colors
if(TRUE){
  open3d()
  for (k in 1:nSamplesG1) {
    plot3d(SkeletalPDMG1[1:upSpoeksNumber,,k],type="p",col = "yellow",expand = 10,box=FALSE,add = TRUE)
    plot3d(BoundaryPDMG1[1:upSpoeksNumber,,k],type="p",col = "orange",expand = 10,box=FALSE,add = TRUE)
    plot3d(BoundaryPDMG1[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),,k],type="p",col = "green",expand = 10,box=FALSE,add = TRUE)
    plot3d(SkeletalPDMG1[(2*upSpoeksNumber+1):nTotalRadii,,k],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  open3d()
  for (k in 1:nSamplesG2) {
    plot3d(SkeletalPDMG2[1:upSpoeksNumber,,k],type="p",col = "yellow",expand = 10,box=FALSE,add = TRUE)
    plot3d(BoundaryPDMG2[1:upSpoeksNumber,,k],type="p",col = "orange",expand = 10,box=FALSE,add = TRUE)
    plot3d(BoundaryPDMG2[(upSpoeksNumber+1):(upSpoeksNumber+downSpoeksNumber),,k],type="p",col = "green",expand = 10,box=FALSE,add = TRUE)
    plot3d(SkeletalPDMG2[(2*upSpoeksNumber+1):nTotalRadii,,k],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  rglwidget()
}

#####################################################################################################
#####################################################################################################
# Extract DSRep features

# Calculate spokes' lengths (radii)
radii_G1<-array(NA, dim=c(nTotalRadii,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:nTotalRadii) {
    radii_G1[i,k]<-norm(SkeletalPDMG1[i,,k]-
                          BoundaryPDMG1[i,,k],type = "2")
  }
}
radii_G2<-array(NA, dim=c(nTotalRadii,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:nTotalRadii) {
    radii_G2[i,k]<-norm(SkeletalPDMG2[i,,k]-
                          BoundaryPDMG2[i,,k],type = "2")
  }
}

# Calsulate spokes' directions (in global coordinate system)
spokeDirections_G1<-array(NA, dim=c(nTotalRadii,3,nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:nTotalRadii) {
    spokeDirections_G1[i,,k]<-convertVec2unitVec(BoundaryPDMG1[i,,k]-SkeletalPDMG1[i,,k])
  }
}
spokeDirections_G2<-array(NA, dim=c(nTotalRadii,3,nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:nTotalRadii) {
    spokeDirections_G2[i,,k]<-convertVec2unitVec(BoundaryPDMG2[i,,k]-SkeletalPDMG2[i,,k])
  }
}


# Plot all directions of an specific spoke

spokeNumber<-10 # pick a number in 1,...,122 as we have 122 spokes
if(TRUE){
  open3d()
  spheres3d(c(0,0,0),col="white",alpha=0.2)
  plot3d(t(spokeDirections_G1[spokeNumber,,]),type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
  plot3d(t(spokeDirections_G2[spokeNumber,,]),type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
  rglwidget()
}

# plot all the spokes' directions
if(TRUE){
  open3d()
  spheres3d(c(0,0,0),col="white",alpha=0.2)
  for (i in 1:dim(spokeDirections_G1)[1]) {
    plot3d(t(spokeDirections_G1[i,,]),type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
    plot3d(t(spokeDirections_G2[i,,]),type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  rglwidget()
}

# plot extrinsic and pns mean of an specific observation
spokeNumber<-10 # pick a number in 1,...,122 as we have 122 spokes
smallPns<-pns(spokeDirections_G1[spokeNumber,,],sphere.type = "small")
greatPns<-pns(spokeDirections_G1[spokeNumber,,],sphere.type = "great")
extrinsicMean<-colMeans(t(spokeDirections_G1[spokeNumber,,]))/norm(colMeans(t(spokeDirections_G1[spokeNumber,,])),type = "2")
open3d()
plot3d(t(spokeDirections_G1[spokeNumber,,]),type="p",col = "blue",expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(smallPns$PNS$mean,smallPns$PNS$mean),type="s",radius = 0.01,col = "orange",expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(greatPns$PNS$mean,greatPns$PNS$mean),type="s",radius = 0.01,col = "green",expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(extrinsicMean,extrinsicMean),type="s",radius = 0.01,col = "grey",expand = 10,box=TRUE,add = TRUE)
drawCircleS2(center = smallPns$PNS$orthaxis[[1]],theta = asin(smallPns$PNS$radii[[2]]))
drawCircleS2(center = greatPns$PNS$orthaxis[[1]],theta = asin(greatPns$PNS$radii[[2]]))
rglwidget()

# To choose between small and great circle we can use 
sphereType<-kurtosisTestFunction(spokeDirections_G1[spokeNumber,,])
sphereType

# plot the euclideanized data by pns small circle
plot(t(smallPns$resmat),xlim = c(-2,2),ylim = c(-0.5,0.5))
# plot the euclideanized data by pns great circle
plot(t(greatPns$resmat),xlim = c(-2,2),ylim = c(-0.5,0.5))



################################################################################
############################# LPDSRep ##########################################
################################################################################
# Fitting LPDSRep

# Define labels of the frames for the grid of the skeletal sheet
# NB frame 16 is its own parent
framesCenters   <-c(16,13,10,7 ,4 ,1 ,2 ,3 ,19,22,25,28,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,29,30,26,27,23,24,20,21,17,18,14,15,11,12,8 ,9 ,5 ,6 ,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
framesParents   <-c(16,16,13,10,7 ,4 ,1 ,2 ,16,19,22,25,28,31,32,28,34,25,36,22,38,19,40,16,42,13,44,10,46,7 ,48,4 ,50,28,29,25,26,22,23,19,20,16,17,13,14,10,11,7 ,8 ,4 ,5 ,3 ,6 ,9 ,12,15,18,21,24,27,30,33,35,37,39,41,43,45,47,49,51)
framesBackPoints<-c(13,16,13,10,7 ,4 ,1 ,2 ,16,19,22,25,28,31,32,28,34,25,36,22,38,19,40,16,42,13,44,10,46,7 ,48,4 ,50,28,29,25,26,22,23,19,20,16,17,13,14,10,11,7 ,8 ,4 ,5 ,3 ,6 ,9 ,12,15,18,21,24,27,30,33,35,37,39,41,43,45,47,49,51)
framesFronts    <-c(19,10,7 ,4 ,1 ,2 ,3 ,52,22,25,28,31,32,33,62,35,63,37,64,39,65,41,66,43,67,45,68,47,69,49,70,51,71,30,61,27,60,24,59,21,58,18,57,15,56,12,55,9 ,54,6 ,53,rep(Inf,20)) #NB! crest frames don't have front point

# number of frames
numberOfFrames<-length(framesCenters)

#####################################################################################################
#####################################################################################################
# Calculate normal vectors of the skeletal sheet


# Generate normal vectors based on the quadrilateral structure of the skeletal sheet.
if(TRUE){
  
  skeletalSheet_G1<-SkeletalPDMG1[skelRange,,]
  skeletalSheet_G2<-SkeletalPDMG2[skelRange,,]
  medialNormals_G1<-array(NA,dim = dim(skeletalSheet_G1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    setTxtProgressBar(pb, i/nSamplesG1)
    medialNormals_G1[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = skeletalSheet_G1[,,i])
  }
  print("Group 1 is done!")
  medialNormals_G2<-array(NA,dim = dim(skeletalSheet_G2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    setTxtProgressBar(pb, i/nSamplesG2)
    medialNormals_G2[,,i]<-normalsOfSkeletalSheetByTriangles(skeletalPDM = skeletalSheet_G2[,,i])
  }
  print("Group 2 is done!")
  
}

#####################################################################################################
#####################################################################################################
# Calculate frames vectors 

# frames in global coordinate system
if(TRUE){
  framesFirstVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesSecondVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  framesThirdVectors_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG1) {
    setTxtProgressBar(pb, i/nSamplesG1)
    temp<-frameGenerator(centeredSkel = skeletalSheet_G1[,,i],medialNormals = medialNormals_G1[,,i],
                         framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectors_G1[,,i]<-temp$framesFirstVec
    framesSecondVectors_G1[,,i]<-temp$framesSecondVec
    framesThirdVectors_G1[,,i]<-temp$framesThirdVec
  }
  close(pb)
  print("Group 1 is done!")
  framesFirstVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesSecondVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  framesThirdVectors_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  pb <- txtProgressBar(style = 3)
  for (i in 1:nSamplesG2) {
    setTxtProgressBar(pb, i/nSamplesG2)
    temp<-frameGenerator(centeredSkel = skeletalSheet_G2[,,i],medialNormals = medialNormals_G2[,,i],
                         framesCenters = framesCenters,framesBackPoints = framesBackPoints,framesFronts = framesFronts)
    
    framesFirstVectors_G2[,,i]<-temp$framesFirstVec
    framesSecondVectors_G2[,,i]<-temp$framesSecondVec
    framesThirdVectors_G2[,,i]<-temp$framesThirdVec
  }
  close(pb)
  print("Group 2 is done!")
}

#####################################################################################################
#####################################################################################################
# Plot an LPDSRep from each groups

# Choose sampleNo between 1 to nSamplesG1 to see other LPDSReps
sampleNo<-1
#plot
if(TRUE){
  open3d()
  for (i in 2:numberOfFrames) {
    vectors3d(skeletalSheet_G1[framesCenters[i],,sampleNo],origin = skeletalSheet_G1[framesParents[i],,sampleNo],
              headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in framesCenters) {
    vectors3d(skeletalSheet_G1[i,,sampleNo]+framesFirstVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(skeletalSheet_G1[i,,sampleNo]+framesSecondVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(skeletalSheet_G1[i,,sampleNo]+framesThirdVectors_G1[i,,sampleNo],origin = skeletalSheet_G1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  }
  for (i in 1:nTotalRadii) {
    vectors3d(SkeletalPDMG1[i,,sampleNo]+(spokeDirections_G1[i,,sampleNo]*radii_G1[i,sampleNo]),origin = SkeletalPDMG1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="red", lwd=1)
  }
  for (i in 1:nTotalRadii) {
    plot3d(rbind(SkeletalPDMG1[i,,sampleNo]+(spokeDirections_G1[i,,sampleNo]*radii_G1[i,sampleNo]),
                 SkeletalPDMG1[i,,sampleNo]),type="l",lwd = 2,col = "grey",expand = 10,box=FALSE,add = TRUE)
  }
  verts <- rbind(t(as.matrix(aligned_SpharmPDM_G1[,,sampleNo])),1)
  trgls <- as.matrix(data_obj1$it)
  tmesh <- tmesh3d(verts, trgls)
  # wire3d(tmesh, col=sample(1:10,size = 1),alpha=1)  #wire mesh
  shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  # decorate3d(xlab = "x", ylab = "y", zlab = "z",
  #            box = F, axes = TRUE, main = NULL, sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  # 
  rglwidget()
}

# Choose sampleNo between 1 to nSamplesG2 to see other LPDSReps
sampleNo<-1
#plot
if(TRUE){
  open3d()
  for (i in 2:numberOfFrames) {
    vectors3d(skeletalSheet_G2[framesCenters[i],,sampleNo],origin = skeletalSheet_G2[framesParents[i],,sampleNo],
              headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in framesCenters) {
    vectors3d(skeletalSheet_G2[i,,sampleNo]+framesFirstVectors_G2[i,,sampleNo],origin = skeletalSheet_G2[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(skeletalSheet_G2[i,,sampleNo]+framesSecondVectors_G2[i,,sampleNo],origin = skeletalSheet_G2[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
    vectors3d(skeletalSheet_G2[i,,sampleNo]+framesThirdVectors_G2[i,,sampleNo],origin = skeletalSheet_G2[i,,sampleNo],headlength = 0.1,radius = 1/10, col="orange", lwd=1)
  }
  for (i in 1:nTotalRadii) {
    vectors3d(SkeletalPDMG2[i,,sampleNo]+(spokeDirections_G2[i,,sampleNo]*radii_G2[i,sampleNo]),origin = SkeletalPDMG2[i,,sampleNo],headlength = 0.1,radius = 1/10, col="red", lwd=1)
  }
  for (i in 1:nTotalRadii) {
    plot3d(rbind(SkeletalPDMG2[i,,sampleNo]+(spokeDirections_G2[i,,sampleNo]*radii_G2[i,sampleNo]),
                 SkeletalPDMG2[i,,sampleNo]),type="l",lwd = 2,col = "grey",expand = 10,box=FALSE,add = TRUE)
  }
  verts <- rbind(t(as.matrix(aligned_SpharmPDM_G2[,,sampleNo])),1)
  trgls <- as.matrix(data_obj1$it)
  tmesh <- tmesh3d(verts, trgls)
  # wire3d(tmesh, col=sample(1:10,size = 1),alpha=1)  #wire mesh
  shade3d(tmesh, col="white",alpha=0.2)  #surface mech
  # decorate3d(xlab = "x", ylab = "y", zlab = "z",
  #            box = F, axes = TRUE, main = NULL, sub = NULL,
  #            top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
}

#####################################################################################################
#####################################################################################################
# Combine frames vectors to make SO(3) frames in global coordinate system

if(TRUE){
  frames_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in framesCenters) {
      frames_G1[,,i,k]<- rbind(framesFirstVectors_G1[i,,k],
                               framesSecondVectors_G1[i,,k],
                               framesThirdVectors_G1[i,,k])
    }
  }
  frames_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in framesCenters) {
      frames_G2[,,i,k]<- rbind(framesFirstVectors_G2[i,,k],
                               framesSecondVectors_G2[i,,k],
                               framesThirdVectors_G2[i,,k])
    }
  }
}

# plot a distribution of an specific frame in global coordinate system
frameNumber<-10 #choose a number in 1,...,71 as we have 71 frames
open3d()
spheres3d(c(0,0,0),col="white",alpha=0.2)
vectors3d(t(frames_G1[1,,frameNumber,]),headlength = 0.1,radius = 1/10, col="red", lwd=1)
vectors3d(t(frames_G1[2,,frameNumber,]),headlength = 0.1,radius = 1/10, col="blue", lwd=1)
vectors3d(t(frames_G1[3,,frameNumber,]),headlength = 0.1,radius = 1/10, col="green", lwd=1)
vectors3d(diag(3),headlength = 0.1,radius = 1/10, col="black", lwd=1)
rglwidget()

# Calculate the mean of an specific frame
frameNumber<-10 #choose a number in 1,...,71 as we have 71 frames
# Extrinsic mean from "rotations" library
framesVectorized<-array(NA,dim = c(9,nSamplesG1))
for (i in 1:nSamplesG1) {
  framesVectorized[,i]<-as.vector(t(frames_G1[,,frameNumber,i]))
}
tempVec<-mean(as.SO3(t(framesVectorized)),type = 'geometric')
meanFrame<-matrix(tempVec,nrow = 3,byrow = TRUE)

open3d()
vectors3d(t(frames_G1[1,,frameNumber,]),headlength = 0.1,radius = 1/100, col="yellow", lwd=1)
vectors3d(t(frames_G1[2,,frameNumber,]),headlength = 0.1,radius = 1/100, col="yellow", lwd=1)
vectors3d(t(frames_G1[3,,frameNumber,]),headlength = 0.1,radius = 1/100, col="yellow", lwd=1)
vectors3d(meanFrame,headlength = 0.1,radius = 1/20, col="blue", lwd=1)
vectors3d(diag(3),headlength = 0.1,radius = 1/20, col="black", lwd=1)
rglwidget()

# Exercise:
# A coordinate frame can be represented as a unit quaternion,
# i.e., a point on the 3-dimensional unit sphere S^3.
# Use Principal Nested Spheres (PNS) to compute the mean frame.

#####################################################################################################
#####################################################################################################
# Calculate children frames coordinates based on their parents frames

if(TRUE){
  framesBasedOnParents_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      framesBasedOnParents_G1[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k2,k],
                                                             vectors2rotate = frames_G1[,,k1,k])
    } 
  }
  print("Group 1 is done!")
  framesBasedOnParents_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      framesBasedOnParents_G2[,,k1,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k2,k],
                                                             vectors2rotate = frames_G2[,,k1,k])
    } 
  }
  print("Group 2 is done!")
}


# plot a distribution of an specific frame in local coordinate system based on the parent frame
frameNumber<-10 #choose a number in 1,...,71 as we have 71 frames
open3d()
spheres3d(c(0,0,0),col="white",alpha=0.2)
vectors3d(t(framesBasedOnParents_G1[1,,frameNumber,]),headlength = 0.1,radius = 1/10, col="red", lwd=1)
vectors3d(t(framesBasedOnParents_G1[2,,frameNumber,]),headlength = 0.1,radius = 1/10, col="blue", lwd=1)
vectors3d(t(framesBasedOnParents_G1[3,,frameNumber,]),headlength = 0.1,radius = 1/10, col="green", lwd=1)
vectors3d(diag(3),headlength = 0.1,radius = 1/10, col="black", lwd=1)
rglwidget()


#####################################################################################################
#####################################################################################################
# Calculate spokes directions based on their local frames

if(TRUE){
  spokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(nTotalRadii,3,nSamplesG1))
  for (k in 1:nSamplesG1) {
    for (i in 1:nTotalRadii) {
      spokeNo<-i # 1<=spokeNo<=nTotalRadii
      frameOfSpokeNo<-NA
      if(spokeNo<=upSpoeksNumber){
        frameOfSpokeNo<-spokeNo
      }else if(spokeNo<=2*upSpoeksNumber){
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }
      
      spokesDirectionsBasedOnFrames_G1[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,frameOfSpokeNo,k],
                                                                    vectors2rotate = spokeDirections_G1[i,,k])
      
    }
  }
  print("Group 1 is done!")
  spokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(nTotalRadii,3,nSamplesG2))
  for (k in 1:nSamplesG2) {
    for (i in 1:nTotalRadii) {
      spokeNo<-i # 1<=spokeNo<=nTotalRadii
      frameOfSpokeNo<-NA
      if(spokeNo<=upSpoeksNumber){
        frameOfSpokeNo<-spokeNo
      }else if(spokeNo<=2*upSpoeksNumber){
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }
      
      spokesDirectionsBasedOnFrames_G2[i,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,frameOfSpokeNo,k],
                                                                    vectors2rotate = spokeDirections_G2[i,,k])
      
      
    }
  }
  print("Group 2 is done!")
}

#####################################################################################################
#####################################################################################################
# Calculate connection lengths and directions based on their local frames

if(TRUE){
  
  connections_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  connectionsLengths_G1<-array(NA,dim=c(numberOfFrames,nSamplesG1))
  connectionsBasedOnParentFrames_G1<-array(NA,dim = dim(connections_G1))
  for (k in 1:nSamplesG1) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      tempVec<-skeletalSheet_G1[k1,,k]- skeletalSheet_G1[k2,,k]
      
      connectionsLengths_G1[k1,k]<-norm(tempVec,type = "2")
      
      if(norm(tempVec,type = "2")==0){
        connections_G1[k1,,k]<-c(0,0,0)
      }else{
        connections_G1[k1,,k]<-convertVec2unitVec(tempVec)
      }
      
      connectionsBasedOnParentFrames_G1[k1,,k]<-rotateFrameToMainAxes(myFrame = frames_G1[,,k2,k],
                                                                      vectors2rotate = connections_G1[k1,,k])
    }
  }
  print("Group 1 is done!")
  connections_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  connectionsLengths_G2<-array(NA,dim=c(numberOfFrames,nSamplesG2))
  connectionsBasedOnParentFrames_G2<-array(NA,dim = dim(connections_G2))
  for (k in 1:nSamplesG2) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      tempVec<-skeletalSheet_G2[k1,,k]- skeletalSheet_G2[k2,,k]
      
      connectionsLengths_G2[k1,k]<-norm(tempVec,type = "2")
      
      if(norm(tempVec,type = "2")==0){
        connections_G2[k1,,k]<-c(0,0,0)
      }else{
        connections_G2[k1,,k]<-convertVec2unitVec(tempVec)
      }
      
      connectionsBasedOnParentFrames_G2[k1,,k]<-rotateFrameToMainAxes(myFrame = frames_G2[,,k2,k],
                                                                      vectors2rotate = connections_G2[k1,,k])
    }
  }
  print("Group 2 is done!")
}

#####################################################################################################
#####################################################################################################
# Plot aligned skeletal sheets

framesGlobalCoordinate_G1<-array(NA,dim = dim(framesBasedOnParents_G1))
framesGlobalCoordinate_G1[,,16,]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
for (j in 1:nSamplesG1) {
  for (k in 2:numberOfFrames) {
    parent_Index<-framesParents[k]
    child_Index<-framesCenters[k]
    updatedParent<-framesGlobalCoordinate_G1[,,parent_Index,j]
    framesGlobalCoordinate_G1[,,child_Index,j]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                         vectorsInMainAxes = framesBasedOnParents_G1[,,child_Index,j])
  }
}
framesGlobalCoordinate_G2<-array(NA,dim = dim(framesBasedOnParents_G2))
framesGlobalCoordinate_G2[,,16,]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
for (j in 1:nSamplesG2) {
  for (k in 2:numberOfFrames) {
    parent_Index<-framesParents[k]
    child_Index<-framesCenters[k]
    updatedParent<-framesGlobalCoordinate_G2[,,parent_Index,j]
    framesGlobalCoordinate_G2[,,child_Index,j]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                         vectorsInMainAxes = framesBasedOnParents_G2[,,child_Index,j])
  }
}

connectionsGlobalCoordinate_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
for (j in 1:nSamplesG1) {
  for (i in 1:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    connectionsGlobalCoordinate_G1[k1,,j]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = framesGlobalCoordinate_G1[,,k2,j],
                                         vectorsInMainAxes = connectionsBasedOnParentFrames_G1[k1,,j])
  }  
}
connectionsGlobalCoordinate_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
for (j in 1:nSamplesG2) {
  for (i in 1:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    connectionsGlobalCoordinate_G2[k1,,j]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = framesGlobalCoordinate_G2[,,k2,j],
                                         vectorsInMainAxes = connectionsBasedOnParentFrames_G2[k1,,j])
  }  
}


positions_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
positions_G1[16,,]<-c(0,0,0)
for (j in 1:nSamplesG1) {
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    positions_G1[k1,,j]<-positions_G1[k2,,j]+
      connectionsLengths_G1[k1,j]*connectionsGlobalCoordinate_G1[k1,,j]
  }
}
positions_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
positions_G2[16,,]<-c(0,0,0)
for (j in 1:nSamplesG2) {
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    positions_G2[k1,,j]<-positions_G2[k2,,j]+
      connectionsLengths_G2[k1,j]*connectionsGlobalCoordinate_G2[k1,,j]
  }
}

spine_indices<-c(c(16, 13, 10, 7, 4, 1, 2, 3, 52, 19, 22, 25, 28, 31, 32, 33, 62))

open3d()
for (j in 1:nSamplesG1) {
  skelTemp<-positions_G1[,,j]
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    plot3d(rbind(skelTemp[k1,],skelTemp[k2,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
}
rglwidget()
for (j in 1:nSamplesG2) {
  skelTemp<-positions_G2[,,j]
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    plot3d(rbind(skelTemp[k1,],skelTemp[k2,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
  }
}
rglwidget()

# plot spines
open3d()
for (j in 1:nSamplesG1) {
  skelTemp<-positions_G1[,,j]
  for (k in spine_indices) {
    if(k==16){next}
    i<-which(framesCenters==k)
    k1<-framesCenters[i]
    k2<-framesParents[i]
    plot3d(rbind(skelTemp[k1,],skelTemp[k2,]),type="l",col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
}
rglwidget()
for (j in 1:nSamplesG2) {
  skelTemp<-positions_G2[,,j]
  for (k in spine_indices) {
    if(k==16){next}
    i<-which(framesCenters==k)
    k1<-framesCenters[i]
    k2<-framesParents[i]
    plot3d(rbind(skelTemp[k1,],skelTemp[k2,]),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
  }
}
rglwidget()

# calculate connections lengths of the spines
connectionsLengths_spine_G1<-array(NA,dim = c(nSamplesG1,17-1))
for (j in 1:nSamplesG1) {
  t<-1
  skelTemp<-positions_G1[,,j]
  for (k in spine_indices) {
    if(k==16){next}
    i<-which(framesCenters==k)
    k1<-framesCenters[i]
    k2<-framesParents[i]
    #cat(norm(skelTemp[k1,]-skelTemp[k2,],type = "2"),"\n")
    connectionsLengths_spine_G1[j,t]<-norm(c(skelTemp[k1,]-skelTemp[k2,]),type = "2")
    t<-t+1
  }
}
connectionsLengths_spine_G2<-array(NA,dim = c(nSamplesG2,17-1))
for (j in 1:nSamplesG2) {
  t<-1
  skelTemp<-positions_G2[,,j]
  for (k in spine_indices) {
    if(k==16){next}
    i<-which(framesCenters==k)
    k1<-framesCenters[i]
    k2<-framesParents[i]
    #cat(norm(skelTemp[k1,]-skelTemp[k2,],type = "2"),"\n")
    connectionsLengths_spine_G2[j,t]<-norm(c(skelTemp[k1,]-skelTemp[k2,]),type = "2")
    t<-t+1
  }
}

spineLengths_G1<-rowSums(connectionsLengths_spine_G1)
spineLengths_G2<-rowSums(connectionsLengths_spine_G2)

#test the spine lengths
boxplot(spineLengths_G1,spineLengths_G2)
t.test(spineLengths_G1,spineLengths_G2)

#####################################################################################################
# calculate LP sizes

#NB! length of the center frame is 0 and must be excluded
LP_sizes_G1<-rep(NA,nSamplesG1)
for (i in 1:nSamplesG1) {
  LP_sizes_G1[i]<-exp(mean(c(log(connectionsLengths_G1[-16,i]),log(radii_G1[,i]))))
}
LP_sizes_G2<-rep(NA,nSamplesG2)
for (i in 1:nSamplesG2) {
  LP_sizes_G2[i]<-exp(mean(c(log(connectionsLengths_G2[-16,i]),log(radii_G2[,i]))))
}

#####################################################################################################
#####################################################################################################
#Removing or preserving the scale by LP-size

if(typeOfStudy=="sizeAndShapeAnalysis"){
  
  # sizes_G1 and sizes_G2 are LP size but we use them here for plot and scaling
  
  sizes_G1<-rep(1,nSamplesG1)
  sizes_G2<-rep(1,nSamplesG2)
  
  radiiScaled_G1<-radii_G1 #we don't have scaling in size-and-shape analysis
  radiiScaled_G2<-radii_G2
  
  connectionsLengthsScaled_G1<-connectionsLengths_G1
  connectionsLengthsScaled_G2<-connectionsLengths_G2
  
}else if(typeOfStudy=="shapeAnalysis"){
  
  #NB! length of the center frame is 0 and must be excluded
  sizes_G1<-rep(NA,nSamplesG1)
  for (i in 1:nSamplesG1) {
    sizes_G1[i]<-exp(mean(c(log(connectionsLengths_G1[-16,i]),log(radii_G1[,i]))))
  }
  sizes_G2<-rep(NA,nSamplesG2)
  for (i in 1:nSamplesG2) {
    sizes_G2[i]<-exp(mean(c(log(connectionsLengths_G2[-16,i]),log(radii_G2[,i]))))
  }
  
  radiiScaled_G1<-sweep(radii_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  radiiScaled_G2<-sweep(radii_G2, 2, sizes_G2, "/") 
  
  connectionsLengthsScaled_G1<-sweep(connectionsLengths_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
  connectionsLengthsScaled_G2<-sweep(connectionsLengths_G2, 2, sizes_G2, "/") 
  
}

#####################################################################################################
#####################################################################################################
# Calculate mean LPDSRep

#calculate mean frames in local and global coordinate systems
if(TRUE){
  
  framesBasedOnParentsVectorized_G1<-array(NA,dim = c(numberOfFrames,9,nSamplesG1))
  for (i in 1:nSamplesG1) {
    for (k in 1:numberOfFrames) {
      framesBasedOnParentsVectorized_G1[k,,i]<-as.vector(t(framesBasedOnParents_G1[,,k,i]))
    }
  }
  framesBasedOnParentsVectorized_G2<-array(NA,dim = c(numberOfFrames,9,nSamplesG2))
  for (i in 1:nSamplesG2) {
    for (k in 1:numberOfFrames) {
      framesBasedOnParentsVectorized_G2[k,,i]<-as.vector(t(framesBasedOnParents_G2[,,k,i]))
    }
  }
  
  meanFramesBasedOnParents_G1<-array(NA, dim = c(3,3,numberOfFrames))
  for (k in framesCenters) {
    if(k==16){
      meanFramesBasedOnParents_G1[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
    }else{
      tempVec<-mean(as.SO3(t(framesBasedOnParentsVectorized_G1[k,,])),type = 'geometric')
      meanFramesBasedOnParents_G1[,,k]<-matrix(tempVec,nrow = 3,byrow = TRUE)
    }
  }
  meanFramesBasedOnParents_G2<-array(NA, dim = c(3,3,numberOfFrames))
  for (k in framesCenters) {
    if(k==16){
      meanFramesBasedOnParents_G2[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
    }else{
      tempVec<-mean(as.SO3(t(framesBasedOnParentsVectorized_G2[k,,])),type = 'geometric')
      meanFramesBasedOnParents_G2[,,k]<-matrix(tempVec,nrow = 3,byrow = TRUE)
    }
  }
  
  meanFramesGlobalCoordinate_G1<-array(NA,dim = dim(meanFramesBasedOnParents_G1))
  meanFramesGlobalCoordinate_G1[,,framesCenters[1]]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  for (k in 2:numberOfFrames) {
    parent_Index<-framesParents[k]
    child_Index<-framesCenters[k]
    updatedParent<-meanFramesGlobalCoordinate_G1[,,parent_Index]
    meanFramesGlobalCoordinate_G1[,,child_Index]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                         vectorsInMainAxes = meanFramesBasedOnParents_G1[,,child_Index])
  }
  meanFramesGlobalCoordinate_G2<-array(NA,dim = dim(meanFramesBasedOnParents_G2))
  meanFramesGlobalCoordinate_G2[,,framesCenters[1]]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  for (k in 2:numberOfFrames) {
    parent_Index<-framesParents[k]
    child_Index<-framesCenters[k]
    updatedParent<-meanFramesGlobalCoordinate_G2[,,parent_Index]
    meanFramesGlobalCoordinate_G2[,,child_Index]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                         vectorsInMainAxes = meanFramesBasedOnParents_G2[,,child_Index])
  }
  print("Done!")
}

# Calculate mean spokes' directions based on frames
if(TRUE){
  meanSpokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    # For extremely concentrated data we use Mardia mean direction 
    pcaTemp<-prcomp(t(spokesDirectionsBasedOnFrames_G1[i,,]))
    if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-convertVec2unitVec(colMeans(t(spokesDirectionsBasedOnFrames_G1[i,,])))
    }else if(typeOfMeanDirection=="Frechet"){
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-frechetMean(spokesDirectionsBasedOnFrames_G1[i,,]) 
    }else if(typeOfMeanDirection=="PNS"){
      sphereType<-kurtosisTestFunction(spokesDirectionsBasedOnFrames_G1[i,,])
      meanSpokesDirectionsBasedOnFrames_G1[i,]<-pns(spokesDirectionsBasedOnFrames_G1[i,,],sphere.type = sphereType)$PNS$mean 
    }else{
      stop("Please specify the typeOfMeanDirection by PNS or Frechet!")
    }
  }
  meanSpokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    # For extremely concentrated data we use Mardia mean direction 
    pcaTemp<-prcomp(t(spokesDirectionsBasedOnFrames_G2[i,,]))
    if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-convertVec2unitVec(colMeans(t(spokesDirectionsBasedOnFrames_G2[i,,])))
    }else if(typeOfMeanDirection=="Frechet"){
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-frechetMean(spokesDirectionsBasedOnFrames_G2[i,,])
    }else if(typeOfMeanDirection=="PNS"){
      sphereType<-kurtosisTestFunction(spokesDirectionsBasedOnFrames_G2[i,,])
      meanSpokesDirectionsBasedOnFrames_G2[i,]<-pns(spokesDirectionsBasedOnFrames_G2[i,,],sphere.type = sphereType)$PNS$mean 
    }else{
      stop("Please specify the typeOfMeanDirection by PNS or Frechet!")
    }
  }
  print("Done!")
}

# Calculate mean spokes' directions based on global coordinate (using mean frames in global coordinate)
if(TRUE){
  meanSpokesDirectionsGlobalCoordinate_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    spokeNo<-i # 1<=spokeNo<=nTotalRadii
    frameOfSpokeNo<-NA
    if(spokeNo<=upSpoeksNumber){
      frameOfSpokeNo<-spokeNo
    }else if(spokeNo<=2*upSpoeksNumber){
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }
    meanSpokesDirectionsGlobalCoordinate_G1[i,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,frameOfSpokeNo],
                                         vectorsInMainAxes = meanSpokesDirectionsBasedOnFrames_G1[i,])
    
  }
  meanSpokesDirectionsGlobalCoordinate_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    spokeNo<-i # 1<=spokeNo<=nTotalRadii
    frameOfSpokeNo<-NA
    if(spokeNo<=upSpoeksNumber){
      frameOfSpokeNo<-spokeNo
    }else if(spokeNo<=2*upSpoeksNumber){
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
    }
    meanSpokesDirectionsGlobalCoordinate_G2[i,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,frameOfSpokeNo],
                                         vectorsInMainAxes = meanSpokesDirectionsBasedOnFrames_G2[i,])
    
  }
}

# Calculate geometric mean of spokes' lengths
if(TRUE){
  radiiMean_G1<-exp(rowMeans(log(radiiScaled_G1)))
  radiiMean_G2<-exp(rowMeans(log(radiiScaled_G2)))
}

# Calculate geometric mean of connections' lengths
if(TRUE){
  meanConnectionsLengths_G1<-exp(rowMeans(log(connectionsLengthsScaled_G1)))
  meanConnectionsLengths_G2<-exp(rowMeans(log(connectionsLengthsScaled_G2)))
}

# Calculate mean connection directions based on frames
if(TRUE){
  meanConnectionsBasedOnParentFrames_G1<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    if(i==framesCenters[1]){
      meanConnectionsBasedOnParentFrames_G1[i,]<-c(0,0,0)
    }else{
      # For extremely concentrated data we use Mardia mean direction 
      pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G1[i,,]))
      if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
        meanConnectionsBasedOnParentFrames_G1[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G1[i,,])))
      }else if(typeOfMeanDirection=="Frechet"){
        meanConnectionsBasedOnParentFrames_G1[i,]<-frechetMean(connectionsBasedOnParentFrames_G1[i,,]) 
      }else if(typeOfMeanDirection=="PNS"){
        sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G1[i,,])
        meanConnectionsBasedOnParentFrames_G1[i,]<-pns(connectionsBasedOnParentFrames_G1[i,,],sphere.type = sphereType)$PNS$mean  
      }else{
        stop("Please specify the typeOfMeanDirection by PNS or Frechet")
      }
    }
  }
  meanConnectionsBasedOnParentFrames_G2<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    if(i==framesCenters[1]){
      meanConnectionsBasedOnParentFrames_G2[i,]<-c(0,0,0)
    }else{
      # For extremely concentrated data we use Mardia mean direction 
      pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G2[i,,]))
      if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
        meanConnectionsBasedOnParentFrames_G2[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G2[i,,])))
      }else if(typeOfMeanDirection=="Frechet"){
        meanConnectionsBasedOnParentFrames_G2[i,]<-frechetMean(connectionsBasedOnParentFrames_G2[i,,]) 
      }else if(typeOfMeanDirection=="PNS"){
        sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G2[i,,])
        meanConnectionsBasedOnParentFrames_G2[i,]<-pns(connectionsBasedOnParentFrames_G2[i,,],sphere.type = sphereType)$PNS$mean  
      }else{
        stop("Please specify the typeOfMeanDirection by PNS or Frechet")
      }
    }
  }
  print("Done!")
}

# Calculate mean connection based on global coordinate (using mean frames in global coordinate)
if(TRUE){
  meanConnectionsGlobalCoordinate_G1<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    meanConnectionsGlobalCoordinate_G1[k1,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,k2],
                                         vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G1[k1,])
  }
  meanConnectionsGlobalCoordinate_G2<-array(NA,dim = c(numberOfFrames,3))
  for (i in 1:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    meanConnectionsGlobalCoordinate_G2[k1,]<-
      rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,k2],
                                         vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G2[k1,])
  }
}

#####################################################################################################
#####################################################################################################
# Convert mean LPDSRep to a GP-ds-rep

if(TRUE){
  meanPositions_G1<-array(NA,dim = c(numberOfFrames,3))
  meanPositions_G1[16,]<-c(0,0,0)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    meanPositions_G1[k1,]<-meanPositions_G1[k2,]+
      meanConnectionsLengths_G1[k1]*meanConnectionsGlobalCoordinate_G1[k1,]
    
  }
  meanPositions_G2<-array(NA,dim = c(numberOfFrames,3))
  meanPositions_G2[16,]<-c(0,0,0)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    meanPositions_G2[k1,]<-meanPositions_G2[k2,]+
      meanConnectionsLengths_G2[k1]*meanConnectionsGlobalCoordinate_G2[k1,]
    
  }
  meanSpokesTails_G1<-array(NA,dim = c(nTotalRadii,3))
  meanSpokesTips_G1<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    frameOfSpokeNo<-NA
    if(i<=upSpoeksNumber){
      frameOfSpokeNo<-i
    }else if(i<=2*upSpoeksNumber){
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }
    
    meanSpokesTails_G1[i,]<-meanPositions_G1[frameOfSpokeNo,]
    meanSpokesTips_G1[i,]<-meanPositions_G1[frameOfSpokeNo,]+meanSpokesDirectionsGlobalCoordinate_G1[i,]*radiiMean_G1[i]
  }
  meanSpokesTails_G2<-array(NA,dim = c(nTotalRadii,3))
  meanSpokesTips_G2<-array(NA,dim = c(nTotalRadii,3))
  for (i in 1:nTotalRadii) {
    frameOfSpokeNo<-NA
    if(i<=upSpoeksNumber){
      frameOfSpokeNo<-i
    }else if(i<=2*upSpoeksNumber){
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }else{ #crest
      frameOfSpokeNo<-(i-upSpoeksNumber)
    }
    
    meanSpokesTails_G2[i,]<-meanPositions_G2[frameOfSpokeNo,]
    meanSpokesTips_G2[i,]<-meanPositions_G2[frameOfSpokeNo,]+meanSpokesDirectionsGlobalCoordinate_G2[i,]*radiiMean_G2[i]
  }
  print("Done!")
}

# plot mean spines
open3d()
skelTemp<-meanPositions_G1
for (k in spine_indices) {
  if(k==16){next}
  i<-which(framesCenters==k)
  k1<-framesCenters[i]
  k2<-framesParents[i]
  vectors3d(skelTemp[k1,],origin = skelTemp[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
}
skelTemp<-meanPositions_G2
for (k in spine_indices) {
  if(k==16){next}
  i<-which(framesCenters==k)
  k1<-framesCenters[i]
  k2<-framesParents[i]
  vectors3d(skelTemp[k1,],origin = skelTemp[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=1)
}
rglwidget()


# plot mean spines
skelTemp<-meanPositions_G1
for (k in framesCenters) {
  if(k==16){next}
  i<-which(framesCenters==k)
  k1<-framesCenters[i]
  k2<-framesParents[i]
  #vectors3d(skelTemp[k1,],origin = skelTemp[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  plot3d(rbind(skelTemp[k1,],skelTemp[k2,]),type="l",lwd=1,col = "blue",expand = 10,box=FALSE,add = TRUE)
}
skelTemp<-meanPositions_G2
for (k in framesCenters) {
  if(k==16){next}
  i<-which(framesCenters==k)
  k1<-framesCenters[i]
  k2<-framesParents[i]
  plot3d(rbind(skelTemp[k1,],skelTemp[k2,]),type="l",lwd=1,col = "red",expand = 10,box=FALSE,add = TRUE)
}
decorate3d()
rglwidget()

#####################################################################################################
#####################################################################################################
# Plot overlaid LPDSRep means of PD and CG

if(TRUE){
  open3d()
  srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1) #we scale back to the original size by *mean(sizes_G1)
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  srep2<-rbind(meanSpokesTails_G2,meanSpokesTips_G2)*mean(sizes_G2)
  for (i in 1:nTotalRadii) {
    plot3d(srep2[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  # legend3d("topright", legend = paste(c('Mean G1', 'Mean G2')), pch = 16, col = c("blue","red"), cex=1, inset=c(0.02))
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             # xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
             box = TRUE, axes = TRUE, main = "Overlaid mean shapes", sub = NULL,
             top = TRUE, aspect = FALSE, expand = 1.1)
  rglwidget()
}

#####################################################################################################
#####################################################################################################
# Hypothesis testing

# hypothesis test on LP size
pValues_LP_sizes<-meanDifferenceTest1D(log(LP_sizes_G1),log(LP_sizes_G2),type = typeOfTest) 
cat("pValue of LP sizes is:",pValues_LP_sizes,"\n")
boxplot(LP_sizes_G1, LP_sizes_G2, names = c("PG","CG"),main="LP-size")
cat("sd LP size G1:",sd(LP_sizes_G1),"sd LP size G2:",sd(LP_sizes_G2),"\n")
cat("Mean LP size G1:",mean(LP_sizes_G1),"mean LP size G2:",mean(LP_sizes_G2),"\n")

volume_G1<-rep(0,nSamplesG1)
for (i in 1:nSamplesG1) {
  tempMesh<-vcgImport(paste0(G1_path,"/",G1_list[i]))
  volume_G1[i]<-Rvcg::vcgVolume(tempMesh)
}
volume_G2<-rep(0,nSamplesG2)
for (i in 1:nSamplesG2) {
  tempMesh<-vcgImport(paste0(G2_path,"/",G2_list[i]))
  volume_G2[i]<-Rvcg::vcgVolume(tempMesh)
}
boxplot(volume_G1, volume_G2, names = c("PG","CG"),main="SPHARM-PDM Volume")

# hypothesis test on spokes' lengths
pValues_TtestRadii<-rep(NA,nTotalRadii)
pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
for (i in 1:nTotalRadii) {
  setTxtProgressBar(pb, i) #create progress bar
  
  pValues_TtestRadii[i]<-meanDifferenceTest1D(log(radiiScaled_G1[i,]),
                                              log(radiiScaled_G2[i,]),
                                              type = typeOfTest) 
  
}
# which(pValues_TtestRadii<=0.05)


# hypothesis test on connections' length
pValues_TtestConnectionsLengths<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for (i in 1:numberOfFrames) {
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValues_TtestConnectionsLengths[i]<-1
  }else{
    pValues_TtestConnectionsLengths[i]<-meanDifferenceTest1D(log(connectionsLengthsScaled_G1[i,]),
                                                             log(connectionsLengthsScaled_G2[i,]),
                                                             type = typeOfTest)
  }
}
# which(pValues_TtestConnectionsLengths<=0.05)


# hypothesis test on spokes' directions based on local frames
euclideanizedSpokesDirBasedOnFramesG1<-array(NA,dim = c(nSamplesG1,2,nTotalRadii))
euclideanizedSpokesDirBasedOnFramesG2<-array(NA,dim = c(nSamplesG2,2,nTotalRadii))
pValspokesDirectionsBasedOnFrames<-rep(NA,nTotalRadii)
pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
for(i in 1:nTotalRadii){
  setTxtProgressBar(pb, i) #create progress bar
  #NB! euclideanization must contain two groups because it uses the pooled mean 
  euclideanizedTemp<-euclideanization(spokesDirectionsBasedOnFrames_G1[i,,],
                                      spokesDirectionsBasedOnFrames_G2[i,,],
                                      type = typeOfStudy4directions)
  
  pValspokesDirectionsBasedOnFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                       euclideanizedTemp$euclideanG2,
                                                                       type=typeOfTest)

  euclideanizedSpokesDirBasedOnFramesG1[,,i]<-euclideanizedTemp$euclideanG1
  euclideanizedSpokesDirBasedOnFramesG2[,,i]<-euclideanizedTemp$euclideanG2
  
}
# which(pValspokesDirectionsBasedOnFrames<=0.05)


# hypothesis test on connections' directions based on local frames
euclideanizedConnectionsBasedOnParentFramesG1<-array(0,dim = c(nSamplesG1,2,numberOfFrames))
euclideanizedConnectionsBasedOnParentFramesG2<-array(0,dim = c(nSamplesG2,2,numberOfFrames))
pValConnectionsBasedOnParentFrames<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for(i in 1:numberOfFrames){
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValConnectionsBasedOnParentFrames[i]<-1
    next
  }
  
  euclideanizedTemp<-euclideanization(connectionsBasedOnParentFrames_G1[i,,],
                                      connectionsBasedOnParentFrames_G2[i,,],
                                      type = typeOfStudy4directions)
  
  pValConnectionsBasedOnParentFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                        euclideanizedTemp$euclideanG2,
                                                                        type=typeOfTest)
  
  
  euclideanizedConnectionsBasedOnParentFramesG1[,,i]<-euclideanizedTemp$euclideanG1
  euclideanizedConnectionsBasedOnParentFramesG2[,,i]<-euclideanizedTemp$euclideanG2
  
}
# which(pValConnectionsBasedOnParentFrames<=0.05)

framesBasedOnParentsVectorized_G1<-array(NA,dim = c(numberOfFrames,9,nSamplesG1))
for (i in 1:nSamplesG1) {
  for (k in 1:numberOfFrames) {
    framesBasedOnParentsVectorized_G1[k,,i]<-as.vector(t(framesBasedOnParents_G1[,,k,i]))
  }
}
framesBasedOnParentsVectorized_G2<-array(NA,dim = c(numberOfFrames,9,nSamplesG2))
for (i in 1:nSamplesG2) {
  for (k in 1:numberOfFrames) {
    framesBasedOnParentsVectorized_G2[k,,i]<-as.vector(t(framesBasedOnParents_G2[,,k,i]))
  }
}

# hypothesis test on frames' normal directions based on parent frames
euclideanizedFrameBasedOnParentG1<-array(0,dim = c(nSamplesG1,3,numberOfFrames))
euclideanizedFrameBasedOnParentG2<-array(0,dim = c(nSamplesG2,3,numberOfFrames))
pValFramesBasedOnParent<-rep(NA,numberOfFrames)
pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
for(i in 1:numberOfFrames){
  setTxtProgressBar(pb, i) #create progress bar
  if(i==16){
    pValFramesBasedOnParent[i]<-1
    next
  }
  
  Q4Temp<-as.Q4(as.SO3(t(framesBasedOnParentsVectorized_G1[i,,])))
  Q4Temp2<-matrix(as.numeric(t(Q4Temp)),ncol = 4,byrow = TRUE)
  for (j in 1:dim(Q4Temp2)[1]) {
    Q4Temp2[j,]<-Q4Temp2[j,]/norm(Q4Temp2[j,],type = '2')
  }
  Q4Temp3<-as.Q4(as.SO3(t(framesBasedOnParentsVectorized_G2[i,,])))
  Q4Temp4<-matrix(as.numeric(t(Q4Temp3)),ncol = 4,byrow = TRUE)
  for (j in 1:dim(Q4Temp4)[1]) {
    Q4Temp4[j,]<-Q4Temp4[j,]/norm(Q4Temp4[j,],type = '2')
  }
  
  euclideanizedTemp<-euclideanization(t(Q4Temp2),t(Q4Temp4),type = typeOfStudy4directions)
  
  
  pValFramesBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                             euclideanizedTemp$euclideanG2,
                                                             type=typeOfTest) 
  
  euclideanizedFrameBasedOnParentG1[,,i]<-euclideanizedTemp$euclideanG1
  euclideanizedFrameBasedOnParentG2[,,i]<-euclideanizedTemp$euclideanG2
  
}
# which(pValFramesBasedOnParent<=0.05)


#####################################################################################################
#####################################################################################################
# plot significant GOPs

pvalues_LP_ds_rep <- c(pValues_TtestRadii,                 
                       pValues_TtestConnectionsLengths,    
                       pValspokesDirectionsBasedOnFrames,  
                       pValConnectionsBasedOnParentFrames, 
                       pValFramesBasedOnParent,            
                       pValues_LP_sizes)

n_s<-nTotalRadii
n_f<-numberOfFrames

length(pvalues_LP_ds_rep)
alpha<-0.05
significantPvalues<-which(pvalues_LP_ds_rep<=alpha)
significantPvalues

#adjust p-values by Benjamini-Hochberg
FDR<-0.15
pvalues_LP_ds_rep_BH<-p.adjust(pvalues_LP_ds_rep,method = "BH")
pvalues_LP_ds_rep_Bonferroni<-p.adjust(pvalues_LP_ds_rep,method = "bonferroni")
significantPvalues_BH<-which(pvalues_LP_ds_rep_BH<=FDR)
significantPvalues_BH
significantPvalues_Bonferroni<-which(pvalues_LP_ds_rep_Bonferroni<=FDR)
significantPvalues_Bonferroni

cat("\n","Percentage of sig raw p-value is:",length(significantPvalues)/length(pvalues_LP_ds_rep),"\n")
cat("\n","Percentage of BH adjusted p-value is:",length(significantPvalues_BH)/length(pvalues_LP_ds_rep),"\n")



# plot by ggplot
df_LP <- data.frame(Type=c(rep("Raw p-value",length(pvalues_LP_ds_rep)),
                           rep("Bonferroni",length(pvalues_LP_ds_rep)),
                           rep("BH",length(pvalues_LP_ds_rep))),
                    ordereOfPvalues=1:length(pvalues_LP_ds_rep),
                    Values=c(sort(pvalues_LP_ds_rep),sort(pvalues_LP_ds_rep_Bonferroni),sort(pvalues_LP_ds_rep_BH)))
p<-ggplot(df_LP, aes(x=ordereOfPvalues, y=Values, group=Type))
p + geom_line(aes(linetype=Type),size=1)+
  # geom_line(aes(linetype=Type, col=Type),size=1)+
  # geom_point(aes(shape=Type),alpha=0.7)+
  geom_hline(yintercept=alpha,linetype="solid", color = "red")+
  geom_hline(yintercept=FDR,linetype="solid", color = "blue")+
  scale_linetype_manual(values=c("solid","dotdash", "dotted")) +
  theme_bw()+
  theme(plot.title = element_text(size = 17, hjust = 0.5),
        legend.text=element_text(size=17),
        legend.title=element_blank(),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # legend.title=element_text(size=17),
        legend.position="bottom",
        axis.title=element_text(size=17),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  guides(colour = guide_legend(title.hjust = 0.5))+
  xlab("Ranking of p-values") + ylab("p-values")


#1
significantRadii<-significantPvalues[which(significantPvalues<=n_s)]
significantRadii
significantRadii_BH<-significantPvalues_BH[which(significantPvalues_BH<=n_s)]
significantRadii_BH
#2
significantConnectionsLengths<-significantPvalues[which(n_s+1<=significantPvalues
                                                        & significantPvalues<=(n_s+n_f))]-n_s
significantConnectionsLengths
significantConnectionsLengths_BH<-significantPvalues_BH[which((n_s+1)<=significantPvalues_BH
                                                              & significantPvalues_BH<=(n_s+n_f))]-n_s
significantConnectionsLengths_BH
#3
significantspokesDirections<-significantPvalues[which((n_s+n_f+1)<=significantPvalues
                                                      & significantPvalues<=(2*n_s+n_f))]-(n_s+n_f)
significantspokesDirections
significantspokesDirections_BH<-significantPvalues_BH[which((n_s+n_f+1)<=significantPvalues_BH
                                                            & significantPvalues_BH<=(2*n_s+n_f))]-(n_s+n_f)
significantspokesDirections_BH
#4
significantConnectionsDirections<-significantPvalues[which((2*n_s+n_f+1)<=significantPvalues
                                                           & significantPvalues<=(2*n_s+2*n_f))]-(2*n_s+n_f)
significantConnectionsDirections
significantConnectionsDirections_BH<-significantPvalues_BH[which((2*n_s+n_f+1)<=significantPvalues_BH
                                                                 & significantPvalues_BH<=(2*n_s+2*n_f))]-(2*n_s+n_f)
significantConnectionsDirections_BH
#5
significantFrame<-significantPvalues[which((2*n_s+2*n_f+1)<=significantPvalues
                                           & significantPvalues<=(2*n_s+3*n_f))]-(2*n_s+2*n_f)
significantFrame
significantFrame_BH<-significantPvalues_BH[which((2*n_s+2*n_f+1)<=significantPvalues_BH
                                                 & significantPvalues_BH<=(2*n_s+3*n_f))]-(2*n_s+2*n_f)
significantFrame_BH


#plot significant GOPs before and after the BH adjustment
if(TRUE){
  
  #1 plot
  srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1)
  open3d()
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant spokes' lengths",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #with correction
  open3d()
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantRadii_BH) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant spokes' lengths after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #2 plot
  open3d()
  skelG1_1<-meanSpokesTails_G1[skelRange,]*mean(sizes_G1)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in significantConnectionsLengths) {
    k1<-i
    k2<-framesParents[which(framesCenters==i)]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant connections' lengths",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #with correction
  open3d()
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in significantConnectionsLengths_BH) {
    k1<-i
    k2<-framesParents[which(framesCenters==i)]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant connections' lengths after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #3 plot
  srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1)
  open3d()
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantspokesDirections) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant spokes' directions",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #with correction
  open3d()
  for (i in 1:nTotalRadii) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
  }
  for (i in significantspokesDirections_BH) {
    plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 7,col = "red",expand = 10,box=FALSE,add = TRUE)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant spokes' directions after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #4 plot
  open3d()
  skelG1_1<-meanSpokesTails_G1[skelRange,]*mean(sizes_G1)
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in significantConnectionsDirections) {
    k1<-i
    k2<-framesParents[which(framesCenters==i)]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant connections' directions",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #with correction
  open3d()
  for (i in 2:numberOfFrames) {
    k1<-framesCenters[i]
    k2<-framesParents[i]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
  }
  for (i in significantConnectionsDirections_BH) {
    k1<-i
    k2<-framesParents[which(framesCenters==i)]
    vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant connections' directions after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #5 plot
  open3d()
  skeletalSheet<-meanPositions_G1*mean(sizes_G1)
  for (i in 2:numberOfFrames) {
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="darkblue", lwd=2)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="blue", lwd=2)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="lightblue", lwd=2)
  }
  for (i in significantFrame) {
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant frames",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
  
  #with correction
  open3d()
  for (i in 2:numberOfFrames) {
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="darkblue", lwd=2)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="blue", lwd=2)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="lightblue", lwd=2)
  }
  for (i in significantFrame_BH) {
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
    vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=7)
  }
  decorate3d(xlab = "x", ylab = "y", zlab = "z",
             box = F, axes = TRUE,
             main = "Significant frames after BH adjustment",
             sub = NULL, top = T, aspect = FALSE, expand = 1.1)
  rglwidget()
}


####################################################################################################
####################################################################################################
#################################### Classification ################################################
####################################################################################################
####################################################################################################
# classification and cross-validation

library(caret)
library(e1071)
library("kdensity")

significantRadii
significantConnectionsLengths
significantspokesDirections
significantConnectionsDirections
significantFrame


# basedOnSignificantGOPs<-TRUE
basedOnSignificantGOPs<-FALSE
if(basedOnSignificantGOPs==TRUE){
  classG1<-cbind(t(radiiScaled_G1[significantRadii,]),
                 t(connectionsLengthsScaled_G1[significantConnectionsLengths,]),
                 euclideanizedSpokesDirBasedOnFramesG1[,1,significantspokesDirections],
                 euclideanizedSpokesDirBasedOnFramesG1[,2,significantspokesDirections],
                 euclideanizedConnectionsBasedOnParentFramesG1[,1,significantConnectionsDirections],
                 euclideanizedConnectionsBasedOnParentFramesG1[,2,significantConnectionsDirections],
                 euclideanizedFrameBasedOnParentG1[,1,significantFrame],
                 euclideanizedFrameBasedOnParentG1[,2,significantFrame],
                 euclideanizedFrameBasedOnParentG1[,3,significantFrame],
                 LP_sizes_G1)
  classG2<-cbind(t(radiiScaled_G2[significantRadii,]),
                 t(connectionsLengthsScaled_G2[significantConnectionsLengths,]),
                 euclideanizedSpokesDirBasedOnFramesG2[,1,significantspokesDirections],
                 euclideanizedSpokesDirBasedOnFramesG2[,2,significantspokesDirections],
                 euclideanizedConnectionsBasedOnParentFramesG2[,1,significantConnectionsDirections],
                 euclideanizedConnectionsBasedOnParentFramesG2[,2,significantConnectionsDirections],
                 euclideanizedFrameBasedOnParentG2[,1,significantFrame],
                 euclideanizedFrameBasedOnParentG2[,2,significantFrame],
                 euclideanizedFrameBasedOnParentG2[,3,significantFrame],
                 LP_sizes_G2)
  
}else{
  classG1<-cbind(t(radiiScaled_G1),
                 t(connectionsLengthsScaled_G1[-16,]),
                 euclideanizedSpokesDirBasedOnFramesG1[,1,],
                 euclideanizedSpokesDirBasedOnFramesG1[,2,],
                 euclideanizedConnectionsBasedOnParentFramesG1[,1,-16],
                 euclideanizedConnectionsBasedOnParentFramesG1[,2,-16],
                 euclideanizedFrameBasedOnParentG1[,1,-16],
                 euclideanizedFrameBasedOnParentG1[,2,-16],
                 euclideanizedFrameBasedOnParentG1[,3,-16],
                 LP_sizes_G1)
  classG2<-cbind(t(radiiScaled_G2),
                 t(connectionsLengthsScaled_G2[-16,]),
                 euclideanizedSpokesDirBasedOnFramesG2[,1,],
                 euclideanizedSpokesDirBasedOnFramesG2[,2,],
                 euclideanizedConnectionsBasedOnParentFramesG2[,1,-16],
                 euclideanizedConnectionsBasedOnParentFramesG2[,2,-16],
                 euclideanizedFrameBasedOnParentG2[,1,-16],
                 euclideanizedFrameBasedOnParentG2[,2,-16],
                 euclideanizedFrameBasedOnParentG2[,3,-16],
                 LP_sizes_G2)
}


classes<-as.factor(c(rep('CG',nSamplesG1),rep('PD',nSamplesG2)))

LP_srep_Data<-data.frame(rbind(classG1,classG2),classes)

crossValidation <- train(classes~., data=LP_srep_Data,
                         "svmLinear",
                         # "svmRadial",
                         # "svmPoly",
                         # "lda",
                         # "rpart",
                         # "nnet",
                         tuneLength = 10,
                         trControl = trainControl(method = "cv"))
crossValidation

mean(crossValidation$results$Accuracy)
mean(crossValidation$results$Kappa)

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

crossValidation <- train(classes~., data=LP_srep_Data,
                         "svmLinear",
                         # "knn",
                         tuneLength = 10,
                         trControl = fitControl,
                         metric = "ROC")

crossValidation