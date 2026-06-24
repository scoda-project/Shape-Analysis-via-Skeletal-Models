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
# Hotelling T2 Test for SPHARM-PDM shape analysis

numberOfPoints<-1002

# Hotelling test
pValues_Hotelling_SPHARM_PDM<-c()
for (i in 1:numberOfPoints) {
  pValues_Hotelling_SPHARM_PDM<-c(pValues_Hotelling_SPHARM_PDM,
                                  HotellingT2(t(aligned_SpharmPDM_G1[i,,]),t(aligned_SpharmPDM_G2[i,,])))
}

alpha<-0.05
significantpValues_HotellingSpharmPDM<-which(pValues_Hotelling_SPHARM_PDM<=alpha)
significantpValues_HotellingSpharmPDM

FDR<-0.05
pValues_Hotelling_SPHARM_PDM_BH<-p.adjust(pValues_Hotelling_SPHARM_PDM,method = "BH")
significantpValues_HotellingSpharmPDM_BH<-which(pValues_Hotelling_SPHARM_PDM_BH<=FDR)
significantpValues_HotellingSpharmPDM_BH 

#plot significant without correction
open3d()
plot3d(proc$mshape,type="s", radius = 0.1,col = "blue",expand = 10, box=FALSE,add = TRUE)
plot3d(proc$mshape[significantpValues_HotellingSpharmPDM,],type="s",
       radius = 0.3 ,col = "red",expand = 10,box=FALSE,add = TRUE)
verts <- rbind(t(as.matrix(proc$mshape)),1)
trgls <- as.matrix(data_obj1$it)
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="white",alpha=0.2)  #surface mesh
wire3d(tmesh, col="lightgrey")  #surface mesh
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)
rglwidget()

#plot significant by FDR correction
open3d()
plot3d(proc$mshape,type="s", radius = 0.1,col = "blue",expand = 10, box=FALSE,add = TRUE)
plot3d(proc$mshape[significantpValues_HotellingSpharmPDM_BH,],type="s",
       radius = 0.3 ,col = "red",expand = 10,box=FALSE,add = TRUE)
verts <- rbind(t(as.matrix(proc$mshape)),1)
trgls <- as.matrix(data_obj1$it)
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="white",alpha=0.2)  #surface mesh
wire3d(tmesh, col="lightgrey")  #surface mesh
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)
rglwidget()

