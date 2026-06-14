#####################################################################################################
#####################################################################################################
# libraries

library(reticulate)
library(ETRep)
library(rgl)
library(Rvcg)
library(Morpho)
library(ptinpoly)
library(keras)
library(tensorflow)
library(plotly)
library(keras3)

#####################################################################################################
#####################################################################################################

# Clear the environment
remove(list=ls())

# Set working directory to file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#####################################################################################################
#####################################################################################################
# Functions 

uniformScalingETube <- function(tube, scaleFactor = 1) {
  scaledTube<- ETRep::create_Elliptical_Tube(
    numberOfFrames=nrow(tube$spinalPoints3D),
    method = "basedOnMaterialFrames",
    materialFramesBasedOnParents = tube$materialFramesBasedOnParents,
    initialFrame = diag(3),
    initialPoint = c(0, 0, 0),
    ellipseRadii_a = tube$ellipseRadii_a*scaleFactor,
    ellipseRadii_b = tube$ellipseRadii_b*scaleFactor,
    connectionsLengths = tube$connectionsLengths*scaleFactor,
    plotting = FALSE)
  
  return(scaledTube)
}

#####################################################################################################
#####################################################################################################
# Simulation 

# Load artificial samples
load("../Files/hippoAsTube_1.RData")
load("../Files/hippoAsTube_2.RData")

# plot the samples
tubeMesh1<-ETRep::tube_Surface_Mesh(tube = hippoAsTube_1,
                                   meshType = "triangular", plotMesh = FALSE)
open3d()
shade3d(tubeMesh1, col="blue", alpha=0.5)
wire3d(tubeMesh1, col="black")
rglwidget()
tubeMesh2<-ETRep::tube_Surface_Mesh(tube = hippoAsTube_2,
                                   meshType = "triangular", plotMesh = FALSE)
open3d()
shade3d(tubeMesh2, col="red", alpha=0.5)
wire3d(tubeMesh2, col="black")
rglwidget()


# Using the samples to generate two groups of simulated objects as G1 and G2
if(TRUE){
  set.seed(4)
  dir.create("../Files/G1_Obj")
  dir.create("../Files/G2_Obj")
  nSim<-20
  # Generate random samples for G1
  for (i in 1:nSim) {
    tubeTemp <- ETRep::simulate_etube(
      referenceTube = hippoAsTube_1,
      sd_v = 0.1, sd_psi = 0.01, sd_x = 0.0001,
      sd_a = 0.000001, sd_b = 0.000001,
      rangeSdScale = c(1,1.1),
      numberOfSimulation = 1,
      plotting = FALSE
    )[[1]]
    
    scaledTube<- uniformScalingETube(tube = tubeTemp, 
                                     scaleFactor = 400)
    
    tubeMesh<-ETRep::tube_Surface_Mesh(tube = scaledTube,
                                       meshType = "triangular", plotMesh = FALSE)
    
    tubeMesh<-Rvcg::vcgUniformRemesh(tubeMesh)
    tubeMesh_reduced <- Rvcg::vcgQEdecim(tubeMesh, percent = 0.1)
    
    open3d()
    shade3d(tubeMesh_reduced, col="blue", alpha=0.5)
    wire3d(tubeMesh_reduced, col="black")
    rglwidget()
    
    writeOBJ(paste0("../Files/G1_Obj/sample_",i,".obj"), tubeMesh_reduced)
  }
  
  set.seed(4)
  # Generate random samples for G2
  for (i in 1:nSim) {
    tubeTemp <- ETRep::simulate_etube(
      referenceTube = hippoAsTube_2,
      sd_v = 0.1, sd_psi = 0.01, sd_x = 0.0001,
      sd_a = 0.000001, sd_b = 0.000001,
      rangeSdScale = c(1,1.1),
      numberOfSimulation = 1,
      plotting = FALSE
    )[[1]]
    
    scaledTube<- uniformScalingETube(tube = tubeTemp,
                                     scaleFactor = 400)
    
    tubeMesh<-ETRep::tube_Surface_Mesh(tube = scaledTube,
                                       meshType = "triangular", plotMesh = FALSE)
    
    tubeMesh<-Rvcg::vcgUniformRemesh(tubeMesh)
    tubeMesh_reduced <- Rvcg::vcgQEdecim(tubeMesh, percent = 0.1)
    
    open3d()
    shade3d(tubeMesh_reduced, col="blue", alpha=0.5)
    wire3d(tubeMesh_reduced, col="black")
    rglwidget()
    
    writeOBJ(paste0("../Files/G2_Obj/sample_",i,".obj"), tubeMesh_reduced)
  }

}