#include "geoEff.h"

#include <iostream>
#include <vector>
#include <random>
#include <math.h>

#include <Eigen/Dense>

geoEff::geoEff(int seed, bool verbose){

  verbosity = verbose;

  if (verbosity) {
    std::cout << "geoEff constructor called" << std::endl;
  }

  if (verbosity) {
    std::cout << "geoEff set seed to " << seed << std::endl;
  }
  
  prnGenerator = std::mt19937_64(seed);
  if (verbosity) {
    std::cout << "geoEff set random number generator to mt19937_64" << std::endl;
  }

  uniform = std::uniform_real_distribution<>(0., 1.);
  if (verbosity){
    std::cout << "geoEff set uniform distribution in [0, 1]" << std::endl;
  }

  translations[0].reserve(N_THROWS);
  translations[1].reserve(N_THROWS);
  translations[2].reserve(N_THROWS);
  rotations.reserve(N_THROWS);
  transforms.reserve(N_THROWS);
  if (verbosity){
    std::cout << "geoEff allocated memory for transformation vectors" << std::endl;
  }
  
  vetoSize = std::vector<float>(1,30.);
  vetoEnergy = std::vector<float>(1,30.);
  if (verbosity){
    std::cout << "geoEff set veto size to 30 cm and energy threshold to 30 MeV" << std::endl;
  }

}

void geoEff::setVertex(float x, float y, float z){
  vertex[0] = x;
  vertex[1] = y;
  vertex[2] = z;
}

void geoEff::setHitSegEdeps(std::vector<float> thishitSegEdeps){
  hitSegEdeps = thishitSegEdeps;
}

void geoEff::setHitSegPoss(std::vector<float> thishitSegPoss){

  // Set the vector
  hitSegPoss = thishitSegPoss;
}

void geoEff::setRangeX(float xmin, float xmax){
  range[0][0] = xmin;
  range[0][1] = xmax;
}
void geoEff::setRangeY(float ymin, float ymax){
  range[1][0] = ymin;
  range[1][1] = ymax;
}
void geoEff::setRangeZ(float zmin, float zmax){
  range[2][0] = zmin;
  range[2][1] = zmax;
}

void geoEff::setActiveX(float xmin, float xmax){
  active[0][0] = xmin;
  active[0][1] = xmax;
}
void geoEff::setActiveY(float ymin, float ymax){
  active[1][0] = ymin;
  active[1][1] = ymax;
}
void geoEff::setActiveZ(float zmin, float zmax){
  active[2][0] = zmin;
  active[2][1] = zmax;
}

void geoEff::setOffsetX(float x){
  offset[0] = x;
}

void geoEff::setOffsetY(float y){
  offset[1] = y;
}

void geoEff::setOffsetZ(float z){
  offset[2] = z;
}
  
void geoEff::setBeamDir(float xdir, float ydir, float zdir){
  beamdir[0] = xdir;
  beamdir[1] = ydir;
  beamdir[2] = zdir;
}

void geoEff::setVetoSizes(std::vector< float > vSizes){
  vetoSize = vSizes;
}
void geoEff::setVetoEnergyThresholds(std::vector< float > vThresholds){
  vetoEnergy = vThresholds;
}

void geoEff::throwTransforms(){

  // Clear vectors
  translations[0].clear();
  translations[1].clear();
  translations[2].clear();
  rotations.clear();
  transforms.clear();

  for (int dim = 0; dim < 3; dim++){
    if ((range[dim][0] < 0) or (range[dim][1] < 0)){
      translations[dim].resize(N_THROWS, vertex[dim]);
    } else {
      translations[dim].clear();
      for (int i = 0; i < N_THROWS; i++){
        translations[dim].emplace_back(uniform(prnGenerator)*(range[dim][1]-range[dim][0])+range[dim][0]+ offset[dim]);
      }
    }
  }

  rotations.clear();
  for (int i = 0; i < N_THROWS; i++){
    rotations.emplace_back((uniform(prnGenerator)-0.5)*2*M_PI);
  }

  // Tranformations that do not depend on the throws:
  // Move vertex to coordinate system origin to apply rotation
  Eigen::Affine3f tThere(Eigen::Translation3f(Eigen::Vector3f(-vertex[0], -vertex[1], -vertex[2])));
  // Move vertex back
  Eigen::Affine3f tBack(Eigen::Translation3f(Eigen::Vector3f(vertex[0], vertex[1], vertex[2])));

  // Store Eigen transforms
  transforms.clear();
  for (int i = 0; i < N_THROWS; i++){
    // Vertex displacement:
    Eigen::Affine3f tThrow(Eigen::Translation3f(Eigen::Vector3f(translations[0][i]-vertex[0],
                                                                translations[1][i]-vertex[1],
                                                                translations[2][i]-vertex[2])));

    // Rotation
    Eigen::Affine3f rThrow(Eigen::Affine3f(Eigen::AngleAxisf(rotations[i], Eigen::Vector3f(beamdir[0], beamdir[1], beamdir[2]))));

    // Put everything together in single transform and store.
    transforms.emplace_back(tThrow * tBack * rThrow * tThere);
  }
  
}

std::vector<float> geoEff::getCurrentThrowTranslationsX(){
  return translations[0];
}
std::vector<float> geoEff::getCurrentThrowTranslationsY(){
  return translations[1];
}
std::vector<float> geoEff::getCurrentThrowTranslationsZ(){
  return translations[2];
}
std::vector<float> geoEff::getCurrentThrowRotations(){
  return rotations;
}

std::vector< std::vector< uint32_t > > geoEff::getHadronContainment(){
  // Pass/fail for each set of vetoSize and vetoEnergy
  std::vector< std::vector< uint32_t > > hadronContainment(vetoSize.size(), std::vector< uint32_t >(vetoEnergy.size(), 0));

  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,hitSegPoss.size()/3,Eigen::OuterStride<>(3));
  
  // Check if event is contained by any of the existing conditions
  int origContained = 0;
  for (int i = 0; i < vetoSize.size(); i++){
    for (int j = 0; j < vetoEnergy.size(); j++){
      if (isContained(hitSegPosOrig, hitSegEdeps, vetoSize[i], vetoEnergy[i])) origContained++;
    }
  }
  
  // If not, then return
  if (origContained == 0) return hadronContainment;

  // Else, loop through set of rotation translations
  for (int t = 0; t < N_THROWS; t++){
    // Apply transformation to energy deposit positions
    Eigen::Matrix3Xf transformedEdeps = transforms[t] * hitSegPosOrig;
    // Loop through conditions
    for (int i = 0; i < vetoSize.size(); i++){
      for (int j = 0; j < vetoEnergy.size(); j++){
        // Check containment and set bit
        if (isContained(hitSegPosOrig, hitSegEdeps, vetoSize[i], vetoEnergy[i])) hadronContainment[i][j] += 1<<t;
      }
    }
  }
  
  return hadronContainment;
}

void geoEff::setSeed(int seed){
  prnGenerator = std::mt19937_64(seed);
}

bool geoEff::isContained( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits, float vSize, float vetoEnergyThreshold ){

  float vetoEnergy = 0.;

  for (int i = 0; i < energyDeposits.size(); i++){
    for (int dim = 0; dim < 3; dim++){
      // low
      if ( (hitSegments(dim, i)-offset[dim] < active[dim][0]+vSize) and
           (hitSegments(dim, i)-offset[dim] > active[dim][0]) ) {
        vetoEnergy += energyDeposits[i];
        break; // Only count each energy deposit once
      }
      // high
      if ( (hitSegments(dim, i)-offset[dim] > active[dim][1]-vSize) and
           (hitSegments(dim, i)-offset[dim] < active[dim][1]) ) {
        vetoEnergy += energyDeposits[i];
        break; // Only count each energy deposit once
      }
    }
  }

  return vetoEnergy < vetoEnergyThreshold;
}
