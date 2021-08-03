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

  N_THROWS = 64*64;
  if (verbosity){
    std::cout << "Number of throws set to " << N_THROWS << std::endl;
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

  if (verbosity){
    std::cout << "geoEff allocated memory for transformation vectors" << std::endl;
  }
  
  vetoSize = std::vector<float>(1,30.);
  vetoEnergy = std::vector<float>(1,30.);
  if (verbosity){
    std::cout << "geoEff set veto size to 30 cm and energy threshold to 30 MeV" << std::endl;
  }

  useFixedBeamDir = false;    

  // Initialize to all dimensions randomized
  for (int dim = 0; dim < 3; dim++) randomizeVertex[dim] = true;
}

void geoEff::setNthrows(unsigned long n){
  N_THROWS = n;
  if (verbosity){
    std::cout << "geoEff set number of throws to " << N_THROWS << std::endl;
  }

  if (N_THROWS%64) std::cout << "geoEff warning: number of throws should be multiple of 64 for optimal use of output format."  << std::endl;
}

void geoEff::setVertex(float x, float y, float z){
  vertex[0] = x;
  vertex[1] = y;
  vertex[2] = z;
  if(verbosity){
    std::cout << "geoEff set vertex to " << vertex[0] << " "<< vertex[1] << " "<< vertex[2] << std::endl;
  }
}

void geoEff::setHitSegEdeps(std::vector<float> thishitSegEdeps){
  hitSegEdeps = thishitSegEdeps;
  if (verbosity) {
    std::cout << "geoEff setting hit segment energy deposits to ";
    for (unsigned int i = 0; i < hitSegEdeps.size(); i++) std::cout << hitSegEdeps[i] << " ";
    std::cout << std::endl;
  }
}

void geoEff::setHitSegPoss(std::vector<float> thishitSegPoss){

  // Set the vector
  hitSegPoss = thishitSegPoss;
  if (verbosity) {
    std::cout << "geoEff setting hit segment positions to ";
    for (unsigned int i = 0; i < hitSegPoss.size(); i++) std::cout << hitSegPoss[i] << " ";
    std::cout << std::endl;
  }

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

void geoEff::setRandomizeX(bool r){
  randomizeVertex[0] = r;
}

void geoEff::setRandomizeY(bool r){
  randomizeVertex[1] = r;
}

void geoEff::setRandomizeZ(bool r){
  randomizeVertex[2] = r;
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

void geoEff::setDecayPos(float x, float y, float z){
  decaypos[0] = x;
  decaypos[1] = y;
  decaypos[2] = z;
}

void geoEff::setUseFixedBeamDir(bool use){
  useFixedBeamDir = use;
}

void geoEff::setVetoSizes(std::vector< float > vSizes){
  vetoSize = vSizes;
  if(verbosity){
    std::cout << "geoEff set veto sizes to ";
    for (unsigned int i = 0; i < vetoSize.size(); i++) std::cout << vetoSize[i] << " ";
  }
  std::cout << std::endl;
}
void geoEff::setVetoEnergyThresholds(std::vector< float > vThresholds){
  vetoEnergy = vThresholds;
  if(verbosity){
    std::cout << "geoEff set veto energy thresholds to ";
    for (unsigned int i = 0; i < vetoEnergy.size(); i++) std::cout << vetoEnergy[i] << " ";
  }
  std::cout << std::endl;
}

void geoEff::throwTransforms(){

  // Clear vectors
  translations[0].clear();
  translations[1].clear();
  translations[2].clear();
  rotations.clear();

  for (int dim = 0; dim < 3; dim++){
    if (not randomizeVertex[dim]){
      translations[dim].resize(0, 0);
    } else {
      translations[dim].clear();
      for (unsigned int i = 0; i < N_THROWS; i++){
        translations[dim].emplace_back(uniform(prnGenerator)*(range[dim][1]-range[dim][0])+range[dim][0]+offset[dim]);
      }
    }
  }

  rotations.clear();
  for (unsigned int i = 0; i < N_THROWS; i++){
    rotations.emplace_back((uniform(prnGenerator)-0.5)*2*M_PI);
  }

}

std::vector< Eigen::Transform<float,3,Eigen::Affine> > geoEff::getTransforms(unsigned int iStart, int iEnd){
  
  unsigned int thisEnd;
  if (iEnd < 0) thisEnd = N_THROWS;
  else if (iEnd >= 0){
    thisEnd = iEnd;
    if (thisEnd > N_THROWS) thisEnd = N_THROWS;
  }

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms;

  // Tranformations that do not depend on the throws:
  // Move vertex to coordinate system origin to apply rotation
  Eigen::Affine3f tThere(Eigen::Translation3f(Eigen::Vector3f(-vertex[0], -vertex[1], -vertex[2])));
  // Move vertex back
  Eigen::Affine3f tBack(Eigen::Translation3f(Eigen::Vector3f(vertex[0], vertex[1], vertex[2])));
  
  for (unsigned int iThrow = iStart; iThrow < thisEnd; iThrow++){
    // Vertex displacement:
    Eigen::Affine3f tThrow(Eigen::Translation3f(Eigen::Vector3f(randomizeVertex[0] ? translations[0][iThrow]-vertex[0] : 0.,
								randomizeVertex[1] ? translations[1][iThrow]-vertex[1] : 0.,
								randomizeVertex[2] ? translations[2][iThrow]-vertex[2] : 0.)));

    // Rotation
    Eigen::Affine3f rThrow;
    if (useFixedBeamDir){
      rThrow = Eigen::Affine3f(Eigen::AngleAxisf(rotations[iThrow], Eigen::Vector3f(beamdir[0], beamdir[1], beamdir[2])));
    } else {
      // Calculate rotation due to translation
      // Calculate rotation angle
      float decayToVertex[3] = {0};
      float decayToTranslated[3] = {0};
      float translationAngle = 0, magDecayToVertex = 0, magDecayToTranslated = 0;
      for (int dim = 0; dim < 3; dim++) {
        decayToVertex[dim] = vertex[dim]-decaypos[dim];
        decayToTranslated[dim] = translations[dim][iThrow]-decaypos[dim];
        
        translationAngle += (decayToVertex[dim])*(decayToTranslated[dim]);
        magDecayToVertex += pow(decayToVertex[dim], 2);
        magDecayToTranslated += pow(decayToTranslated[dim], 2);
      }
      magDecayToVertex = sqrt(magDecayToVertex);
      magDecayToTranslated = sqrt(magDecayToTranslated);
      translationAngle /= (magDecayToVertex*magDecayToTranslated);
      translationAngle = acos(translationAngle);
      
      // Calculate rotation axis
      // Cross-product
      float translationAxis[3] = {0};
      translationAxis[0] = decayToVertex[1]*decayToTranslated[2] - decayToVertex[2]*decayToTranslated[1];
      translationAxis[1] = decayToVertex[2]*decayToTranslated[0] - decayToVertex[0]*decayToTranslated[2];
      translationAxis[2] = decayToVertex[0]*decayToTranslated[1] - decayToVertex[1]*decayToTranslated[0];
      float magTranslationAxis = 0.;
      for (int dim = 0; dim < 3; dim++) magTranslationAxis += pow(translationAxis[dim], 2);
      magTranslationAxis = sqrt(magTranslationAxis);
      for (int dim = 0; dim < 3; dim++) translationAxis[dim] /= magTranslationAxis;
      
      Eigen::Affine3f rTranslation(Eigen::Affine3f(Eigen::AngleAxisf(translationAngle, Eigen::Vector3f(translationAxis[0], translationAxis[1], translationAxis[2]))));
      
      // Calculate rotation due to thrown angle
      Eigen::Affine3f rPhiThrow(Eigen::Affine3f(Eigen::AngleAxisf(rotations[iThrow], Eigen::Vector3f(decayToTranslated[0]/magDecayToTranslated, decayToTranslated[1]/magDecayToTranslated, decayToTranslated[2]/magDecayToTranslated))));
      
      // Combine
      rThrow = rPhiThrow * rTranslation;
    }
    
    // Put everything together in single transform and store.
    transforms.emplace_back(tThrow * tBack * rThrow * tThere);
  }
  
  return transforms;
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

std::vector< float > geoEff::getCurrentThrowDeps(int i, int dim){
  
  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,hitSegPoss.size()/3,Eigen::OuterStride<>(3));

  Eigen::Matrix3Xf transformedEdeps = getTransforms(i, i+1)[0] * hitSegPosOrig;
  
  int nEdeps = hitSegEdeps.size();

  std::vector< float > ret(nEdeps);
  
  for (int iDep = 0; iDep < nEdeps; iDep++){
    ret[iDep] = transformedEdeps(dim, iDep);
  }

  return ret;
}

std::vector< float > geoEff::getCurrentThrowDepsX(int i){
  return getCurrentThrowDeps(i, 0);
}
std::vector< float > geoEff::getCurrentThrowDepsY(int i){
  return getCurrentThrowDeps(i, 1);
}
std::vector< float > geoEff::getCurrentThrowDepsZ(int i){
  return getCurrentThrowDeps(i, 2);
}


std::vector< std::vector< std::vector< uint64_t > > > geoEff::getHadronContainmentThrows(bool ignore uncontained){

  // Figure out how many multiples of 64 bits needed to store output
  int n_longs = N_THROWS / 64;
  if (N_THROWS % 64) n_longs++;
  
  // Pass/fail for each set of vetoSize and vetoEnergy 
  std::vector< std::vector< std::vector< uint64_t > > > hadronContainment(vetoSize.size(), std::vector< std::vector< uint64_t > >(vetoEnergy.size(), std::vector < uint64_t >(n_longs, 0)));
  
  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,hitSegPoss.size()/3,Eigen::OuterStride<>(3));
  
  // Check if event is contained by any of the existing conditions
  if (ignore_uncontained) {
    int origContained = 0;
    std::vector< std::vector< bool > > vecOrigContained = getHadronContainmentOrigin();
    for (unsigned int i = 0; i < vetoSize.size(); i++){
      for (unsigned int j = 0; j < vetoEnergy.size(); j++){
	if (vecOrigContained[i][j]) origContained++;
      }
    }
    
    // If not, then return
    if (origContained == 0) return hadronContainment;
  }
  
  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms = getTransforms();
  // Else, loop through set of rotation translations
  for (unsigned int t = 0; t < N_THROWS; t++){
    // Apply transformation to energy deposit positions
    Eigen::Matrix3Xf transformedEdeps = transforms[t] * hitSegPosOrig;
    // Loop through conditions
    for (unsigned int i = 0; i < vetoSize.size(); i++){
      for (unsigned int j = 0; j < vetoEnergy.size(); j++){
        // Check containment and set bit
        if (isContained(transformedEdeps, hitSegEdeps, vetoSize[i], vetoEnergy[j])) {
	  hadronContainment[i][j][t/64] |= ((uint64_t)1)<<(t%64);
	}
      }
    }
  }
  
  return hadronContainment;
}

void geoEff::setSeed(int seed){
  prnGenerator = std::mt19937_64(seed);
}

std::vector< std::vector< bool > > geoEff::getHadronContainmentOrigin(){
  // Initialize return vector
  std::vector< std::vector< bool > > hadronContainment(vetoSize.size(), std::vector< bool >(vetoEnergy.size(), false));

  // Set Eigen Map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,hitSegPoss.size()/3,Eigen::OuterStride<>(3));
  
  for (unsigned int i = 0; i < vetoSize.size(); i++){
    for (unsigned int j = 0; j < vetoEnergy.size(); j++){
      if (isContained(hitSegPosOrig, hitSegEdeps, vetoSize[i], vetoEnergy[j])) hadronContainment[i][j] = true;
    }
  }

  return hadronContainment;
}

bool geoEff::isContained( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits, float vSize, float vetoEnergyThreshold ){

  float vetoEnergy = 0.;

  for (unsigned int i = 0; i < energyDeposits.size(); i++){
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

