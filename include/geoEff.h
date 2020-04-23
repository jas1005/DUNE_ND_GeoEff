#ifndef GEOEFF_H
#define GEOEFF_H

#include <vector>
#include <random>

#include <Eigen/Dense>

class geoEff
{

 private:

  // Event vertex
  float vertex[3];
  // Vector to store energy deposits corresponding to hit segments
  std::vector<float> hitSegEdeps;
  // Vector to store positions of hit segments 
  std::vector<float> hitSegPoss;
  
  // Range of coordinates to randomize over. Set negative for no randomization
  float range[3][2];

  // Active volume:
  float active[3][2];

  // Detector coordinates offset:
  float offset[3];
  
  // Beam direction:
  float beamdir[3];
  
  // Hard-coded to optimize output size (N_THROWS matches hadThrowResult's uint32_t)
  const unsigned long N_THROWS = sizeof(uint32_t)*8;

  // Veto size. Each element defines a veto region in cm from the end of the active volume in all directions.
  std::vector < float > vetoSize;
  // Veto energy. Each element defines a veto energy threshold.
  std::vector < float > vetoEnergy;
  
  // Current throws
  std::vector< float > translations[3];
  std::vector< float > rotations;

  // Corresponding vector of Eigen transforms
  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms;

  bool verbosity;

  std::mt19937_64 prnGenerator;
  std::uniform_real_distribution<> uniform;

  bool isContained( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits, float vSize, float vetoEnergyThreshold );
  
 public:
  geoEff(int seed, bool verbose = false);
  ~geoEff(){;}

  void setVertex(float x, float y, float z);
  void setHitSegEdeps(std::vector<float> thishitSegEdeps);
  void setHitSegPoss(std::vector<float> thishitSegPoss);
  
  void setRangeX(float xmin, float xmax);
  void setRangeY(float ymin, float ymax);
  void setRangeZ(float zmin, float zmax);

  void setActiveX(float xmin, float xmax);
  void setActiveY(float ymin, float ymax);
  void setActiveZ(float zmin, float zmax);

  void setOffsetX(float x);
  void setOffsetY(float y);
  void setOffsetZ(float z);
  
  void setBeamDir(float xdir, float ydir, float zdir);

  void setVetoSizes(std::vector< float > vSizes);
  void setVetoEnergyThresholds(std::vector< float > vThresholds);
  
  void throwTransforms();

  std::vector< float > getCurrentThrowTranslationsX();
  std::vector< float > getCurrentThrowTranslationsY();
  std::vector< float > getCurrentThrowTranslationsZ();
  std::vector< float > getCurrentThrowRotations();

  // Pass/fail for each set of vetoSize and vetoEnergy
  std::vector< std::vector< uint32_t > > getHadronContainment();

  void setSeed(int seed);
  
};

#endif
