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
  
  // Range of coordinates to randomize over. Set negative for no randomization.
  float range[3][2];

  // Active volume:
  float active[3][2];

  // Detector coordinates offset:
  float offset[3];
  
  // Beam direction:
  float beamdir[3];
  
  // Decay position used to calculate neutrino direction.
  float decaypos[3];

  // Flag to determine whether to use fixed beam direction (appropriate for FD) or calculate direction from vertex position and decay point.
  bool useFixedBeamDir;
  
  // Number of throws. Should be a multiple of 64 to optimize output efficiency
  unsigned long N_THROWS;

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

  void setNthrows(unsigned long n);
  
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
  void setDecayPos(float x, float y, float z);
  void setUseFixedBeamDir(bool use);

  void setVetoSizes(std::vector< float > vSizes);
  void setVetoEnergyThresholds(std::vector< float > vThresholds);
  
  void throwTransforms();

  std::vector< float > getCurrentThrowTranslationsX();
  std::vector< float > getCurrentThrowTranslationsY();
  std::vector< float > getCurrentThrowTranslationsZ();
  std::vector< float > getCurrentThrowRotations();

  // Pass/fail for each set of vetoSize and vetoEnergy. Storing in TTree as uint64_t seems to take ~half the space of the equivalent vector< bool >.
  std::vector< std::vector< std::vector< uint64_t > > > getHadronContainmentThrows();

  // Get pass/fail containment criterion for original event
  std::vector< std::vector< bool > > getHadronContainmentOrigin();
  
  void setSeed(int seed);
  
};

#endif
