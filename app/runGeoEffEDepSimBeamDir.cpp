#include "geoEff.h"

#include <vector>
#include <string>

#include "TFile.h"
#include "TChain.h"

// TG4Event.h assumes this... blame ROOT's MakeProject
using namespace std;

#include "TG4Event.h"

float offset[] = { 0., 305., 5. };
float collarLo[] = {-320., -120., 30. };
float collarHi[] = { 320.,  120., 470.};
float detLo[] = {-350., -150, 0.};
float detHi[] = {+350., +150, 500.};

// ND fiducial volume.
bool isFVfun ( float * vtx){
  bool inDeadRegion = false;
  for( int i = -3; i <= 3; ++i ) {
    // 0.5cm cathode in the middle of each module, plus 0.5cm buffer
    double cathode_center = i*102.1;
    if( vtx[0] - offset[0] > cathode_center - 0.75 && vtx[0] - offset[0] < cathode_center + 0.75 ) inDeadRegion = true;
    
    // 1.6cm dead region between modules (0.5cm module wall and 0.3cm pixel plane, x2)
    // don't worry about outer boundary because events are only generated in active Ar + insides
    double module_boundary = i*102.1 + 51.05;
    if( i <= 2 && vtx[0] - offset[0] > module_boundary - 1.3 && vtx[0] - offset[0] < module_boundary + 1.3 ) inDeadRegion = true;
  }
  for( int i = 1; i <= 4; ++i ) {
    // module boundaries in z are 1.8cm (0.4cm ArCLight plane + 0.5cm module wall, x2)
    // module is 102.1cm wide, but only 101.8cm long due to cathode (0.5cm) being absent in length
    // but ArCLight is 0.1cm thicker than pixel plane so it's 0.3cm difference
    // positions are off-set by 0.6 because I defined 0 to be the upstream edge based on the active volume
    // by inspecting a plot, and aparently missed by 3 mm, but whatever
    // add 8mm = 2 pad buffer due to worse position resolution in spatial dimension z compared to timing direction x
    // so total FV gap will be 1.8 + 2*0.8 = 3.4cm
    double module_boundary = i*101.8 - 0.6;
    if( vtx[2] - offset[2] > module_boundary - 1.7 && vtx[2] - offset[2] < module_boundary + 1.7 ) inDeadRegion = true;
  }
  
  return (
	  //	  abs(vtx[0] - offset[0]) < 300 &&
	  abs(vtx[0] - offset[0]) < 200 && // 1 meter vertex desert
	  abs(vtx[1] - offset[1]) < 100 &&
	  vtx[2] - offset[2] > 50 &&
	  //	  vtx[2] - offset[2] < 350 &&
	  vtx[2] - offset[2] < 450 && // No z-asymmetry in FV
	  !inDeadRegion
	  );
}


// Takes EDepsim file and produces friend TTree with geometric efficiency
// using 30 MeV threshold at 30 cm from edges of active volume using 
// fixed neutrino beam direction. Mainly for validating new code.
int main(int argc, char** argv){

  string inFname, outFname;
  
  if (argc != 3){
    cout << "Usage: runGeoEffEDepSimBeamDir [EDepSim filename] [output filename]" << endl;
  } else {
    inFname = string(argv[1]);
    outFname = string(argv[2]);
  }

  string detName = "ArgonCube";

  TChain * inChain  = new TChain ("EDepSimEvents");
  TG4Event * ev = new TG4Event();
  inChain->SetBranchAddress("Event", &ev);
  inChain->Add(inFname.c_str());
  
  TTree * effTree = new TTree("geoEff", "Geometric efficiency");
  float efficiency;
  double lepE;
  int lepPDG;
  bool isFV;
  bool isSelected;
  
  float vtxArr[3];
  double hadEdep;
  effTree->Branch("efficiency", &efficiency, "efficiency/F");
  effTree->Branch("lepE", &lepE, "lepE/D");
  effTree->Branch("lepPDG", &lepPDG, "lepPDG/I");
  effTree->Branch("hadEdep", &hadEdep, "hadEdep/D");
  effTree->Branch("vtxArr", &vtxArr, "vtxArr[3]/F");
  effTree->Branch("isFV", &isFV, "isFV/O");
  effTree->Branch("isSelected", &isSelected, "isSelected/O");
    
  geoEff * eff = new geoEff(314, true);

  eff->setNthrows(1000);

  eff->setUseFixedBeamDir(true);
  double avgCosDip = -0.09737;
  double avgBeamDir[] = {0., avgCosDip, sqrt(1. - pow(avgCosDip,2) )};

  eff->setBeamDir(avgBeamDir[0], avgBeamDir[1], avgBeamDir[2]);

  // 30 cm
  eff->setVetoSizes(vector<float>(1, 30.));
  // 30 MeV
  eff->setVetoEnergyThresholds(vector<float>(1, 30.));

  // Detector active volume
  eff->setActiveX(detLo[0], detHi[0]);
  eff->setActiveY(detLo[1], detHi[1]);
  eff->setActiveZ(detLo[2], detHi[2]);

  // Range for translation throws. Use full active volume but fix X.
  eff->setRangeX(-1, -1);
  eff->setRangeY(detLo[1], detHi[1]);
  eff->setRangeZ(detLo[2], detHi[2]);

  eff->setOffsetX(offset[0]);
  eff->setOffsetY(offset[1]);
  eff->setOffsetZ(offset[2]);
  
  vector<float> hitSegEdeps;
  vector<float> hitSegPoss;

  unsigned long int nEvents = inChain->GetEntries();
  for (unsigned long int entry = 0; entry < nEvents; entry++){
    inChain->GetEntry(entry);

    hitSegEdeps.clear();
    hitSegPoss.clear();

    TLorentzVector thisPos;
    int lepId = -999;
    lepE = -999.;
    lepPDG = -999;
    efficiency = -1.;
    isFV = false;
    isSelected = false;

    // Find lepton
    for (TG4PrimaryVertex vtx : ev->Primaries ){
      thisPos = vtx.Position; // PILE UP NOT CONSIDERED !!!
      vtxArr[0] = thisPos.X();
      vtxArr[1] = thisPos.Y();
      vtxArr[2] = thisPos.Z();
      for (TG4PrimaryParticle part : vtx.Particles){
	if ( (abs(part.PDGCode) == 11) || (abs(part.PDGCode) == 13) ) {
	  //	if ( (abs(part.PDGCode) >= 11) and (abs(part.PDGCode) <= 16) ){ // Look for any primary lepton
          lepPDG = part.PDGCode;
          lepE = part.Momentum.E();
	  lepId = part.TrackId;
	  break;
	}
      }
    }

    // Get hit segments
    vector<TG4HitSegment> hitSegs;
    for ( pair<string,vector<TG4HitSegment> > pair : ev->SegmentDetectors ){
      if (not pair.first.compare(detName) ) {
	hitSegs = pair.second;
	break;
      }
    }

    // Set Energy deposit and position arrays
    int nHitSegs = hitSegs.size();
    hitSegEdeps.reserve(nHitSegs);
    hitSegPoss.reserve(nHitSegs*3);
    
    hadEdep = 0.;
    
    for ( TG4HitSegment hitSeg : hitSegs ){
      if (hitSeg.PrimaryId == lepId) continue;
      hadEdep += hitSeg.EnergyDeposit;
      hitSegEdeps.emplace_back(hitSeg.EnergyDeposit);
      hitSegPoss.emplace_back(hitSeg.Start.X()/10.); // Change to average between start and end position!
      hitSegPoss.emplace_back(hitSeg.Start.Y()/10.); // Change to average between start and end position!
      hitSegPoss.emplace_back(hitSeg.Start.Z()/10.); // Change to average between start and end position!
    }

    eff->setVertex(vtxArr[0], vtxArr[1], vtxArr[2]);
    eff->setHitSegEdeps(hitSegEdeps);
    eff->setHitSegPoss(hitSegPoss);
    
    vector< vector< bool > > originContained = eff->getHadronContainmentOrigin();

    if (originContained[0][0]) isSelected = true;
    
    if (isFVfun (vtxArr) ) isFV = true;

    if (isSelected and isFV){
      eff->throwTransforms();

      vector< float > throwX = eff->getCurrentThrowTranslationsX();
      vector< float > throwY = eff->getCurrentThrowTranslationsY();
      vector< float > throwZ = eff->getCurrentThrowTranslationsZ();
      
      vector< vector< vector< uint64_t > > > throwsContained = eff->getHadronContainmentThrows();

      unsigned long throwInFV = 0;
      unsigned long throwSelected = 0;

      for (int iThrow = 0; iThrow < throwX.size(); iThrow++){
        float thisThrowVtx[] = {throwX[iThrow], throwY[iThrow], throwZ[iThrow]};

        if (isFVfun ( thisThrowVtx) ) throwInFV++;
        else continue;

        if ( (throwsContained[0][0][iThrow/64] >> iThrow) & 1) throwSelected++;
      }

      efficiency = ((float) throwSelected)/throwInFV;
    }

    effTree->Fill();
  }

  TFile * outFile = new TFile(outFname.c_str(), "RECREATE");
  effTree->Write();
  outFile->Close();
}
