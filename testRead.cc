#include <iostream>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TH2F.h>
#include <TH3F.h>

using namespace std;

#include <TG4Event.h>
#include <TG4HitSegment.h>
#include <TG4PrimaryParticle.h>
#include <TG4PrimaryVertex.h>
#include <TG4Trajectory.h>
#include <TG4TrajectoryPoint.h>

#include <GENIE/FluxDrivers/GSimpleNtpFlux.h>

#include <Eigen/Dense>

float offset[] = { 0., 305., 5. };
float collarLo[] = {-320., -120., 30. };
float collarHi[] = { 320.,  120., 470.};
float detLo[] = {-350., -150, 0.};
float detHi[] = {+350., +150, 500.};

float containmentCut = 30;

int N = 500; // 10000 random rotations per event

bool isFVfun ( float * vtx){
  bool inDeadRegion = false;
  for( int i = -3; i <= 3; ++i ) {
    // 0.5cm cathode in the middle of each module, plus 0.5cm buffer
    double cathode_center = i*102.1;
    if( vtx[0]/10. - offset[0] > cathode_center - 0.75 && vtx[0]/10. - offset[0] < cathode_center + 0.75 ) inDeadRegion = true;
    
    // 1.6cm dead region between modules (0.5cm module wall and 0.3cm pixel plane, x2)
    // don't worry about outer boundary because events are only generated in active Ar + insides
    double module_boundary = i*102.1 + 51.05;
    if( i <= 2 && vtx[0]/10. - offset[0] > module_boundary - 1.3 && vtx[0]/10. - offset[0] < module_boundary + 1.3 ) inDeadRegion = true;
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
    if( vtx[2]/10. - offset[2] > module_boundary - 1.7 && vtx[2]/10. - offset[2] < module_boundary + 1.7 ) inDeadRegion = true;
  }
  
  return (
	  //	  abs(vtx[0]/10. - offset[0]) < 300 &&
	  abs(vtx[0]/10. - offset[0]) < 200 && // 1 meter vertex desert
	  abs(vtx[1]/10. - offset[1]) < 100 &&
	  vtx[2]/10. - offset[2] > 50 &&
	  //	  vtx[2]/10. - offset[2] < 350 &&
	  vtx[2]/10. - offset[2] < 450 && // No z-asymmetry in FV
	  !inDeadRegion
	  );
}


bool isContained( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits ){
  
  float collarEnergy = 0.;

  for (int i = 0; i < energyDeposits.size(); i++){
    for (int dim = 0; dim < 3; dim++){
      if ( hitSegments(dim, i)/10.-offset[dim] < collarLo[dim] and hitSegments(dim, i)/10.-offset[dim] > detLo[dim]) 
	{
	  collarEnergy += energyDeposits[i];
	  break;
	}
      else if ( hitSegments(dim, i)/10.-offset[dim] > collarHi[dim] and hitSegments(dim, i)/10.-offset[dim] < detHi[dim])
	{
	  collarEnergy += energyDeposits[i];
	  break;
	}
    }
  }

  return collarEnergy < containmentCut;
	   
}

std::pair<TH2F*, TH2F*> displayTwoViews(Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits,  std::string histName) {

  TH2F * hXY = new TH2F((histName+"_XY").c_str(), "; x [cm]; y [cm]", 200, -500., 500., 200, -500., 500.);
  TH2F * hZY = new TH2F((histName+"_ZY").c_str(), "; z [cm]; y [cm]", 200, -250., 750., 200, -500., 500.);

  for (int i = 0; i < energyDeposits.size(); i++){
    hXY->Fill( hitSegments(0, i)/10.-offset[0], hitSegments(1, i)/10.-offset[1], energyDeposits[i] );
    hZY->Fill( hitSegments(2, i)/10.-offset[2], hitSegments(1, i)/10.-offset[1], energyDeposits[i] );
  }

  return std::pair<TH2F*, TH2F*> (hXY, hZY);
}

TH3F* display3D(Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits,  std::string histName) {

  TH3F * hXYZ = new TH3F((histName+"_XYZ").c_str(), "; x [cm]; y [cm]; z [cm]", 200, -500., 500., 200, -500., 500., 200, -250., 750.);

  for (int i = 0; i < energyDeposits.size(); i++){
    hXYZ->Fill( hitSegments(0, i)/10.-offset[0], hitSegments(1, i)/10.-offset[1], hitSegments(2, i)/10.-offset[2], energyDeposits[i]);
  }

  return hXYZ;
}

int main (int argc, char **argv){

  // For debugging purposes
  std::vector<TH2F*> evDisplaysXY;
  std::vector<TH2F*> evDisplaysZY;
  std::vector<TH3F*> evDisplaysXYZ;
  int nEventToDisplay = 0; // What event will be displayed -- this is index of selected event
  int NumberOfSamplesToDisplay = 20; // Number of samples saved
  // end debug

  double avgCosDip = -0.09737;
  double avgBeamDir[] = {0., avgCosDip, sqrt(1. - pow(avgCosDip,2) )};

  std::string inFname = std::string(argv[1]);
  std::string inGenieFname = std::string(argv[2]);
  std::string outFname = std::string(argv[3]);
  
  //  std::string inFname = "/pnfs/dune/persistent/users/marshalc/CAF/edepNewFluxv2/LAr/FHC/00/LAr.neutrino.[0-9].edepsim.root";
  std::string detName = "ArgonCube";
  
  TRandom3 * ran = new TRandom3();

  //  /pnfs/dune/persistent/users/marshalc/CAF/genieNewFluxv2/LAr/FHC/00/LAr.neutrino.0.ghep.root
  TChain * inGenieChain = new TChain("gtree");
  genie::flux::GSimpleNtpEntry * simple = new genie::flux::GSimpleNtpEntry();

  double Ev;
  inGenieChain->SetBranchAddress("simple", &simple);
  inGenieChain->Add(inGenieFname.c_str());

  TChain * inChain  = new TChain ("EDepSimEvents");
  TG4Event * ev = new TG4Event();
  inChain->SetBranchAddress("Event", &ev);
  inChain->Add(inFname.c_str());
  
  inChain->AddFriend(inGenieChain);
  
  //  inFname = "/pnfs/dune/persistent/users/marshalc/CAF/edepNewFluxv2/LAr/FHC/00/LAr.neutrino.[0-9][0-9].edepsim.root";
  //  inChain->Add(inFname.c_str());
  //  inFname = "/pnfs/dune/persistent/users/marshalc/CAF/edepNewFluxv2/LAr/FHC/00/LAr.neutrino.[0-9][0-9][0-9].edepsim.root";
  //  inChain->Add(inFname.c_str());

  TTree * effTree = new TTree("geoEff", "Geometric efficiency");
  float efficiency;
  double lepE;
  int lepPDG;
  bool isFV;
  bool isSelected;
  
  double vtxArr[3];
  double hadEdep;
  effTree->Branch("efficiency", &efficiency, "efficiency/F");
  effTree->Branch("lepE", &lepE, "lepE/D");
  effTree->Branch("lepPDG", &lepPDG, "lepPDG/I");
  effTree->Branch("Ev", &Ev, "Ev/D");
  effTree->Branch("hadEdep", &hadEdep, "hadEdep/D");
  effTree->Branch("vtxArr", &vtxArr, "vtxArr[3]/D");
  effTree->Branch("isFV", &isFV, "isFV/O");
  effTree->Branch("isSelected", &isSelected, "isSelected/O");

  std::vector<float> hitSegEdeps;
  std::vector<float> hitSegPoss;

  unsigned long int nEvents = inChain->GetEntries();
  
  int nSelected = 0;
  
  nEvents = 100;

  for (unsigned long int entry = 0; entry < nEvents; entry++){
    
    std::cout << "Starting event " << entry << std::endl;
    
    inChain->GetEntry(entry);
    inGenieChain->GetEntry(entry);
    hitSegEdeps.clear();
    hitSegPoss.clear();
    
    TLorentzVector thisPos;
    int lepId = -999;
    lepE = -999.;
    lepPDG = -999;

    Ev = simple->E;
    
    std::cout << Ev << std::endl;

    for (TG4PrimaryVertex vtx : ev->Primaries ){
      thisPos = vtx.Position; // PILE UP NOT CONSIDERED !!!
      vtxArr[0] = thisPos.X();
      vtxArr[1] = thisPos.Y();
      vtxArr[2] = thisPos.Z();
      for (TG4PrimaryParticle part : vtx.Particles){
	if ( (abs(part.PDGCode) == 11) || (abs(part.PDGCode) == 13) ) {
	  //	if ( (abs(part.PDGCode) >= 11) and (abs(part.PDGCode) <= 16) ){ // Look for any primary lepton
	  lepId = part.TrackId;
	  lepPDG = part.PDGCode;
	  lepE = part.Momentum.E();
	  break;
	}
      }
    }
    
    if (lepId == -999) continue;

    //    std::cout << "PRIMARY LEPTON TRACK ID " << lepId << std::endl;
    
    //    TH2F* hBeforeRot = new TH2F("before", "before;x;y", 600, -6000, 6000,  600, -6000, 6000);
    //    TH2F* hAfterRot = new TH2F("after", "after;x;y", 600, -6000, 6000,  600, -6000, 6000);
    
    std::vector<TG4HitSegment> hitSegs;
    for ( std::pair<std::string,std::vector<TG4HitSegment> > pair : ev->SegmentDetectors ){
      if (not pair.first.compare(detName) ) {
	hitSegs = pair.second;
	break;
      }
    }
    
    int nHitSegs = hitSegs.size();
    //    std::cout << "NUMBER OF HIT SEGMENTS " << nHitSegs << std::endl;
    hitSegEdeps.reserve(nHitSegs);
    hitSegPoss.reserve(nHitSegs*3);
    
    hadEdep = 0.;
    
    for ( TG4HitSegment hitSeg : hitSegs ){
      if (hitSeg.PrimaryId == lepId) continue;
      hadEdep += hitSeg.EnergyDeposit;
      hitSegEdeps.emplace_back(hitSeg.EnergyDeposit);
      hitSegPoss.emplace_back(hitSeg.Start.X());
      hitSegPoss.emplace_back(hitSeg.Start.Y());
      hitSegPoss.emplace_back(hitSeg.Start.Z());
    }
    
    const int nNonLepHitSegs = hitSegEdeps.size();
    
    Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,nNonLepHitSegs,Eigen::OuterStride<>(3)); 
    
    float thisPosArr[] = {thisPos.X(),thisPos.Y(),thisPos.Z()}; // clean this up;....
    if (isContained(hitSegPosOrig, hitSegEdeps)) isSelected = true;
    else isSelected = false;
    
    if (isFVfun( thisPosArr )) isFV = true;
    else isFV = false;
    

    if ( (not isSelected) or (not isFV ) )  {
      efficiency = -1;
    } else {
      efficiency = 0.;
      //    }  
      if ( nSelected == nEventToDisplay ) {
	std::pair<TH2F*, TH2F*> hPair = displayTwoViews(hitSegPosOrig, hitSegEdeps, std::string("original"));
	evDisplaysXY.push_back(hPair.first);
	evDisplaysZY.push_back(hPair.second);
	evDisplaysXYZ.push_back(display3D(hitSegPosOrig, hitSegEdeps, std::string("original")));
      }
    
      //      std::cout << "Starting loop " << std::endl;
      for (int rep = 0; rep < N; rep++){
	
	// Randomize direction
	double phi = (ran->Rndm()-0.5)*2*TMath::Pi();
	
	// Randomize position Y, Z (keep X)
	float randomVtx[] = {-1e10, -1e10, -1e10};
	while(not isFVfun(randomVtx) ){
	  randomVtx[0] = thisPos.X();
	  randomVtx[1] = 10*((ran->Rndm()-0.5)*2*100 + offset[1]);
	  randomVtx[2] = 10*((ran->Rndm()-0.5)*2*(450-50)/2. + 200. + offset[2]);
	  //	  std::cout << "Trying random vertex " << randomVtx[0] << " "<< randomVtx[2] << " "<< randomVtx[1] << std::endl;
	}
	
	//	std::cout << "Got a random vertex in the fiducial volume " << std::endl;
	
	Eigen::Affine3f randomDisplacement(Eigen::Translation3f(Eigen::Vector3f(randomVtx[0]-thisPos.X(),randomVtx[1]-thisPos.Y(),randomVtx[2]-thisPos.Z())));
	
	
	Eigen::Affine3f tThere(Eigen::Translation3f(Eigen::Vector3f(-thisPos.X(),-thisPos.Y(),-thisPos.Z())));
	Eigen::Affine3f r = Eigen::Affine3f(Eigen::AngleAxisf(phi, Eigen::Vector3f(avgBeamDir[0], avgBeamDir[1], avgBeamDir[2])));
	Eigen::Affine3f tBack(Eigen::Translation3f(Eigen::Vector3f(thisPos.X(),thisPos.Y(),thisPos.Z())));
	
	Eigen::Transform<float,3,Eigen::Affine> m = randomDisplacement * tBack * r * tThere;
	
	bool throwSelected = false;
	if (isContained(m * hitSegPosOrig, hitSegEdeps)) throwSelected = true;
 
	if ( (nSelected == nEventToDisplay) && (rep < NumberOfSamplesToDisplay) ) {
	  if (throwSelected) {
	    std::pair<TH2F*, TH2F*> hPair = displayTwoViews(m*hitSegPosOrig, hitSegEdeps, ("sample_"+std::to_string(rep)+"_sel1_Ev"+std::to_string(Ev)));
	    evDisplaysXY.push_back(hPair.first);
	    evDisplaysZY.push_back(hPair.second);
	    evDisplaysXYZ.push_back(display3D(hitSegPosOrig, hitSegEdeps,("sample_"+std::to_string(rep)+"_sel1_Ev"+std::to_string(Ev))));
	  } else {
	    std::pair<TH2F*, TH2F*> hPair = displayTwoViews(m*hitSegPosOrig, hitSegEdeps, ("sample_"+std::to_string(rep)+"_sel0_Ev"+std::to_string(Ev)));
	    evDisplaysXY.push_back(hPair.first);
	    evDisplaysZY.push_back(hPair.second);
	    evDisplaysXYZ.push_back(display3D(hitSegPosOrig, hitSegEdeps,("sample_"+std::to_string(rep)+"_sel0_Ev"+std::to_string(Ev))));
	  }
	}
	
	if (throwSelected) efficiency += 1.;
	
      }
      
      efficiency /= N;
      std::cout << "Loop ends. Efficiency: " << efficiency << std::endl;
      nSelected++;
    }
    effTree->Fill();
    
  }
  
  TFile *out = new TFile(outFname.c_str(), "RECREATE");
  effTree->Write();

  
  //  TFile *evDisplaysOut = new TFile("evDisplaysOut.root", "RECREATE");
  for (int i = 0; i < evDisplaysXY.size(); i++){
    std::cout << "Writing histogram " << i << " " << evDisplaysXY[i]->GetEntries() << std::endl;
    evDisplaysXY[i]->Write();
    evDisplaysZY[i]->Write();
    evDisplaysXYZ[i]->Write();
  }
  //  evDisplaysOut->Close();
  out->Close();
}


