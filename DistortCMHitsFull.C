// step 1 for both sides of tpc
#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include "/sphenix/u/skurdi/CMCalibration/PHG4TpcCentralMembrane.h"
R__LOAD_LIBRARY(build/.libs/libg4tpccentralmembrane)

//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
R__LOAD_LIBRARY(libphg4hit.so)


using namespace std;

class Shifter {
public:
  Shifter(TString truefilename, TString averagefilename);
  TVector3 Shift(TVector3 position);
  TVector3 ShiftForwardPos(TVector3 position); //only shift with forward histogram
  TVector3 ShiftForwardNeg(TVector3 position); 
  TVector3 ShiftBack(TVector3 position); //
  TFile *forward, *back, *average;
  TH3F *hPosX, *hPosY, *hPosZ, *hPosR, *hPosPhi, *hPosXave, *hPosYave, *hPosZave, *hPosRave, *hPosPhiave;
  TH3F *hNegX, *hNegY, *hNegZ, *hNegR, *hNegPhi, *hNegXave, *hNegYave, *hNegZave, *hNegRave, *hNegPhiave;
  TH3F *hXBack, *hYBack, *hZBack;  
};

Shifter::Shifter(TString truefilename, TString averagefilename){
  //single event distortion file
  forward=TFile::Open(truefilename,"READ"); 

  //positive side (z = 0 to 105.5)
  hPosX=(TH3F*)forward->Get("hIntDistortionPosX");
  hPosY=(TH3F*)forward->Get("hIntDistortionPosY");
  hPosZ=(TH3F*)forward->Get("hIntDistortionPosZ");

  hPosR=(TH3F*)forward->Get("hIntDistortionPosR");
  hPosPhi=(TH3F*)forward->Get("hIntDistortionPosP");

  hPosX=(TH3F*)forward->Get("hIntDistortionPosX");
  hPosY=(TH3F*)forward->Get("hIntDistortionPosY");
  hPosZ=(TH3F*)forward->Get("hIntDistortionPosZ");

  hPosR=(TH3F*)forward->Get("hIntDistortionPosR");
  hPosPhi=(TH3F*)forward->Get("hIntDistortionPosP");

  //negative side (z = -105.5 to 0)
  hNegX=(TH3F*)forward->Get("hIntDistortionNegX");
  hNegY=(TH3F*)forward->Get("hIntDistortionNegY");
  hNegZ=(TH3F*)forward->Get("hIntDistortionNegZ");

  hNegR=(TH3F*)forward->Get("hIntDistortionNegR");
  hNegPhi=(TH3F*)forward->Get("hIntDistortionNegP");

  hNegX=(TH3F*)forward->Get("hIntDistortionNegX");
  hNegY=(TH3F*)forward->Get("hIntDistortionNegY");
  hNegZ=(TH3F*)forward->Get("hIntDistortionNegZ");

  hNegR=(TH3F*)forward->Get("hIntDistortionNegR");
  hNegPhi=(TH3F*)forward->Get("hIntDistortionNegP");

  //average distortion file
  average=TFile::Open(averagefilename,"READ"); 
  
  hPosXave=(TH3F*)average->Get("hIntDistortionPosX");
  hPosYave=(TH3F*)average->Get("hIntDistortionPosY");
  hPosZave=(TH3F*)average->Get("hIntDistortionPosZ");
  
  hPosRave=(TH3F*)average->Get("hIntDistortionPosR");
  hPosPhiave=(TH3F*)average->Get("hIntDistortionPosP");

  hNegXave=(TH3F*)average->Get("hIntDistortionNegX");
  hNegYave=(TH3F*)average->Get("hIntDistortionNegY");
  hNegZave=(TH3F*)average->Get("hIntDistortionNegZ");
  
  hNegRave=(TH3F*)average->Get("hIntDistortionNegR");
  hNegPhiave=(TH3F*)average->Get("hIntDistortionNegP");

  
  back=TFile::Open("/sphenix/user/rcorliss/distortion_maps/averages/empty.2sides.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); 
   
  hXBack=(TH3F*)back->Get("hIntDistortionX");
  hYBack=(TH3F*)back->Get("hIntDistortionY");
  hZBack=(TH3F*)back->Get("hIntDistortionZ");
}

TVector3 Shifter::ShiftForwardPos(TVector3 position){
  double x, y, z, xshift, yshift, zshift, raveshift, phiaveshift;
  double cosphi, sinphi;
  const double mm = 1.0;
  const double cm = 10.0;
  TVector3 shiftposition;

  x= position.X();
  y= position.Y();
  z= position.Z();

  double r=position.Perp();
  double phi=position.Phi();
  if(position.Phi() < 0.0){
    phi = position.Phi() + 2.0*TMath::Pi(); 
  }

  //distort coordinate of stripe 
  xshift=hPosX->Interpolate(phi,r,z);
  yshift=hPosY->Interpolate(phi,r,z);
  zshift=hPosZ->Interpolate(phi,r,z);
  
  //subtract average from total distortions
  xshift -= hPosXave->Interpolate(phi,r,z);
  yshift -= hPosYave->Interpolate(phi,r,z);
  zshift -= hPosZave->Interpolate(phi,r,z);
  
  TVector3 forwardshiftpos(x+xshift,y+yshift,z+zshift);

  return forwardshiftpos;
}

TVector3 Shifter::ShiftForwardNeg(TVector3 position){
  double x, y, z, xshift, yshift, zshift, raveshift, phiaveshift;
  double cosphi, sinphi;
  const double mm = 1.0;
  const double cm = 10.0;
  TVector3 shiftposition;

  x= position.X();
  y= position.Y();
  z= position.Z();

  double r=position.Perp();
  double phi=position.Phi();
  if(position.Phi() < 0.0){
    phi = position.Phi() + 2.0*TMath::Pi(); 
  }

  //distort coordinate of stripe 
  xshift=hNegX->Interpolate(phi,r,z);
  yshift=hNegY->Interpolate(phi,r,z);
  zshift=hNegZ->Interpolate(phi,r,z);
  
  //subtract average from total distortions
  xshift -= hNegXave->Interpolate(phi,r,z);
  yshift -= hNegYave->Interpolate(phi,r,z);
  zshift -= hNegZave->Interpolate(phi,r,z);
  
  TVector3 forwardshiftneg(x+xshift,y+yshift,z+zshift);

  return forwardshiftneg;
}

TVector3 Shifter::ShiftBack(TVector3 forwardshift){
double x, y, z, xshift, yshift, zshift;
  const double mm = 1.0;
  const double cm = 10.0;
  TVector3 shiftposition;

  x= forwardshift.X();
  y= forwardshift.Y();
  z= forwardshift.Z();

  double rforward=forwardshift.Perp();
  double phiforward=forwardshift.Phi();
  if(forwardshift.Phi() < 0.0){
    phiforward += 2.0*TMath::Pi();
  }

  double xshiftback=-1*hXBack->Interpolate(phiforward,rforward,z);
  double yshiftback=-1*hYBack->Interpolate(phiforward,rforward,z);
  double zshiftback=-1*hZBack->Interpolate(phiforward,rforward,z);
  
  shiftposition.SetXYZ(x+xshiftback,y+yshiftback,z+zshiftback);

  return shiftposition;
}

TVector3 Shifter::Shift(TVector3 position){
  
  return ShiftBack(ShiftForwardPos(position));
}

void TestPlots(int nEvents);

int DistortCMHitsFull(int nMaxEvents = -1) {
  Shifter *shifter;
  PHG4TpcCentralMembrane stripes;
  vector<PHG4Hitv1*> Hits = stripes.PHG4Hits;
  double x, y, z;
  TVector3 *position, *newpositionPos, *newpositionNeg;

  position = new TVector3(1.,1.,1.);
  newpositionPos = new TVector3(1.,1.,1.);
  newpositionNeg = new TVector3(1.,1.,1.);

  //take in events
  const char * inputpattern="/sphenix/user/rcorliss/distortion_maps/2021.04/*h_Charge_*.root"; 
  
  //find all files that match the input string (includes wildcards)
  TFileCollection *filelist=new TFileCollection(); 
  filelist->Add(inputpattern);
  TString truefilename;

  TString averagefilename = "/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";

  int nEvents; 
  //how many events
  if (nMaxEvents<0){
    nEvents=filelist->GetNFiles();
  } else if(nMaxEvents<filelist->GetNFiles()){
    nEvents=nMaxEvents;
  } else {
    nEvents= filelist->GetNFiles();
  }
    
  for (int ifile=0;ifile < nEvents;ifile++){
    //for each file, find all histograms in that file
    truefilename=((TFileInfo*)(filelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();

    //create shifter
    shifter = new Shifter(truefilename, averagefilename);

    //set up TFile for TTree
    TFile *output = TFile::Open(Form("cmDistHitsFullTree_Event%d.root", ifile),"RECREATE");
    
    //set up TTree to store position and newposition
    TTree *cmHitsTree = new TTree("tree","cmDistHitsTree");
    cmHitsTree->Branch("position","TVector3",&position);
    cmHitsTree->Branch("newpositionPos","TVector3",&newpositionPos);
    cmHitsTree->Branch("newpositionNeg","TVector3",&newpositionNeg);
  
    for (int i = 0; i < Hits.size(); i++){
      //store each stripe center's coordinates in position vector
      x = (Hits[i]->get_x(0) + Hits[i]->get_x(1))/2; 
      y = (Hits[i]->get_y(0) + Hits[i]->get_y(1))/2;
      z = 0.0;
      position->SetXYZ(x,y,z);
      
      //shift position
      *newpositionPos = shifter->ShiftForwardPos(*position);
      *newpositionNeg = shifter->ShiftForwardNeg(*position);

      //store position and newposition in tree
      cmHitsTree->Fill();  
    }

    //save tree
    cmHitsTree->Write();
    output->Close();
  }

  //use test function to check if data opens properly
  //TestPlots(nEvents); 

  return 0;
}


void TestPlots(int nEvents){
  //setup for models
  int nbins = 35;
  double low = -80.0;
  double high = 80.0;
  double deltaX, deltaY, deltaZ, deltaR, deltaPhi; 

  TH2F *hStripesPerBin = new TH2F("hStripesPerBin","CM Stripes Per Bin (z in stripes); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    
  TH2F *hCartesianForward[3];
  hCartesianForward[0] = new TH2F("hForwardX","X Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianForward[1] = new TH2F("hForwardY","Y Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCartesianForward[2] = new TH2F("hForwardZ","Z Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    
  TH2F *hCylindricalForward[2];
  hCylindricalForward[0] = new TH2F("hForwardR","Radial Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
  hCylindricalForward[1] = new TH2F("hForwardPhi","Phi Shift Forward of Stripe Centers (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high);

  TCanvas *canvas=new TCanvas("canvas","DistortCMHitsTest",1200,800);
  canvas->Divide(3,2);
  
  TVector3 *positionT, *newpositionT;

  positionT = new TVector3(1.,1.,1.);
  newpositionT = new TVector3(1.,1.,1.);
  
  for (int ifile=0;ifile < nEvents;ifile++){ 

    //Clear out histograms
    hStripesPerBin->Reset();
    hCartesianForward[0]->Reset();
    hCartesianForward[1]->Reset();
    hCartesianForward[2]->Reset();
    hCylindricalForward[0]->Reset();
    hCylindricalForward[1]->Reset();
    
    //Get data from TTree  
    char const *treename="cmDistHitsTree";
    TFile *input=TFile::Open(Form("cmDistHitsFullTree_Event%d.root", ifile),"READ");
    TTree *inTree=(TTree*)input->Get("tree");
    
    inTree->SetBranchAddress("position",&positionT);
    inTree->SetBranchAddress("newpositionNeg",&newpositionT);
    
    for (int i=0;i<inTree->GetEntries();i++){
      inTree->GetEntry(i);

      double r = positionT->Perp();
    
      double phi = positionT->Phi();
      if(positionT->Phi() < 0.0){
	phi = positionT->Phi() + TMath::TwoPi(); 
      }

      hStripesPerBin->Fill(positionT->X(),positionT->Y(),1);
      
      deltaX = (newpositionT->X() - positionT->X())*(1e4); //convert from cm to micron 
      deltaY = (newpositionT->Y() - positionT->Y())*(1e4);
      deltaZ = (newpositionT->Z() - positionT->Z())*(1e4);

      deltaR = (newpositionT->Perp() - positionT->Perp())*(1e4);
      deltaPhi = newpositionT->DeltaPhi(*positionT);

      hCartesianForward[0]->Fill(positionT->X(),positionT->Y(),deltaX);
      hCartesianForward[1]->Fill(positionT->X(),positionT->Y(),deltaY);
      hCartesianForward[2]->Fill(positionT->X(),positionT->Y(),deltaZ);

      hCylindricalForward[0]->Fill(positionT->X(),positionT->Y(),deltaR);
      hCylindricalForward[1]->Fill(positionT->X(),positionT->Y(),deltaPhi);
    }
    
    input->Close();
    
    canvas->cd(1);
    hCartesianForward[0]->Draw("colz");
    canvas->cd(2);
    hCartesianForward[1]->Draw("colz");
    canvas->cd(3);
    hCartesianForward[2]->Draw("colz");
    canvas->cd(4);
    hCylindricalForward[0]->Draw("colz");
    canvas->cd(5);
    hCylindricalForward[1]->Draw("colz");
    canvas->cd(6)->Clear();
    
  
    if(ifile == 0){ 
      canvas->Print("DistortCMHitsTest.pdf(","pdf");
    } else if (ifile == nEvents - 1){
      canvas->Print("DistortCMHitsTest.pdf)","pdf");
    } else{
      canvas->Print("DistortCMHitsTest.pdf","pdf");
    }
  }
  
}
