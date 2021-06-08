// step 2
#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include <TTime.h>
#include "/sphenix/u/skurdi/CMCalibration/PHG4TpcCentralMembrane.h"
R__LOAD_LIBRARY(build/.libs/libg4tpccentralmembrane)

//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
R__LOAD_LIBRARY(libphg4hit.so)


using namespace std;

class Shifter {
public:
Shifter(TString sourcefilename);
  TFile *forward, *average;
  TH3F *hX, *hY, *hZ, *hR, *hXave, *hYave, *hZave, *hRave;  
};

Shifter::Shifter(TString sourcefilename){
  //single event distortion file
  forward=TFile::Open(sourcefilename,"READ"); 

  hX=(TH3F*)forward->Get("hIntDistortionPosX");
  hY=(TH3F*)forward->Get("hIntDistortionPosY");
  hZ=(TH3F*)forward->Get("hIntDistortionPosZ");

  hR=(TH3F*)forward->Get("hIntDistortionPosR");

  //average distortion file
  average=TFile::Open("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); 
  
  hXave=(TH3F*)average->Get("hIntDistortionX");
  hYave=(TH3F*)average->Get("hIntDistortionY");
  hZave=(TH3F*)average->Get("hIntDistortionZ");
  
  hRave=(TH3F*)average->Get("hIntDistortionR");

  //subtract average from total distortions to study fluctuations
  hX->Add(hXave,-1);
  hY->Add(hYave,-1);
  hZ->Add(hZave,-1);
  
  hR->Add(hRave,-1);
}

int CMDistortionReco() {
  Shifter *shifter;
  int nbins = 35; 
  double low = -80.0;
  double high = 80.0;
  double deltaX, deltaY, deltaZ, deltaR, deltaPhi;

  int nEvents = 3; //change based on number of tree files available in source directory
    
  //take in events
  const char * inputpattern="/sphenix/user/rcorliss/distortion_maps/2021.04/*h_Charge_*.root"; 
  
  //find all files that match the input string (includes wildcards)
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  TString sourcefilename;
  
  /*TCanvas *canvas=new TCanvas("canvas","CMDistortionReco1",1200,800);
    canvas->Divide(3,2);*/

  //canvas for time plot
  TCanvas *canvas=new TCanvas("canvas","CMDistortionReco2",400,400);
  
  TVector3 *position, *newposition;
  position = new TVector3(1.,1.,1.);
  newposition = new TVector3(1.,1.,1.);

  //histogram to compare times
    TH1F *hTimePerEvent = new TH1F("hTimePerEvent","Time Per Event; time (ms)",20,0,10000);
    
  for (int ifile=0;ifile < nEvents;ifile++){
    //for each file, find all histograms in that file
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();

    //create shifter
    shifter = new Shifter(sourcefilename);

    //call to TTime before opening ttree
    TTime now;
    now=gSystem->Now();
    //printf("the time is %lu\n",(unsigned long)before);
    unsigned long before = now;
    
    //get data from ttree
    char const *treename="cmDistHitsTree";
    TFile *input=TFile::Open(Form("cmDistHitsTree_Event%d.root", ifile), "READ");
    TTree *inTree=(TTree*)input->Get("tree");
    
    inTree->SetBranchAddress("position",&position);
    inTree->SetBranchAddress("newposition",&newposition);

    // 0 to 2pi for phi, 0 to 90 for r
    
    //for forward only
   
    TH2F *hStripesPerBin = new TH2F("hStripesPerBin","CM Stripes Per Bin (z in stripes); x (cm); y (cm)",nbins,low,high,nbins,low,high); 

    TH2F *hCartesianForward[3];
    hCartesianForward[0] = new TH2F("hForwardX","X Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[1] = new TH2F("hForwardY","Y Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianForward[2] = new TH2F("hForwardZ","Z Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);

    TH2F *hCylindricalForward[2];
    hCylindricalForward[0] = new TH2F("hForwardR","Radial Shift Forward of Stripe Centers (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalForward[1] = new TH2F("hForwardPhi","Phi Shift Forward of Stripe Centers (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    
    for (int i=0;i<inTree->GetEntries();i++){
      inTree->GetEntry(i);

      hStripesPerBin->Fill(position->X(),position->Y(),1);
      
      deltaX = (newposition->X() - position->X())*(1e4); //convert from cm to micron 
      deltaY = (newposition->Y() - position->Y())*(1e4);
      deltaZ = (newposition->Z() - position->Z())*(1e4);

      deltaR = (newposition->Perp() - position->Perp())*(1e4);
      deltaPhi = newposition->DeltaPhi(*position);

      hCartesianForward[0]->Fill(position->X(),position->Y(),deltaX);
      hCartesianForward[1]->Fill(position->X(),position->Y(),deltaY);
      hCartesianForward[2]->Fill(position->X(),position->Y(),deltaZ);

      hCylindricalForward[0]->Fill(position->X(),position->Y(),deltaR);
      hCylindricalForward[1]->Fill(position->X(),position->Y(),deltaPhi);
    }

    TH2F *hCartesianAveShift[3];
    hCartesianAveShift[0] = new TH2F("AveShiftX","Average of CM Model X over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[1] = new TH2F("AveShiftY","Average of CM Model Y over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCartesianAveShift[2] = new TH2F("AveShiftZ","Average of CM Model Z over Stripes per Bin (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 

    TH2F *hCylindricalAveShift[2];
     hCylindricalAveShift[0] = new TH2F("AveShiftRCart","Average of CM Model R over Stripes per Bin from Cartesian (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    hCylindricalAveShift[1] = new TH2F("AveShiftPhiCart","Average of CM Model Phi over Stripes per Bin from Cartesian (rad); x (cm); y (cm)",nbins,low,high,nbins,low,high); 
    
    for (int i = 0; i < 3; i ++){
      hCartesianAveShift[i]->Divide(hCartesianForward[i],hStripesPerBin);
    }

      for(int i = 0; i < nbins; i++){
      double x = low + ((high - low)/(1.0*nbins))*(i+0.5); //center of bin
      for(int j = 0; j < nbins; j++){
	double y = low + ((high - low)/(1.0*nbins))*(j+0.5); //center of bin
      
	int xbin = hCartesianAveShift[0]->FindBin(x,y);
	int ybin = hCartesianAveShift[1]->FindBin(x,y);
	double xaveshift = (hCartesianAveShift[0]->GetBinContent(xbin))*(1e-4); // converts  microns to cm 
	double yaveshift = (hCartesianAveShift[1]->GetBinContent(ybin))*(1e-4);
	
	TVector3 shifted, original;
	original.SetX(x);
	original.SetY(y);
	shifted.SetX(x+xaveshift);
	shifted.SetY(y+yaveshift);
	// have x n y above for orig
	//shifted is orig + ave shift
	
	double raveshift = (shifted.Perp() - original.Perp())*(1e4);
	double phiaveshift = shifted.DeltaPhi(original);

	//fill with r from x n y
	hCylindricalAveShift[0]->Fill(x,y,raveshift);
	hCylindricalAveShift[1]->Fill(x,y,phiaveshift);
      } 
    } 
    
    //same range and bins for each coordinate, can use hR for all, binned in cm
    int nphi = shifter->hR->GetXaxis()->GetNbins();
    int nr = shifter->hR->GetYaxis()->GetNbins();
    int nz = (shifter->hR->GetZaxis()->GetNbins())/2;
    
    double minphi = shifter->hR->GetXaxis()->GetXmin();
    double minr = shifter->hR->GetYaxis()->GetXmin();
    double minz = 5.0;
    
    double maxphi = shifter->hR->GetXaxis()->GetXmax();
    double maxr = shifter->hR->GetYaxis()->GetXmax();
    double maxz = shifter->hR->GetZaxis()->GetXmax();
    
    TH3F *hCartesianCMModel[3];
    hCartesianCMModel[0]=new TH3F("hCMModelX", "CM Model: X Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz); //rad, cm, cm
    hCartesianCMModel[1]=new TH3F("hCMModelY", "CM Model: Y Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCartesianCMModel[2]=new TH3F("hCMModelZ", "CM Model: Z Shift Forward of Stripe Centers", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);

    TH3F *hCylindricalCMModel[2];
    hCylindricalCMModel[0]=new TH3F("hCMModelRCart", "CM Model: Radial Shift Forward of Stripe Centers from Cartesian", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);
    hCylindricalCMModel[1]=new TH3F("hCMModelPhiCart", "CM Model: Phi Shift Forward of Stripe Centers from Cartesian", nphi,minphi,maxphi, nr,minr,maxr, nz,minz,maxz);  
      
    double xshift, yshift, zshift, rshiftcart, phishiftcart;
  
    for(int i = 0; i < nphi; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin

      for(int j = 0; j < nr; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin

	double x = r*cos(phi); //cm
	double y = r*sin(phi);

	for(int k = 0; k < nz; k++){
	  double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	  xshift=(hCartesianAveShift[0]->Interpolate(x,y))*(1e-4);//coordinate of your stripe
	  yshift=(hCartesianAveShift[1]->Interpolate(x,y))*(1e-4);//convert micron to cm
	  zshift=(hCartesianAveShift[2]->Interpolate(x,y))*(1e-4);

	  rshiftcart=(hCylindricalAveShift[0]->Interpolate(x,y))*(1e-4);
	  phishiftcart=hCylindricalAveShift[1]->Interpolate(x,y);
	  
	  hCartesianCMModel[0]->Fill(phi,r,z,xshift*(1-z/105.5));
	  hCartesianCMModel[1]->Fill(phi,r,z,yshift*(1-z/105.5));
	  hCartesianCMModel[2]->Fill(phi,r,z,zshift*(1-z/105.5));

	  hCylindricalCMModel[0]->Fill(phi,r,z,rshiftcart*(1-z/105.5));
	  hCylindricalCMModel[1]->Fill(phi,r,z,phishiftcart*(1-z/105.5));
	}
      }
    }
 
   
    TFile *plots;

    plots=TFile::Open(Form("CMModels_Event%d.root",ifile),"RECREATE");
    hStripesPerBin->Write(); 

    for(int i = 0; i < 3; i++){
      hCartesianForward[i]->Write();
      hCartesianAveShift[i]->Write();
      hCartesianCMModel[i]->Write();
    }

    for(int i = 0; i < 2; i++){
      hCylindricalForward[i]->Write();
      hCylindricalAveShift[i]->Write();
      hCylindricalCMModel[i]->Write();
    }
    
    plots->Close();

    //call to TTime after outputting TH3Fs
    //TTime after;
    now=gSystem->Now();
    //printf("the time is %lu\n",(unsigned long)after);
    unsigned long after = now;
    
    hTimePerEvent->Fill(after-before);
    
    /*
      for (int i = 0; i < 3; i++){
      hCartesianForward[i]->SetStats(0);
      hCartesianAveShift[i]->SetStats(0);
      }

      hStripesPerBin->SetStats(0);
      canvas->cd(1);
      hCartesianForward[0]->Draw("colz");
      canvas->cd(2);
      hCartesianForward[1]->Draw("colz");
      canvas->cd(3);
      hCartesianForward[2]->Draw("colz");
      canvas->cd(4)->Clear();
      canvas->cd(5)->Clear();
      canvas->cd(6)->Clear();
    
  
      if(ifile == 0){ 
      canvas->Print("CMDistortionReco1.pdf(","pdf");
      } else if (ifile == nEvents - 1){
      canvas->Print("CMDistortionReco1.pdf)","pdf");
      } else{
      canvas->Print("CMDistortionReco1.pdf","pdf");
      }
    */
  }

  canvas->cd();
  hTimePerEvent->Draw();

  return 0;
}
