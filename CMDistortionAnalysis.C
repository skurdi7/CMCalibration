//step 3
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
Shifter(TString sourcefilename);
  TFile *forward, *average;
  TH3F *hX, *hY, *hZ, *hR, *hPhi, *hXave, *hYave, *hZave, *hRave, *hPhiave;  
};

Shifter::Shifter(TString sourcefilename){
  //single event distortion file
  forward=TFile::Open(sourcefilename,"READ"); 

  hX=(TH3F*)forward->Get("hIntDistortionPosX");
  hY=(TH3F*)forward->Get("hIntDistortionPosY");
  hZ=(TH3F*)forward->Get("hIntDistortionPosZ");

  hR=(TH3F*)forward->Get("hIntDistortionPosR");
  hPhi=(TH3F*)forward->Get("hIntDistortionPosP");

  //average distortion file
  average=TFile::Open("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); 
  
  hXave=(TH3F*)average->Get("hIntDistortionX");
  hYave=(TH3F*)average->Get("hIntDistortionY");
  hZave=(TH3F*)average->Get("hIntDistortionZ");
  
  hRave=(TH3F*)average->Get("hIntDistortionR");
  hPhiave=(TH3F*)average->Get("hIntDistortionP");
 
  //subtract average from total distortions to study fluctuations
  hX->Add(hXave,-1);
  hY->Add(hYave,-1);
  hZ->Add(hZave,-1);
  
  hR->Add(hRave,-1);
  hPhi->Add(hPhiave,-1);
}

int CMDistortionAnalysis() {
  Shifter *shifter;
  int nbins = 35; 
  double x, y, z;
  double low = -80.0;
  double high = 80.0;
  double deltaX, deltaY, deltaZ, deltaR, deltaPhi;
  int nEvents = 3;
  
  TCanvas *canvas=new TCanvas("canvas","CMDistortionAnalysis",2000,3000);

  int nsumbins = 20;
  int minsum = -10;
  int maxsum = 10;
  
  //set up summary plots
  TH1F *hDifferenceMeanR = new TH1F("hDifferenceMeanR", "Average Difference between R Model and True of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    TH1F *hDifferenceStdDevR = new TH1F("hDifferenceStdDevR", "Std Dev of Difference between R Model and True of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hTrueMeanR = new TH1F("hTrueMeanR", "Mean True R Distortion Model of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    TH1F *hTrueStdDevR = new TH1F("hTrueStdDevR", "Std Dev of True R Distortion Model of All Events (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hDifferenceMeanPhi = new TH1F("hDifferenceMeanPhi", "Average Difference between Phi Model and True of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    TH1F *hDifferenceStdDevPhi = new TH1F("hDifferenceStdDevPhi", "Std Dev of Difference between Phi Model and True of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    
    TH1F *hTrueMeanPhi = new TH1F("hTrueMeanPhi", "Mean True Phi Distortion Model of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    TH1F *hTrueStdDevPhi = new TH1F("hTrueStdDevPhi", "Std Dev of True Phi Distortion Model of All Events (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);

    const char * inputpattern="/sphenix/user/rcorliss/distortion_maps/2021.04/*h_Charge_*.root"; //updated
    
  //find all files that match the input string (includes wildcards)
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  TString sourcefilename;

  for (int ifile=0;ifile < nEvents;ifile++){
    //for each file, find all histograms in that file.
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();

    //create shifter
    shifter = new Shifter(sourcefilename);
    
    TFile *plots;

    plots=TFile::Open(Form("CMModels_Event%d.root",ifile),"READ");
   
    TH3F *hCartCMModel[3];
    hCartCMModel[0]=(TH3F*)plots->Get("hCMModelX");
    hCartCMModel[1]=(TH3F*)plots->Get("hCMModelY");
    hCartCMModel[2]=(TH3F*)plots->Get("hCMModelZ");

    TH3F *hCylCMModel[2];
    hCylCMModel[0]=(TH3F*)plots->Get("hCMModelRCart");
    hCylCMModel[1]=(TH3F*)plots->Get("hCMModelPhiCart");
   
    //for forward only

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

    double rshiftcart, phishiftcart;

    int ndiff = 300;
    int mindiff = -20;
    int maxdiff = 20;
  
    TH1F *hCartesianShiftDifference[3];
    hCartesianShiftDifference[0] = new TH1F("hShiftDifferenceX", "Difference between CM Model X and True (R > 30); #Delta X (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifference[1] = new TH1F("hShiftDifferenceY", "Difference between CM Model Y and True (R > 30); #Delta Y (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifference[2] = new TH1F("hShiftDifferenceZ", "Difference between CM Model Z and True (R > 30); #Delta Z (#mum)", ndiff, mindiff, maxdiff);
  
    TH1F *hCylindricalShiftDifference[2];
    hCylindricalShiftDifference[0] = new TH1F("hShiftDifferenceRCart", "Difference between CM Model R from Cartesian and True (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    hCylindricalShiftDifference[1] = new TH1F("hShiftDifferencePhiCart", "Difference between CM Model Phi from Cartesian and True (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);

    TH1F *hRShiftTrue = new TH1F("hRShiftTrue", "True R Distortion Model (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    TH1F *hPhiShiftTrue = new TH1F("hPhiShiftTrue", "True Phi Distortion Model (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);
  
    TH2F *hCartesianDiff[6];
    hCartesianDiff[0] = new TH2F("hDiffXYX", "Difference in XY for CM Model X; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianDiff[1] = new TH2F("hDiffRZX", "Difference in RZ for CM Model X; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianDiff[2] = new TH2F("hDiffXYY", "Difference in XY for CM Model Y; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianDiff[3] = new TH2F("hDiffRZY", "Difference in RZ for CM Model Y; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianDiff[4] = new TH2F("hDiffXYZ", "Difference in XY for CM Model Z; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianDiff[5] = new TH2F("hDiffRZZ", "Difference in RZ for CM Model Z; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hCylindricalDiff[4];
    hCylindricalDiff[0] = new TH2F("hDiffXYRCart", "Difference in XY for CM Model R from Cartesian; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalDiff[1] = new TH2F("hDiffRZRCart", "Difference in RZ for CM Model R from Cartesian; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);
    hCylindricalDiff[2] = new TH2F("hDiffXYPhiCart", "Difference in XY for CM Model Phi from Cartesian; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalDiff[3] = new TH2F("hDiffRZPhiCart", "Difference in RZ for CM Model Phi from Cartesian; z (cm); r (cm)",nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hCartesianAveDiff[6];
    hCartesianAveDiff[0] = new TH2F("hAveDiffXYX", "X Model - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianAveDiff[1] = new TH2F("hAveDiffRZX", "X Model - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianAveDiff[2] = new TH2F("hAveDiffXYY", "Y Model - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianAveDiff[3] = new TH2F("hAveDiffRZY", "Y Model - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCartesianAveDiff[4] = new TH2F("hAveDiffXYZ", "Z Model - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCartesianAveDiff[5] = new TH2F("hAveDiffRZZ", "Z Model - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hCylindricalAveDiff[4];
    hCylindricalAveDiff[0] = new TH2F("hAveDiffXYRCart", "R Model from Cartesian - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalAveDiff[1] = new TH2F("hAveDiffRZRCart", "R Model from Cartesian - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
    hCylindricalAveDiff[2] = new TH2F("hAveDiffXYPhiCart", "Phi Model from Cartesian - Truth Averaged Over z (#mum); x (cm); y (cm)",nbins,low,high,nbins,low,high);
    hCylindricalAveDiff[3] = new TH2F("hAveDiffRZPhiCart", "Phi Model from Cartesian - Truth Averaged Over phi (#mum); z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);
  
    TH2F *hSamplePerBinXY = new TH2F("hSamplePerBinXY", "Filling each xy bin; x (cm); y (cm)",nbins,low,high,nbins,low,high);
    TH2F *hSamplePerBinRZ = new TH2F("hSamplePerBinRZ", "Filling each rz bin; z (cm); r (cm)", nz,minz,maxz,nr,minr,maxr);

    TH2F *hCompareRTrue = new TH2F("hCompareRTrue", "Compare Difference from R Model and True (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);
    TH2F *hComparePhiTrue = new TH2F("hComparePhiTrue", "Compare Difference from Phi Model and True (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);

    TH2F *hRDiffvR = new TH2F("hRDiffvR", "Difference between R Model and True vs. r (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvZ = new TH2F("hRDiffvZ", "Difference between R Model and True vs. z (R > 30); z (cm); shift difference (#mum)",nz,minz,maxz,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvPhi = new TH2F("hRDiffvPhi", "Difference between R Model and True vs. phi (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    TH2F *hPhiDiffvR = new TH2F("hPhiDiffvR", "Difference between Phi Model and True vs. r (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvZ = new TH2F("hPhiDiffvZ", "Difference between Phi Model and True vs. z (R > 30); z (cm); shift difference (#mum)",nz,minz,maxz,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvPhi = new TH2F("hPhiDiffvPhi", "Difference between Phi Model and True vs. phi (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    for(int i = 1; i < nphi - 1; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
      for(int j = 1; j < nr - 1; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
	for(int k = 1; k < nz - 1; k++){
	  double z = minz + ((maxz - minz)/(1.0*nz))*(k+0.5); //center of bin

	  double shiftrecoCart[3];
	  double shifttrueCart[3];
	  double differenceCart[3];
	
	  double shiftrecoCyl[2];
	  double shifttrueCyl[2];
	  double differenceCyl[2];

	  double differenceR, differencePhi;

	  int bin = hCartCMModel[0]->FindBin(phi,r,z); //same for all

	  if((r > 30.0) && (r < 76.0)){
	    shifttrueCart[0] = (shifter->hX->Interpolate(phi,r,z))*(1e4); //convert from cm to micron
	    shifttrueCart[1] = (shifter->hY->Interpolate(phi,r,z))*(1e4); //convert from cm to micron 
	    shifttrueCart[2] = (shifter->hZ->Interpolate(phi,r,z))*(1e4); //convert from cm to micron 

	    //x y and z
	    for(int l = 0; l < 3; l ++){
	      shiftrecoCart[l] =  (hCartCMModel[l]->GetBinContent(bin))*(1e4);
	  
	      differenceCart[l] = shiftrecoCart[l] - shifttrueCart[l]; 

	      hCartesianShiftDifference[l]->Fill(differenceCart[l]);
	    }

	    //r from cart
	    shiftrecoCyl[0] =  (hCylCMModel[0]->GetBinContent(bin))*(1e4);
	    shifttrueCyl[0] = (shifter->hR->Interpolate(phi,r,z))*(1e4); //convert from cm to micron 
	    differenceCyl[0] = shiftrecoCyl[0] - shifttrueCyl[0]; 
	    hCylindricalShiftDifference[0]->Fill(differenceCyl[0]);
	      
	    //phi from cart
	    shiftrecoCyl[1] = r*(1e4)*(hCylCMModel[1]->GetBinContent(bin));
	    shifttrueCyl[1] = (shifter->hPhi->Interpolate(phi,r,z))*(1e4); 
	    differenceCyl[1] = (shiftrecoCyl[1] - shifttrueCyl[1]); 
	    hCylindricalShiftDifference[1]->Fill(differenceCyl[1]);

	    hRShiftTrue->Fill(shifttrueCyl[0]);
	    hPhiShiftTrue->Fill(shifttrueCyl[1]);
	
	    double x = r*cos(phi);
	    double y = r*sin(phi);
       	
	    //x
	    hCartesianDiff[0]->Fill(x,y, differenceCart[0]);
	    hCartesianDiff[1]->Fill(z,r, differenceCart[0]);
	    //y
	    hCartesianDiff[2]->Fill(x,y, differenceCart[1]);	  
	    hCartesianDiff[3]->Fill(z,r, differenceCart[1]);
	    //z
	    hCartesianDiff[4]->Fill(x,y, differenceCart[2]);
	    hCartesianDiff[5]->Fill(z,r, differenceCart[2]);

	    //r cart
	    hCylindricalDiff[0]->Fill(x,y, differenceCyl[0]);
	    hCylindricalDiff[1]->Fill(z,r, differenceCyl[0]);
	    //phi cart
	    hCylindricalDiff[2]->Fill(x,y, differenceCyl[1]);
	    hCylindricalDiff[3]->Fill(z,r, differenceCyl[1]);
	    
	    hCompareRTrue->Fill(shiftrecoCyl[0],shifttrueCyl[0]);
	    hComparePhiTrue->Fill(shiftrecoCyl[1],shifttrueCyl[1]);

	    hRDiffvR->Fill(r,differenceCyl[0],1);
	    hRDiffvPhi->Fill(phi,differenceCyl[0],1);
	    hRDiffvZ->Fill(z,differenceCyl[0],1);
	    
	    hPhiDiffvR->Fill(r,differenceCyl[1],1);
	    hPhiDiffvPhi->Fill(phi,differenceCyl[1],1);
	    hPhiDiffvZ->Fill(z,differenceCyl[1],1);
	    
	    hSamplePerBinXY->Fill(x,y,1);

	    hSamplePerBinRZ->Fill(z,r,1);
	    
	  }
	}
      }
    }

    //average over z
    for (int m = 0; m < 6; m = m+2){
      hCartesianAveDiff[m]->Divide(hCartesianDiff[m],hSamplePerBinXY);
    }
    for (int m = 0; m < 4; m = m+2){
      hCylindricalAveDiff[m]->Divide(hCylindricalDiff[m],hSamplePerBinXY);
    }
    
    //average over phi
    for (int m = 1; m < 6; m = m+2){
      hCartesianAveDiff[m]->Divide(hCartesianDiff[m],hSamplePerBinRZ);
    }
    for (int m = 1; m < 4; m = m+2){
      hCylindricalAveDiff[m]->Divide(hCylindricalDiff[m],hSamplePerBinRZ);
    }

    //summary plots
    hDifferenceMeanR->Fill(hCylindricalShiftDifference[2]->GetMean(1));
    hDifferenceStdDevR->Fill(hCylindricalShiftDifference[2]->GetStdDev(1));

    hTrueMeanR->Fill(hRShiftTrue->GetMean(1));
    hTrueStdDevR->Fill(hRShiftTrue->GetStdDev(1));
    
    hDifferenceMeanPhi->Fill(hCylindricalShiftDifference[3]->GetMean(1));
    hDifferenceStdDevPhi->Fill(hCylindricalShiftDifference[3]->GetStdDev(1));

    hTrueMeanPhi->Fill(hPhiShiftTrue->GetMean(1));
    hTrueStdDevPhi->Fill(hPhiShiftTrue->GetStdDev(1));

    for (int m = 0; m < 6; m++){
      hCartesianAveDiff[m]->SetStats(0);
    }
    for (int m = 0; m < 4; m++){
      hCylindricalAveDiff[m]->SetStats(0);
    }
  
    hCompareRTrue->SetStats(0);
    hComparePhiTrue->SetStats(0);

    hRDiffvR->SetStats(0);
    hRDiffvZ->SetStats(0);
    hRDiffvPhi->SetStats(0);
  
    hPhiDiffvR->SetStats(0);
    hPhiDiffvZ->SetStats(0);
    hPhiDiffvPhi->SetStats(0);
    
    TPad *c1=new TPad("c1","",0.0,0.8,1.0,0.93); //can i do an array of pads?
    TPad *c2=new TPad("c2","",0.0,0.64,1.0,0.77);
    TPad *c3=new TPad("c3","",0.0,0.48,1.0,0.61);
    TPad *c4=new TPad("c4","",0.0,0.32,1.0,0.45);
    TPad *c5=new TPad("c5","",0.0,0.16,1.0,0.29);
    TPad *c6=new TPad("c6","",0.0,0.0,1.0,0.13);
    
    TPad *titlepad=new TPad("titlepad","",0.0,0.96,1.0,1.0);

    TPad *stitlepad1=new TPad("stitlepad1","",0.0,0.93,1.0,0.96);
    TPad *stitlepad2=new TPad("stitlepad2","",0.0,0.77,1.0,0.8);
    TPad *stitlepad3=new TPad("stitlepad3","",0.0,0.61,1.0,0.64);
    TPad *stitlepad4=new TPad("stitlepad4","",0.0,0.45,1.0,0.48);
    TPad *stitlepad5=new TPad("stitlepad5","",0.0,0.29,1.0,0.32);
    TPad *stitlepad6=new TPad("stitlepad6","",0.0,0.13,1.0,0.16);
    
    TLatex * title = new TLatex(0.0,0.0,"");

    TLatex * stitle1 = new TLatex(0.0,0.0,""); //array?
    TLatex * stitle2 = new TLatex(0.0,0.0,"");
    TLatex * stitle3 = new TLatex(0.0,0.0,"");
    TLatex * stitle4 = new TLatex(0.0,0.0,"");
    TLatex * stitle5 = new TLatex(0.0,0.0,"");
    TLatex * stitle6 = new TLatex(0.0,0.0,"");
    
    title->SetNDC();
    stitle1->SetNDC();
    stitle2->SetNDC();
    stitle3->SetNDC();
    stitle4->SetNDC();
    stitle5->SetNDC();
    stitle6->SetNDC();
    
    title->SetTextSize(0.32);
    stitle1->SetTextSize(0.35);
    stitle2->SetTextSize(0.35);
    stitle3->SetTextSize(0.35);
    stitle4->SetTextSize(0.35);
    stitle5->SetTextSize(0.35);
    stitle6->SetTextSize(0.35);
    
    canvas->cd();
    c1->Draw();
    stitlepad1->Draw();
    c2->Draw();
    stitlepad2->Draw();
    c3->Draw();
    stitlepad3->Draw();
    c4->Draw();
    stitlepad4->Draw();
    c5->Draw();
    stitlepad5->Draw();
    c6->Draw();
    stitlepad6->Draw();
    titlepad->Draw();

    //x plots
    c1->Divide(4,1);
    c1->cd(1);
    hCartesianAveDiff[0]->Draw("colz");
    c1->cd(2);
    hCartesianAveDiff[1]->Draw("colz");
    c1->cd(3);
    hCartesianShiftDifference[0]->Draw();
    //c1->cd(4)->Clear();  
    c1->cd(4);
    //hCMmodelSliceRvTrue->Draw("colz");
    hSamplePerBinRZ->Draw("colz");
    
    //y plots
    c2->Divide(4,1);
    c2->cd(1);
    hCartesianAveDiff[2]->Draw("colz");
    c2->cd(2);
    hCartesianAveDiff[3]->Draw("colz");
    c2->cd(3);
    hCartesianShiftDifference[1]->Draw();
    //c2->cd(4)->Clear();
    c2->cd(4);
    //hStripesPerBin->Draw("colz");
    hSamplePerBinXY->Draw("colz");
    
    //r cart
    c3->Divide(4,1);
    c3->cd(1);
    hCylindricalAveDiff[0]->Draw("colz");
    c3->cd(2);
    hCylindricalAveDiff[1]->Draw("colz");
    c3->cd(3);
    hCylindricalShiftDifference[0]->Draw();
    c3->cd(4);
    hRShiftTrue->Draw();
    
    //phi cart
    c4->Divide(4,1);
    c4->cd(1);
    hCylindricalAveDiff[2]->Draw("colz");
    c4->cd(2);
    hCylindricalAveDiff[3]->Draw("colz");
    c4->cd(3);
    hCylindricalShiftDifference[1]->Draw();
    c4->cd(4);
    hPhiShiftTrue->Draw();

    //r to true comparison
    c5->Divide(4,1);
    c5->cd(1);
    hCompareRTrue->Draw("colz");
    c5->cd(2);
    hRDiffvR->Draw("colz");
    c5->cd(3);
    hRDiffvZ->Draw("colz");
    c5->cd(4);
    hRDiffvPhi->Draw("colz");

    //phi to true comparison
    c6->Divide(4,1);
    c6->cd(1);
    hComparePhiTrue->Draw("colz");
    c6->cd(2);
    hPhiDiffvR->Draw("colz");
    c6->cd(3);
    hPhiDiffvZ->Draw("colz");
    c6->cd(4);
    hPhiDiffvPhi->Draw("colz");

    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.01,0.4,Form("Event %d; %s", ifile, sourcefilename.Data())); 
    title->Draw();
    
    stitlepad1->cd();
    stitlepad1->Clear();
    stitle1->DrawLatex(0.45,0.2,"X Model"); 
    stitle1->Draw();
     
    stitlepad2->cd();
    stitlepad2->Clear();
    stitle2->DrawLatex(0.45,0.2,"Y Model"); 
    stitle2->Draw();

    stitlepad3->cd();
    stitlepad3->Clear();
    stitle3->DrawLatex(0.45,0.2,"R Model"); 
    stitle3->Draw();

    stitlepad4->cd();
    stitlepad4->Clear();
    stitle4->DrawLatex(0.45,0.2,"Phi Model"); 
    stitle4->Draw();

    stitlepad5->cd();
    stitlepad5->Clear();
    stitle5->DrawLatex(0.4,0.2,"Comparing R Model to True"); 
    stitle5->Draw();

    stitlepad6->cd();
    stitlepad6->Clear();
    stitle6->DrawLatex(0.4,0.2,"Comparing Phi Model to True"); 
    stitle6->Draw();

    if(ifile == 0){ 
      //if(ifile == 1){
      canvas->Print("CMDistortionAnalysis.pdf(","pdf");
    } else{
      canvas->Print("CMDistortionAnalysis.pdf","pdf");
    }
  }

  TCanvas *summary = new TCanvas("summary","ShiftPlotsSummary",2000,3000);

  TPad *sumtitlepad = new TPad("sumtitlepad","",0.0,0.96,1.0,1.0);
  TPad *sumplots = new TPad("sumplotspad","",0.0,0.0,1.0,0.96);

  TLatex *sumtitle = new TLatex(0.0,0.0,"");

  sumtitle->SetNDC();
  sumtitle->SetTextSize(0.4);

  summary->cd();
  sumplots->Draw();
  sumtitlepad->Draw();

  sumplots->Divide(4,6);
  sumplots->cd(1);
  hDifferenceMeanR->Draw();
  sumplots->cd(2);
  hDifferenceStdDevR->Draw();
  sumplots->cd(3);
  hTrueMeanR->Draw();
  sumplots->cd(4);
  hTrueStdDevR->Draw();
  sumplots->cd(5);
  hDifferenceMeanPhi->Draw();
  sumplots->cd(6);
  hDifferenceStdDevPhi->Draw();
  sumplots->cd(7);
  hTrueMeanPhi->Draw();
  sumplots->cd(8);
  hTrueStdDevPhi->Draw();
  sumplots->cd(9);
  sumplots->cd(10)->Clear();
  sumplots->cd(11)->Clear();
  sumplots->cd(12)->Clear();
  sumplots->cd(13)->Clear();
  sumplots->cd(14)->Clear();
  sumplots->cd(15)->Clear();
  sumplots->cd(16)->Clear();
  sumplots->cd(17)->Clear();
  sumplots->cd(18)->Clear();
  sumplots->cd(19)->Clear();
  sumplots->cd(20)->Clear();
  sumplots->cd(21)->Clear();
  sumplots->cd(22)->Clear();
  sumplots->cd(23)->Clear();
  sumplots->cd(24)->Clear();

  sumtitlepad->cd();
  sumtitlepad->Clear();
  sumtitle->DrawLatex(0.4,0.4,"Summary of Events"); 
  summary->Print("CMDistortionAnalysis.pdf)","pdf");


  
    //canvas->Print(Form("ShiftPlotsEvent%03i\n.gif", ifile),"gif"); //left padding for up to 3 digits
    
    /*   TCanvas *canvas=new TCanvas("canvas","ShiftPlots",1500,1000);
    TPad *c=new TPad("c","",0.0,0.0,1.0,0.9);
    TPad *titlepad=new TPad("titlepad","",0.0,0.9,1.0,1.0);
    TLatex * title = new TLatex(0.0,0.0,"");
    title->SetNDC();
    title->SetTextSize(0.5);
    canvas->cd();
    c->Draw();
    titlepad->Draw();
  
    // x plots
    c->Divide(3,2);
    c->cd(1);
    hCartesianForward[0]->Draw("colz");
    c->cd(2);
    hStripesPerBin->Draw("colz");
    c->cd(3);
    hCartesianAveShift[0]->Draw("colz");
    c->cd(4);
    hCartesianAveDiff[0]->Draw("colz");
    c->cd(5);
    hCartesianAveDiff[1]->Draw("colz");
    c->cd(6);
    hCartesianShiftDifference[0]->Draw();
    //c->cd();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"X Shift Model");
    title->Draw();
    canvas->Print("ShiftPlots.pdf(","pdf");
  
    // y plots
    //c->Divide(3,2);
    c->cd(1);
    hCartesianForward[1]->Draw("colz");
    c->cd(2);
    hStripesPerBin->Draw("colz");
    c->cd(3);
    hCartesianAveShift[1]->Draw("colz");
    c->cd(4);
    hCartesianAveDiff[2]->Draw("colz");
    c->cd(5);
    hCartesianAveDiff[3]->Draw("colz");
    c->cd(6);
    hCartesianShiftDifference[1]->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Y Shift Model");
    canvas->Print("ShiftPlots.pdf","pdf");
  
    // z plots
    // c->Divide(3,2);
    c->cd(1);
    hCartesianForward[2]->Draw("colz");
    c->cd(2);
    hStripesPerBin->Draw("colz");
    c->cd(3);
    hCartesianAveShift[2]->Draw("colz");
    c->cd(4);
    hCartesianAveDiff[4]->Draw("colz");
    c->cd(5);
    hCartesianAveDiff[5]->Draw("colz");
    c->cd(6);
    hCartesianShiftDifference[2]->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Z Shift Model");
    canvas->Print("ShiftPlots.pdf","pdf");
  
    // r plots
    //c->Divide(3,2);
    c->cd(1);
    hCylindricalForward[0]->Draw("colz");
    c->cd(2);
    hStripesPerBin->Draw("colz");
    c->cd(3);
    hCylindricalAveShift[0]->Draw("colz");
    c->cd(4);
    //PhiCheck->Draw();
    hCylindricalAveDiff[0]->Draw("colz");
    c->cd(5);
    // hPhiCheck2d->Draw("colz");
    hCylindricalAveDiff[1]->Draw("colz");
    c->cd(6);
    hCylindricalShiftDifference[0]->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"R Shift Model");
    canvas->Print("ShiftPlots.pdf","pdf");

    // r plots from cart
    //c->Divide(3,2);
    c->cd(1)->Clear();
    //hCylindricalForward[2]->Draw("colz");
    c->cd(2)->Clear();
    //hStripesPerBin->Draw("colz");
    c->cd(3);
    hCylindricalAveShift[2]->Draw("colz");
    c->cd(4);
    hCylindricalAveDiff[4]->Draw("colz");
    c->cd(5);
    hCylindricalAveDiff[5]->Draw("colz");
    c->cd(6);
    hCylindricalShiftDifference[2]->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.3,0.2,"R Shift from Cartesian Model");
    canvas->Print("ShiftPlots.pdf","pdf");

    // compare the two R models
    c->cd(1);
    hCylindricalAveShift[0]->Draw("colz");
    c->cd(2);
    hCylindricalAveShift[2]->Draw("colz");
    c->cd(3)->Clear();
    c->cd(4);
    hRAveDiff[0]->Draw("colz");
    c->cd(5);
    hRAveDiff[1]->Draw("colz");
    c->cd(6);
    hRShiftDifference->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Comparing R Models");
    canvas->Print("ShiftPlots.pdf","pdf");

    //compare R cart to true
    c->cd(1);
    hCylindricalShiftDifference[2]->Draw();
    c->cd(2);
    hRShiftTrue->Draw();
    c->cd(3);
    hCompareRTrue->Draw("colz");
    c->cd(4);
    hRDiffvR->Draw("colz");
    c->cd(5);
    hRDiffvZ->Draw("colz");
    c->cd(6);
    hRDiffvPhi->Draw("colz");
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.3,0.2,"Comparing R from Cart Model to True");
    canvas->Print("ShiftPlots.pdf","pdf");
  
    // phi plots
    //c->Divide(3,2);
    c->cd(1);
    hCylindricalForward[1]->Draw("colz");
    c->cd(2);
    hStripesPerBin->Draw("colz");
    c->cd(3);
    hCylindricalAveShift[1]->Draw("colz");
    c->cd(4);
    hCylindricalAveDiff[2]->Draw("colz");
    c->cd(5);
    hCylindricalAveDiff[3]->Draw("colz");
    c->cd(6);
    hCylindricalShiftDifference[1]->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Phi Shift Model");
    canvas->Print("ShiftPlots.pdf","pdf");

    // phi plots from cart
    c->cd(1)->Clear();
    c->cd(2)->Clear();
    c->cd(3);
    hCylindricalAveShift[3]->Draw("colz");
    c->cd(4);
    hCylindricalAveDiff[6]->Draw("colz");
    c->cd(5);
    hCylindricalAveDiff[7]->Draw("colz");
    c->cd(6);
    hCylindricalShiftDifference[3]->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.3,0.2,"Phi Shift from Cartesian Model");
    canvas->Print("ShiftPlots.pdf","pdf");

    // compare the two Phi models
    c->cd(1);
    hCylindricalAveShift[1]->Draw("colz");
    c->cd(2);
    hCylindricalAveShift[3]->Draw("colz");
    c->cd(3)->Clear();
    c->cd(4);
    hPhiAveDiff[0]->Draw("colz");
    c->cd(5);
    hPhiAveDiff[1]->Draw("colz");
    c->cd(6);
    hPhiShiftDifference->Draw();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Comparing Phi Models");
    canvas->Print("ShiftPlots.pdf","pdf");
  
    //compare Phi cart to true
    c->cd(1);
    hCylindricalShiftDifference[3]->Draw();
    c->cd(2);
    hPhiShiftTrue->Draw();
    c->cd(3);
    hComparePhiTrue->Draw("colz");
    c->cd(4);
    hPhiDiffvR->Draw("colz");
    c->cd(5);
    hPhiDiffvZ->Draw("colz");
    c->cd(6);
    hPhiDiffvPhi->Draw("colz");
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.3,0.2,"Comparing Phi from Cart Model to True");
    canvas->Print("ShiftPlots.pdf","pdf");
  
    //Compare X and Y diff
    c->cd(1);
    hCompareXY->Draw("colz");
    c->cd(2)->Clear();
    c->cd(3)->Clear();
    c->cd(4)->Clear();
    c->cd(5)->Clear();
    c->cd(6)->Clear();
    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.4,0.2,"Comparing X and Y Models");
    canvas->Print("ShiftPlots.pdf)","pdf"); */
  
    // c->SaveAs("RShift.pdf"); // replace w print
  

  
  return 0;
}

void ScanHist(int nbins, double low, double high, double x, double y){
  PHG4TpcCentralMembrane stripes;
  int stripeID; 
  //histogram from search
  TH2F *Pattern1 = new TH2F("Pattern1","X,Y Scan if in Stripe;X (mm);Y (mm)",nbins,low,high,nbins,low,high); // min n max just beyond extent of CM so it's easier to see
  
  //TLatex *tex=new TLatex(x,y,"Stripe");
  //tex->SetTextSize(0.005);
  // for (r = stripes.begin_CM; r < stripes.end_CM; r = r + rstepsize){ // radii spanning full CM
  // for (phi = 0.0; phi < 2*TMath::Pi(); phi = phi + phistepsize){ // angles spanning full CM
  for (int i=0;i<nbins;i++){
    x=low+(high-low)/nbins*(i+0.5);
    for (int j=0;j<nbins;j++){
      y=low+(high-low)/nbins*(j+0.5);
 
    // x = r*cos(phi);
    // y = r*sin(phi);
      // cout << x << endl;
      //cout << y << endl; 
      
      stripeID = stripes.getStripeID(x, y);
      /* TLatex *tex=new TLatex(x,y,"StripeID");
	 tex->SetTextSize(0.005);
	 tex->DrawLatex(x,y,Form("%d",stripeID));
      */
      
      if(stripeID == -1){
        Pattern1->Fill(x,y,0);
      } else{
        Pattern1->Fill(x,y,1);
      }
      
    }
  }

  TCanvas *c=new TCanvas("a","CheckStripeID.cpp",500,500); 
  Pattern1->Draw("colz"); //check if theres still white spaces
  c->SaveAs("cmScan.pdf");
  //
}

void IDLabels(){
  PHG4TpcCentralMembrane stripes;
  int stripeID;
  vector<PHG4Hitv1*> Hits = stripes.PHG4Hits;
  const double mm = 1.0;
  const double cm = 10.0;
  
  vector<double> xhit;
  vector<double> yhit;
  
  //build tgraph from dummy hits
  for (int i = 0; i < Hits.size(); i++){
    xhit.push_back(Hits[i]->get_x(0)*cm/mm); 
    yhit.push_back(Hits[i]->get_y(0)*cm/mm);
    xhit.push_back(Hits[i]->get_x(1)*cm/mm);
    yhit.push_back(Hits[i]->get_y(1)*cm/mm);
  }
  
  int npts = 2*Hits.size();
  TGraph *gDummyHits = new TGraph(npts, &xhit[0], &yhit[0]);
  gDummyHits->SetMarkerColor(2);
  gDummyHits->SetMarkerSize(0.5);
  
  gStyle->SetOptStat(0);
  TCanvas *c=new TCanvas("a","CheckStripeID.cpp",500,500);
  gDummyHits->Draw("AP");

  //loop thru hits again
  double xav, yav, xa, xb, ya, yb;

  for (int i = 0; i < Hits.size(); i++){
    //avg x0 n x1, y0 n y1 and use to draw stripeID
    xav = (Hits[i]->get_x(0)*cm/mm + Hits[i]->get_x(1)*cm/mm)/2;
    yav = (Hits[i]->get_y(0)*cm/mm + Hits[i]->get_y(1)*cm/mm)/2;
    //cout << "i: " << i << endl;
    xa = Hits[i]->get_x(0)*cm/mm ;
    xb= Hits[i]->get_x(1)*cm/mm ;
    ya= Hits[i]->get_y(0)*cm/mm ;
    yb = Hits[i]->get_y(1)*cm/mm ;
    //cout << "xav: " << xav << endl;
    //cout << "yav: " << yav << endl; 
    stripeID = stripes.getStripeID(xav, yav);
    TLatex *tex=new TLatex(xav,yav,"StripeID");
    tex->SetTextSize(0.005);
    tex->DrawLatex(xav,yav,Form("%d",stripeID));
  }

  c->SaveAs("cmStripeID.pdf");
}  
  //loop thru hits again
  /*
  double xav, yav, xa, xb, ya, yb;
 for (int i = 0; i < Hits.size(); i++){
  
    //avg x0 n x1, y0 n y1 and use to draw stripeID
    xav = (Hits[i]->get_x(0)*cm/mm + Hits[i]->get_x(1)*cm/mm)/2;
    yav = (Hits[i]->get_y(0)*cm/mm + Hits[i]->get_y(1)*cm/mm)/2;
    //cout << "i: " << i << endl;
    xa = Hits[i]->get_x(0)*cm/mm ;
    xb= Hits[i]->get_x(1)*cm/mm ;
    ya= Hits[i]->get_y(0)*cm/mm ;
    yb = Hits[i]->get_y(1)*cm/mm ;
    //cout << "xav: " << xav << endl;
    //cout << "yav: " << yav << endl; 
    stripeID = stripes.getStripeID(xav, yav);
    TLatex *tex=new TLatex(xav,yav,"StripeID");
    tex->SetTextSize(0.005);
    tex->DrawLatex(xav,yav,Form("%d",stripeID));
    

    TLatex *tex=new TLatex(xav,yav,"Stripe");
    tex->SetTextSize(0.005);
    if(stripeID == -1){
      tex->DrawLatex(xav,yav,Form("%d",0));
    } else{
      tex->DrawLatex(xav,yav,Form("%d",1));
    }
}
  */
      
    //TLine *line=new TLine;
    //line->DrawLine(Hits[i]->get_x(0)*cm/mm ,Hits[i]->get_y(0)*cm/mm,Hits[i]->get_x(1)*cm/mm, Hits[i]->get_y(1)*cm/mm);
    //line->DrawLine(0, 0, 200, 200*((yb- ya)/(xb-xa)));
    
  
  
  /*for (r = stripes.begin_CM; r < stripes.end_CM; r = r + rstepsize){ // radii spanning full CM
    for (phi = 0.0; phi < 2.0*TMath::Pi()/9.0; phi = phi + phistepsize){ // angles spanning full CM
      
      x = r*cos(phi);
      y = r*sin(phi);
      
      
      
      //TLatex tex;
      
      //if(stripeID == 1)
      //Pattern1->Fill(x,y);
      cout << stripeID << endl;
    }
  }
  */
  
 
