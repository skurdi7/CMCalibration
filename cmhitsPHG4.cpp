#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "PHG4TpcCentralMembrane.h"
#include "TTree.h"
R_LOAD_LIBRARY(.libs/libg4tpccentralmembrane)

//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
R__LOAD_LIBRARY(libphg4hit.so)


// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps

using namespace std;


int cmhitsPHG4() {
  StripesClass stripes;

  int result, nbins, rsteps, phisteps; 
  double r, phi, x, y, xmod, ymod, phimod, rstepsize, phistepsize;
  
  nbins = 100;
  rsteps = 100;
  phisteps = 100;
  
  rstepsize = (stripes.end_CM - stripes.begin_CM)/rsteps;
  phistepsize = 2*TMath::Pi()/phisteps;
  
  //histogram from search
  TH2F *Pattern1 = new TH2F("Pattern1","Pattern1",nbins,-770.0,770.0,nbins,-770.0,770.0); // min n max just beyond extent of CM so it's easier to see
  
  for (r = stripes.begin_CM; r < stripes.end_CM; r = r + rstepsize){ // radii spanning full CM
    for (phi = 0.0; phi < 2*TMath::Pi(); phi = phi + phistepsize){ // angles spanning full CM
      
      x = r*cos(phi);
      y = r*sin(phi);

      result = stripes.getSearchResult(x, y);

      if(result == 1)
	Pattern1->Fill(x,y);
    }
  }

  vector<PHG4Hitv1*> Hits = stripes.PHG4Hits;

  const double mm = 1.0;
  const double cm = 10.0;
  
  vector<double> xhitS;
  vector<double> yhitS;
  
  //build tgraph from dummy hits
  for (int i = 0; i < Hits.size(); i++){
    xhitS.push_back(Hits[i]->get_x(0)*cm/mm); 
    yhitS.push_back(Hits[i]->get_y(0)*cm/mm);
    xhitS.push_back(Hits[i]->get_x(1)*cm/mm);
    yhitS.push_back(Hits[i]->get_y(1)*cm/mm);
    
  }
  /*
  int npts = 2*Hits.size();
  TGraph *gDummyHits = new TGraph(npts, &xhit[0], &yhit[0]);
  gDummyHits->SetMarkerColor(2);

  gStyle->SetOptStat(0);
  TCanvas *c=new TCanvas("a","cmhitsPHG4.cpp",500,500);
  Pattern1->Draw();
  gDummyHits->Draw("P");
  c->SaveAs("cmhitsPHG4.pdf");
  */



  TTree *sTree=new TTree("tree","phg4hits");
  double xhitfortree, yhitfortree;
  sTree->Branch("xhit",&xhitfortree);
  sTree->Branch("yhit",&yhitfortree);

  for (int i=0;i<xhitS.size();i++){
    xhitfortree=xhitS[i];
    yhitfortree=yhitS[i];
   
    sTree->Fill();
  }
  
  sTree->SaveAs("phg4hitsTree.root");

  vector<double> xhit;
  vector<double> yhit;

  //this is how to get the hits
  char const *treename="sTree";
  TFile *input=TFile::Open("phg4hitsTree.root");
  TTree *inTree=(TTree*)input->Get("tree");
  inTree->SetBranchAddress("xhit",&xhitfortree);
  inTree->SetBranchAddress("yhit",&yhitfortree);
  for (int i=0;i<inTree->GetEntries();i++){
    inTree->GetEntry(i);
    xhit.push_back(xhitfortree);
    yhit.push_back(yhitfortree);   
  }
  input->Close();

  int npts = 2*Hits.size();
  TGraph *gDummyHits = new TGraph(npts, &xhit[0], &yhit[0]);
  gDummyHits->SetMarkerColor(2);
  
  gStyle->SetOptStat(0);
  TCanvas *c=new TCanvas("a","cmhitsPHG4.cpp",500,500);
  Pattern1->Draw();
  gDummyHits->Draw("P");
  c->SaveAs("cmhitsPHG4.pdf");
    
  return 0;
}
