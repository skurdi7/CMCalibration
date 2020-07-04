#include <iostream>
#include <cmath>
#include "TVector3.h"

using namespace std;

// generates histogram of stripes and dummy hit coordinates 
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps
// all distances in mm, all angles in rad

const int nRadii = 8; 
const double str_width = 1.0; // width of a stripe
const double arc_r = 0.5; // radius of arc on end of a stripe

struct Container{
  double x0,y0,z0;
  double x1,y1,z1;
  float px0,py0,pz0;
  float px1,py1,pz1;
  float t0,t1;
  int trkid;
  float edep0,eion0;
  float edep1,eion1;
};

void Vertices(int nStripes, double R[], double spacing[], double x1a[][nRadii], double y1a[][nRadii], double x1b[][nRadii], double y1b[][nRadii], double x2a[][nRadii], double y2a[][nRadii], double x2b[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double padfrac, int nGoodStripes[], int keepThisAndAfter[], int keepUntil[]);

int search(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double x, double y, int nGoodStripes[], int keepThisAndAfter[], int keepUntil[]);

void GetG4Hit(Container hit[][nRadii], int nStripes, double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], int nGoodStripes[], int keepThisAndAfter[], int keepUntil[], int nElectrons, float time, TVector3 pos0[][nRadii], TVector3 pos1[][nRadii]);

int cmhits() {
  
  // Radii of centers (taken from Inventor measurements)
  //R1 extension:
  double R1_e[nRadii] =
    {227.0902789, 238.4100043, 249.7297296, 261.049455, 272.3691804, 283.6889058, 295.0086312, 306.3283566};
  //R1:
  double R1[nRadii] =
    {317.648082, 328.9678074, 340.2875328, 351.6072582, 362.9269836, 374.246709,  385.5664344, 396.8861597};
  //R2:
  double R2[nRadii] =
    {421.705532, 442.119258, 462.532984, 482.9467608, 503.36069, 523.774416, 544.188015, 564.601868};
  //R3:
  double R3[nRadii] =
    {594.6048725, 616.545823, 638.4867738, 660.4277246, 682.3686754, 704.3096262, 726.250577, 748.1915277};
  

  // Angle calculation
  const double pi = 3.14159265358979323846;
  const double phi_module = pi/6.0; // angle span of a module
  const int nPads_R1 = 6*16; //pads per row (256 pads per card in module, 16 rows in module), same for R1_e
  const int nPads_R2 = 8*16;
  const int nPads_R3 = 12*16;
  const int pr_mult = 3; // multiples of intrinsic resolution of pads
  const int dw_mult = 8; // multiples of diffusion width
  const double diffwidth = 0.6; // diffusion width
  double spacing_R1_e[nRadii], spacing_R1[nRadii], spacing_R2[nRadii], spacing_R3[nRadii];

  for (int i=0; i<nRadii; i++){
    spacing_R1_e[i] = 2.0*((dw_mult*diffwidth/R1_e[i]) + (pr_mult*phi_module/nPads_R1));
    spacing_R1[i] = 2.0*((dw_mult*diffwidth/R1[i]) + (pr_mult*phi_module/nPads_R1));
    spacing_R2[i] = 2.0*((dw_mult*diffwidth/R2[i]) + (pr_mult*phi_module/nPads_R2));
    spacing_R3[i] = 2.0*((dw_mult*diffwidth/R3[i]) + (pr_mult*phi_module/nPads_R3));
  }

  
  // vertex coordinates calculation
  const double padfrac_R1 = 0.5*5.59106385; // 1/2 of R1 pad length (sectors 2-15), same for R1_e
  const double padfrac_R2 = 0.5*10.13836283; // 1/2 of R2 pad length (sectors 2-15)
  const double padfrac_R3 = 0.5*10.90189537; // 1/2 of R3 pad length (sectors 2-15)
  const int nStripes_R1 = 6; // Same for R1_e
  const int nStripes_R2 = 8;
  const int nStripes_R3 = 12;
  
  //bottom left - 1a
  double x1a_R1_e[nStripes_R1][nRadii], y1a_R1_e[nStripes_R1][nRadii];
  double x1a_R1[nStripes_R1][nRadii], y1a_R1[nStripes_R1][nRadii];
  double x1a_R2[nStripes_R2][nRadii], y1a_R2[nStripes_R2][nRadii];
  double x1a_R3[nStripes_R3][nRadii], y1a_R3[nStripes_R3][nRadii];

  //bottom right - 1b
  double x1b_R1_e[nStripes_R1][nRadii], y1b_R1_e[nStripes_R1][nRadii];
  double x1b_R1[nStripes_R1][nRadii], y1b_R1[nStripes_R1][nRadii];
  double x1b_R2[nStripes_R2][nRadii], y1b_R2[nStripes_R2][nRadii];
  double x1b_R3[nStripes_R3][nRadii], y1b_R3[nStripes_R3][nRadii];

  //top left - 2a
  double x2a_R1_e[nStripes_R1][nRadii], y2a_R1_e[nStripes_R1][nRadii];
  double x2a_R1[nStripes_R1][nRadii], y2a_R1[nStripes_R1][nRadii];
  double x2a_R2[nStripes_R2][nRadii], y2a_R2[nStripes_R2][nRadii];
  double x2a_R3[nStripes_R3][nRadii], y2a_R3[nStripes_R3][nRadii];

  //top right - 2b
  double x2b_R1_e[nStripes_R1][nRadii], y2b_R1_e[nStripes_R1][nRadii];
  double x2b_R1[nStripes_R1][nRadii], y2b_R1[nStripes_R1][nRadii];
  double x2b_R2[nStripes_R2][nRadii], y2b_R2[nStripes_R2][nRadii];
  double x2b_R3[nStripes_R3][nRadii], y2b_R3[nStripes_R3][nRadii];

  //left midpoint - 3a
  double x3a_R1_e[nStripes_R1][nRadii], y3a_R1_e[nStripes_R1][nRadii];
  double x3a_R1[nStripes_R1][nRadii], y3a_R1[nStripes_R1][nRadii];
  double x3a_R2[nStripes_R2][nRadii], y3a_R2[nStripes_R2][nRadii];
  double x3a_R3[nStripes_R3][nRadii], y3a_R3[nStripes_R3][nRadii];
  
  //right midpoint - 3b
  double x3b_R1_e[nStripes_R1][nRadii], y3b_R1_e[nStripes_R1][nRadii];
  double x3b_R1[nStripes_R1][nRadii], y3b_R1[nStripes_R1][nRadii];
  double x3b_R2[nStripes_R2][nRadii], y3b_R2[nStripes_R2][nRadii];
  double x3b_R3[nStripes_R3][nRadii], y3b_R3[nStripes_R3][nRadii];

  //Check which stripes get removed
  int nGoodStripes_R1_e[nRadii];
  int nGoodStripes_R1[nRadii];
  int nGoodStripes_R2[nRadii];
  int nGoodStripes_R3[nRadii];
  int keepThisAndAfter[nRadii] = {1,0,1,0,1,0,1,0}; //min stripe index
  int keepUntil_R1_e[nRadii] = {4,4,5,4,5,5,5,5}; //max stripe index
  int keepUntil_R1[nRadii] = {5,5,6,5,6,5,6,5};
  int keepUntil_R2[nRadii] = {7,7,8,7,8,8,8,8};
  int keepUntil_R3[nRadii] = {11,10,11,11,11,11,12,11};
  
  Vertices(nStripes_R1, R1_e, spacing_R1_e, x1a_R1_e, y1a_R1_e, x1b_R1_e, y1b_R1_e, x2a_R1_e, y2a_R1_e, x2b_R1_e, y2b_R1_e, x3a_R1_e, y3a_R1_e, x3b_R1_e, y3b_R1_e, padfrac_R1, nGoodStripes_R1_e, keepThisAndAfter, keepUntil_R1_e);
  Vertices(nStripes_R1, R1, spacing_R1, x1a_R1, y1a_R1, x1b_R1, y1b_R1, x2a_R1, y2a_R1, x2b_R1, y2b_R1, x3a_R1, y3a_R1, x3b_R1, y3b_R1, padfrac_R1, nGoodStripes_R1, keepThisAndAfter, keepUntil_R1);
  Vertices(nStripes_R2, R2, spacing_R2, x1a_R2, y1a_R2, x1b_R2, y1b_R2, x2a_R2, y2a_R2, x2b_R2, y2b_R2, x3a_R2, y3a_R2, x3b_R2, y3b_R2, padfrac_R2, nGoodStripes_R2, keepThisAndAfter, keepUntil_R2);
  Vertices(nStripes_R3, R3, spacing_R3, x1a_R3, y1a_R3, x1b_R3, y1b_R3, x2a_R3, y2a_R3, x2b_R3, y2b_R3, x3a_R3, y3a_R3, x3b_R3, y3b_R3, padfrac_R3, nGoodStripes_R3, keepThisAndAfter, keepUntil_R3);  

  // Are you in a stripe?
  const double begin_CM = 221.4019814; // inner radius of CM
  const double end_CM = 759.2138; // outer radius of CM
  const double end_R1_e = 312.0; // arbitrary radius between R1_e and R1
  const double end_R1 = 408.0; // arbitrary radius between R1 and R2
  const double end_R2 = 580.0; // arbitrary radius between R2 and R3
  const double phi_petal = pi/9.0; // angle span of one petal
  int result, nbins;
  double x, y, r, phi;
  double xmod, ymod, phimod;

  nbins = 1000;
  
  //histogram from search
  TH2F *Pattern1 = new TH2F("Pattern1","Pattern1",nbins,-770.0,770.0,nbins,-770.0,770.0); // min n max just beyond extent of CM so it's easier to see
  for (r = begin_CM; r < end_CM; r = r + 0.5){ // radii spanning full CM
    for (phi = 0.0; phi < 2.0*pi; phi = phi + 0.00005){ // angles spanning full CM
      
      x = r*cos(phi);
      y = r*sin(phi);
      
      phimod = fmod(phi,phi_petal);
      xmod = r*cos(phimod);
      ymod = r*sin(phimod);
      
      if (r <= end_R1_e){
	result = search(nStripes_R1, x1a_R1_e, x1b_R1_e, x2a_R1_e, x2b_R1_e, y1a_R1_e, y1b_R1_e, y2a_R1_e, y2b_R1_e, x3a_R1_e, y3a_R1_e, x3b_R1_e, y3b_R1_e, xmod, ymod, nGoodStripes_R1_e, keepThisAndAfter, keepUntil_R1_e); 
      } else if ((r > end_R1_e) && (r <= end_R1)){
	result = search(nStripes_R1, x1a_R1, x1b_R1, x2a_R1, x2b_R1, y1a_R1, y1b_R1, y2a_R1, y2b_R1, x3a_R1, y3a_R1, x3b_R1, y3b_R1, xmod, ymod, nGoodStripes_R1, keepThisAndAfter, keepUntil_R1); 
      } else if ((r > end_R1) && (r <= end_R2)){
	result = search(nStripes_R2, x1a_R2, x1b_R2, x2a_R2, x2b_R2, y1a_R2, y1b_R2, y2a_R2, y2b_R2, x3a_R2, y3a_R2, x3b_R2, y3b_R2, xmod, ymod, nGoodStripes_R2, keepThisAndAfter, keepUntil_R2);
      } else if ((r > end_R2) && (r <= end_CM)){
	result = search(nStripes_R3, x1a_R3, x1b_R3, x2a_R3, x2b_R3, y1a_R3, y1b_R3, y2a_R3, y2b_R3, x3a_R3, y3a_R3, x3b_R3, y3b_R3, xmod, ymod, nGoodStripes_R3, keepThisAndAfter, keepUntil_R3);
      }
      
      if (result == 1)
	Pattern1->Fill(x,y);
    }
  }

  
   // build hits
  Container hit_R1_e[nStripes_R1][nRadii];
  Container hit_R1[nStripes_R1][nRadii];
  Container hit_R2[nStripes_R2][nRadii];
  Container hit_R3[nStripes_R3][nRadii];
  
  TVector3 pos0_R1_e[nStripes_R1][nRadii];
  TVector3 pos0_R1[nStripes_R1][nRadii];
  TVector3 pos0_R2[nStripes_R2][nRadii];
  TVector3 pos0_R3[nStripes_R3][nRadii];

  TVector3 pos1_R1_e[nStripes_R1][nRadii];
  TVector3 pos1_R1[nStripes_R1][nRadii];
  TVector3 pos1_R2[nStripes_R2][nRadii];
  TVector3 pos1_R3[nStripes_R3][nRadii];

  int nElectrons = 100;
  float time = 0.5;

  GetG4Hit(hit_R1_e, nStripes_R1, x3a_R1_e, y3a_R1_e, x3b_R1_e, y3b_R1_e, nGoodStripes_R1_e, keepThisAndAfter, keepUntil_R1_e, nElectrons, time, pos0_R1_e, pos1_R1_e);
  GetG4Hit(hit_R1, nStripes_R1, x3a_R1, y3a_R1, x3b_R1, y3b_R1, nGoodStripes_R1, keepThisAndAfter, keepUntil_R1, nElectrons, time, pos0_R1, pos1_R1);
  GetG4Hit(hit_R2, nStripes_R2, x3a_R2, y3a_R2, x3b_R2, y3b_R2, nGoodStripes_R2, keepThisAndAfter, keepUntil_R2, nElectrons, time, pos0_R2, pos1_R2);
  GetG4Hit(hit_R3, nStripes_R3, x3a_R3, y3a_R3, x3b_R3, y3b_R3, nGoodStripes_R3, keepThisAndAfter, keepUntil_R3, nElectrons, time, pos0_R3, pos1_R3);

  //histogram of start- and end-points
  TH2F *Pattern2 = new TH2F("Pattern2","Pattern2",nbins,-770.0,770.0,nbins,-770.0,770.0); // min n max just beyond extent of CM so it's easier to see
  Pattern2->SetMarkerColor(kRed);
  
  for(int j=0; j<nRadii; j++){
    for (int i=keepThisAndAfter[j]; i<nGoodStripes_R1_e[j]; i++){		
      Pattern2->Fill(hit_R1_e[i][j].x0, hit_R1_e[i][j].y0);
      Pattern2->Fill(hit_R1_e[i][j].x1, hit_R1_e[i][j].y1);
    }
    for (int i=keepThisAndAfter[j]; i<nGoodStripes_R1[j]; i++){
      Pattern2->Fill(hit_R1[i][j].x0, hit_R1[i][j].y0);
      Pattern2->Fill(hit_R1[i][j].x1, hit_R1[i][j].y1);
    }
    for (int i=keepThisAndAfter[j]; i<nGoodStripes_R2[j]; i++){
      Pattern2->Fill(hit_R2[i][j].x0, hit_R2[i][j].y0);
      Pattern2->Fill(hit_R2[i][j].x1, hit_R2[i][j].y1);
    }
    for (int i=keepThisAndAfter[j]; i<nGoodStripes_R3[j]; i++){
      Pattern2->Fill(hit_R3[i][j].x0, hit_R3[i][j].y0);
      Pattern2->Fill(hit_R3[i][j].x1, hit_R3[i][j].y1);
    }
  }

  //rotate to the rest of the petals
  for(int k=1; k<18; k++){
    for(int j=0; j<nRadii; j++){
      for (int i=keepThisAndAfter[j]; i<nGoodStripes_R1_e[j]; i++){
	pos0_R1_e[i][j].RotateZ(phi_petal);
	pos1_R1_e[i][j].RotateZ(phi_petal);
	Pattern2->Fill(pos0_R1_e[i][j].X(), pos0_R1_e[i][j].Y());
	Pattern2->Fill(pos1_R1_e[i][j].X(), pos1_R1_e[i][j].Y());	
      }
      for (int i=keepThisAndAfter[j]; i<nGoodStripes_R1[j]; i++){
	pos0_R1[i][j].RotateZ(phi_petal);
	pos1_R1[i][j].RotateZ(phi_petal);
	Pattern2->Fill(pos0_R1[i][j].X(), pos0_R1[i][j].Y());
	Pattern2->Fill(pos1_R1[i][j].X(), pos1_R1[i][j].Y());	
      }
      for (int i=keepThisAndAfter[j]; i<nGoodStripes_R2[j]; i++){
	pos0_R2[i][j].RotateZ(phi_petal);
	pos1_R2[i][j].RotateZ(phi_petal);
	Pattern2->Fill(pos0_R2[i][j].X(), pos0_R2[i][j].Y());
	Pattern2->Fill(pos1_R2[i][j].X(), pos1_R2[i][j].Y());	
      }
      for (int i=keepThisAndAfter[j]; i<nGoodStripes_R3[j]; i++){
	pos0_R3[i][j].RotateZ(phi_petal);
	pos1_R3[i][j].RotateZ(phi_petal);
	Pattern2->Fill(pos0_R3[i][j].X(), pos0_R3[i][j].Y());
	Pattern2->Fill(pos1_R3[i][j].X(), pos1_R3[i][j].Y());	
      }
    }
  }
  
  gStyle->SetOptStat(0);
  TCanvas *c=new TCanvas("a","cmhits.cpp",500,500);
  Pattern1->Draw();
  Pattern2->Draw("same");
  c->SaveAs("cmhits.pdf");
  
  return 0;
}


//vertex calculation function
void Vertices(int nStripes, double R[], double spacing[], double x1a[][nRadii], double y1a[][nRadii], double x1b[][nRadii], double y1b[][nRadii], double x2a[][nRadii], double y2a[][nRadii], double x2b[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double padfrac, int nGoodStripes[], int keepThisAndAfter[], int keepUntil[]){
  const double adjust = 0.015; //arbitrary angle to center the pattern in a petal
  double theta = 0.0;
  //center coords
  double cx[nStripes][nRadii], cy[nStripes][nRadii];
  //corner coords
  double tempX1a[nStripes][nRadii], tempY1a[nStripes][nRadii];
  double tempX1b[nStripes][nRadii], tempY1b[nStripes][nRadii];
  double tempX2a[nStripes][nRadii], tempY2a[nStripes][nRadii];
  double tempX2b[nStripes][nRadii], tempY2b[nStripes][nRadii];
  double rotatedX1a[nStripes][nRadii], rotatedY1a[nStripes][nRadii];
  double rotatedX1b[nStripes][nRadii], rotatedY1b[nStripes][nRadii];
  double rotatedX2a[nStripes][nRadii], rotatedY2a[nStripes][nRadii];
  double rotatedX2b[nStripes][nRadii], rotatedY2b[nStripes][nRadii];

  for (int j=0; j<(nRadii - 1); j=j+2){
    int i_out = keepThisAndAfter[j];
    for (int i=keepThisAndAfter[j]; i<nStripes; i++){
      theta = i_out*spacing[j];

      cx[i_out][j] = R[j]*cos(theta + (spacing[j]/2) -adjust);
      cy[i_out][j] = R[j]*sin(theta + (spacing[j]/2) -adjust);

      x1a[i_out][j] = cx[i_out][j] - padfrac + arc_r;
      y1a[i_out][j] = cy[i_out][j] - str_width/2;
      x1b[i_out][j] = cx[i_out][j] + padfrac - arc_r;
      y1b[i_out][j] = cy[i_out][j] - str_width/2;
      x2a[i_out][j] = cx[i_out][j] - padfrac + arc_r;
      y2a[i_out][j] = cy[i_out][j] + str_width/2;
      x2b[i_out][j] = cx[i_out][j] + padfrac - arc_r;
      y2b[i_out][j] = cy[i_out][j] + str_width/2;
      
      tempX1a[i_out][j] = x1a[i_out][j] - cx[i_out][j];
      tempY1a[i_out][j] = y1a[i_out][j] - cy[i_out][j];
      tempX1b[i_out][j] = x1b[i_out][j] - cx[i_out][j];
      tempY1b[i_out][j] = y1b[i_out][j] - cy[i_out][j];
      tempX2a[i_out][j] = x2a[i_out][j] - cx[i_out][j];
      tempY2a[i_out][j] = y2a[i_out][j] - cy[i_out][j];
      tempX2b[i_out][j] = x2b[i_out][j] - cx[i_out][j];
      tempY2b[i_out][j] = y2b[i_out][j] - cy[i_out][j];

      rotatedX1a[i_out][j] = tempX1a[i_out][j]*cos(theta) - tempY1a[i_out][j]*sin(theta);
      rotatedY1a[i_out][j] = tempX1a[i_out][j]*sin(theta) + tempY1a[i_out][j]*cos(theta);
      rotatedX1b[i_out][j] = tempX1b[i_out][j]*cos(theta) - tempY1b[i_out][j]*sin(theta);
      rotatedY1b[i_out][j] = tempX1b[i_out][j]*sin(theta) + tempY1b[i_out][j]*cos(theta);
      rotatedX2a[i_out][j] = tempX2a[i_out][j]*cos(theta) - tempY2a[i_out][j]*sin(theta);
      rotatedY2a[i_out][j] = tempX2a[i_out][j]*sin(theta) + tempY2a[i_out][j]*cos(theta);
      rotatedX2b[i_out][j] = tempX2b[i_out][j]*cos(theta) - tempY2b[i_out][j]*sin(theta);
      rotatedY2b[i_out][j] = tempX2b[i_out][j]*sin(theta) + tempY2b[i_out][j]*cos(theta);

      x1a[i_out][j] = rotatedX1a[i_out][j] + cx[i_out][j];
      y1a[i_out][j] = rotatedY1a[i_out][j] + cy[i_out][j];
      x1b[i_out][j] = rotatedX1b[i_out][j] + cx[i_out][j];
      y1b[i_out][j] = rotatedY1b[i_out][j] + cy[i_out][j];
      x2a[i_out][j] = rotatedX2a[i_out][j] + cx[i_out][j];
      y2a[i_out][j] = rotatedY2a[i_out][j] + cy[i_out][j];
      x2b[i_out][j] = rotatedX2b[i_out][j] + cx[i_out][j];
      y2b[i_out][j] = rotatedY2b[i_out][j] + cy[i_out][j];

      x3a[i_out][j] = (x1a[i_out][j] +  x2a[i_out][j])/ 2.0;
      y3a[i_out][j] = (y1a[i_out][j] +  y2a[i_out][j])/ 2.0;
      x3b[i_out][j] = (x1b[i_out][j] +  x2b[i_out][j])/ 2.0;
      y3b[i_out][j] = (y1b[i_out][j] +  y2b[i_out][j])/ 2.0;

      if(i<keepUntil[j]) i_out++;
    }
    nGoodStripes[j]=i_out;
  }

  for (int j=1; j<nRadii; j=j+2){
    int i_out = keepThisAndAfter[j];
    for (int i=keepThisAndAfter[j]; i<nStripes; i++){
      theta = (i_out+1)*spacing[j];

      cx[i_out][j] = R[j]*cos(theta-adjust);
      cy[i_out][j] = R[j]*sin(theta-adjust);
      
      x1a[i_out][j] = cx[i_out][j] - padfrac + arc_r;
      y1a[i_out][j] = cy[i_out][j] - str_width/2;
      x1b[i_out][j] = cx[i_out][j] + padfrac - arc_r;
      y1b[i_out][j] = cy[i_out][j] - str_width/2;
      x2a[i_out][j] = cx[i_out][j] - padfrac + arc_r;
      y2a[i_out][j] = cy[i_out][j] + str_width/2;
      x2b[i_out][j] = cx[i_out][j] + padfrac - arc_r;
      y2b[i_out][j] = cy[i_out][j] + str_width/2;
      
      tempX1a[i_out][j] = x1a[i_out][j] - cx[i_out][j];
      tempY1a[i_out][j] = y1a[i_out][j] - cy[i_out][j];
      tempX1b[i_out][j] = x1b[i_out][j] - cx[i_out][j];
      tempY1b[i_out][j] = y1b[i_out][j] - cy[i_out][j];
      tempX2a[i_out][j] = x2a[i_out][j] - cx[i_out][j];
      tempY2a[i_out][j] = y2a[i_out][j] - cy[i_out][j];
      tempX2b[i_out][j] = x2b[i_out][j] - cx[i_out][j];
      tempY2b[i_out][j] = y2b[i_out][j] - cy[i_out][j];

      rotatedX1a[i_out][j] = tempX1a[i_out][j]*cos(theta) - tempY1a[i_out][j]*sin(theta);
      rotatedY1a[i_out][j] = tempX1a[i_out][j]*sin(theta) + tempY1a[i_out][j]*cos(theta);
      rotatedX1b[i_out][j] = tempX1b[i_out][j]*cos(theta) - tempY1b[i_out][j]*sin(theta);
      rotatedY1b[i_out][j] = tempX1b[i_out][j]*sin(theta) + tempY1b[i_out][j]*cos(theta);
      rotatedX2a[i_out][j] = tempX2a[i_out][j]*cos(theta) - tempY2a[i_out][j]*sin(theta);
      rotatedY2a[i_out][j] = tempX2a[i_out][j]*sin(theta) + tempY2a[i_out][j]*cos(theta);
      rotatedX2b[i_out][j] = tempX2b[i_out][j]*cos(theta) - tempY2b[i_out][j]*sin(theta);
      rotatedY2b[i_out][j] = tempX2b[i_out][j]*sin(theta) + tempY2b[i_out][j]*cos(theta);

      x1a[i_out][j] = rotatedX1a[i_out][j] + cx[i_out][j];
      y1a[i_out][j] = rotatedY1a[i_out][j] + cy[i_out][j];
      x1b[i_out][j] = rotatedX1b[i_out][j] + cx[i_out][j];
      y1b[i_out][j] = rotatedY1b[i_out][j] + cy[i_out][j];
      x2a[i_out][j] = rotatedX2a[i_out][j] + cx[i_out][j];
      y2a[i_out][j] = rotatedY2a[i_out][j] + cy[i_out][j];
      x2b[i_out][j] = rotatedX2b[i_out][j] + cx[i_out][j];
      y2b[i_out][j] = rotatedY2b[i_out][j] + cy[i_out][j];

      x3a[i_out][j] = (x1a[i_out][j] +  x2a[i_out][j])/ 2.0;
      y3a[i_out][j] = (y1a[i_out][j] +  y2a[i_out][j])/ 2.0;
      x3b[i_out][j] = (x1b[i_out][j] +  x2b[i_out][j])/ 2.0;
      y3b[i_out][j] = (y1b[i_out][j] +  y2b[i_out][j])/ 2.0;

      if(i<keepUntil[j]) i_out++;
    }
    nGoodStripes[j]=i_out;
  }
}

//search function based on PnPoly problem by WRF:
//for (i = 0, j = nvert-1; i < nvert; j = i++) {
//if( ((verty[i] >= y) != (verty[j] >= y)) &&
//	(x <= (vertx[j]-vertx[i]) * (y-verty[i]) / (verty[j]-verty[i]) + vertx[i])) )
//c = !c; 
//}
//nvert is number of vertices of a polygon
//vertx,verty are x,y coordinates of a polygon
// our code has the vertices of a single polygon as the i,j entries of four different arrays
//the loop is:
// i=0 --> 1a, j=nvert-1 --> 2a
// i=1 --> 1b, j=0 --> 1a
// i=2 --> 2b, j=1 --> 1b
// i=3 --> 2a, j=2 --> 2b
int search(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double x, double y, int nGoodStripes[], int keepThisAndAfter[], int keepUntil[]){
  int c = 0;
  
  for(int j=0; j<nRadii; j++){
    for (int i=keepThisAndAfter[j]; i<nGoodStripes[j]; i++){
      if( ((y1a[i][j]>y) != (y2a[i][j]>y) && (x<(x2a[i][j]-x1a[i][j])*(y-y1a[i][j])/(y2a[i][j]-y1a[i][j])+x1a[i][j])))
	c = !c;
      if( ((y1b[i][j]>y) != (y1a[i][j]>y) && (x<(x1a[i][j]-x1b[i][j])*(y-y1b[i][j])/(y1a[i][j]-y1b[i][j])+x1b[i][j])))
	c = !c;
      if( ((y2b[i][j]>y) != (y1b[i][j]>y) && (x<(x1b[i][j]-x2b[i][j])*(y-y2b[i][j])/(y1b[i][j]-y2b[i][j])+x2b[i][j])))
	c = !c;
      if( ((y2a[i][j]>y) != (y2b[i][j]>y) && (x<(x2b[i][j]-x2a[i][j])*(y-y2a[i][j])/(y2b[i][j]-y2a[i][j])+x2a[i][j])))
	c = !c;

      //check inside arcs
      if (c==0){
	if (((x - x3a[i][j])*(x-x3a[i][j]) + (y-y3a[i][j])*(y-y3a[i][j])) <= arc_r*arc_r){
	  c =!c;
	} else if (((x - x3b[i][j])*(x-x3b[i][j]) + (y-y3b[i][j])*(y-y3b[i][j])) <= arc_r*arc_r){
	  c =!c;
	}
      }
    }
  }
  return c;
}

void GetG4Hit(Container hit[][nRadii], int nStripes, double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], int nGoodStripes[], int keepThisAndAfter[], int keepUntil[], int nElectrons, float time, TVector3 pos0[][nRadii], TVector3 pos1[][nRadii]){
  int savetrkid[nStripes][nRadii]; //?
  
  for(int j=0; j<nRadii; j++){
    for (int i=keepThisAndAfter[j]; i<nGoodStripes[j]; i++){
      //here we set the entrance values in mm
      hit[i][j].x0 = x3a[i][j] ; 
      hit[i][j].y0 = y3a[i][j] ;
      hit[i][j].z0 = 0;     

      pos0[i][j].SetXYZ(hit[i][j].x0, hit[i][j].y0, hit[i][j].z0);
      	
      // momentum
      hit[i][j].px0 = 0;
      hit[i][j].py0 = 0;
      hit[i][j].pz0 = 0;
            
      // time in ns
      hit[i][j].t0 = 0;

      //set and save the track ID
      hit[i][j].trkid = 0;
      savetrkid[i][j] = 0; //?

      //set the initial energy deposit
      hit[i][j].edep0 = 0;

      hit[i][j].eion0 = 0;

      // here we just update the exit values, it will be overwritten
      // for every step until we leave the volume or the particle
      // ceases to exist      
      hit[i][j].x1 = x3b[i][j] ; 
      hit[i][j].y1 = y3b[i][j] ;
      hit[i][j].z1 = 0;

      pos1[i][j].SetXYZ(hit[i][j].x1, hit[i][j].y1, hit[i][j].z1);

      hit[i][j].px1 = 500; // think abt what # would make sense
      hit[i][j].py1 = 500;
      hit[i][j].pz1 = 500;

      hit[i][j].t1 = time;

      //sum up the energy to get total deposited
      hit[i][j].edep1 = hit[i][j].edep0; // still need to make smth for get_edep 

      hit[i][j].eion1 = hit[i][j].eion0;
    }
  }
}

