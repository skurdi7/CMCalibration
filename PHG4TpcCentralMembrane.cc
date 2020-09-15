#include "PHG4TpcCentralMembrane.h"

#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TH2F.h"


//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
//R__LOAD_LIBRARY(libphg4hit.so)

// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps

using namespace std;

StripesClass::StripesClass()
  : R1_e {227.0902789 * mm, 238.4100043 * mm, 249.7297296 * mm, 261.049455 * mm, 272.3691804 * mm, 283.6889058 * mm, 295.0086312 * mm, 306.3283566 * mm},
    R1 {317.648082 * mm, 328.9678074 * mm, 340.2875328 * mm, 351.6072582 * mm, 362.9269836 * mm, 374.246709 * mm,  385.5664344 * mm, 396.8861597 * mm},
    R2 {421.705532 * mm, 442.119258 * mm, 462.532984 * mm, 482.9467608 * mm, 503.36069 * mm, 523.774416 * mm, 544.188015 * mm, 564.601868 * mm},
    R3 {594.6048725 * mm, 616.545823 * mm, 638.4867738 * mm, 660.4277246 * mm, 682.3686754 * mm, 704.3096262 * mm, 726.250577 * mm, 748.1915277 * mm},
    keepThisAndAfter {1,0,1,0,1,0,1,0},
    keepUntil_R1_e {4,4,5,4,5,5,5,5},
    keepUntil_R1 {5,5,6,5,6,5,6,5},
    keepUntil_R2 {7,7,8,7,8,8,8,8},
    keepUntil_R3 {11,10,11,11,11,11,12,11}
{
  begin_CM = 221.4019814 * mm; // inner radius of CM
  end_CM = 759.2138 * mm; // outer radius of CM
  
  nPads_R1 = 6*16;
  nPads_R2 = 8*16;
  nPads_R3 = 12*16;
  
  padfrac_R1 = 0.5*5.59106385 * mm;
  padfrac_R2 = 0.5*10.13836283 * mm;
  padfrac_R3 = 0.5*10.90189537 * mm;
  
  str_width = 1.0 * mm;
  arc_r = 0.5 * mm;

  nStripesPerPetal = 213;
  nPetals = 18;

  nTotStripes = nStripesPerPetal * nPetals;
  
  nElectrons = 100;

  CalculateVertices(nStripes_R1, nPads_R1, R1_e, spacing_R1_e, x1a_R1_e, y1a_R1_e, x1b_R1_e, y1b_R1_e, x2a_R1_e, y2a_R1_e, x2b_R1_e, y2b_R1_e, x3a_R1_e, y3a_R1_e, x3b_R1_e, y3b_R1_e, padfrac_R1, nGoodStripes_R1_e, keepUntil_R1_e);
  CalculateVertices(nStripes_R1, nPads_R1, R1, spacing_R1, x1a_R1, y1a_R1, x1b_R1, y1b_R1, x2a_R1, y2a_R1, x2b_R1, y2b_R1, x3a_R1, y3a_R1, x3b_R1, y3b_R1, padfrac_R1, nGoodStripes_R1, keepUntil_R1);
  CalculateVertices(nStripes_R2, nPads_R2, R2, spacing_R2, x1a_R2, y1a_R2, x1b_R2, y1b_R2, x2a_R2, y2a_R2, x2b_R2, y2b_R2, x3a_R2, y3a_R2, x3b_R2, y3b_R2, padfrac_R2, nGoodStripes_R2, keepUntil_R2);
  CalculateVertices(nStripes_R3, nPads_R3, R3, spacing_R3, x1a_R3, y1a_R3, x1b_R3, y1b_R3, x2a_R3, y2a_R3, x2b_R3, y2b_R3, x3a_R3, y3a_R3, x3b_R3, y3b_R3, padfrac_R3, nGoodStripes_R3, keepUntil_R3);
   
  for (int i = 0; i < 18; i++){ // loop over petalID
    for (int j = 0; j < 8; j++){ // loop over radiusID
      for (int k = 0; k < nGoodStripes_R1_e[j]; k++){ // loop over stripeID
	PHG4Hits.push_back(GetPHG4HitFromStripe(i, 0, j, k, nElectrons));
      }
      for (int k = 0; k < nGoodStripes_R1[j]; k++){ // loop over stripeID
	PHG4Hits.push_back(GetPHG4HitFromStripe(i, 1, j, k, nElectrons));
      }
      for (int k = 0; k < nGoodStripes_R2[j]; k++){ // loop over stripeID
	PHG4Hits.push_back(GetPHG4HitFromStripe(i, 2, j, k, nElectrons));
      }
      for (int k = 0; k < nGoodStripes_R3[j]; k++){ // loop over stripeID
	PHG4Hits.push_back(GetPHG4HitFromStripe(i, 3, j, k, nElectrons));
      }
    }
  }
  
 return;
}


void StripesClass::CalculateVertices(int nStripes, int nPads, double R[], double spacing[], double x1a[][nRadii], double y1a[][nRadii], double x1b[][nRadii], double y1b[][nRadii], double x2a[][nRadii], double y2a[][nRadii], double x2b[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double padfrac, int nGoodStripes[],  int keepUntil[]) {
  const double phi_module = TMath::Pi()/6.0; // angle span of a module
  const int pr_mult = 3; // multiples of intrinsic resolution of pads
  const int dw_mult = 8; // multiples of diffusion width
  const double diffwidth = 0.6 * mm; // diffusion width
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

  //calculate spacing first:
  for (int i=0; i<nRadii; i++){
    spacing[i] = 2.0*((dw_mult*diffwidth/R[i]) + (pr_mult*phi_module/nPads));
  }
  
  //vertex calculation
  for (int j=0; j<nRadii; j++){
    int i_out = 0;
    for (int i=keepThisAndAfter[j]; i<keepUntil[j]; i++){
      if (j % 2 == 0){
	theta = i*spacing[j];
	cx[i_out][j] = R[j]*cos(theta + (spacing[j]/2) -adjust);
	cy[i_out][j] = R[j]*sin(theta + (spacing[j]/2) -adjust);
      } else {
	theta = (i+1)*spacing[j];
	cx[i_out][j] = R[j]*cos(theta-adjust);
	cy[i_out][j] = R[j]*sin(theta-adjust);
      }

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

      i_out++;
    }
    nGoodStripes[j]=i_out;
  }
}


int StripesClass::SearchModule(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x3a[][nRadii], double y3a[][nRadii], double x3b[][nRadii], double y3b[][nRadii], double x, double y, int nGoodStripes[]){
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

int StripesClass::getSearchResult(double xcheck, double ycheck){
  const double phi_petal = TMath::Pi()/9.0; // angle span of one petal
  const double end_R1_e = 312.0 * mm; // arbitrary radius between R1_e and R1
  const double end_R1 = 408.0 * mm; // arbitrary radius between R1 and R2
  const double end_R2 = 580.0 * mm; // arbitrary radius between R2 and R3

  double r, phi, phimod, xmod, ymod;
  
  r = sqrt(xcheck*xcheck + ycheck*ycheck);
  phi = atan(ycheck/xcheck);
  if((xcheck < 0.0) && (ycheck > 0.0)){
    phi = phi + TMath::Pi();
  } else if ((xcheck > 0.0) && (ycheck < 0.0)){
    phi = phi + 2.0*TMath::Pi();
  }
  
  phimod = fmod(phi,phi_petal);
  xmod = r*cos(phimod);
  ymod = r*sin(phimod); 
  
  if (r <= end_R1_e){ 
    result = SearchModule(nStripes_R1, x1a_R1_e, x1b_R1_e, x2a_R1_e, x2b_R1_e, y1a_R1_e, y1b_R1_e, y2a_R1_e, y2b_R1_e, x3a_R1_e, y3a_R1_e, x3b_R1_e, y3b_R1_e, xmod, ymod, nGoodStripes_R1_e);
  } else if ((r > end_R1_e) && (r <= end_R1)){
    result = SearchModule(nStripes_R1, x1a_R1, x1b_R1, x2a_R1, x2b_R1, y1a_R1, y1b_R1, y2a_R1, y2b_R1, x3a_R1, y3a_R1, x3b_R1, y3b_R1, xmod, ymod, nGoodStripes_R1);
  } else if ((r > end_R1) && (r <= end_R2)){
    result = SearchModule(nStripes_R2, x1a_R2, x1b_R2, x2a_R2, x2b_R2, y1a_R2, y1b_R2, y2a_R2, y2b_R2, x3a_R2, y3a_R2, x3b_R2, y3b_R2, xmod, ymod, nGoodStripes_R2);
  } else if ((r > end_R2) && (r <= end_CM)){
    result = SearchModule(nStripes_R3, x1a_R3, x1b_R3, x2a_R3, x2b_R3, y1a_R3, y1b_R3, y2a_R3, y2b_R3, x3a_R3, y3a_R3, x3b_R3, y3b_R3, xmod, ymod, nGoodStripes_R3);
  }
  
  return result;
}

PHG4Hitv1* StripesClass::GetPHG4HitFromStripe(int petalID, int moduleID, int radiusID, int stripeID, int nElectrons)
{ //this function generates a PHG4 hit using coordinates from a stripe
  const double phi_petal = TMath::Pi()/9.0; // angle span of one petal
  PHG4Hitv1 *hit;
  TVector3 dummyPos0, dummyPos1;
  
  //could put in some sanity checks here but probably not necessary since this is only really used within the class
  //petalID ranges 0-17, module ID 0-3, stripeID varies - nGoodStripes for each module
  //radiusID ranges 0-7

  //from phg4tpcsteppingaction.cc
  hit = new PHG4Hitv1();
  hit->set_layer(-1); // dummy number 
  //here we set the entrance values in cm
  if (moduleID == 0){
    hit->set_x(0, x3a_R1_e[stripeID][radiusID] / cm);
    hit->set_y(0, y3a_R1_e[stripeID][radiusID] / cm);
  } else if (moduleID == 1){
    hit->set_x(0, x3a_R1[stripeID][radiusID] / cm);
    hit->set_y(0, y3a_R1[stripeID][radiusID] / cm);
  } else if (moduleID == 2){
    hit->set_x(0, x3a_R2[stripeID][radiusID] / cm);
    hit->set_y(0, y3a_R2[stripeID][radiusID] / cm);
  } else if (moduleID == 3){
    hit->set_x(0, x3a_R3[stripeID][radiusID] / cm);
    hit->set_y(0, y3a_R3[stripeID][radiusID] / cm);
  }
  hit->set_z(0, 0.0 / cm);

  // check if you need to rotate coords to another petal
  if(petalID > 0){
    dummyPos0.SetXYZ(hit->get_x(0), hit->get_y(0), hit->get_z(0));
    dummyPos0.RotateZ(petalID * phi_petal);
    hit->set_x(0, dummyPos0.X());
    hit->set_y(0, dummyPos0.Y());
  }
  
  // momentum
  hit->set_px(0, 0.0); // GeV
  hit->set_py(0, 0.0);
  hit->set_pz(0, 0.0);
  
  // time in ns
  hit->set_t(0, 0.0); // nanosecond
  //set and save the track ID
  hit->set_trkid(-1); // dummy number

  // here we just update the exit values, it will be overwritten
  // for every step until we leave the volume or the particle
  // ceases to exist
  if (moduleID == 0){
    hit->set_x(1, x3b_R1_e[stripeID][radiusID] / cm);
    hit->set_y(1, y3b_R1_e[stripeID][radiusID] / cm);
  } else if (moduleID == 1){
    hit->set_x(1, x3b_R1[stripeID][radiusID] / cm);
    hit->set_y(1, y3b_R1[stripeID][radiusID] / cm);
  } else if (moduleID == 2){
    hit->set_x(1, x3b_R2[stripeID][radiusID] / cm);
    hit->set_y(1, y3b_R2[stripeID][radiusID] / cm);
  } else if (moduleID == 3){
    hit->set_x(1, x3b_R3[stripeID][radiusID] / cm);
    hit->set_y(1, y3b_R3[stripeID][radiusID] / cm);
  }
  hit->set_z(1, 0.0 / cm);

  // check if you need to rotate coords to another petal
  if(petalID > 0){
    dummyPos1.SetXYZ(hit->get_x(1), hit->get_y(1), hit->get_z(1));
    dummyPos1.RotateZ(petalID * phi_petal);
    hit->set_x(1, dummyPos1.X());
    hit->set_y(1, dummyPos1.Y());
  }
  
  hit->set_px(1, 500.0); // dummy large number, in GeV
  hit->set_py(1, 500.0);
  hit->set_pz(1, 500.0);
  
  hit->set_t(1, 1.0); // dummy number, nanosecond

  //calculate the total energy deposited
  
  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm

  //double Ne_NTotal = 43;    // Number/cm
  //double CF4_NTotal = 100;  // Number/cm
  //double Tpc_NTot = 0.90 * Ne_NTotal + 0.10 * CF4_NTotal;

  double Tpc_NTot = nElectrons;
  double Tpc_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;

  //double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  //double Tpc_ElectronsPerGeV = Tpc_NTot / Tpc_dEdx*1e6; //electrons per gev.

  double edep = Tpc_dEdx*1e6 / Tpc_NTot; // GeV dep per electron
  hit->set_edep(edep); // dont need get edep
  //calculate eion - make same as edep
  hit->set_eion(edep);// dont need get eion

  /*
  if (hit->get_edep()){ //print out hits
    double rin = sqrt(hit->get_x(0) * hit->get_x(0) + hit->get_y(0) * hit->get_y(0));
    double rout = sqrt(hit->get_x(1) * hit->get_x(1) + hit->get_y(1) * hit->get_y(1));
    cout << "Added Tpc g4hit with rin, rout = " << rin << "  " << rout
	 << " g4hitid " << hit->get_hit_id() << endl;
    cout << " xin " << hit->get_x(0)
	 << " yin " << hit->get_y(0)
	 << " zin " << hit->get_z(0)
	 << " rin " << rin
	 << endl;
    cout << " xout " << hit->get_x(1)
	 << " yout " << hit->get_y(1)
	 << " zout " << hit->get_z(1)
	 << " rout " << rout
	 << endl;
    cout << " xav " << (hit->get_x(1) + hit->get_x(0)) / 2.0
	 << " yav " << (hit->get_y(1) + hit->get_y(0)) / 2.0
	 << " zav " << (hit->get_z(1) + hit->get_z(0)) / 2.0
	 << " rav " << (rout + rin) / 2.0
	 << endl;
  }
  */
  
  return hit;
}

int StripesClass::getStripeID(double xcheck, double ycheck){
  //check if point came from stripe then see which stripe it is
  //213 stripes in a petal, 18 petals, ntotstripes = 3834
  int result, rID, petalID, nStripesPerR;
  int phiID = 0;
  int fullID = -1;
  //double theta, spacing[nRadii], angle, m, dist;
  double m, dist;
  //const double adjust = 0.015; //arbitrary angle to center the pattern in a petal
  const double phi_petal = TMath::Pi()/9.0; // angle span of one petal

  double r, phi;
  //phimod, xmod, ymod;

  // check if in a stripe
  result = getSearchResult(xcheck, ycheck);
  
  // find which stripe
  if(result == 1){
    
    //convert coords to radius n angle
    r = sqrt(xcheck*xcheck + ycheck*ycheck);
    phi = atan(ycheck/xcheck);
    if((xcheck < 0.0) && (ycheck > 0.0)){
      phi = phi + TMath::Pi();
    } else if ((xcheck > 0.0) && (ycheck < 0.0)){
      phi = phi + 2.0*TMath::Pi();
    }
    //get angle within first petal
    //phimod = fmod(phi,phi_petal);
    //xmod = r*cos(phimod);
    //ymod = r*sin(phimod);

    petalID = phi/phi_petal; 
    
    for(int j=0; j<nRadii; j++){
      if(((R1_e[j] - padfrac_R1) < r) && (r < (R1_e[j] + padfrac_R1))){ // check if radius is in stripe 
	rID = j; 
	//cout << "R1_E" << endl;
	//'angle' is to the center of a stripe
	for (int i=keepThisAndAfter[j]; i<keepUntil_R1_e[j]; i++){
	  //if (j % 2 == 0){
	  //theta = i*spacing[j];
	  //angle = theta + (spacing[j]/2) - adjust;
	    // look at distance from center line of stripe
	    // if distance from x,y to center line < str_width
	    //dist = fabs((y3b[i][j] - y3a[i][j])*xcheck - (x3b[i][j] - x3a[i][j])*ycheck + x3b[i][j]*y3a[i][j] - y3b[i][j]*x3a[i][j])/sqrt((y3b[i][j]-y3a[i][j])*(y3b[i][j]-y3a[i][j]) + (x3b[i][j]-x3a[i][j])*(x3b[i][j]-x3a[i][j]));
	    // or calculate slope n then do dist

	  m = (y3b_R1_e[i][j] - y3a_R1_e[i][j])/(x3b_R1_e[i][j] - x3a_R1_e[i][j]);
	  //cout << "m: " << m << endl;
	  cout << fabs((-m)*xcheck + ycheck) << endl;
	  dist = fabs((-m)*xcheck + ycheck)/sqrt(1 + m*m);
	  //cout << "dist: " << dist << endl;
	  if(dist < (str_width/2.0)){ 
	    phiID = i;
	    cout << "phiID: " << phiID << endl;
	  }
	    //} else {
	    //theta = (i+1)*spacing[j];
	    //angle = theta-adjust;
	    //m = (y3b_R1_e[i][j] - y3a_R1_e[i][j])/(y3b_R1_e[i][j] - y3a_R1_e[i][j]);
	    //dist = fabs(m*xcheck - ycheck)/sqrt(1 + m*m);
	    //if(dist < str_width){ 
	    // phiID = i;
	    //}
		  
	}
	nStripesPerR = keepUntil_R1_e[j] - keepThisAndAfter[j];
	fullID = petalID*nStripesPerPetal + rID*nStripesPerR + phiID;
	//cout << "fullID: " << fullID << endl;
      } else if (((R1[j]- padfrac_R1) < r) && (r < (R1[j]+ padfrac_R1))){
	rID = j+nRadii;
	//cout << "R1" << endl;
	for (int i=keepThisAndAfter[j]; i<keepUntil_R1[j]; i++){
	  // look at distance from center line of stripe
	  m = (y3b_R1[i][j] - y3a_R1[i][j])/(x3b_R1[i][j] - x3a_R1[i][j]);
	  dist = fabs(m*xcheck - ycheck)/sqrt(1 + m*m);
	  if(dist < (str_width/2.0)){ 
	    phiID = i;
	  }
	}
	nStripesPerR = keepUntil_R1[j] - keepThisAndAfter[j];
	fullID = petalID*nStripesPerPetal + rID*nStripesPerR + phiID;
	
      } else if (((R2[j]- padfrac_R2) < r) && (r < (R2[j]+ padfrac_R2))){
	rID = j+(2*nRadii);
	//cout << "R2" << endl;
	for (int i=keepThisAndAfter[j]; i<keepUntil_R2[j]; i++){
	  // look at distance from center line of stripe
	  m = (y3b_R2[i][j] - y3a_R2[i][j])/(x3b_R2[i][j] - x3a_R2[i][j]);
	  dist = fabs(m*xcheck - ycheck)/sqrt(1 + m*m);
	  if(dist < (str_width/2.0)){ 
	    phiID = i;
	  }
	}	  
	
	nStripesPerR = keepUntil_R2[j] - keepThisAndAfter[j];
	fullID = petalID*nStripesPerPetal + rID*nStripesPerR + phiID;
	
      } else if (((R3[j]- padfrac_R3) < r) && (r < (R3[j]+ padfrac_R3))){
	rID = j+(3*nRadii);
	//cout << "R3" << endl;
	for (int i=keepThisAndAfter[j]; i<keepUntil_R3[j]; i++){
	  // look at distance from center line of stripe
	  m = (y3b_R3[i][j] - y3a_R3[i][j])/(x3b_R3[i][j] - x3a_R3[i][j]);
	  dist = fabs(m*xcheck - ycheck)/sqrt(1 + m*m);
	  if(dist < (str_width/2.0)){ 
	    phiID = i;
	  }
	}
	nStripesPerR = keepUntil_R3[j] - keepThisAndAfter[j];
	fullID = petalID*nStripesPerPetal + rID*nStripesPerR + phiID;
      }
    }
  } else {
    fullID = -1;
  }
  
  return fullID;
      
}
