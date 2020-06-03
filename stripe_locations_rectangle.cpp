#include <iostream>
#include <cmath>
using namespace std;

// positions of stripes (rectangles for now)
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps
// all distances in mm, all angles in rad

const int nRadii = 8;

void corners(int nStripes, double R[], double spacing[], double x1a[][nRadii], double y1a[][nRadii], double x1b[][nRadii], double y1b[][nRadii], double x2a[][nRadii], double y2a[][nRadii], double x2b[][nRadii], double y2b[][nRadii], double padfrac);

int searchR1_e(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x, double y);

int searchR1(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x, double y);

int searchR2(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x, double y);

int searchR3(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x, double y);

int main() {
  
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
  const int nPads_R1 = 6*16; //pads per row (256 pads per card in module, 16 rows in module), same for R1_e
  const int nPads_R2 = 8*16;
  const int nPads_R3 = 12*16;
  double spacing_R1_e[nRadii], spacing_R1[nRadii], spacing_R2[nRadii], spacing_R3[nRadii];

  for (int i=0; i<nRadii; i++){
    spacing_R1_e[i] = 2.0*((8.0*0.6/R1_e[i]) + (3.0*pi/6.0/nPads_R1));
    spacing_R1[i] = 2.0*((8.0*0.6/R1[i]) + (3.0*pi/6.0/nPads_R1));
    spacing_R2[i] = 2.0*((8.0*0.6/R2[i]) + (3.0*pi/6.0/nPads_R2));
    spacing_R3[i] = 2.0*((8.0*0.6/R3[i]) + (3.0*pi/6.0/nPads_R3));
  }

  
  // Corner coordinates calculation
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

  corners(nStripes_R1, R1_e, spacing_R1_e, x1a_R1_e, y1a_R1_e, x1b_R1_e, y1b_R1_e, x2a_R1_e, y2a_R1_e, x2b_R1_e, y2b_R1_e, padfrac_R1);
  corners(nStripes_R1, R1, spacing_R1, x1a_R1, y1a_R1, x1b_R1, y1b_R1, x2a_R1, y2a_R1, x2b_R1, y2b_R1, padfrac_R1);
  corners(nStripes_R2, R2, spacing_R2, x1a_R2, y1a_R2, x1b_R2, y1b_R2, x2a_R2, y2a_R2, x2b_R2, y2b_R2, padfrac_R2);
  corners(nStripes_R3, R3, spacing_R3, x1a_R3, y1a_R3, x1b_R3, y1b_R3, x2a_R3, y2a_R3, x2b_R3, y2b_R3, padfrac_R3);  

  
  // Are you in a stripe?
  char ans;
  int resultR1_e;
  int resultR1;
  int resultR2;
  int resultR3;
  double x, y;
  
  do {
    cout << "Enter x,y coordinates: ";
    cin >> x;
    cin >> y;

    //I did this part as functions just to keep them out of the way, but they could all go here.
    resultR1_e = searchR1_e(nStripes_R1, x1a_R1_e, x1b_R1_e, x2a_R1_e, x2b_R1_e, y1a_R1_e, y1b_R1_e, y2a_R1_e, y2b_R1_e, x, y);
    resultR1 = searchR1(nStripes_R1, x1a_R1, x1b_R1, x2a_R1, x2b_R1, y1a_R1, y1b_R1, y2a_R1, y2b_R1, x, y);
    resultR2 = searchR2(nStripes_R2, x1a_R2, x1b_R2, x2a_R2, x2b_R2, y1a_R2, y1b_R2, y2a_R2, y2b_R2, x, y);
    resultR3 = searchR3(nStripes_R3, x1a_R3, x1b_R3, x2a_R3, x2b_R3, y1a_R3, y1b_R3, y2a_R3, y2b_R3, x, y);
    
    if ( (resultR1_e == 1) || (resultR1 == 1) || (resultR2 == 1) || (resultR3 == 1) )
      cout << "The point is in a stripe." << endl;
    else
      cout << "The point is not in a stripe." << endl;
    cout << "Search again? (y/n): ";
    cin >> ans;
  } while ((ans != 'n'));
  
  return 0;
}

//corners calculation function
void corners(int nStripes, double R[], double spacing[], double x1a[][nRadii], double y1a[][nRadii], double x1b[][nRadii], double y1b[][nRadii], double x2a[][nRadii], double y2a[][nRadii], double x2b[][nRadii], double y2b[][nRadii], double padfrac){
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
    for (int i=0; i<nStripes; i++){
      theta = i*spacing[j];

      cx[i][j] = R[j]*cos(theta + (spacing[j]/2) -adjust);
      cy[i][j] = R[j]*sin(theta + (spacing[j]/2) -adjust);

      x1a[i][j] = cx[i][j] - padfrac + 0.5;
      y1a[i][j] = cy[i][j] - 0.5;
      x1b[i][j] = cx[i][j] + padfrac - 0.5;
      y1b[i][j] = cy[i][j] - 0.5;
      x2a[i][j] = cx[i][j] - padfrac + 0.5;
      y2a[i][j] = cy[i][j] + 0.5;
      x2b[i][j] = cx[i][j] + padfrac - 0.5;
      y2b[i][j] = cy[i][j] + 0.5;
      
      tempX1a[i][j] = x1a[i][j] - cx[i][j];
      tempY1a[i][j] = y1a[i][j] - cy[i][j];
      tempX1b[i][j] = x1b[i][j] - cx[i][j];
      tempY1b[i][j] = y1b[i][j] - cy[i][j];
      tempX2a[i][j] = x2a[i][j] - cx[i][j];
      tempY2a[i][j] = y2a[i][j] - cy[i][j];
      tempX2b[i][j] = x2b[i][j] - cx[i][j];
      tempY2b[i][j] = y2b[i][j] - cy[i][j];

      rotatedX1a[i][j] = tempX1a[i][j]*cos(theta) - tempY1a[i][j]*sin(theta);
      rotatedY1a[i][j] = tempX1a[i][j]*sin(theta) + tempY1a[i][j]*cos(theta);
      rotatedX1b[i][j] = tempX1b[i][j]*cos(theta) - tempY1b[i][j]*sin(theta);
      rotatedY1b[i][j] = tempX1b[i][j]*sin(theta) + tempY1b[i][j]*cos(theta);
      rotatedX2a[i][j] = tempX2a[i][j]*cos(theta) - tempY2a[i][j]*sin(theta);
      rotatedY2a[i][j] = tempX2a[i][j]*sin(theta) + tempY2a[i][j]*cos(theta);
      rotatedX2b[i][j] = tempX2b[i][j]*cos(theta) - tempY2b[i][j]*sin(theta);
      rotatedY2b[i][j] = tempX2b[i][j]*sin(theta) + tempY2b[i][j]*cos(theta);

      x1a[i][j] = rotatedX1a[i][j] + cx[i][j];
      y1a[i][j] = rotatedY1a[i][j] + cy[i][j];
      x1b[i][j] = rotatedX1b[i][j] + cx[i][j];
      y1b[i][j] = rotatedY1b[i][j] + cy[i][j];
      x2a[i][j] = rotatedX2a[i][j] + cx[i][j];
      y2a[i][j] = rotatedY2a[i][j] + cy[i][j];
      x2b[i][j] = rotatedX2b[i][j] + cx[i][j];
      y2b[i][j] = rotatedY2b[i][j] + cy[i][j];
    }
  }

  for (int j=1; j<nRadii; j=j+2){
    for (int i=0; i<nStripes; i++){
      theta = (i+1)*spacing[j];

      cx[i][j] = R[j]*cos(theta-adjust);
      cy[i][j] = R[j]*sin(theta-adjust);

      x1a[i][j] = cx[i][j] - padfrac + 0.5;
      y1a[i][j] = cy[i][j] - 0.5;
      x1b[i][j] = cx[i][j] + padfrac - 0.5;
      y1b[i][j] = cy[i][j] - 0.5;
      x2a[i][j] = cx[i][j] - padfrac + 0.5;
      y2a[i][j] = cy[i][j] + 0.5;
      x2b[i][j] = cx[i][j] + padfrac - 0.5;
      y2b[i][j] = cy[i][j] + 0.5;
      
      tempX1a[i][j] = x1a[i][j] - cx[i][j];
      tempY1a[i][j] = y1a[i][j] - cy[i][j];
      tempX1b[i][j] = x1b[i][j] - cx[i][j];
      tempY1b[i][j] = y1b[i][j] - cy[i][j];
      tempX2a[i][j] = x2a[i][j] - cx[i][j];
      tempY2a[i][j] = y2a[i][j] - cy[i][j];
      tempX2b[i][j] = x2b[i][j] - cx[i][j];
      tempY2b[i][j] = y2b[i][j] - cy[i][j];

      rotatedX1a[i][j] = tempX1a[i][j]*cos(theta) - tempY1a[i][j]*sin(theta);
      rotatedY1a[i][j] = tempX1a[i][j]*sin(theta) + tempY1a[i][j]*cos(theta);
      rotatedX1b[i][j] = tempX1b[i][j]*cos(theta) - tempY1b[i][j]*sin(theta);
      rotatedY1b[i][j] = tempX1b[i][j]*sin(theta) + tempY1b[i][j]*cos(theta);
      rotatedX2a[i][j] = tempX2a[i][j]*cos(theta) - tempY2a[i][j]*sin(theta);
      rotatedY2a[i][j] = tempX2a[i][j]*sin(theta) + tempY2a[i][j]*cos(theta);
      rotatedX2b[i][j] = tempX2b[i][j]*cos(theta) - tempY2b[i][j]*sin(theta);
      rotatedY2b[i][j] = tempX2b[i][j]*sin(theta) + tempY2b[i][j]*cos(theta);

      x1a[i][j] = rotatedX1a[i][j] + cx[i][j];
      y1a[i][j] = rotatedY1a[i][j] + cy[i][j];
      x1b[i][j] = rotatedX1b[i][j] + cx[i][j];
      y1b[i][j] = rotatedY1b[i][j] + cy[i][j];
      x2a[i][j] = rotatedX2a[i][j] + cx[i][j];
      y2a[i][j] = rotatedY2a[i][j] + cy[i][j];
      x2b[i][j] = rotatedX2b[i][j] + cx[i][j];
      y2b[i][j] = rotatedY2b[i][j] + cy[i][j];
    }
  }
}


//search functions based on PnPoly problem by WRF:
//for (i = 0, j = nvert-1; i < nvert; j = i++) {
//if( ((verty[i] >= y) != (verty[j] >= y)) &&
//	(x <= (vertx[j]-vertx[i]) * (y-verty[i]) / (verty[j]-verty[i]) + vertx[i])) )
//c = !c; 
//}
//nvert is number of vertices of a polygon
//vertx,verty are x,y coordinates of a polygon
// our code has the vertices a single polygon as the i,j entries of four different arrays
//the loop is:
// i=0 --> 1a, j=nvert-1 --> 2a
// i=1 --> 1b, j=0 --> 1a
// i=2 --> 2b, j=1 --> 1b
// i=3 --> 2a, j=2 --> 2b
int searchR1_e(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x, double y){
  int c = 0;
  
  for(int j=0; j<nRadii; j++){
    for(int i=0; i<nStripes; i++){
      if( ((y1a[i][j]>y) != (y2a[i][j]>y) && (x<(x2a[i][j]-x1a[i][j])*(y-y1a[i][j])/(y2a[i][j]-y1a[i][j])+x1a[i][j])))
	c = !c;
      if( ((y1b[i][j]>y) != (y1a[i][j]>y) && (x<(x1a[i][j]-x1b[i][j])*(y-y1b[i][j])/(y1a[i][j]-y1b[i][j])+x1b[i][j])))
	c = !c;
      if( ((y2b[i][j]>y) != (y1b[i][j]>y) && (x<(x1b[i][j]-x2b[i][j])*(y-y2b[i][j])/(y1b[i][j]-y2b[i][j])+x2b[i][j])))
	c = !c;
      if( ((y2a[i][j]>y) != (y2b[i][j]>y) && (x<(x2b[i][j]-x2a[i][j])*(y-y2a[i][j])/(y2b[i][j]-y2a[i][j])+x2a[i][j])))
	c = !c;
    }
  }

  //special case: stripe 0 in radii 0, 2, 4, 6 removed from all module sections
  for(int k=0; k<nRadii; k=k+2){
    if( ((y1a[0][k]>y) != (y2a[0][k]>y) && (x<(x2a[0][k]-x1a[0][k])*(y-y1a[0][k])/(y2a[0][k]-y1a[0][k])+x1a[0][k])))
      c = !c;
    if( ((y1b[0][k]>y) != (y1a[0][k]>y) && (x<(x1a[0][k]-x1b[0][k])*(y-y1b[0][k])/(y1a[0][k]-y1b[0][k])+x1b[0][k])))
      c = !c;
    if( ((y2b[0][k]>y) != (y1b[0][k]>y) && (x<(x1b[0][k]-x2b[0][k])*(y-y2b[0][k])/(y1b[0][k]-y2b[0][k])+x2b[0][k])))
      c = !c;
    if( ((y2a[0][k]>y) != (y2b[0][k]>y) && (x<(x2b[0][k]-x2a[0][k])*(y-y2a[0][k])/(y2b[0][k]-y2a[0][k])+x2a[0][k])))
      c = !c;
  }

  //special case: stripe 4 in radii 0, 1, 3 removed from R1_e
  if( ((y1a[4][0]>y) != (y2a[4][0]>y) && (x<(x2a[4][0]-x1a[4][0])*(y-y1a[4][0])/(y2a[4][0]-y1a[4][0])+x1a[4][0])))
    c = !c;
  if( ((y1b[4][0]>y) != (y1a[4][0]>y) && (x<(x1a[4][0]-x1b[4][0])*(y-y1b[4][0])/(y1a[4][0]-y1b[4][0])+x1b[4][0])))
    c = !c;
  if( ((y2b[4][0]>y) != (y1b[4][0]>y) && (x<(x1b[4][0]-x2b[4][0])*(y-y2b[4][0])/(y1b[4][0]-y2b[4][0])+x2b[4][0])))
    c = !c;
  if( ((y2a[4][0]>y) != (y2b[4][0]>y) && (x<(x2b[4][0]-x2a[4][0])*(y-y2a[4][0])/(y2b[4][0]-y2a[4][0])+x2a[4][0])))
    c = !c;
  
  for(int a=1; a<4; a=a+2){
    if( ((y1a[4][a]>y) != (y2a[4][a]>y) && (x<(x2a[4][a]-x1a[4][a])*(y-y1a[4][a])/(y2a[4][a]-y1a[4][a])+x1a[4][a])))
      c = !c;
    if( ((y1b[4][a]>y) != (y1a[4][a]>y) && (x<(x1a[4][a]-x1b[4][a])*(y-y1b[4][a])/(y1a[4][a]-y1b[4][a])+x1b[4][a])))
      c = !c;
    if( ((y2b[4][a]>y) != (y1b[4][a]>y) && (x<(x1b[4][a]-x2b[4][a])*(y-y2b[4][a])/(y1b[4][a]-y2b[4][a])+x2b[4][a])))
      c = !c;
    if( ((y2a[4][a]>y) != (y2b[4][a]>y) && (x<(x2b[4][a]-x2a[4][a])*(y-y2a[4][a])/(y2b[4][a]-y2a[4][a])+x2a[4][a])))
      c = !c;
  }
  
  //special case: stripe 5 in radii 1, 2, 3, 4, 5, 6, 7 removed from R1_e
  for(int b=1; b<nRadii; b++){
    if( ((y1a[5][b]>y) != (y2a[5][b]>y) && (x<(x2a[5][b]-x1a[5][b])*(y-y1a[5][b])/(y2a[5][b]-y1a[5][b])+x1a[5][b])))
      c = !c;
    if( ((y1b[5][b]>y) != (y1a[5][b]>y) && (x<(x1a[5][b]-x1b[5][b])*(y-y1b[5][b])/(y1a[5][b]-y1b[5][b])+x1b[5][b])))
      c = !c;
    if( ((y2b[5][b]>y) != (y1b[5][b]>y) && (x<(x1b[5][b]-x2b[5][b])*(y-y2b[5][b])/(y1b[5][b]-y2b[5][b])+x2b[5][b])))
      c = !c;
    if( ((y2a[5][b]>y) != (y2b[5][b]>y) && (x<(x2b[5][b]-x2a[5][b])*(y-y2a[5][b])/(y2b[5][b]-y2a[5][b])+x2a[5][b])))
      c = !c;
  } 

  return c;
}

int searchR1(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x, double y){
  int c = 0;
  
  for(int j=0; j<nRadii; j++){
    for(int i=0; i<nStripes; i++){
      if( ((y1a[i][j]>y) != (y2a[i][j]>y) && (x<(x2a[i][j]-x1a[i][j])*(y-y1a[i][j])/(y2a[i][j]-y1a[i][j])+x1a[i][j])))
	c = !c;
      if( ((y1b[i][j]>y) != (y1a[i][j]>y) && (x<(x1a[i][j]-x1b[i][j])*(y-y1b[i][j])/(y1a[i][j]-y1b[i][j])+x1b[i][j])))
	c = !c;
      if( ((y2b[i][j]>y) != (y1b[i][j]>y) && (x<(x1b[i][j]-x2b[i][j])*(y-y2b[i][j])/(y1b[i][j]-y2b[i][j])+x2b[i][j])))
	c = !c;
      if( ((y2a[i][j]>y) != (y2b[i][j]>y) && (x<(x2b[i][j]-x2a[i][j])*(y-y2a[i][j])/(y2b[i][j]-y2a[i][j])+x2a[i][j])))
	c = !c;
    }
  }

  //special case: stripe 0 in radii 0, 2, 4, 6 removed from all module sections
  for(int k=0; k<nRadii; k=k+2){
    if( ((y1a[0][k]>y) != (y2a[0][k]>y) && (x<(x2a[0][k]-x1a[0][k])*(y-y1a[0][k])/(y2a[0][k]-y1a[0][k])+x1a[0][k])))
      c = !c;
    if( ((y1b[0][k]>y) != (y1a[0][k]>y) && (x<(x1a[0][k]-x1b[0][k])*(y-y1b[0][k])/(y1a[0][k]-y1b[0][k])+x1b[0][k])))
      c = !c;
    if( ((y2b[0][k]>y) != (y1b[0][k]>y) && (x<(x1b[0][k]-x2b[0][k])*(y-y2b[0][k])/(y1b[0][k]-y2b[0][k])+x2b[0][k])))
      c = !c;
    if( ((y2a[0][k]>y) != (y2b[0][k]>y) && (x<(x2b[0][k]-x2a[0][k])*(y-y2a[0][k])/(y2b[0][k]-y2a[0][k])+x2a[0][k])))
      c = !c;
  }

    //special case: stripe 5 in radii 0, 1, 3, 5, 7 removed from R1
  if( ((y1a[5][0]>y) != (y2a[5][0]>y) && (x<(x2a[5][0]-x1a[5][0])*(y-y1a[5][0])/(y2a[5][0]-y1a[5][0])+x1a[5][0])))
    c = !c;
  if( ((y1b[5][0]>y) != (y1a[5][0]>y) && (x<(x1a[5][0]-x1b[5][0])*(y-y1b[5][0])/(y1a[5][0]-y1b[5][0])+x1b[5][0])))
    c = !c;
  if( ((y2b[5][0]>y) != (y1b[5][0]>y) && (x<(x1b[5][0]-x2b[5][0])*(y-y2b[5][0])/(y1b[5][0]-y2b[5][0])+x2b[5][0])))
    c = !c;
  if( ((y2a[5][0]>y) != (y2b[5][0]>y) && (x<(x2b[5][0]-x2a[5][0])*(y-y2a[5][0])/(y2b[5][0]-y2a[5][0])+x2a[5][0])))
    c = !c;
  
  for(int b=1; b<nRadii; b=b+2){
    if( ((y1a[5][b]>y) != (y2a[5][b]>y) && (x<(x2a[5][b]-x1a[5][b])*(y-y1a[5][b])/(y2a[5][b]-y1a[5][b])+x1a[5][b])))
      c = !c;
    if( ((y1b[5][b]>y) != (y1a[5][b]>y) && (x<(x1a[5][b]-x1b[5][b])*(y-y1b[5][b])/(y1a[5][b]-y1b[5][b])+x1b[5][b])))
      c = !c;
    if( ((y2b[5][b]>y) != (y1b[5][b]>y) && (x<(x1b[5][b]-x2b[5][b])*(y-y2b[5][b])/(y1b[5][b]-y2b[5][b])+x2b[5][b])))
      c = !c;
    if( ((y2a[5][b]>y) != (y2b[5][b]>y) && (x<(x2b[5][b]-x2a[5][b])*(y-y2a[5][b])/(y2b[5][b]-y2a[5][b])+x2a[5][b])))
      c = !c;
  } 

  return c;
}

int searchR2(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x, double y){
  int c = 0;
  
  for(int j=0; j<nRadii; j++){
    for(int i=0; i<nStripes; i++){
      if( ((y1a[i][j]>y) != (y2a[i][j]>y) && (x<(x2a[i][j]-x1a[i][j])*(y-y1a[i][j])/(y2a[i][j]-y1a[i][j])+x1a[i][j])))
	c = !c;
      if( ((y1b[i][j]>y) != (y1a[i][j]>y) && (x<(x1a[i][j]-x1b[i][j])*(y-y1b[i][j])/(y1a[i][j]-y1b[i][j])+x1b[i][j])))
	c = !c;
      if( ((y2b[i][j]>y) != (y1b[i][j]>y) && (x<(x1b[i][j]-x2b[i][j])*(y-y2b[i][j])/(y1b[i][j]-y2b[i][j])+x2b[i][j])))
	c = !c;
      if( ((y2a[i][j]>y) != (y2b[i][j]>y) && (x<(x2b[i][j]-x2a[i][j])*(y-y2a[i][j])/(y2b[i][j]-y2a[i][j])+x2a[i][j])))
	c = !c;
    }
  }

  //special case: stripe 0 in radii 0, 2, 4, 6 removed from all module sections
  for(int k=0; k<nRadii; k=k+2){
    if( ((y1a[0][k]>y) != (y2a[0][k]>y) && (x<(x2a[0][k]-x1a[0][k])*(y-y1a[0][k])/(y2a[0][k]-y1a[0][k])+x1a[0][k])))
      c = !c;
    if( ((y1b[0][k]>y) != (y1a[0][k]>y) && (x<(x1a[0][k]-x1b[0][k])*(y-y1b[0][k])/(y1a[0][k]-y1b[0][k])+x1b[0][k])))
      c = !c;
    if( ((y2b[0][k]>y) != (y1b[0][k]>y) && (x<(x1b[0][k]-x2b[0][k])*(y-y2b[0][k])/(y1b[0][k]-y2b[0][k])+x2b[0][k])))
      c = !c;
    if( ((y2a[0][k]>y) != (y2b[0][k]>y) && (x<(x2b[0][k]-x2a[0][k])*(y-y2a[0][k])/(y2b[0][k]-y2a[0][k])+x2a[0][k])))
      c = !c;
  }

  //special case: stripe 7 in radii 0, 1, 3, 6 removed from R2
  if( ((y1a[7][0]>y) != (y2a[7][0]>y) && (x<(x2a[7][0]-x1a[7][0])*(y-y1a[7][0])/(y2a[7][0]-y1a[7][0])+x1a[7][0])))
    c = !c;
  if( ((y1b[7][0]>y) != (y1a[7][0]>y) && (x<(x1a[7][0]-x1b[7][0])*(y-y1b[7][0])/(y1a[7][0]-y1b[7][0])+x1b[7][0])))
    c = !c;
  if( ((y2b[7][0]>y) != (y1b[7][0]>y) && (x<(x1b[7][0]-x2b[7][0])*(y-y2b[7][0])/(y1b[7][0]-y2b[7][0])+x2b[7][0])))
    c = !c;
  if( ((y2a[7][0]>y) != (y2b[7][0]>y) && (x<(x2b[7][0]-x2a[7][0])*(y-y2a[7][0])/(y2b[7][0]-y2a[7][0])+x2a[7][0])))
    c = !c;
  
  for(int a=1; a<4; a=a+2){
    if( ((y1a[7][a]>y) != (y2a[7][a]>y) && (x<(x2a[7][a]-x1a[7][a])*(y-y1a[7][a])/(y2a[7][a]-y1a[7][a])+x1a[7][a])))
      c = !c;
    if( ((y1b[7][a]>y) != (y1a[7][a]>y) && (x<(x1a[7][a]-x1b[7][a])*(y-y1b[7][a])/(y1a[7][a]-y1b[7][a])+x1b[7][a])))
      c = !c;
    if( ((y2b[7][a]>y) != (y1b[7][a]>y) && (x<(x1b[7][a]-x2b[7][a])*(y-y2b[7][a])/(y1b[7][a]-y2b[7][a])+x2b[7][a])))
      c = !c;
    if( ((y2a[7][a]>y) != (y2b[7][a]>y) && (x<(x2b[7][a]-x2a[7][a])*(y-y2a[7][a])/(y2b[7][a]-y2a[7][a])+x2a[7][a])))
      c = !c;
  }
  
  if( ((y1a[7][6]>y) != (y2a[7][6]>y) && (x<(x2a[7][6]-x1a[7][6])*(y-y1a[7][6])/(y2a[7][6]-y1a[7][6])+x1a[7][6])))
    c = !c;
  if( ((y1b[7][6]>y) != (y1a[7][6]>y) && (x<(x1a[7][6]-x1b[7][6])*(y-y1b[7][6])/(y1a[7][6]-y1b[7][6])+x1b[7][6])))
    c = !c;
  if( ((y2b[7][6]>y) != (y1b[7][6]>y) && (x<(x1b[7][6]-x2b[7][6])*(y-y2b[7][6])/(y1b[7][6]-y2b[7][6])+x2b[7][6])))
    c = !c;
  if( ((y2a[7][6]>y) != (y2b[7][6]>y) && (x<(x2b[7][6]-x2a[7][6])*(y-y2a[7][6])/(y2b[7][6]-y2a[7][6])+x2a[7][6])))
    c = !c;
  
  return c;
}

int searchR3(int nStripes, double x1a[][nRadii], double x1b[][nRadii], double x2a[][nRadii], double x2b[][nRadii], double y1a[][nRadii], double y1b[][nRadii], double y2a[][nRadii], double y2b[][nRadii], double x, double y){
  int c = 0;
  
  for(int j=0; j<nRadii; j++){
    for(int i=0; i<nStripes; i++){
      if( ((y1a[i][j]>y) != (y2a[i][j]>y) && (x<(x2a[i][j]-x1a[i][j])*(y-y1a[i][j])/(y2a[i][j]-y1a[i][j])+x1a[i][j])))
	c = !c;
      if( ((y1b[i][j]>y) != (y1a[i][j]>y) && (x<(x1a[i][j]-x1b[i][j])*(y-y1b[i][j])/(y1a[i][j]-y1b[i][j])+x1b[i][j])))
	c = !c;
      if( ((y2b[i][j]>y) != (y1b[i][j]>y) && (x<(x1b[i][j]-x2b[i][j])*(y-y2b[i][j])/(y1b[i][j]-y2b[i][j])+x2b[i][j])))
	c = !c;
      if( ((y2a[i][j]>y) != (y2b[i][j]>y) && (x<(x2b[i][j]-x2a[i][j])*(y-y2a[i][j])/(y2b[i][j]-y2a[i][j])+x2a[i][j])))
	c = !c;
    }
  }

  //special case: stripe 0 in radii 0, 2, 4, 6 removed from all module sections
  for(int k=0; k<nRadii; k=k+2){
    if( ((y1a[0][k]>y) != (y2a[0][k]>y) && (x<(x2a[0][k]-x1a[0][k])*(y-y1a[0][k])/(y2a[0][k]-y1a[0][k])+x1a[0][k])))
      c = !c;
    if( ((y1b[0][k]>y) != (y1a[0][k]>y) && (x<(x1a[0][k]-x1b[0][k])*(y-y1b[0][k])/(y1a[0][k]-y1b[0][k])+x1b[0][k])))
      c = !c;
    if( ((y2b[0][k]>y) != (y1b[0][k]>y) && (x<(x1b[0][k]-x2b[0][k])*(y-y2b[0][k])/(y1b[0][k]-y2b[0][k])+x2b[0][k])))
      c = !c;
    if( ((y2a[0][k]>y) != (y2b[0][k]>y) && (x<(x2b[0][k]-x2a[0][k])*(y-y2a[0][k])/(y2b[0][k]-y2a[0][k])+x2a[0][k])))
      c = !c;
  }
  
  //special case: stripe 10 in radius 1 removed from R3
  if( ((y1a[10][1]>y) != (y2a[10][1]>y) && (x<(x2a[10][1]-x1a[10][1])*(y-y1a[10][1])/(y2a[10][1]-y1a[10][1])+x1a[10][1])))
    c = !c;
  if( ((y1b[10][1]>y) != (y1a[10][1]>y) && (x<(x1a[10][1]-x1b[10][1])*(y-y1b[10][1])/(y1a[10][1]-y1b[10][1])+x1b[10][1])))
    c = !c;
  if( ((y2b[10][1]>y) != (y1b[10][1]>y) && (x<(x1b[10][1]-x2b[10][1])*(y-y2b[10][1])/(y1b[10][1]-y2b[10][1])+x2b[10][1])))
    c = !c;
  if( ((y2a[10][1]>y) != (y2b[10][1]>y) && (x<(x2b[10][1]-x2a[10][1])*(y-y2a[10][1])/(y2b[10][1]-y2a[10][1])+x2a[10][1])))
    c = !c;

  //special case: stripe 11 in radii 0, 1, 2, 3, 4, 5, 7 removed from R3
  for(int a=0; a<6; a++){
    if( ((y1a[11][a]>y) != (y2a[11][a]>y) && (x<(x2a[11][a]-x1a[11][a])*(y-y1a[11][a])/(y2a[11][a]-y1a[11][a])+x1a[11][a])))
      c = !c;
    if( ((y1b[11][a]>y) != (y1a[11][a]>y) && (x<(x1a[11][a]-x1b[11][a])*(y-y1b[11][a])/(y1a[11][a]-y1b[11][a])+x1b[11][a])))
      c = !c;
    if( ((y2b[11][a]>y) != (y1b[11][a]>y) && (x<(x1b[11][a]-x2b[11][a])*(y-y2b[11][a])/(y1b[11][a]-y2b[11][a])+x2b[11][a])))
      c = !c;
    if( ((y2a[11][a]>y) != (y2b[11][a]>y) && (x<(x2b[11][a]-x2a[11][a])*(y-y2a[11][a])/(y2b[11][a]-y2a[11][a])+x2a[11][a])))
      c = !c;
  }

  if( ((y1a[11][7]>y) != (y2a[11][7]>y) && (x<(x2a[11][7]-x1a[11][7])*(y-y1a[11][7])/(y2a[11][7]-y1a[11][7])+x1a[11][7])))
    c = !c;
  if( ((y1b[11][7]>y) != (y1a[11][7]>y) && (x<(x1a[11][7]-x1b[11][7])*(y-y1b[11][7])/(y1a[11][7]-y1b[11][7])+x1b[11][7])))
    c = !c;
  if( ((y2b[11][7]>y) != (y1b[11][7]>y) && (x<(x1b[11][7]-x2b[11][7])*(y-y2b[11][7])/(y1b[11][7]-y2b[11][7])+x2b[11][7])))
    c = !c;
  if( ((y2a[11][7]>y) != (y2b[11][7]>y) && (x<(x2b[11][7]-x2a[11][7])*(y-y2a[11][7])/(y2b[11][7]-y2a[11][7])+x2a[11][7])))
    c = !c;
    
  return c;
}
