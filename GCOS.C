// Author: Lukas Scherne, December 2022
/*******************************************************************************
 * Simulation of a MINI MACHETE by ROBAST. - Proposal a)
 ******************************************************************************/
#include "TCanvas.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TH2.h"
#include "TH1.h"
#include "TThread.h"

#include "AFocalSurface.h"
#include "AGeoUtil.h"
#include "AMirror.h"
#include "AObscuration.h"
#include "AOpticsManager.h"
#include "ARayShooter.h"
#include "AGeoAsphericDisk.h"
#include "ALens.h"
#include "ASchottFormula.h"
#include "AGlassCatalog.h"

#include <iostream>
#include <fstream>

// define useful units
const double cm = AOpticsManager::cm();
const double mm = AOpticsManager::mm();
const double um = AOpticsManager::um();
const double nm = AOpticsManager::nm();
const double m = AOpticsManager::m();

//global variables - mirror and camera parameters
double lambda = 380 * nm;
const double kCameraDx = 4.73*m;     //camera diameter in y
const double kCameraDy = 0.93*m;        //camera diameter in x

const double kCameraOR =  4.634*m;        //camera outer Radius
const double kCameraIR = kCameraOR-0.001*m;          //camera inner Radius

const double kMirrorR = 9*m;    // the radius of curvature
const double kMirrorDx = 14*m;    //  diameter of the Mirror in x
const double kMirrorDy = 6.5 * m;    //  diameter of the Mirror in y
const double kMirrorT = 0.001 * mm;  // mirror thickness

const double pixelwidth = 0.6; //deg

//define some functions
void AddMirror(AOpticalComponent* opt);
void AddSegMirror(AOpticalComponent* opt);
void AddCamera(AOpticalComponent* opt);
void DrawAxes(AOpticalComponent* opt);
void RayTrace(AOpticsManager* manager, TCanvas* can3D);

void GCOS() { //Mainfunction - builds the world
  TThread::Initialize();  // call this first when you use the multi-thread mode

  AOpticsManager* manager = new AOpticsManager("manager", "MACHETE-FD");
  manager->DisableFresnelReflection(kTRUE); //

  // Make the OpenGL objects more smooth
  manager->SetNsegments(100);

  // Make the world of 40-m cube
  TGeoBBox* boxWorld = new TGeoBBox("boxWorld", 50 * m, 50 * m, 50 * m);
  AOpticalComponent* world = new AOpticalComponent("world", boxWorld);
  manager->SetTopVolume(world);

  TCanvas* can = new TCanvas("can3D", "can3D", 1000, 1000);

  //Fill the world
  AddCamera(world);
  AddMirror(world);
  //DrawAxes(world);
  manager->CloseGeometry();  // finalize the geometry construction

  #if ROOT_VERSION_CODE < ROOT_VERSION(6, 2, 0)
  manager->SetMultiThread(kTRUE);  // enable multi threading
  #endif
  manager->SetMaxThreads(4);

  world->Draw("ogl"); //Display Simulation

  RayTrace(manager, can);
}

void RayTrace(AOpticsManager* manager, TCanvas* can3D) { // raytracing, pointspread histogramms and cummulativ plots

  ofstream file;
  ofstream file2;
  ofstream file3;

  //result files
  file.open("photondata.txt");  //contains postion and direction information
  file2.open("resultingdata.txt"); //contains analyses results as r95, Aeff and so on
  file3.open("cummulativdata.txt"); //contains data to plot cumulative distributions


  file << "#Number of simulated photons: 50000" << "\n";
  file << "#Simulationdata i, || x, y, z, t, dx, dy, dz Postion of all registered photons - dx dy dz is the direction" << "\n";
  file2 << "#Number of simulated photons: 50000" << "\n";
  file2<< "#Aeff, incident angle [deg], acceptance angle[deg], acceptance radius [m], spotszize, r95 [mm], number focused photons, number of stopped photons" << "\n";
  file3 <<"#number of bins to be filled, binnumber, bincontent" << "\n";


  const int kNdeg = 31; //number of loops, 6 is minimum
  const int histsinx = kNdeg/2 ; // # of histogramms in x
  const int histsiny = 3 ;// # of histogramms in y

  //varibales to find max of hists
  double maxVal[kNdeg] ;
  int xBinMax[kNdeg]  ;
  int yBinMax[kNdeg] ;

  const double delta = 3*cm;
  TH2D* h2[kNdeg];
  TH1D* h1[kNdeg];
  gStyle->SetPalette(kBird);
  TColor::InvertPalette();
  TGraph* graph = new TGraph();
  TCanvas* can = new TCanvas("can", "can", 900, 600); //canvas for PSF
  TCanvas* can1D = new TCanvas("can1D", "can1D", 900, 600); //canvas for 1D hists

  can->Divide(histsinx, histsiny, 1e-10, 1e-10);
  can1D->Divide(histsinx, histsiny, 1e-10, 1e-10);

  double sum_focused_array[kNdeg];
  double Alpha[kNdeg];
  double Aeff[kNdeg];

  for (int i = 0 ; i < kNdeg; i++) {//for loop to fill histogramms
    double r_inc = (2.45)*m; // radius of incident ray field - acceptance angle of one pmt pixel
    double deg = i;

    double deg2 = 5;

    double accang = 2*TMath::ASin((r_inc)/(2*kCameraOR))*TMath::RadToDeg(); // opening angle of the light cone each pixel accepts - given trough the aperture


    //define direction of incident photon-field
    TGeoTranslation raytr("raytr", 0, 0, 0);

    TGeoRotation *rayrot = new TGeoRotation();
    rayrot->SetAngles(0,deg2,0);

    TVector3 dir;
    dir.SetMagThetaPhi(1, (deg)*TMath::DegToRad(), (0)*TMath::DegToRad());


    const int nBinsX = 5000; // number of bin in 2D Hists
    const int nBinsY = 5000; // number of bin in 2D Hists
    const double N_inc = 50000; // number of simulated rays

    h2[i] = new TH2D(Form("h%d", i),Form("#it{#theta} = %3.1f#circ;x (mm); y (mm)", deg), nBinsX, (-kCameraDx/2-delta)/mm ,  (kCameraDx/2+delta)/mm , nBinsY, -(kCameraDy/2-delta)/mm , (kCameraDy/2+delta)/mm); // use size of camera, numbers in mm
    h1[i] = new TH1D(Form("H%d",i),Form("#it{#theta} = %3.1f#circ; Distance from spot centre (mm); fraction of encircled Energy", deg), 1000, 0, 30);     // hier noch anpassen

    ARayArray* array = ARayShooter::RandomCircle(lambda, r_inc, N_inc, rayrot, &raytr, &dir);//shoots rays in given dir with given trans raytr
    // r = 6 m because that is the effectiv apperture

    manager->TraceNonSequential(*array);

    TObjArray* focused = array -> GetFocused();
    TObjArray* stopped = array -> GetStopped();

    //ratio of incident rays that are focused to the stopped ones
    double sum_focused = 0;
    double sum_stopped = 0;


//GetFocused
    for (Int_t k = 0; k <= focused->GetLast(); k++) { //loop over all rays that hitted the camera
      ARay* ray = (ARay*)(*focused)[k];

      Double_t p[4];
      ray->GetLastPoint(p);

      Double_t d[3];
      ray -> GetDirection(d);
      //x,y,z of direction vector to calculate angular distribution of the incoming rays


      h2[i]->Fill((p[0]/mm) , (p[1])/mm);
      //save data in txt file
      file <<  deg << ",";
      file << d[0] << ",";
      file << d[1] << ",";
      file << d[2] << ",";

      file << p[0]/mm << ",";
      file << p[1]/mm << ",";
      file << p[2]/mm << ",";
      file << p[3]/mm;
      file << "\n";
      //füllt das histogramm

      sum_focused++;
      if (deg == 0 && k <= 100) {
        TPolyLine3D* pol = ray->MakePolyLine3D();
        pol->SetLineColor(3);
        pol->SetLineWidth(2);
        can3D->cd();
        pol->Draw("p");
      }  // if
      if ( i  == kNdeg-1 && k <= 100) {
        TPolyLine3D* pol = ray->MakePolyLine3D();
        pol->SetLineColor(4);
        pol->SetLineWidth(2);
        can3D->cd();
        pol->Draw("p");
      }
    }    // k

    sum_focused_array[i] =  sum_focused;
//Get number of Stopped arrays
    for (Int_t k = 0 ; k <= stopped -> GetLast(); k++) { //loop over all rays where stopped by the camera back
      ARay* ray = (ARay*)(*stopped)[k];
      sum_stopped++;
    }

	//find max in 2D hists
	for(int iX = 0; iX < h2[i]->GetXaxis()->GetNbins(); ++iX){
		for(int iY = 0; iY < h2[i]->GetYaxis()->GetNbins(); ++iY){
			const double cont = h2[i]->GetBinContent(iX+1, iY+1);
				if(cont > maxVal[i]){
					maxVal[i] = cont;
					xBinMax[i] = iX+1;
					yBinMax[i] = iY+1;
				}
		}
	}

    //for debugging
    //cout << "max from 1D hist " << i << " x " << xBinMax[i] << " y " << yBinMax[i] << endl;

    //define centre and plot
    double xCenter[kNdeg];
    xCenter[i] = h2[i]->GetXaxis()->GetBinCenter(xBinMax[i]);

    double yCenter[kNdeg];
    yCenter[i] = h2[i]->GetYaxis()->GetBinCenter(yBinMax[i]);

    //some variables
    int dimx = h2[i]->GetXaxis()->GetNbins();
    int dimy = h2[i]->GetYaxis()->GetNbins();
    double x[dimx];
    double y[dimy];
    double r[dimx];


    double cont[kNdeg] ;
    double sum[kNdeg] ;



    //calculate radial distribution for each hist
  for (int iX = 0; iX < h2[i]->GetXaxis()->GetNbins(); ++iX) {
    x[i] = h2[i]->GetXaxis()->GetBinCenter(iX+1);
    for (int iY = 0; iY < h2[i]->GetYaxis()->GetNbins(); ++iY) {
      y[i] = h2[i]->GetYaxis()->GetBinCenter(iY+1);
      r[i] = sqrt(pow(x[i]-xCenter[i], 2) + pow(y[i]-yCenter[i], 2));
      // cont : contains bin content of each bin
      cont[i] = h2[i]->GetBinContent(iX+1, iY+1);
      //sum : totale number of entrys
      sum[i] += cont[i];
      h1[i]->Fill(r[i], cont[i]); // filling new 1D hist, radial distance with content of each bin of 2D Hist
    }
  }

  //calculate cumulative plot
  double cumu[kNdeg];

  double r95[kNdeg];


  for(int iX = 0; iX < h1[i]->GetXaxis()->GetNbins(); iX++){
   //cont : contains bin content of each bin in the 1D Hist
   cont[i] = h1[i]->GetBinContent(iX+1);
   // cumu : total number of bin entrys
   cumu[i] += cont[i];
   //find r95
   if(cumu[i]/sum[i] < 0.95 and cumu[i]/sum[i] > 0.93){
        //iX : binnumber * bins/mm in cummulativ plot
        // saves the last value smaller than 0.95
       r95[i] = iX * 30./1000.;
   }

   //fill each bin with #entrys in h1D / #entrys is h2D
   h1[i] -> SetBinContent(iX+1, cumu[i]/sum[i]);
   //write cummulativ data to the files
   file3 << h1[i]->GetXaxis()->GetNbins() << ",";
   file3 << iX+1 << ",";
   file3 << cumu[i]/sum[i];
   file3 << "\n";


 }

    can->cd(i + 1); // histogram
    h2[i]->Draw("COLZ");
    const double ddelta = delta/mm;
    h2[i]->GetXaxis()->SetRangeUser(xCenter[i]-ddelta,xCenter[i]+ddelta);
    h2[i]->GetYaxis()->SetRangeUser(yCenter[i]-ddelta,yCenter[i]+ddelta);
    can->Update();

    can1D->cd(i+1);
    h1[i]->Draw("CHIST");
    //draw lines where 95% of the light is encircled, as in other paper
    TLine *xa = new TLine(0,0.95,18,0.95);
    TLine *ya = new TLine(r95[i],0, r95[i], 1);
    ya -> Draw();
    xa -> Draw();
    can1D->Update();

    delete array;
    double aeff = r_inc/m*r_inc/m*TMath::Pi() * (sum_focused/50000);
    Aeff[i] = aeff;
    Alpha[i] = accang;


    double r95indeg  = 2*TMath::ASin((r95[i]*mm)/(2*kCameraOR))*TMath::RadToDeg();
    cout << "----------------------------------------------------------------------------------------------" << endl;
    cout << "iteration " << i << endl;
    cout << "Spotcentre, in mm on the camera (x,y) " << xCenter[i] << "," << yCenter[i] << endl;
    cout << "Radius of generated photonfield " << r_inc/m << " m " << endl;
    cout << "Aperture in deg " << accang << " deg " << endl;
    cout << "ratio of incident to focused photons at " << deg << " deg " << (sum_focused/50000)*100 << "%" << endl;
    cout << "ratio of incident to stopped photons at " << deg << " deg " << (sum_stopped/50000)*100 << "%"<< endl;
    cout << "ratio of focused to stopped photons at " << deg << " deg " << sum_focused/sum_stopped << endl;
    cout << "r95 : "<< r95[i] << " mm" << " in deg: " << r95indeg << endl;
    cout << "Effective imaging area at " << deg << " deg" << Aeff[i] << "m²" << endl;
    cout << "Propotionality factor of Signal to Noise Ration at " << deg << " deg: " << TMath::Sqrt((aeff/(pixelwidth*pixelwidth))) << endl;
    cout << "----------------------------------------------------------------------------------------------" << endl;

    file2 << aeff << ","; //effective aperture in m²
    file2 << deg << ","; // angle of incidence in deg
    file2 << accang << ","; // angle of acceptance in deg
    file2 << r_inc/m << ","; // radius of generated photon field
    file2 << r95[i] << ","; // radius where 95% of light is encircled in mm
    file2 << sum_focused << ","; //number of focused photons
    file2 << sum_stopped; // number of stopped photons
    file2 << "\n";
  }  // i

}

void DrawAxes(AOpticalComponent* opt){ //Draws a coordinate system and a aperture, for debugging
  //Display koordinate-system and a "wall"
  TGeoBBox* xaxis = new TGeoBBox("xa", 2.5*m,0.01*mm, 0.01*mm);
  TGeoBBox* yaxis = new TGeoBBox("ya", 0.01*mm,2.5*m, 0.01*mm);
  TGeoBBox* zaxis = new TGeoBBox("za", 0.01*mm,0.01*mm, 2.5*m);

  //Aperture
  const double z = 0;
  TGeoTubeSeg* Aperture = new TGeoTubeSeg("T1", 1.1*m, 2*m, 0.01*mm, 0, 360);
  TGeoTubeSeg* Aperture2 = new TGeoTubeSeg("T2", 0*m, 2.19/2*m, 0.01*mm, 0,360);
  TGeoTubeSeg* Aperture3 = new TGeoTubeSeg("T3", 1.15*m, 1.1*m, 1.15*mm, 82.5, 82.5+15);

  AObscuration* Wall = new AObscuration("Wall", Aperture);
  AObscuration* Wall2 = new AObscuration("Wall2", Aperture2);
  AObscuration* Wall3 = new AObscuration("Wall3", Aperture3);

  //coordinate axis, you can see them in wire mode
  AObscuration* xa = new AObscuration("X", xaxis);
  AObscuration* ya = new AObscuration("Y", yaxis);
  AObscuration* za = new AObscuration("Z", zaxis);


  //opt -> AddNode(Wall, 1, new TGeoTranslation(0,0,0*mm));
  //opt -> AddNode(Wall2, 1, new TGeoTranslation(0,0,0*mm));
  //opt -> AddNode(Wall3, 1, new TGeoTranslation(0,0,-5.1/2*mm));


  opt -> AddNode(xa,1, new TGeoTranslation(0,0,0));
  opt -> AddNode(ya, 1, new TGeoTranslation(0,0,0));
  opt -> AddNode(za, 1, new TGeoTranslation(0,0,0));

}

void AddCamera(AOpticalComponent* opt) { // adds focal plane (camera) and camera shadow to the world
  // Make a focal plane and display camera
  //camera paramters

  const double alpha = 30;

  const double dz = 0*mm; //moves focal surface along the optical axis (z axis)

  TGeoSphere* sphericalCamera = new TGeoSphere("S2", kCameraIR , kCameraOR, 0, 2*alpha, -180, 180);
  TGeoBBox* rect = new TGeoBBox("CB",kCameraDx/2, kCameraDy/2, 20*m);

  TGeoTranslation* transZ2 = new TGeoTranslation("transZ2", 0,0, 0*m);
  transZ2->RegisterYourself();

  TGeoCompositeShape *cs2 = new TGeoCompositeShape("Cam" , "CB:transZ2*S2");

  AFocalSurface* focalPlane = new AFocalSurface("focalPlane", cs2);
  opt->AddNode(focalPlane, 1, new TGeoTranslation(0*m, 0*m, dz)); //

  //make a camera box as obscuration
  //Back of the Camera
  TGeoSphere* SphericalCameraBox = new TGeoSphere("SCB", kCameraIR+0.1*mm,kCameraOR+0.1*mm,0,2*alpha,-180,180);
  TGeoBBox* rectforCam = new TGeoBBox("RFC", kCameraDx/2+1*mm, kCameraDy/2+1*mm, 20*m);

  TGeoTranslation* transZ3 = new TGeoTranslation("transZ3", 0,0, 0*m);
  transZ3->RegisterYourself();

  TGeoCompositeShape *camBox = new TGeoCompositeShape("Camera Box", "RFC:transZ3*SCB");
  const double delta = kCameraOR-kCameraIR;

  AObscuration* camBoxObs = new AObscuration("Camera Box", camBox);
  opt -> AddNode(camBoxObs, 1, new TGeoTranslation(0,0, dz-0.005*m));

}

void AddMirror(AOpticalComponent* opt){ //adds segmentation to the world

  //Add Mirror as Composite shapes

 	TGeoSphere* mirSphereid = new TGeoSphere("Sid", kMirrorR , (kMirrorR+kMirrorT), 0, 2*30, 0, 360);

  TGeoBBox* boxid = new TGeoBBox("Bid", (kMirrorDx)/2, (kMirrorDy)/2, (kMirrorDy)/2); // x,y,z sind halbe längen
  TGeoTranslation* transZ1 = new TGeoTranslation("transZ1", 0,0,kMirrorDy);
  transZ1->RegisterYourself();

  //Create composite, as mirror is a rectangle with a spherical surface
  TGeoCompositeShape *recmirrid = new TGeoCompositeShape("recmirrid" , "Bid:transZ1*Sid");
  AMirror* mirrorid = new AMirror("mirrorid", recmirrid);

  mirrorid->SetReflectance(1);
  opt->AddNode(mirrorid, 1, new TGeoTranslation(0*m, 0*m, 0*mm));

}
