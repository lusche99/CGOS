// Author: Lukas Scherne, November 2022
/*******************************************************************************
 * Simulation of the Pierre Auger-FD with ROBAST.
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

//global variables - mirror and camera
double lambda = 400 * nm;
const double kCameraDy = 0.93*m;     //camera diameter in x
const double kCameraDx = 0.86*m;     //camera diameter in y

const double kCameraIR = 1742.4199*mm;    //camera inner Radius
const double kCameraOR =  1742.42*mm;        //camera outer Radius

const double kMirrorR = 3.4 * m;    // the radius of curvature
const double kMirrorD = 3.6 * m;    //  diameter of the Mirror
const double kMirrorT = 0.001 * mm;  // mirror thickness

const double pixelwidth = 1.5; //deg

//get refractic index for the corrector ring
AGlassCatalog schott("../misc/schottzemax-20180601.agf");
auto bk7 = schott.GetRefractiveIndex("N-BK7");

//define some functions
void AddMirror(AOpticalComponent* opt);
void AddCamera(AOpticalComponent* opt);
void AddLenses(AOpticalComponent* opt);
void DrawAxes(AOpticalComponent* opt);
void RayTrace(AOpticsManager* manager, TCanvas* can3D);

void PierreAuger() { //Mainfunction - builds the world
  TThread::Initialize();  // call this first when you use the multi-thread mode

  AOpticsManager* manager = new AOpticsManager("manager", "PierreAuger-FD");
  manager->DisableFresnelReflection(kTRUE); //

  // Make the OpenGL objects more smooth
  manager->SetNsegments(100);

  // Make the world of 40-m cube
  TGeoBBox* boxWorld = new TGeoBBox("boxWorld", 20 * m, 20 * m, 20 * m);
  AOpticalComponent* world = new AOpticalComponent("world", boxWorld);
  manager->SetTopVolume(world);

  TCanvas* can = new TCanvas("can3D", "can3D", 1000, 1000);

  //Fill the world
  AddCamera(world);
  //AddLenses(world);
  AddMirror(world);
  DrawAxes(world);
  manager->CloseGeometry();  // finalize the geometry construction

  #if ROOT_VERSION_CODE < ROOT_VERSION(6, 2, 0)
  manager->SetMultiThread(kTRUE);  // enable multi threading
  #endif
  manager->SetMaxThreads(4);

  world->Draw("ogl"); //Display Simulation

  RayTrace(manager, can);
}

void RayTrace(AOpticsManager* manager, TCanvas* can3D) { // raytracing and pointspread hist as well as cummulativ plots

  ofstream file;
  ofstream file2;
  ofstream file3;

  //save information to files
  file.open("photondata_PA_ncr.txt"); //contains postion and direction information
  file2.open("resultingdata_PA_ncr.txt"); //contains analyses results as r95, Aeff and so on
  file3.open("cummulativdata_PA_ncr.txt"); //contains data to plot cumulative distributions

  file << "#Number of simulated photons: 50000" << "\n";
  file << "#Simulationdata i, || x, y, z, t, dx, dy, dz Postion of all registered photons" << "\n";
  file2 << "#Number of simulated photons: 50000" << "\n";
  file2<< "#Aeff, incident angle [deg], acceptance angle[deg], acceptance radius [m], spotszize, r95 [mm], number focused photons, number of stopped photons" << "\n";
  file3 <<"#number of bins to be filled, binnumber, bincontent" << "\n";

  const int kNdeg = 21;
  const int histsinx = kNdeg/2; // # of histogramms in x
  const int histsiny = 3; // # of histogramms in y

  //varibales to find max of hists
  double maxVal[kNdeg] ;
  int xBinMax[kNdeg] ;
  int yBinMax[kNdeg] ;

  const double delta = 3*cm;
  TH2D* h2[kNdeg];
  TH1D* h1[kNdeg]; 
  gStyle->SetPalette(kBird);
  TColor::InvertPalette();

  TCanvas* can = new TCanvas("can", "can", 900, 600); //canvas for PSF
  TCanvas* can1D = new TCanvas("can1D", "can1D", 900, 600); //canvas for 1D hists

  can->Divide(histsinx, histsiny, 1e-10, 1e-10);
  can1D->Divide(histsinx, histsiny, 1e-10, 1e-10);


  for (int i = 0 ; i < kNdeg; i++) {//main for loop to fill histogramms
    double deg = i;
    if(i == 5){deg = 20;}

    double r_inc = .85*m;

    TGeoTranslation raytr("raytr", 0, 0, -1*m);
    double theta = 0;
    TGeoRotation rayrot("rayrot", 0, 0, 0);
    TVector3 dir;

    dir.SetMagThetaPhi(1, (deg)*TMath::DegToRad(), (90-45)*TMath::DegToRad());


    const int nBinsX = 3000; // number of bin in 2D Hists
    const int nBinsY = 3000; // number of bin in 2D Hists

    h2[i] = new TH2D(Form("h%d", i),Form("#it{#theta} = %3.1f#circ;x (mm); y (mm)", deg), nBinsX, (-kCameraDx/2-delta)/mm ,  (kCameraDx/2+delta)/mm , nBinsY, -(kCameraDy/2-delta)/mm , (kCameraDy/2+delta)/mm); // use size of camera, numbers in mm
    h1[i] = new TH1D(Form("H%d",i),Form("#it{#theta} = %3.1f#circ; Distance from spot centre (mm); fraction of encircled light", deg), 3000, 0, 20);

    ARayArray* array = ARayShooter::RandomCircle(lambda, r_inc, 50000, &rayrot, &raytr, &dir);//shoots rays in given dir with gicen trans raytr

    manager->TraceNonSequential(*array);

    TObjArray* focused = array -> GetFocused();
    TObjArray* stopped = array -> GetStopped();

    double sum_focused = 0;
    double sum_stopped = 0;


    //GetFocused
    for (Int_t k = 0; k <= focused->GetLast(); k++) { //loop over all rays that hitted the camera
      ARay* ray = (ARay*)(*focused)[k];

      Double_t p[4];
      ray->GetLastPoint(p);

      Double_t d[3];
      ray -> GetDirection(d);

      sum_focused++;
      //fill histogram
      h2[i]->Fill((p[0]/mm) , (p[1]/mm));

      file << deg << ",";
      file << d[0] << ",";
      file << d[1] << ",";
      file << d[2] << ",";

      file << p[0]/mm << ",";
      file << p[1]/mm << ",";
      file << p[2]/mm << ",";
      file << p[3]/mm;
      file << "\n";

      if (i == 0 && k <= 50 ) {
        TPolyLine3D* pol = ray->MakePolyLine3D();
        pol->SetLineColor(4);
        pol->SetLineWidth(2);
        can3D->cd();
        pol->Draw("p");
      }  // if
      if (i == kNdeg-1 && k <= 50 ) {
        TPolyLine3D* pol = ray->MakePolyLine3D();
        pol->SetLineColor(3);
        pol->SetLineWidth(2);
        can3D->cd();
        pol->Draw("p");
      }  // if
    }    // k

    //Get number of Stopped rays
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
    double sum[kNdeg]  ;

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
    double r95;

    for(int iX = 0; iX < h1[i]->GetXaxis()->GetNbins(); iX++){
    	//cont : contains bin content of each bin in the 1D Hist
    	cont[i] = h1[i]->GetBinContent(iX+1);
    	// cumu : totale number of bin entrys
    	cumu[i] += cont[i];
    	//fill each bin with #entrys in h1D / #entrys is h2D
    	h1[i] -> SetBinContent(iX+1, cumu[i]/sum[i]);

      //write cummulativ data to the files
      file3 << h1[i]->GetXaxis()->GetNbins() << ",";
      file3 << iX+1 << ",";
      file3 << cumu[i]/sum[i];
      file3 << "\n";

      if(cumu[i]/sum[i] < .95 and cumu[i]/sum[i] > 0.90){
           //iX : binnumber * bins/mm in cummulativ plot
           // saves the last value smaller than 0.95
          r95 = iX * 0.006666667;
          //cout << r95 << endl;
      }
    }

    //Draw hisogramms
    can->cd(i + 1);
    h2[i]->Draw(""); //change draw option of PSF here
    const double ddelta = (delta/3)/mm;
    h2[i]->GetXaxis()->SetRangeUser(xCenter[i]-ddelta,xCenter[i]+ddelta);
    h2[i]->GetYaxis()->SetRangeUser(yCenter[i]-ddelta,yCenter[i]+ddelta);
    can->Update();


    //display cummulativ plots
    can1D->cd(i+1);
    h1[i]->Draw("CHIST");
    //draw lines where 95% of the light is encircled
    TLine *xa = new TLine(0,0.95,10,0.95);
    TLine *ya = new TLine(7.5,0, 7.5, 1);
    ya -> Draw();
    xa -> Draw();
    can1D->Update();


  //calculate and print out useful information
    double aeff = r_inc/m*r_inc/m*TMath::Pi() * (sum_focused/50000);
    double accang = 2*TMath::ASin((r_inc)/(2*kCameraOR))*TMath::RadToDeg();
    double r95indeg  = 2*TMath::ASin((r95*mm)/(2*kCameraOR))*TMath::RadToDeg();

    cout << "----------------------------------------------------------------------------------------------" << endl;
    cout << "iteration " << i << endl;
    cout << "Radius of generated photonfield " << r_inc/m << " m " << endl;
    cout << "Angular realm of the mirror accepted by the pixel" << accang << " deg " << endl;
    cout << "ratio of incident to focused photons at " << deg << " deg " << (sum_focused/50000)*100 << "%" << endl;
    cout << "ratio of incident to stopped photons at " << deg << " deg " << (sum_stopped/50000)*100 << "%"<< endl;
    cout << "ratio of focused to stopped photons at " << deg << " deg " << sum_focused/sum_stopped << endl;
    cout << "r95 : "<< r95 << " mm " << " in deg: " << r95indeg << endl;
    cout << "Effective imaging area at " << deg << " deg incidence: " << aeff << "m²"<< endl;
    cout << "Propotionality factor of >Signal to Noise Ration< at " << deg << " deg: " << TMath::Sqrt((aeff/(pixelwidth*pixelwidth))) << endl;
    cout << "----------------------------------------------------------------------------------------------" << endl;

    file2 << aeff << ",";
    file2 << deg << ",";
    file2 << accang << ","; // angle of acceptance in deg
    file2 << r_inc/m << ","; // radius of generated photon field
    file2 << r95 << ",";
    file2 << sum_focused<<",";
    file2 << sum_stopped;
    file2 <<"\n";


    //can->cd();
    delete array;
  }  // i
}

void DrawAxes(AOpticalComponent* opt){ //Draws a coordinate system and a aperture, for debugging
  //Display koordinate-system and a "wall"
  TGeoBBox* xaxis = new TGeoBBox("xa", 2.5*m,0.01*mm, 0.01*mm);
  TGeoBBox* yaxis = new TGeoBBox("ya", 0.01*mm,2.5*m, 0.01*mm);
  TGeoBBox* zaxis = new TGeoBBox("za", 0.01*mm,0.01*mm, 2.5*m);

  //Aperture
  const double z = 0;
  TGeoTubeSeg* Aperture = new TGeoTubeSeg("T1", 1.1*m, 2.2*m, 0.01*mm, 0, 360);
  TGeoTubeSeg* Aperture2 = new TGeoTubeSeg("T2", 0.85*m, 2.2*m, 0.01*mm, 0,360);
  TGeoTubeSeg* Aperture3 = new TGeoTubeSeg("T3", 1.15*m, 1.1*m, 1.15*mm, 82.5, 82.5+15);

  AObscuration* Wall = new AObscuration("Wall", Aperture);
  AObscuration* Wall2 = new AObscuration("Wall2", Aperture2);
  AObscuration* Wall3 = new AObscuration("Wall3", Aperture3);

  //coordinate axis, you can see them in wire mode
  AObscuration* xa = new AObscuration("X", xaxis);
  AObscuration* ya = new AObscuration("Y", yaxis);
  AObscuration* za = new AObscuration("Z", zaxis);


  //opt -> AddNode(Wall, 1, new TGeoTranslation(0,0,0*mm));
  opt -> AddNode(Wall2, 1, new TGeoTranslation(0,0,0*mm));
  //opt -> AddNode(Wall3, 1, new TGeoTranslation(0,0,-5.1/2*mm));

  //opt -> AddNode(xa,1, new TGeoTranslation(0,0,0));
  //opt -> AddNode(ya, 1, new TGeoTranslation(0,0,0));
  //opt -> AddNode(za, 1, new TGeoTranslation(0,0,0));

}

void AddCamera(AOpticalComponent* opt) { // adds focal plane (camera) and camera shadow to the world
  //Make a focal plane and display camera

  const double alpha = TMath::ASin((3.6/2)/3.4) * TMath::RadToDeg();
  const double dz = 0*mm; //moves focal surface along the optical axis

  //create shapes to build camera and camera shadow
  TGeoSphere* sphericalCamera = new TGeoSphere("S2", kCameraIR , kCameraOR, 0, 2*alpha, -180, 180);
  TGeoBBox* rect = new TGeoBBox("CB",kCameraDx/2, kCameraDy/2, 1*m);
  //move rectangle
  TGeoTranslation* transZ2 = new TGeoTranslation("transZ2", 0,0, 1.6*m);
  transZ2->RegisterYourself();
  //create camera
  TGeoCompositeShape *cs2 = new TGeoCompositeShape("Cam" , "CB:transZ2*S2");
  AFocalSurface* focalPlane = new AFocalSurface("focalPlane", cs2);
  opt->AddNode(focalPlane, 1, new TGeoTranslation(0*m, 0*m, dz)); //

  //make a camera box as obscuration
  //Back of the Camera
  TGeoSphere* SphericalCameraBox = new TGeoSphere("SCB", kCameraIR+0.1*mm,kCameraOR+0.1*mm,0,2*alpha,-180,180);
  TGeoBBox* rectforCam = new TGeoBBox("RFC", kCameraDx/2+1*mm, kCameraDy/2+1*mm, 1*m);

  TGeoTranslation* transZ3 = new TGeoTranslation("transZ3", 0,0, 1.2*m);
  transZ3->RegisterYourself();

  TGeoCompositeShape *camBox = new TGeoCompositeShape("Camera Box", "RFC:transZ3*SCB");
  const double delta = kCameraOR-kCameraIR;

  AObscuration* camBoxObs = new AObscuration("Camera Box", camBox);
  opt -> AddNode(camBoxObs, 1, new TGeoTranslation(0,0, dz-2*mm));

}

void AddLenses(AOpticalComponent* opt){// adds schmidt corrector ring to the world

  //Adding Schmidt Corrector Ring, spherical approximation
  //parameters
  const double t = 12*mm;
  const double t2 = 6.9*mm;
  const double d1 = 1700*mm;
  const double d2 = 2200*mm;

  const double r = 8380*mm;
  const double y = 795*mm;
  //**********************building segments - tubeSegment**************************
  //build tube segment in the right size, only 1 deg of the whole circle
  TGeoTubeSeg* TubeSegement = new TGeoTubeSeg("TubeSegement",d1/2, d2/2, t/2 , 82.5,97.5);
  //Postion segment
  TGeoTranslation* tr0 = new TGeoTranslation("tr0",0,0,+t/2);
  tr0->RegisterYourself();
  //***********************building segment - with a tube**************************
  TGeoTubeSeg* tube = new TGeoTubeSeg("tube", 7*m, r, 1*m,82.5,97.5);

  //Postion Tube, needs to be rotated by 90 deg and placed on cutting edge
  TGeoRotation* rot = new TGeoRotation("rot",90, 90, 0);
  rot -> RegisterYourself();
  //new Postion of the big tube
  Double_t dx = 0;
  Double_t dy = y;//-(d1/2-y);
  Double_t dz = -r-t2+t;
  //rotate and translate
  TGeoCombiTrans *c1 = new TGeoCombiTrans("c1",dx,dy,dz, rot);
  c1 -> RegisterYourself();
  //building corrector ring segment
  TGeoCompositeShape* seg2 = new TGeoCompositeShape("seg2", "TubeSegement:tr0-tube:c1");
  //*******************************************************************************
  //making lens
  ALens* Segment = new ALens("Segment", seg2); //seg: with sphere, seg2: with tube
  Segment-> SetRefractiveIndex(bk7);

  //flat side facing the mirror
  //spherical side on the left per construction, rotating by 180°-> spherical side facing the mirror

  for(int deg = 0; deg < 360; deg += 15){
    opt -> AddNode(Segment, 1, new TGeoRotation(Form("rot%d",deg),deg,180,0));
  }

}

void AddMirror(AOpticalComponent* opt){ //adds Mirror without segmentation to the world

  //Add Mirror as Composite shapes
  double theta = TMath::ASin((kMirrorD/2)/kMirrorR) * TMath::RadToDeg();
 	TGeoSphere* mirSphere = new TGeoSphere("S", kMirrorR , (kMirrorR+kMirrorT), 0, 2*theta, 0, 360);

  TGeoBBox* box = new TGeoBBox("B", kMirrorD/2, kMirrorD/2, kMirrorD/2); // x,y,z sind halbe längen
  TGeoTranslation* transZ1 = new TGeoTranslation("transZ1", 0,0,kMirrorD);
  transZ1->RegisterYourself();

  //Create composite, as mirror is a rectangle with a spherical surface
  TGeoCompositeShape *recmirr = new TGeoCompositeShape("recmirr" , "B:transZ1*S");
  AMirror* mirror = new AMirror("mirror", recmirr);

  mirror->SetReflectance(1);
  opt->AddNode(mirror, 1, new TGeoTranslation(0*mm, 0*mm, 0*mm));

}
