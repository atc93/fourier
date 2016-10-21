// 	-------------------------------
//	|														  |
//	| G-2 Fast Rotation toy model |
//	|														  |
//	|       Antoine Chapelain		 	|
//	|       atc93@cornell.edu		  |
//	|														 	|
// 	-------------------------------

// ROOT include files
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TTree.h>
#include <time.h>
#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <vector>
#include <TPad.h>
#include <TLatex.h>
#include <TMath.h>

// C++ include files
#include <fstream>    // for reading files
#include <sstream>

// Simulation paramaters
const double N_muon = 100000; 	// Number of muons to simulate
const double N_turn = 5000;			// How many cyclotron revolution to simulate
const double p_spread = 0.112;	// Muon momentum distribution spread in %
const double t_spread = 25;			// Muon beam length spread in nano second
const bool   fillTree = false;  // Fill ROOT tree with time, momentum...

// Physical constants
const double M_mu = 0.1056583715 ; 								// Muon mass in GeV
const double M_e = 0.510998928e-3 ;  							// Electron mass in GeV
const double q_e = 1.602176565e-19 ;							// Charge of the positron in Coulomb
const double c_light= 299792458 ; 								// Speed of light in m/s
const double B = 1.4513; 													// B-field in T
const double R0= 7.112;														// The radius of the ideal orbit in meter
const double p_mu_magic = 3.094; 									// Muon magic momentum, in GeV/c
const double mass_factor = 1.782661845e-36*1e9 ; 	// kg for 1GeV
const double muon_life = 2.1969811e-6; 						// Muon life time in seconds
const double inch = 0.0254;

// Math constants
const double PI = TMath::Pi();

// Debugging mode
const int DEBUG=0;

// Namespaces
using namespace std;

// Ppre-declaration of functions
Double_t ComputeTravelTime(double gamma, int n, double p_mu, double tzero); 
void DrawATLASLabel(TCanvas* _c);


// ---------------
//  Main function
// ---------------

int main() {

	// Initializaion of variables
	Double_t time;
	Double_t p_mu;
	Double_t E_mu;
	Double_t gamma;
	TRandom *r3 = new TRandom3(0);
	gRandom = r3;

	// ROOT objects
	TFile *rootFile = new TFile("root/FastRotation.root","recreate");
	TTree *tr = new TTree("tr","tr");
	tr->Branch("time",&time,"time/D");
	tr->Branch("p_mu",&p_mu,"p_mu/D");

	TH1D *h_fr = new TH1D("h_fr", "Fast Rotation", 7000000, 0, 700);
	h_fr->GetXaxis()->SetTitle("Time [#mus]");
	h_fr->GetXaxis()->CenterTitle();
	h_fr->GetXaxis()->SetTitleOffset(1.1);
	h_fr->GetYaxis()->SetTitle("# entry");
	h_fr->GetYaxis()->CenterTitle();
	h_fr->GetYaxis()->SetTitleOffset(1.3);
	gStyle->SetOptStat(0);
	TCanvas c;c.cd();
  string _level="internal simulation";
  string _ATLAS_desc="#font[72]{g-2}";

  TLatex _g;
	_g.SetTextSize(0.035);

  //draw horizontal label  
  string _sh="";
  _sh+=_ATLAS_desc+" "+_level;
	gPad->Update();

	//DrawATLASLabel(c);
  gPad->SetTicks(1);
  gPad->SetTicks(1);

	// Processing tim book-kepping
	clock_t t1,t2;
	t1 = clock();

	cout << "\n-- Start Fast Rotation routine --\n" << endl;

	cout << " Momentum spread (1 sigma): " << p_spread << "%" << endl;
	cout << "     Beam length (1 sigma): " << t_spread << " ns" << endl;
	cout << " Number of simulated muons: " << N_muon << endl;
	cout << "               Fill length: " << N_turn * 0.150 << " us" << endl;
	//cout << " draw random number for momentum spread" << endl;

	// Draw random numbers
	vector<double> v_p_spread;
	vector<double> v_t_spread;
	for (int i=0; i<N_muon; ++i) 	{
		v_p_spread.push_back(gRandom -> Gaus(0, p_spread/100)); 
		v_t_spread.push_back(gRandom -> Gaus(0, t_spread*1E-9));
	}

	// Muon loop
	//cout << " loop over muons\n" << endl;
	for (int i=0; i<N_muon; ++i) {

		// draw a random muon momentum
		p_mu = (1+v_p_spread.at(i))*p_mu_magic;
		// compute muon energy
		E_mu = sqrt( p_mu * p_mu + M_mu * M_mu);
		// compute gamma factor
		gamma = E_mu / M_mu;

		// turn loop
		for (int turn=0; turn<N_turn; ++turn) {
			// compute travel time
			time= ComputeTravelTime(gamma, turn, p_mu, v_t_spread.at(i))*1E6;
			if (DEBUG) cout << "time = " << time << "   -   turn = " << turn << endl;
			// fill ROOT Tree if enabled
			if (fillTree) tr->Fill();
			// fill histogram
			h_fr->Fill(time);
		} // end turn loop

		// Monitor progress
		int progress = N_muon / 10;
		if(i % progress == 0){
			float percent = i / N_muon;
			cout << "--> " << percent*100 << "\% accomplished        (" << i << "/" << N_muon << "  muons processed)" << endl;
		}

	} // end Muon loop

	// Processing time book-keeping
	t2 = clock();
	float diff ((float)t2-(float)t1);
	cout << "\n  -- Routine DONE --\n" << endl;
	cout << "Time elapsed: " << diff/CLOCKS_PER_SEC/60 << " minutes." << endl;

	// ROOT objects saving
	h_fr->Draw();
  //double _xleft = -0.87;
  //double _ytop = (1+0.135)*h_fr->GetMaximum();
	//_g.DrawLatex(_xleft,_ytop, _sh.c_str());
	h_fr->Print();
	h_fr->Write();
	c.SaveAs("plot/FastRotation.eps");
	c.SaveAs("plot/FastRotation.png");
	tr->Write();
	rootFile->Close();

	return 0;

}

// function to compute travel time for the case:
// - longitudinal initial beam width
// - momentum spread

Double_t ComputeTravelTime(double gamma, int n, double p_mu, double tzero) {

  return (1 + n) * (2 * PI * gamma * M_mu / (c_light * c_light * 1E-9 * B)) - tzero;

}


// draw g-2 label

void DrawATLASLabel(TCanvas* _c) {

  string _level="Internal Simulation";
  string _ATLAS_desc="#font[72]{g-2}";

	TLatex* _g;
	double _xleft = 0.1;
	double _ytop = 1.0;

  _c->cd();
  //draw horizontal label  
    string _sh="";
  _sh+=_ATLAS_desc+" "+_level;
  _g->DrawLatex(_xleft,_ytop, "");
  //_g->DrawLatex(_xleft,_ytop-_g->GetTextSize()*1.15, _sh.c_str());

}
