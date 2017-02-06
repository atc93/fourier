// 	-------------------------------
//	|														  |
//	| g-2 fast rotation toy model |
//	|														  |
//	|       antoine chapelain		 	|
//	|       atc93@cornell.edu		  |
//	|														 	|
// 	-------------------------------

// root include files
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
#include <TLatex.h>
#include <TMath.h>
#include <TF1.h>

// c++ include files
#include <fstream>    // for reading files
#include <sstream>

// simulation paramaters
const double N_muon = 10000000; 	// number of muons to simulate
const double N_turn = 1;			// how many cyclotron revolution to simulate
const double p_spread = 0.;//0.112;	// muon momentum distribution spread in %
const double t_spread = 0.;//25;			// muon beam length spread in nano second
const bool   fillTree = false;  // fill root tree with time, momentum...
const bool   useMuon = true;		// simulate muons
const bool   useProton = false;		// simulate protons

// physical constants
const double M_mu = 0.1056583715 ; 								// muon mass in gev
const double M_p = 0.9382720813 ; 								// proton mass in gev
const double M_e = 0.510998928e-3 ;  							// electron mass in gev
const double q_e = 1.602176565e-19 ;							// charge of the positron in coulomb
const double c_light= 299792458 ; 								// speed of light in m/s
const double B = 1.4513; 													// b-field in t
const double R0= 7.112;														// the radius of the ideal orbit in meter
const double p_ptcl_magic = 3.094; 									// muon magic momentum, in gev/c
const double mass_factor = 1.782661845e-36*1e9 ; 	// kg for 1gev
const double muon_life = 2.1969811e-6; 						// muon life time in seconds
const double inch = 0.0254;

// math constants
const double PI = TMath::Pi();

// debugging mode
const int DEBUG=1;

// Set moomentum distribution
const int Gaus		= 1;
const int Landau	= 0;
const int Unif		= 0;
const int Custom 	= 0;

// namespaces
using namespace std;

// ppre-declaration of functions
Double_t ComputeTravelTimeMuon(double gamma, int n, double p_ptcl, double tzero); 
Double_t ComputeTravelTimeProton(double gamma, int n, double p_ptcl, double tzero); 
void DrawATLASLabel(TCanvas* _c);


// ---------------
//  main function
// ---------------

int main() {

	// initializaion of variables
	Double_t time;
	Double_t p_ptcl;
	Double_t E_ptcl;
	Double_t gamma;
	TRandom *r3 = new TRandom3(0);
	gRandom = r3;

	// root objects
	TFile *rootFile = new TFile("root/FastRotation.root","recreate");
	TTree *tr = new TTree("tr","tr");
	tr->Branch("time",&time,"time/d");
	tr->Branch("p_ptcl",&p_ptcl,"p_ptcl/d");

	TH1D *h_frs = new TH1D("h_frs", "fast rotation", 200000, 0, 20);
	h_frs->GetXaxis()->SetTitle("time [#mus]");
	h_frs->GetXaxis()->CenterTitle();
	h_frs->GetXaxis()->SetTitleOffset(1.1);
	h_frs->GetXaxis()->SetTitle("intensity");
	h_frs->GetXaxis()->CenterTitle();
	h_frs->GetXaxis()->SetTitleOffset(1.3);
	gStyle->SetOptStat(0);
	TCanvas c; c.cd();
  string _level="internal simulation";
  string _atlas_desc="#font[72]{g-2}";

	TH1D *h_freq = new TH1D("h_freq", "freq dist", 500,6.64e6,6.77e6);
	h_freq->GetXaxis()->SetTitle("Frequency [Hz]");
	TH1D *h_p		 = new TH1D("h_p", 		"momentum dist", 500,3.06,3.13);
	h_p->GetXaxis()->SetTitle("Momentum [GeV]");

  TLatex _g;
	_g.SetTextSize(0.035);

  //draw horizontal label  
  string _sh="";
  _sh+=_atlas_desc+" "+_level;
	gPad->Update();

	//drawatlaslabel(c);
  gPad->SetTicks(1);
  gPad->SetTicks(1);

	// processing tim book-kepping
	clock_t t1,t2;
	t1 = clock();

	cout << "\n-- start fast rotation routine --\n" << endl;

	cout << " momentum spread (1 sigma): " << p_spread << "%" << endl;
	cout << "     beam length (1 sigma): " << t_spread << " ns" << endl;
	cout << " number of simulated ptcle: " << N_muon << endl;
	cout << "               fill length: " << N_turn * 0.150 << " us" << endl;
	//cout << " draw random number for momentum spread" << endl;

	// draw random numbers
	vector<double> v_p_spread;
	vector<double> v_t_spread;
	TF1 *f1 = new TF1("f1","((sin(50*(x))/(50*(x))+1)*exp(-0.5*((x)/0.3)**2))",-1,1);
	f1->SetNpx(10000); 
	for (int i=0; i<N_muon; ++i) 	{
		if (Gaus) 	v_p_spread.push_back(gRandom -> Gaus(0, p_spread/100)); 
		if (Landau) v_p_spread.push_back(gRandom -> Landau(0, p_spread/250)); 
		if (Unif) 	v_p_spread.push_back(gRandom -> Uniform(-5*p_spread/100,5*p_spread/100)); 
		if (Custom) v_p_spread.push_back((f1->GetRandom())/150); 
		v_t_spread.push_back(gRandom -> Gaus(0, t_spread*1E-9));
	}

	// Muon loop
	//cout << " loop over muons\n" << endl;
	for (int i=0; i<N_muon; ++i) {

		// draw a random muon momentum
		p_ptcl = (1+v_p_spread.at(i))*p_ptcl_magic;
		h_p->Fill(p_ptcl);
		// compute particle energy
		if (useMuon)				E_ptcl = sqrt( p_ptcl * p_ptcl + M_mu * M_mu);
		else if (useProton)	E_ptcl = sqrt( p_ptcl * p_ptcl + M_p * M_p);
		// compute gamma factor
		if (useMuon)				gamma = E_ptcl / M_mu;
    else if (useProton) gamma = E_ptcl / M_p;

		// turn loop
		double time1 = 0., time2 = 0.;
		for (int turn=0; turn<N_turn; ++turn) {
			if (turn==0) time1=0;
			else time1 = time2;
			// compute travel time
			if (useMuon) 		time= ComputeTravelTimeMuon(gamma, turn, p_ptcl, v_t_spread.at(i))*1E6;
			if (useProton) 	time= ComputeTravelTimeProton(gamma, turn, p_ptcl, v_t_spread.at(i))*1E6;
			time2 = time;
			if (DEBUG) cout << "time = " << time << "   -   turn = " << turn << endl;
			if (DEBUG) cout << "time2 - time1 = " << time2 - time1 << "   -   turn = " << turn << endl;
			// fill ROOT Tree if enabled
			if (fillTree) tr->Fill();
			// fill histogram
			h_frs->Fill(time);
			h_freq->Fill(1E6*1/(time2-time1));
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
  //double _xleft = -0.87;
  //double _ytop = (1+0.135)*h_frs->GetMaximum();
	//_g.DrawLatex(_xleft,_ytop, _sh.c_str());
	h_frs->Print();
	h_freq->Write();
	h_p->Write();

	h_frs->Draw();
	c.SaveAs("plot/FastRotation.eps");
	c.SaveAs("plot/FastRotation.png");

  h_freq->Draw();
  c.SaveAs("plot/FreqDist.eps");
  c.SaveAs("plot/FreqDist.png");

  h_p->Draw();
  c.SaveAs("plot/MomentumDist.eps");
  c.SaveAs("plot/MomentumDist.png");

	tr->Write();
	rootFile->Close();

	return 0;

}

// function to compute travel time for the case:
// - longitudinal initial beam width
// - momentum spread

Double_t ComputeTravelTimeMuon(double gamma, int n, double p_ptcl, double tzero) {

  return (1 + n) * (2 * PI * gamma * M_mu / (c_light * c_light * 1E-9 * B)) - tzero;

}

Double_t ComputeTravelTimeProton(double gamma, int n, double p_ptcl, double tzero) {

  return (1 + n) * (2 * PI * gamma * M_p / (c_light * c_light * 1E-9 * B)) - tzero;

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
