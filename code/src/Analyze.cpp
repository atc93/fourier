/* Code containing several function to analyze the tracking results from 
 * BMAD. Started by Antoine Chapelain (antoine.chapelain@cornell.edu) April 2015 */

#include "../inc/header.hh"

#define Debug 0
const double TWOPI = 2.0 * 3.1415926535897932384626433832795 ;
const double PI    = 3.1415926535897932384626433832795 ;

int main() {

	// Open ROOT file and associated tree
	TString inRootFileName = "/nfs/gm2/data1/achapela/BMAD/Data/StartTrackingInflectorExitPerfectBeamWithRealisticPandTspread/fiber.root" ;
	TFile *inFile = new TFile(inRootFileName, "UPDATE") ;
	TTree *inTree = (TTree*) inFile->Get("tree") ;

	// Define all the required variables
	Double_t turn, mu_x, mu_px, mu_y, mu_py, mu_z, mu_pz,time, number;
	TVector3 VecPolXZ ; // There is a hack: using component X and Y as X and Z for the g-2 geometry.
	Double_t PolOrientation ;

	// Assigned the above variables to the tree branches
	inTree->SetBranchAddress("<x>",  	&mu_x) ;
	inTree->SetBranchAddress("<px>",  &mu_px) ;
	inTree->SetBranchAddress("<y>",  	&mu_y) ;
	inTree->SetBranchAddress("<py>",  &mu_py) ;
	inTree->SetBranchAddress("<z>",  	&mu_z) ;
	inTree->SetBranchAddress("<pz>",  &mu_pz) ;
	inTree->SetBranchAddress("time",  &time) ;
	inTree->SetBranchAddress("number",  &number) ;

	// Define variable for the TGraph to plot polarization evolution
	int nentries = inTree->GetEntries() ;

	//y Define 1D histogram to plot the Muon detection time at the "fiber harp"
	TH1D *h_mu_t = new TH1D("h_fr","Muon detection time",700000,0,700E-6) ;

	// Loop over tree entries 
	for (int i=0; i<nentries; i++) {

		// Get event
		inTree->GetEntry(i) ;

		// Fill histogram
		h_mu_t->Fill(time, number) ;

	}

	// Plot and save histogram
	TCanvas c;
	h_mu_t->Draw();
	c.SaveAs("MuTime.eps") ;
	c.SaveAs("MuTime.png") ;
	TString outRootFileName = "/nfs/gm2/data1/achapela/BMAD/Data/StartTrackingInflectorExitPerfectBeamWithRealisticPandTspread/frs.root" ;
	TFile *outFile = new TFile(outRootFileName, "RECREATE") ;
	h_mu_t->Write();

	// Close Tfile
	outFile->Close();
	inFile->Close();

}
