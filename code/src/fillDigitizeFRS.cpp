#include "../inc/headers.hh"

int main (int argc, char* argv[]) {


	std::string st_inputFile = argv[1];
	TString rootFileName = "frs_fromAnalytical.root";
	TFile *ff = new TFile(rootFileName,"recreate");

	TH1D *h_frs = new TH1D("h_frs", "h_frs", 700000, 0, 700);

	cout << st_inputFile << endl;
	std::ifstream infile(st_inputFile);
	double binID, intensity;
	while (infile >> binID >> intensity)
	{
		//cout << binID << " " << intensity << endl;
		h_frs->SetBinContent(binID, intensity);
	}

	h_frs->Write();
	ff->Close();

}
