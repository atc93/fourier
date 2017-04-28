#include "../inc/headers.hh"

int main (int argc, char* argv[]) {


	std::string st_inputFile = argv[1];
	TString rootFileName = "wShape.root";
	TFile *ff = new TFile(rootFileName,"recreate");

	TH1D *h_wShape = new TH1D("h_wShape", "h_wShape", 121, -60, 60);

	cout << st_inputFile << endl;
	std::ifstream infile(st_inputFile);
	double binID, intensity;
	while (infile >> binID >> intensity)
	{
		//cout << binID << " " << intensity << endl;
		h_wShape->SetBinContent(binID+61, intensity);
	}

	h_wShape->Write();
	ff->Close();

}
