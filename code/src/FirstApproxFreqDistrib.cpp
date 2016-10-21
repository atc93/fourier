#include "../inc/headers.hh"

const double minFreqCutoff = 6.675E6; // lower bound of the "central part"
const double maxFreqCutoff = 6.725E6; // upper bound of the "central part"

int main(){

//access the raw FFT
TFile *f = new TFile("RawFFT.root", "READ");
TH1F *h1 = (TH1F*)f->Get("hfft");
double samplerate=h1->GetBinWidth(1);
double maxFreq = -999;
int size = h1->GetNbinsX();
double fft_input[size];

//set the frequencies corresponding to those outside the vacuum chamber to 0,
//and use those inside the vacuum chamber. This is the central part of the FFT,
//to be used as part of the first approximation for the frequency distribution,
//see Orlov et al
maxFreq = h1->GetBinLowEdge(size)+h1->GetBinWidth(size);
for (int i = 1; i<=size; i++) {
  int j=i-1;
	//identify the turning points of the FFT in the Vacuum chamber frequency range
	if( (h1->GetBinLowEdge(i) + h1->GetBinWidth(i)) < minFreqCutoff || h1->GetBinLowEdge(i) > maxFreqCutoff) {
		fft_input[j]=0;
	}
	else {fft_input[j]=h1->GetBinContent(i);
		std::cout << fft_input[j] << std::endl;}
}
f->Close();

//Plot the first approximation
TH1F *hfft3 = new TH1F("hfft3","First Approx. of Frequency Spectrum within Vacuum Chamber",size,0,maxFreq);
hfft3->GetXaxis()->SetTitle("Frequency [Hz]");
hfft3->GetXaxis()->CenterTitle();

for(int ii=1;ii<=size;ii++){
    hfft3->SetBinContent(ii,fft_input[ii-1]);
}
TCanvas c;
hfft3->Draw();
hfft3->SetAxisRange(6667000, 6745000);
c.SaveAs("firstApprox.eps");
TFile *FFT = new TFile("firstApprox.root","RECREATE");
hfft3->Write();
FFT->Close();
return 0;
}
