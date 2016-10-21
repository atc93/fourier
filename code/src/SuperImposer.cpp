#include "../inc/headers.hh"

int main(){

//Access first approximation and Delta
TFile *f = new TFile("firstApprox.root","READ");
TH1F *h = (TH1F*)f->Get("hfft3");
int size = h->GetNbinsX();
TFile *f1 = new TFile("Delta.root","READ");
TH1F *h1 = (TH1F*)f1->Get("hfft2");
double fft_input[size];

//add first approximation and Delta
h->Add(h1,1);
for (int i = 1; i<=size; i++) {
        fft_input[i-1]=h->GetBinContent(i);
}
TH1F *h3 = new TH1F("hfft","Delta Freq Distribution, Signal Start Time 150ns, t0=75ns",size,0,1E9);
h3->GetXaxis()->SetTitle("Frequency [Hz]");
h3->GetXaxis()->CenterTitle();

for(int ii=1;ii<=size;ii++){
    h3->SetBinContent(ii,fft_input[ii-1]);
}

TCanvas c1;
h3->Draw();
h3->SetAxisRange(6667000, 6745000);
c1.SaveAs("secondApprox.eps");
TFile *FFT = new TFile("secondApprox.root","RECREATE");
h3->Write();
FFT->Close();

f->Close();
f1->Close();
return 0;
}
