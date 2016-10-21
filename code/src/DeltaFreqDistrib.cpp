#include "../inc/headers.hh"

const double PI = 3.141592653589793238460;

int main(){

//Access the first approximation
TFile *f = new TFile("firstApprox.root", "READ");
TH1F *h1 = (TH1F*)f->Get("hfft3");
double samplerate=h1->GetBinWidth(1);
int size = h1->GetNbinsX();
double maxFreq = h1->GetBinLowEdge(size)+h1->GetBinWidth(size);
double fft_input[size];
int bin_i=1;
int bin_f=size;
int minBin=-1;
int maxBin=-1;
bool searchBin=true;

//The derivation of the correction to the first approximation relies on formula 5 in Orlov
//et al, so it is an inverse FFT. We must then introduce the appropriate normalization,
//the range of the frequency spectrum in the raw FFT
for (int i =bin_i; i<=bin_f; i++) {
        int j=i-1;
	fft_input[j]=(maxFreq)*(h1->GetBinContent(i));
	if( (h1->GetBinLowEdge(i) + h1->GetBinWidth(i)) < 6.675E6) minBin = i;
	if (searchBin) if( h1->GetBinLowEdge(i) > 6.725E6) {maxBin=i; searchBin = false;}
}
double delta[size];

//perform the integration in equation 18
for(int n=0;n<size;n++){
	double i=(h1->GetBinLowEdge(bin_i))+n*samplerate;
        for (int m=minBin;m<maxBin;m++){
		double j=(h1->GetBinLowEdge(bin_i))+m*samplerate;
		if (abs(i-j)>1e-10){
			delta[n]+=fft_input[m]*(sin(2*PI*((2.6e-6)-(75e-9))*(i-j))/(i-j));
		}
		else{
			delta[n]+=fft_input[m]*(0.000015865);
		}
	}

	//if (n%1000==0) std::cout<<n<<std::endl;
}
for(int n=0;n<size;n++){
	delta[n]=delta[n]*(0.4)*samplerate/(2*PI*PI*size);
}
f->Close();

//Create a plot of Delta
TH1F *hfft2 = new TH1F("hfft2","Addition to Frequency Spectrum",size,0,1E9);
hfft2->GetXaxis()->SetTitle("Frequency [Hz]");
hfft2->GetXaxis()->CenterTitle();

for(int ii=1;ii<=size;ii++){
    hfft2->SetBinContent(ii,delta[ii-1]);
}

TCanvas c;
hfft2->Draw();
hfft2->SetAxisRange(6667000, 6745000);
c.SaveAs("Delta.eps");

TFile *FFT = new TFile("Delta.root","RECREATE");
hfft2->Write();
FFT->Close();
return 0;
}
