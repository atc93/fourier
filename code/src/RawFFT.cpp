#include "../inc/headers.hh"

const double PI = 3.141592653589793238460;

const double period = 4E-9;
const double freq   = 1/period;
const double t_0    = 74.5E-9;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::string String;
typedef std::valarray<double> DArray;

//Given an interval in the signal, this function finds
//the time corresponding to the signal maximum

int MaxTimeFinder(int bin_i, int size, TH1D *h1){
	double InitialContent= h1->GetBinContent(bin_i);
	double t_max=h1->GetBinLowEdge(bin_i);
	for (int i=bin_i+1; i<=(size+bin_i-1); i++){
		double NewContent=h1->GetBinContent(i);
		if (InitialContent<=NewContent){ 
			t_max=h1->GetBinLowEdge(i);
			InitialContent=NewContent;
		}
	}
	return t_max;
}	
//This function Fourier analyzes an interval to find the
//dominant frequency and then the period
int PeriodFinder(int bin_i, int size, TH1D *h1){
	Complex fft_input[size];
	for (int i = bin_i; i<=(size+bin_i-1); i++) {
		int j=i-bin_i;
		fft_input[j]=h1->GetBinContent(i);
	}
	CArray data(fft_input, size);
	double re[size];
	double im[size];

	for (int i=0;i<size;++i){
		re[i]=std::real(data[i]);
		im[i]=std::imag(data[i]);
	}

	TVirtualFFT *fft_own1=TVirtualFFT::FFT(1, &size, "C2CBACKWARD ES K");
	fft_own1->SetPointsComplex(re,im);
	fft_own1->Transform();
	fft_own1->GetPointsComplex(re,im);

	TH1D *hfft1 = new TH1D("hfft","RE_FFT of Toy Model, Perfect Start Time",size,0,1E9);
	hfft1->GetXaxis()->SetTitle("Frequency [Hz]");
	hfft1->GetXaxis()->CenterTitle();

	for(int ii=1;ii<=size;ii++){
		hfft1->SetBinContent(ii,re[ii-1]);
	}

	double InitialContent= hfft1->GetBinContent(64);
	double f_max=hfft1->GetBinLowEdge(64);
	for (int i=65; i<=72; i++){
		double NewContent=hfft1->GetBinContent(i);
		if (InitialContent<=NewContent){
			f_max=hfft1->GetBinLowEdge(i);
			InitialContent=hfft1->GetBinContent(i);
		}
	}
	double period=1.0/f_max;
	std::cout<<f_max<<std::endl;
	return period;
}

int main(){

	//Access the Fast Rotation Signal File
	TFile *f = new TFile("root/FastRotation10M_0.12%pspread_25nstspread.root", "READ");
	TH1D *h1 = (TH1D*)f->Get("h_fr");
	int nBinsX = h1->GetNbinsX();
	// rebinning to 4 ns bins to get the maximum bins allowed by the FFT
	int size= 700E-6 / period; // compute how many bins there are given 700 us fill and bin width
	int rebin;
	if (0==0) {rebin = nBinsX / size;}
	else {std::cout << "- ERROR - Cannot rebinning with an integer" << std::endl; return 0;}
	h1->Rebin(rebin);
	Complex fft_input[size];

	//Fill an array with the histogram data
	for (int i = 1; i<=size; i++) {
		fft_input[i-1]=h1->GetBinContent(i);
	}
	//f->Close();

	CArray data(fft_input, size);
	double re[size];
	double im[size];

	//to do the C2C FFT, you have to split your data into real and purely imaginary parts
	for (int i=0;i<size;++i){
		re[i]=std::real(data[i]);
		im[i]=std::imag(data[i]);
	}

	TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &size, "C2CBACKWARD ES K");//decalare FFT object
	fft_own->SetPointsComplex(re,im);//transform the data
	fft_own->Transform();
	fft_own->GetPointsComplex(re,im);//overwrite re and im with the transformed data

	//add appropriate normalization; see paper for specifics
	for (int i=0;i<size;++i){
		re[i]=((2e-9)*re[i])/(PI);
		im[i]=((2e-9)*im[i])/(PI);
	}

	//Multiplication by Exp[-i*omega*t_0]
	for(int n=0;n<size;n++){
		re[n]=re[n]*std::real(std::polar(1.0,(-2*PI*n*t_0)/(size*(1e-9))))-im[n]*std::imag(std::polar(1.0,(-2*PI*n*t_0)/(size*(1e-9))));
		//	im[n]=re1[n]*std::imag(std::polar(1.0,(-2*PI*n*t_0)/(size*(1e-9))))+im[n]*std::real(std::polar(1.0,(-2*PI*n*t_0)/(size*(1e-9))));
	}

	//Fill your new histogram with the real output; the maximum of the range is the inverse of the sampling rate of the FRS 
	TH1D *hfft = new TH1D("hfft","",size,0,freq);
	hfft->GetXaxis()->SetTitle("Frequency [Hz]");
	hfft->GetXaxis()->SetTitleOffset(1.3);
	hfft->GetYaxis()->SetTitle("ReFFT");
	hfft->GetXaxis()->CenterTitle();

  int minBin=0;
  int maxBin=0;
  bool searchBin = true;
	for(int ii=1;ii<=size;ii++){
		hfft->SetBinContent(ii,re[ii-1]);
    if (hfft->GetBinCenter(ii)<6630000) minBin =ii-1;
    if (searchBin) if (hfft->GetBinCenter(ii)>6780000) {maxBin =ii+1; searchBin=false;};
	}


	TCanvas c;
	gStyle->SetOptStat(1100);
	gPad->SetTicks(1);
	hfft->Draw();

	// Extract and rescale histo yMax and yMin
	double yMax = -999;
	for (int i=minBin; i<=maxBin; ++i) if (hfft->GetBinContent(i)>yMax) yMax = hfft->GetBinContent(i);
	if (yMax>0) 		 {hfft->SetMaximum(yMax*(1+0.05)); yMax = yMax*(1+0.05);}
	else if (yMax<0) {hfft->SetMaximum(yMax*(1-0.05)); yMax = yMax*(1-0.05);}
	double yMin = 999;
	for (int i=minBin; i<=maxBin; ++i) if (hfft->GetBinContent(i)<yMin) yMin = hfft->GetBinContent(i);
  if (yMin>0)      {hfft->SetMinimum(yMin*(1-0.05)); yMin = yMin*(1-0.05);}
  else if (yMin<0) {hfft->SetMinimum(yMin*(1+0.05)); yMin = yMin*(1+0.05);}

	// Draw vachuum chambers limits 
	double minLine = 0;
	double maxLine = 0;
	if (yMax>0 &&yMin>0) {
		std::cout << yMax-yMin << std::endl;
	}
	else if(yMax<0 && yMin<0) {
		std::cout << fabs(yMin)-fabs(yMax) << std::endl;
	}
	else if (yMax>0 && yMin<0) {
		minLine = yMin+(yMax-yMin)*0.45;
		maxLine = yMax-(yMax-yMin)*0.45;
	}

	TPaveText *pt = new TPaveText(6640000,minLine,6666000,maxLine);
	TPaveText *pt2 = new TPaveText(6747000,minLine,6773000,maxLine);
	pt->AddText("outside");
	pt->AddText("vacuum chamber");
	pt->SetShadowColor(0);
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetLineWidth(0);
	pt->SetLineColor(0);
	pt->SetTextAngle(90);
	pt2->AddText("outside");
	pt2->AddText("vacuum chamber");
	pt2->SetShadowColor(0);
	pt2->SetBorderSize(0);
	pt2->SetFillColor(0);
	pt2->SetLineWidth(0);
	pt2->SetLineColor(0);
	pt2->SetTextAngle(90);
	hfft->SetAxisRange(6667000, 6745000);
	hfft->SetAxisRange(6630000, 6780000);
	TLine *line1 = new TLine(6667000,yMin,6667000,yMax);
	TLine *line2 = new TLine(6745000,yMin,6745000,yMax);
	line1->SetLineWidth(2.5);
	line2->SetLineWidth(2.5);
	line1->SetLineColor(1);
	line2->SetLineColor(1);
	line1->Draw("same");
	line2->Draw("same");
	pt->Draw("same");
	pt2->Draw("same");

	c.SaveAs("RawFFT.eps");
	c.SaveAs("RawFFT.png");
	//c.SetLogy();
	//c.SaveAs("REFFT_0.2pct_momspread_log_ts-100ns_t0-79ns.eps");

	TFile *FFT = new TFile("RawFFT.root","RECREATE");
	hfft->Write();
	h1->Write();
	hfft->SetAxisRange(6400000, 7000000);
	c.SaveAs("RawFFT_zoomOut.eps");
	FFT->Close();
	return 0;
}

