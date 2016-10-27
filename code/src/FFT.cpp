#include "../inc/headers.hh"

const double PI = 3.141592653589793238460;

const double timeBinWidth = 4; // bin width for the F(t) that will be used for the FFTs [ns]
const double freq   = 1/(timeBinWidth*1E-9);
const double t_0    = 152.5E-9;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::string String;
typedef std::valarray<double> DArray;

// === FUNCTION PRE-DECLARATION === //


int MaxTimeFinder(int bin_i, int size, TH1D *h1);
int PeriodFinder(int bin_i, int size, TH1D *h1);


// === MAIN FUNCTION === //


int main(int argc, char* argv[]){

	// Access the Fast Rotation Signal File
	string inRootFileName = argv[1];
	TFile *f = new TFile(inRootFileName.c_str(), "READ");
	TH1D *h1 = (TH1D*)f->Get("h_fr");
	int nBinsX = h1->GetNbinsX();

	// Rebinning to timeBinWidth (4 ns by default to get the almost maximum number of bins allowed by the FFTW package)
	int size= 700000 / timeBinWidth; // compute how many 4 ns bins there are given 700 us fill and bin width. If not integer, size will have the closest interger value
	int rebin;
	if (0==0) {rebin = nBinsX / size;}
	else {std::cout << "- ERROR - Cannot rebinning with an integer" << std::endl; return 0;}
	h1->Rebin(rebin);
	Complex fft_input[size];

	// Fill an array with the histogram data
	for (int i = 1; i<=size; i++) {
		fft_input[i-1]=h1->GetBinContent(i);
	}

	CArray data(fft_input, size);
	double re[size];
	double im[size];
	double mag[size];

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
		im[n]=im[n]*std::imag(std::polar(1.0,(-2*PI*n*t_0)/(size*(1e-9))))+im[n]*std::real(std::polar(1.0,(-2*PI*n*t_0)/(size*(1e-9))));
		mag[n] = re[n]*re[n] + im[n]*im[n];
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
		//hfft->SetBinContent(ii,mag[ii-1]);
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
	std::cout << yMax << std::endl;
	if (yMax>0) 		 {hfft->SetMaximum(yMax*(1+0.05)); yMax = yMax*(1+0.05);}
	else if (yMax<0) {hfft->SetMaximum(yMax*(1-0.05)); yMax = yMax*(1-0.05);}
	std::cout << yMax << std::endl;
	double yMin = 999;
	for (int i=minBin; i<=maxBin; ++i) if (hfft->GetBinContent(i)<yMin) yMin = hfft->GetBinContent(i);
	if (yMin>0)      {hfft->SetMinimum(yMin*(1-0.05)); yMin = yMin*(1-0.05);}
	else if (yMin<0) {hfft->SetMinimum(yMin*(1+0.05)); yMin = yMin*(1+0.05);}
	std::cout << yMin << std::endl;

	// Draw vachuum chambers limits 
	double minLine = 0;
	double maxLine = 0;
	if (yMax>0 &&yMin>0) {
		std::cout << yMax-yMin << std::endl;
		minLine = yMin+(yMax-yMin)*0.45;
		maxLine = yMax-(yMax-yMin)*0.45;
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

	c.SaveAs("plot/FRS_ReFFT.eps");
	c.SaveAs("plot/FRS_ReFFT.png");

	TFile *FFT = new TFile("RawFFT.root","RECREATE");
	hfft->Write();
	h1->Write();
	hfft->SetAxisRange(6400000, 7000000);
	c.SaveAs("analyticREFFT_0.2pct_momspread_ts-0ns_t0-75ns_zoomOut.eps");


	// **** COMPUTE THE MAGNITUDE OF THE FFT **** //


	// For different start times tS: 1 tick = timeBinWidth 
	// FYI: delay = 1 to avoid having to code something specific for delay = 0
	// default tS = 0, 1us, 2us, 3us, 4us, 5us, 6us 
	const int ntS = 7; // number of different tS
	float  delay[ntS] = {1, 	1000/timeBinWidth, 
		2000/timeBinWidth, 
		3000/timeBinWidth, 
		4000/timeBinWidth, 
		5000/timeBinWidth, 
		6000/timeBinWidth}; 

	// Create a legend to label the different tS
	double _xl=.15;
	double _yt=.86;
	double _xr=.3;
	double _yb=.66;
	TLegend *leg = new TLegend(_xl, _yb, _xr, _yt);

	// Loop over the different tS
	for (int i=0; i<ntS; ++i) { 

		cout << delay[i] << endl;

		// For legend labels and histo names
		stringstream ss_name, ss_leg;
		ss_leg << "t_{s} = " << delay[i] * timeBinWidth/1000 << " #mus";
		ss_name << "h_" << delay[i]*4 << "ns";

		// Fill a temporary histogram that allows to shift tS and perform the FFT
		size = size - delay[i];
		TH1 *hm =0;
		TVirtualFFT::SetTransform(0);
		TH1D *h1bis = new TH1D("","",size,0,freq);
		for (int ii=delay[i]; ii<=size+delay[i]; ++ii) h1bis->SetBinContent(ii,h1->GetBinContent(ii));
		hm = h1bis->FFT(hm, "MAG");

		// Crate and fill the output histogram
		TH1F *hmag = new TH1F("hmag","FFT MAG",size,0,freq);
		hmag->GetXaxis()->SetTitle("Frequency [Hz]");
		hmag->GetXaxis()->CenterTitle();
		for(int ii=delay[i];ii<=size;ii++)	hmag->SetBinContent(ii,hm->GetBinContent(ii));
		leg->AddEntry(hmag,ss_leg.str().c_str(), "l");
		hmag->SetLineColor(1+i);

		// Draw histograms
		if (i==0) hmag->Draw();
		else hmag->Draw("same");
		hmag->SetAxisRange(6630000, 6780000);
		hmag->Write(ss_name.str().c_str());
		if (i==ntS-1) leg->Draw("same");

	}


	leg->Draw("same");

	c.SaveAs("plot/FRS_FFT_MAG.eps");
	c.SaveAs("plot/FRS_FFT_MAG.pdf");
	c.SaveAs("plot/FRS_FFT_MAG.png");
	c.Write("MAG");

	FFT->Close();
	f->Close();

	return 0;
}



// === FUNCTION DEFINITION === //


// Given an interval in the signal, this function finds
// the time corresponding to the signal maximum

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


// This function Fourier analyzes an interval to find the
// dominant frequency and then the period
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
