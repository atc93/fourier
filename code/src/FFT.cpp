#include "../inc/headers.hh"

const double PI = TMath::Pi();

const bool muon=true;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::string String;
typedef std::valarray<double> DArray;

//==================================//
//==== FUNCTION PRE-DECLARATION ====//
//==================================//

int MaxTimeFinder(int bin_i, int size, TH1D *h1);
int PeriodFinder(int bin_i, int size, TH1D *h1);
void inlineHelp();


//=======================//
//==== MAIN FUNCTION ====//
//=======================//

int main(int argc, char* argv[]){

	//*************************//
	//**** PARSING IN-LINE ****//
	//*************************//

	// helper for in-line options
	if (argc == 1 ) {
		inlineHelp();
		return 0;

	} else if (argc == 2) {
		string arg = argv[1];

		if (arg.compare("-h") == 0 || arg.compare("-H") == 0) {
			inlineHelp();
			return 0;
		}
	}

	// variable declaration
	string	inRootFileName;
	string	histoName;
	double	t0 = -999; // default non-sensical initialization
	double	tS = -999;
	double	rebinWidth = -999; // default non-sensical initialization
	double	maxTime = -999;
	double	samplingPeriod = -999;
	double	samplingFreq = -999;
	bool		rebin = false;		

	// parse in-line options
	for (int i = 1; i < argc; ++i) {

		string arg = argv[i];

		if ( arg.compare("-f") == 0 ) {          // input file to analyze
			inRootFileName = argv[i+1];
			// check that the input file exists
			if ( !std::ifstream(inRootFileName).good() ) {
				cout << "Error: File does not exist." << endl;
				return -1;
			}

		}	else if (arg.compare("-samplingPeriod") == 0) {
				samplingPeriod 	= atof(argv[i+1])*1E-9;
				samplingFreq		= 1/(samplingPeriod);
				cout << "sampling period: \t" << samplingPeriod << " ns" << endl;
				cout << "sampling frequency: \t" << samplingFreq << " GHz" << endl;
	
		}	else if (arg.compare("-maxTime") == 0) {	
				maxTime = atof(argv[i+1]);
				cout << "max time: " << maxTime << endl;

	
		}	else if (arg.compare("-t0") == 0) {		// t0 (lower limit of Fourier Transform)
			t0 = atof(argv[i+1]);
			cout << "t0: " << t0 << endl;

    } else if (arg.compare("-tS") == 0) {   // t0 (lower limit of Fourier Transform)
      tS = atof(argv[i+1]);
      cout << "tS: " << tS << endl;

		}	else if ( arg.compare("-rebin") == 0 ) {
				arg = argv[i+1];
				if ( arg.compare("true") == 0 ) rebin = true;
				else 														rebin = false; 

		}	else if ( arg.compare("-rebinWidth") == 0 ) {
				rebinWidth = atof(argv[i+1]);

		}	else if ( arg.compare("-histoName") == 0 ) {
				histoName = argv[i+1];

		}

	}

	// check the required options / input
	if ( t0 == -999 ) {
		cout << "Error: t0 was not provided" << endl;
		return -1;
	}
	if ( histoName.empty() ) {
		cout << "Error: FRS histo name not provided" << endl;
		return -1;
	}


	//*****************************//
	//**** ROOT FILE OPERATION ****//
	//*****************************//

	// Access the Fast Rotation Signal File
	TFile *f = new TFile(inRootFileName.c_str(), "READ");
	TH1D *h1 = (TH1D*)f->Get(histoName.c_str());
	int nBinsX = h1->GetNbinsX();

	// Rebinning to samplingPeriod (4 ns by default to get the almost maximum number of bins allowed by the FFTW package)
	int size= ( maxTime -tS ) / (samplingPeriod/1E-9); // compute how many 4 ns bins there are given 700 us fill and bin width. If not integer, size will have the closest interger value
	cout << "number of bins: " << size << endl;
	//int rebin;
	//if (0==0) {rebin = nBinsX / size;}
	//else {std::cout << "- ERROR - Cannot rebinning with an integer" << std::endl; return 0;}
	if (rebin) h1->Rebin(rebinWidth);
	Complex *fft_input = new Complex[size];

	// Fill an array with the histogram data
	for (int i = 1; i<=size; ++i) {
		fft_input[i-1]=h1->GetBinContent(i);
	}

	CArray data(fft_input, size);
	double *re = new double[size];
	double *im = new double[size];
	double *mag = new double[size];

	//need to add empty data at the begining?
	for (int i=0;i<int(tS);++i){
    re[i]=0;
    im[i]=0;
  }

	//to do the C2C FFT, you have to split your data into real and purely imaginary parts
	for (int i=tS;i<size;++i){
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

	//Multiplication by Exp[-i*omega*t0]
	for(int n=0;n<size;++n){
		re[n]=re[n]*std::real(std::polar(1.0,(-2*PI*n*t0)/(size)))-im[n]*std::imag(std::polar(1.0,(-2*PI*n*t0)/(size)));
		im[n]=im[n]*std::imag(std::polar(1.0,(-2*PI*n*t0)/(size)))+im[n]*std::real(std::polar(1.0,(-2*PI*n*t0)/(size)));
		mag[n] = re[n]*re[n] + im[n]*im[n];
	}

	//Fill your new histogram with the real output; the maximum of the range is the inverse of the sampling rate of the FRS 
	TH1D *hfft = new TH1D("hfft","",size/2,0,samplingFreq/2); // divide by 2 cause Nyquist
	hfft->GetXaxis()->SetTitle("Frequency [Hz]");
	hfft->GetXaxis()->SetTitleOffset(1.3);
	hfft->GetYaxis()->SetTitle("ReFFT");
	hfft->GetXaxis()->CenterTitle();

  TH1D *hfftIm = new TH1D("hfftIm","",size/2,0,samplingFreq/2); // divide by 2 cause Nyquist
  hfftIm->GetXaxis()->SetTitle("Frequency [Hz]");
  hfftIm->GetXaxis()->SetTitleOffset(1.3);
  hfftIm->GetYaxis()->SetTitle("ImFFT");
  hfftIm->GetXaxis()->CenterTitle();

  TH1D *hfftMag = new TH1D("hfftMag","",size/2,0,samplingFreq/2); // divide by 2 cause Nyquist
  hfftMag->GetXaxis()->SetTitle("Frequency [Hz]");
  hfftMag->GetXaxis()->SetTitleOffset(1.3);
  hfftMag->GetYaxis()->SetTitle("MagFFT");
  hfftMag->GetXaxis()->CenterTitle();

	int minBin=0;
	int maxBin=0;
	bool searchBin = true;
	for(int ii=1;ii<=size;++ii){
		hfft->SetBinContent(ii,sqrt(re[ii-1]*re[ii-1])); // FORCE POSITIVE AMPLITUDE SPECTRUM
		hfftIm->SetBinContent(ii,sqrt(im[ii-1]*im[ii-1])); // FORCE POSITIVE AMPLITUDE SPECTRUM
		hfftMag->SetBinContent(ii,sqrt(re[ii-1]*re[ii-1] + im[ii-1]*im[ii-1])); // FORCE POSITIVE AMPLITUDE SPECTRUM
		if (hfft->GetBinCenter(ii)<6630000) minBin =ii-1;
		if (searchBin) if (hfft->GetBinCenter(ii)>6780000) {maxBin =ii+1; searchBin=false;};
	}

	TCanvas c;
	gStyle->SetOptStat(1100);
	gPad->SetTicks(1);
	hfft->Draw();
	hfftIm->SetLineColor(2);
	//hfftIm->Draw("same");
	hfftMag->SetLineColor(1);
	//hfftMag->Draw("same");

	if (muon) {
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
	//hfft->SetAxisRange(6667000, 6745000);
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
	}
	else {
	hfft->SetAxisRange(6300000, 6780000);
	}

	c.SaveAs("plot/FRS_ReFFT.eps");
	c.SaveAs("plot/FRS_ReFFT.png");

	TFile *FFT = new TFile("RawFFT.root","RECREATE");
	hfft->Write();
	h1->Write();
	hfft->SetAxisRange(6400000, 7000000);
	c.SaveAs("analyticREFFT_0.2pct_momspread_ts-0ns_t0-75ns_zoomOut.eps");


	// **** COMPUTE THE MAGNITUDE OF THE FFT **** //


	// For different start times tS: 1 tick = samplingPeriod 
	// FYI: delay = 1 to avoid having to code something specific for delay = 0
	// default tS = 0, 1us, 2us, 3us, 4us, 5us, 6us 
	const int ntS = 7; // number of different tS
	float  delay[ntS] = {1, 	1000/samplingPeriod, 
		2000/samplingPeriod, 
		3000/samplingPeriod, 
		4000/samplingPeriod, 
		5000/samplingPeriod, 
		6000/samplingPeriod}; 

	// Create a legend to label the different tS
	double _xl=.15;
	double _yt=.86;
	double _xr=.3;
	double _yb=.66;
	TLegend *leg = new TLegend(_xl, _yb, _xr, _yt);

	// Loop over the different tS
	for (int i=0; i<0; ++i) { 

		cout << delay[i] << endl;

		// For legend labels and histo names
		stringstream ss_name, ss_leg;
		ss_leg << "t_{s} = " << delay[i] * samplingPeriod/1000 << " #mus";
		ss_name << "h_" << delay[i]*4 << "ns";

		// Fill a temporary histogram that allows to shift tS and perform the FFT
		size = size - delay[i];
		TH1 *hm =0;
		TVirtualFFT::SetTransform(0);
		TH1D *h1bis = new TH1D("","",size,0,samplingFreq);
		for (int ii=delay[i]; ii<=size+delay[i]; ++ii) h1bis->SetBinContent(ii,h1->GetBinContent(ii));
		hm = h1bis->FFT(hm, "MAG");

		// Crate and fill the output histogram
		TH1F *hmag = new TH1F("hmag","FFT MAG",size,0,samplingFreq);
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

	delete fft_input, re, im, mag;

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

// ------------------------
// analyze counter data set
// ------------------------
void inlineHelp() {
	cout << "[usage]: ./bin/analyze -F    <input file>" << endl;
	cout << "                       -S    <board number>" << endl;
	cout << "                       -B    <burst count>" << endl;
	cout << "                       -T    <trigger count>" << endl;
	cout << "                       -E    <channel type>" << endl;
	cout << "                       -P    <test pattern #1> <test pattern #2>" << endl;
	cout << "                       -L[F] <output label>" << endl;
	cout << endl;
}

