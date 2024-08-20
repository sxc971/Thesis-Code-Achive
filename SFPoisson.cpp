void prob_2nd_sum(int thread_num,int shave1,int bin,int start,double *array,TH1D *Hist)
{
    int ceiling = bin;
    int shave2 = 0;
    int max_track = -1*thread_num+start;

    for (int a = start; a <= ceiling; a+= thread_num)
    {
    max_track += thread_num;
	double Sum_Track = 0;
        for (int b = shave1; b < max_track + 1 - shave1; b++)
        {
            int z = max_track - b ;
			Sum_Track += array[b]*array[z];
        }
	Hist -> SetBinContent(a,Sum_Track);
    }
}

void prob_mix_sum(int thread_num,int shave1,int bin,int start,double *array0,double *array1,TH1D *Hist)
{
    int ceiling = bin;
    int shave2 = 0;
    int max_track = -1*thread_num+start;

    for (int a = start; a <= ceiling; a+= thread_num)
    {
    max_track += thread_num;
	double Sum_Track = 0;
        for (int b = shave1; b < max_track + 1 - shave1; b++)
        {
            int z = max_track - b ;
			Sum_Track += array0[b]*array1[z];
        }
	Hist -> SetBinContent(a,Sum_Track);
    }
}

int HistToText(TH1D* hdiff1, TH1D* Pomme_Plus, TH1D* k0Integral, TH1D* k1Integral,TH1D* k2Integral,TH1D* k3Integral, TString fname) {
	fname.ReplaceAll("_FittedResults.root","");
	fname.Append("Results.txt");
	TFile f(fname,"new");
	// Open the output text file
	std::ofstream outfile(fname);

	// Write the header line
	outfile << "bin center; hdiff1; modified theory; original theory" << std::endl;

	// Loop over the bins of the histograms
	for (int i=1; i<=hdiff1->GetNbinsX(); i++) {

		// Get the bin center
		double bin_center = hdiff1->GetBinCenter(i);

		// Get the bin content for each histogram
		double hist1_content = hdiff1->GetBinContent(i);
		double hist2_content = Pomme_Plus->GetBinContent(i);
		double hist3_content = k0Integral->GetBinContent(i);
		double hist4_content = k1Integral->GetBinContent(i);
		double hist5_content = k2Integral->GetBinContent(i);
		double hist6_content = k3Integral->GetBinContent(i);

		// Write the data to the output file
		outfile << bin_center << ";" << hist1_content << ";" << hist2_content  << ";" << hist3_content << ";" << hist4_content << ";" << hist5_content << ";" << hist6_content << std::endl;
	}

	// Close the output file and the ROOT file
	outfile.close();
	f.Close();

	return 0;
}

int HistToText(TH1D* hdiff1, TH1D* Pomme_Plus,TString fname) {
	fname.ReplaceAll("_FittedResults.root","");
	fname.Append("Results.txt");
	TFile f(fname,"new");
	// Open the output text file
	std::ofstream outfile(fname);

	// Write the header line
	outfile << "bin center; hdiff1; Pomme_Plus" << std::endl;

	// Loop over the bins of the histograms
	for (int i=1; i<=hdiff1->GetNbinsX(); i++) {

		// Get the bin center
		double bin_center = hdiff1->GetBinCenter(i);

		// Get the bin content for each histogram
		double hdiff1_content = hdiff1->GetBinContent(i);
		double hdiff2_content = Pomme_Plus->GetBinContent(i);

		// Write the data to the output file
		outfile << bin_center << ";" << hdiff1_content << ";" << hdiff2_content << std::endl;
	}

	// Close the output file and the ROOT file
	outfile.close();
	f.Close();

	return 0;
}

double fitU(double_t val)
{
    if (val < 0) {return 0;}
    else {return 1;}
}

double Cin(double i, double n)
{
    double out = TMath::Factorial(n)/(TMath::Factorial(i)*TMath::Factorial(n-i));
    return out;
}

double_t fitPoisson(double_t *x,double_t *par)
{
    double_t t	    = x[0]/par[2];
    double_t Scalar = par[0];
    double_t p	    = par[1]/par[2];
    double_t Tw	    = par[2]/par[2];
    double_t PET    = par[3]/par[2];
    Int_t   iter    = par[4];
    double_t output = 0;
    double_t P1	    = 0;
    double_t P2	    = 0;
    double_t P3	    = 0;
    double_t P4	    = 0;
    double_t P5	    = 0;
    double_t F1	    = 0;
    double_t F2	    = 0;

    for (int j = 1; j <= iter; j++)
    {
	P1 = pow(-1,j+1)*pow(p,j)*exp(-j*p*Tw);
	P2 = fitU(t-j*Tw)*pow(t-j*Tw,j-1)*exp(-p*(t-j*Tw))/TMath::Factorial(j-1);
	for (int a = 1; a <= j; a++)
	{
	    P3 = fitU(t- j*Tw - a*PET)*Cin(a,j)*exp(-a*p*PET);
	    for (int n = 1; n <= a; n++)
	    {
			P4 = P4 + (pow(-1,a-n)*pow(t-j*Tw-a*PET,n-1)/TMath::Factorial(n-1))*(1/pow(p,j-n))*Cin(j-1,a-n+j-1);
	    }
	    	for (int n = 1; n <= j; n++)
	    {
			P5 = P5 +(pow(t-j*Tw-a*PET,n-1)/TMath::Factorial(n-1))*(1/pow(p,j-n))*(Cin(a-1,j-n+a-1)*exp(-p*(t-j*Tw-a*PET)));
	    }
	    F1 = F1 + P3*(P4+pow(-1,a)*P5);
	    P4 = 0;
	    P5 = 0;
	    //cout << F1 << endl;
	}
	F2 = F2 + P1*(P2 + F1);
	F1 = 0;
    }
    output = Scalar*F2;
    return output;
}

TH1D* ConvertToLogScale(TH1D* hist) {
    TH1D* new_hist = dynamic_cast<TH1D*>(hist->Clone("logScaleHist"));
    int numBins = hist->GetNbinsX();

    for (int i = 1; i <= numBins; ++i) {
        double binContent = hist->GetBinContent(i);
        
        // Skip bins with zero entries to avoid log(0) undefined
        if (binContent > 0) {
            double logValue = TMath::Log(binContent);
            new_hist->SetBinContent(i, logValue);
        }
    }
    return new_hist;
}

double_t Func_Angle(double_t *x,double_t *par)
{
    if (x[0]-par[1] < 0) {return 0;}
    else 
    {
        double_t a;
        a = par[0]*x[0]*par[2]+par[0]*par[3];    
        if (x[0] > par[4]){a = par[0]*par[4]*par[2]+par[0]*par[3];}
        return a;
    }
}

double_t HDerivative(TString filename)
{
    TFile *f1	= new TFile(filename.Data());
    TH1D* hist = (TH1D*)f1->Get("hdiff1");
	double_t Deadtime = hist->FindFirstBinAbove(10,1,1,-1)*hist->GetXaxis()->GetBinWidth(1);
    TH1D* hdiff = ConvertToLogScale(hist);
//    hdiff -> Draw();

    double BinNum = hdiff->GetNbinsX();
    double BinWidth = hdiff -> GetXaxis()->GetBinWidth(1);

    double Flat  = Deadtime * 1.5;
    TF1 *Angle = new TF1("Angle",Func_Angle,0,Deadtime*2,5);
    Angle -> SetParameters(10,Deadtime,-1000,5,Deadtime*2);
    Angle -> FixParameter(1,Deadtime);
    Angle -> SetParLimits(4,2*Deadtime,2*Deadtime);
    TFitResultPtr r0 = hdiff -> Fit(Angle,"SRQ0","",Deadtime,Deadtime*1.5);

    for (double i = hdiff->FindFirstBinAbove(2,1,1,-1); i < 2*hist->FindFirstBinAbove(2,1,1,-1); i++)
    {
        double_t a = r0 -> Parameter(0) * (BinWidth * i * r0 -> Parameter(2) + r0 -> Parameter(3));
        double_t diff = a - hdiff -> GetBinContent(i);
        if (diff < -0.1)
        {
            Flat = hdiff -> GetBinCenter(i);
            break;
        }
    }
    Angle -> SetParLimits(2,r0->Parameter(2)*0.98,r0->Parameter(2)*1.02);
    Angle -> SetParLimits(3,r0->Parameter(3)*0.98,r0->Parameter(3)*1.02);
    Angle -> SetParLimits(4,Flat*0.80,Flat*1.2);
    r0 = hdiff -> Fit(Angle,"SRQ","",0,Deadtime*2);
    Flat = r0 -> Parameter(4);
    return Flat;
}

double Chi_Test_HistHist(TH1D* hdiff1, TH1D* HistExpect,double UpTo)
{
	// Perform the Chi-square test weighted		
	double Total = hdiff1 -> GetEntries();
	int FirstBin = hdiff1-> FindFirstBinAbove(10,1,1-1);
	int UpTo98 = 1; 
	double Moving_Total = 0;
	for (int i = 1; i <= hdiff1->GetNbinsX(); i++) 
	{
		Moving_Total += hdiff1->GetBinContent(i);
		if (Moving_Total < UpTo * Total) 
		{
			UpTo98 = i; // Update the last bin index
		}
	}

	double Chi2_W = 0;
	double Total_W = HistExpect->Integral(FirstBin,UpTo98);
	double Avg_W = Total_W/(UpTo98-FirstBin);
	double Hold_A = 0;
	double Hold_B = 0;
	for (int i = FirstBin; i < UpTo98; i++)
	{
		Hold_A = hdiff1->GetBinContent(i);
		Hold_B = HistExpect->GetBinContent(i);
		if (Hold_B < 1)
		{
			// the model histogram didnt extend to this region, so its time to end the cycle a little early, few times have to be calculated
			Total_W = HistExpect->Integral(FirstBin,i);
			Avg_W = Total_W/(i-FirstBin);
			UpTo98 = i;
			break;
		}
		Chi2_W += (TMath::Power(Hold_A - Hold_B,2)/Hold_B)*(Hold_A/Avg_W);
	}
	
	double Chi2_W_NDF = Chi2_W / (UpTo98 - FirstBin);
	cout << fixed << setprecision(2) <<"Non-Normalized Chi-square value :" << Chi2_W_NDF << endl;
	return Chi2_W_NDF;
}


void SFPoisson(int thread_num, int MODA, int nano_time, const char *ext="Hist.root")
{
    TStopwatch t;
    t.Start();
    TString dirname = gSystem -> WorkingDirectory();
	freopen("Terminal.txt", "w", stdout);
    TCanvas *cs = new TCanvas("cs","cs",1920,1080);
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files) {
	TSystemFile *file;
	TString fname;
	TIter next(files);
	while ((file=(TSystemFile*)next())) {
	    fname = file->GetName();
	    if (!file->IsDirectory() && fname.EndsWith(ext)) {
		cout << fname.Data() << endl;

		TFile *f1 = new TFile(fname.Data());
		TH1D* hdiff1 = (TH1D*)f1->Get("hdiff1");
		TH1D* hdiff2 = (TH1D*)f1->Get("hdiff2");

		//rebin until a peak of 40,000 is reached, this is to prevent too much uncertainty
		while (hdiff1->GetMaximum() < 4000) 
		{
			hdiff1->Rebin(2);
			hdiff2->Rebin(2);
		}

		double_t Flat = HDerivative(fname.Data());
//		double_t Flat = 15.25e6;
		TAxis* axis = hdiff1 -> GetXaxis();
		double nbin = hdiff1->GetNbinsX();
		double cwidth = axis->GetBinWidth(1)*(nbin);
    	double BinWidth = axis->GetBinWidth(1);
//		cout << Flat << endl;		
		
		double secsum = 0;
		for (int i = 1; i <= nbin; i++) {
			double freq = hdiff1->GetBinContent(i);
			double width = hdiff1->GetBinWidth(i);
			secsum += freq * width*i;
		}
		double hdiff1integral = hdiff1->Integral(1, -1);
		double MCR;
		if (nano_time == 0){MCR = hdiff1integral / (secsum/1e12);}
		else {MCR = hdiff1integral / (secsum/1e9);}

		double_t Scalar	=   hdiff1->GetMaximum()*3;
		double_t Tw	    =	hdiff1->FindFirstBinAbove(10,1,1,-1)*axis->GetBinWidth(1);
		double_t PET	=	Flat - Tw;
		double_t p	    =	0;
		int	 iter	    =   20;

		double Total = hdiff1 -> GetEntries();
		int FirstBin = hdiff1-> FindFirstBinAbove(10,1,1-1);
		int FinalBin = iter*hdiff1->FindFirstBinAbove(10,1,1-1);
		
		int UpTo9999 = 1; 
		double Moving_Total = 0;
		for (int i = 1; i <= hdiff1->GetNbinsX(); i++) 
		{
			Moving_Total += hdiff1->GetBinContent(i);
			if (Moving_Total < 0.9999 * Total) 
			{
				UpTo9999 = i; // Update the last bin index
			}
		}
		
		TF1 *FitStep = new TF1("FitStep", "[A]",Flat,Tw*2);
		FitStep -> SetParameters(Scalar,10000);
		TFitResultPtr rS = hdiff1 -> Fit(FitStep,"SRQ0");

		double intercept;
		if (nano_time == 0)
		{
			cout << "initial curve [Tw,Tw+PET] fit" << endl;
			TF1 *FitP = new TF1("FitP", "[0]*exp(-x*[1]/1e12)",Tw,Tw*1.5);
			FitP -> SetParameters(Scalar,10000);
			TFitResultPtr r0 = hdiff1 -> Fit(FitP,"SR0");
			p = (r0->Parameter(1)/(1e12/Tw))*Tw;
			intercept = 1e12/r0->Parameter(1) * TMath::Log(r0->Parameter(0)/rS->Parameter(0));	
		}
		else if (nano_time == 1)
		{
			cout << "initial curve [Tw,Tw+PET] fit" << endl;
			TF1 *FitP = new TF1("FitP", "[0]*exp(-x*[1]/1e9)",Tw,Tw*1.5);
			FitP -> SetParameters(Scalar,10000);
			TFitResultPtr r0 = hdiff1 -> Fit(FitP,"SR0");
			p = (r0->Parameter(1)/(1e9/Tw))*Tw;
			intercept = 1e9/r0->Parameter(1) * TMath::Log(r0->Parameter(0)/rS->Parameter(0));	
		}

		TF1 *Fit1 = new TF1("Fit1",fitPoisson,0,iter*Tw,5);
		Fit1 -> SetParameters(Scalar,p,Tw,PET,iter);
		Fit1 -> SetParNames("Scalar","Timing","Resolving Time","Peak Evo. Time","Fitting Iterations");
		Fit1 -> SetNpx(iter*20);
		if (MODA > 0)
		{
			Fit1 -> FixParameter(1,p);
		}
		//Fit1 -> FixParameter(0,412340);
		//Fit1 -> FixParameter(1,12800/1e12*Tw*Tw);
		Fit1 -> FixParameter(2,Tw);
		Fit1 -> FixParameter(3,intercept-Tw);
		//Fit1 -> FixParameter(3,7.71e6);
		Fit1 -> FixParameter(4,iter);

		fname.ReplaceAll("_Flag0_Hist","");
		fname.ReplaceAll(".root","");
		fname.Append("_FittedResults.root");
		TFile f(fname,"new");
		//Always save these  
		hdiff1			-> SetLineColor(kRed);
		hdiff1 			-> Draw();
		hdiff1 			-> Write();
		hdiff2 			-> Write();
		Fit1 			-> Write();

//		TFitResultPtr r1 = hdiff1 -> Fit(Fit1,"MSRW");
		if (MODA == 0)
		{
			TFitResultPtr r1 = hdiff1 -> Fit(Fit1,"MSRW0");
						
			//make a fake array that is based off of the Fit1 function		
			double *hdiffarray = new double[nbin];
			for (int i = 0; i < nbin; i++)
			{
				hdiffarray[i] = 0;
			}

			for (int i = 0; i < FinalBin; i++)
			{
				hdiffarray[i] = Fit1 -> Eval(i*BinWidth);
			}

			//turn that array into a histogram
			TH1D *HextendFit1 = new TH1D("HextendFit1","",nbin,0,cwidth);
			for (int i = 0; i < FinalBin; i++)
			{
				HextendFit1 -> SetBinContent(i,hdiffarray[i]);
			}

			//find the exponential fit from 10*Tw to 20*Tw, and expand the fake array and histogram to be based on that exponential
			TF1 *FextendFit1;
			if (nano_time == 0){	FextendFit1 = new TF1("FextendFit1","[0]*exp(-x*[1]/1e12)", 10*Tw, 20*Tw);	}
			else {FextendFit1 = new TF1("FextendFit1","[0]*exp(-x*[1]/1e9)", 10*Tw, 20*Tw);	}
			FextendFit1->SetParameter(Scalar, 10000);
			TFitResultPtr r2 = HextendFit1 -> Fit(FextendFit1,"SQR0");
			for (int i = FinalBin; i < nbin; i++)
			{
				double Fextender = FextendFit1 -> Eval(i*BinWidth);
				HextendFit1 -> SetBinContent(i,Fextender);
				hdiffarray[i] = HextendFit1->GetBinContent(i);
			}
			HextendFit1		-> SetLineColor(kBlack);
			HextendFit1 	-> Draw("same");
			HextendFit1 	-> Write();
			cs				-> Write();
		
			cout << "Estimated ICR		:" <<  (r1->Parameter(1)/r1->Parameter(2))*(1e12/r1->Parameter(2))  << endl;
			cout << fixed << setprecision(1)<< "Estimated Deadtime :"<< (1 - MCR/ ((r1->Parameter(1)/r1->Parameter(2))*(1e12/r1->Parameter(2))))*100 << "%" << endl;
			Chi_Test_HistHist(hdiff1,HextendFit1,0.98);
		}
		else if (MODA == 1)
		{
			TFitResultPtr r1 = hdiff1 -> Fit(Fit1,"MSW0","L",Tw,Flat);

			//make a fake array that is based off of the Fit1 function		
			double *hdiffarray = new double[nbin];
			for (int i = 0; i < nbin; i++)
			{
				hdiffarray[i] = 0;
			}

			for (int i = 0; i < FinalBin; i++)
			{
				hdiffarray[i] = Fit1 -> Eval(i*BinWidth);
			}

			//turn that array into a histogram
			TH1D *HextendFit1 = new TH1D("HextendFit1","",nbin,0,cwidth);
			for (int i = 0; i < FinalBin; i++)
			{
				HextendFit1 -> SetBinContent(i,hdiffarray[i]);
			}

			//find the exponential fit from 10*Tw to 20*Tw, and expand the fake array and histogram to be based on that exponential
			TF1 *FextendFit1;
			if (nano_time == 0){	FextendFit1 = new TF1("FextendFit1","[0]*exp(-x*[1]/1e12)", 10*Tw, 20*Tw);	}
			else{FextendFit1 = new TF1("FextendFit1","[0]*exp(-x*[1]/1e9)", 10*Tw, 20*Tw);	}
			FextendFit1->SetParameter(Scalar, 10000);
			TFitResultPtr r2 = HextendFit1 -> Fit(FextendFit1,"SQR0");
			for (int i = FinalBin; i < nbin; i++)
			{
				double Fextender = FextendFit1 -> Eval(i*BinWidth);
				HextendFit1 -> SetBinContent(i,Fextender);
				hdiffarray[i] = HextendFit1->GetBinContent(i);
			}

			//perform the convolution between points for hdiff1 using the expanded fake array and create a replication for pomme's theory for second order "HORDER2"
			TH1D *HOrder2 = new TH1D("HOrder2","",nbin,0,cwidth);

			vector<thread> threads(thread_num);
			
			for (int i = 0; i < thread_num; i++)
			{
				threads[i] = thread(prob_2nd_sum,thread_num,FirstBin,UpTo9999,i,hdiffarray,HOrder2);
			}

			// Wait for all the threads to finish
			for (auto& th : threads) 
			{
				th.join();
			}

			//Create a histogram that shows the difference between theory and results
			TH1D *res1 = new TH1D("res1","",nbin,0,cwidth);
			res1 -> Add(hdiff1,1);
			res1 -> Add(HextendFit1,-1);

			res1 -> Scale(-1.0);
			res1 -> SetMinimum(0);
			//The HOrder2 distribution is just a probability distribution right now, a scalar multiplier (t1/t2) is included to make it comparable to the actual hdiff2 distribution
			double t1 = hdiff2   -> Integral(FirstBin*2,FinalBin);
			double t2 = HOrder2  -> Integral(FirstBin*2,FinalBin);
			HOrder2 -> Add(HOrder2,(t1/t2)-1);
			HOrder2 -> SetMinimum(0);
			HextendFit1 -> SetMinimum(0);

			//create a new histogram that is based on HORDER2 with changing sizes to match the residual differences
			TH1D *Fres1 = new TH1D("Fres1","",nbin,0,cwidth);
			t1 = res1       -> Integral(FirstBin*2,FinalBin);
			t2 = HOrder2    -> Integral(FirstBin*2,FinalBin);
			Fres1 -> Add(HOrder2,t1/t2);
			Fres1 -> SetMinimum(0);
			double HOrder2CR = Fres1 -> Integral(FirstBin,FinalBin) / (HextendFit1->Integral(FirstBin,FinalBin) + Fres1->Integral(FirstBin,FinalBin));

			cout << "Convolution Factor  : " << HOrder2CR*100 << "%" << endl;

			//adjust the original theory to the now HORDER2 matched to residual values histogram
			TH1D *ModHist1 = new TH1D("ModHist1","",nbin,0,cwidth);
			ModHist1 -> Add(HextendFit1);
			ModHist1 -> Add(Fres1,-1);
			ModHist1 -> SetMinimum(0);
			
			ModHist1        -> SetLineColor(kBlue);
			ModHist1		-> Draw("same");
			HextendFit1 	-> Write();
			HOrder2 		-> Write();
			res1 			-> Write();
			Fres1 			-> Write();
			ModHist1 		-> Write();
			cs -> Write();

			if (nano_time == 0)
			{
				cout << "Estimated ICR		:" << (r1->Parameter(1)/r1->Parameter(2))*(1e12/r1->Parameter(2)) * (1 + HOrder2CR) << endl;
				cout << fixed << setprecision(1) << "Estimated Deadtime :" << (1 - (MCR/((r1->Parameter(1)/r1->Parameter(2))*(1e12/r1->Parameter(2)) * (1 + HOrder2CR))))*100 << "%" << endl;
			}
			if (nano_time == 1)
			{
				cout << "Estimated ICR		:" << (r1->Parameter(1)/r1->Parameter(2))*(1e9/r1->Parameter(2)) * (1 + HOrder2CR) << endl;
				cout << fixed << setprecision(1) << "Estimated Deadtime :" << (1 - (MCR/((r1->Parameter(1)/r1->Parameter(2))*(1e9/r1->Parameter(2)) * (1 + HOrder2CR))))*100 << "%" << endl;
			}
			HistToText(hdiff1,ModHist1,fname);

			cout << fixed << setprecision(2) << endl;
			Chi_Test_HistHist(hdiff1,ModHist1,0.98);
		}
		else if (MODA == 2)
		{
			TFitResultPtr r1 = hdiff1 -> Fit(Fit1,"MSW0","L",Tw,Flat);

			//make a fake array that is based off of the Fit1 function		
			double *hdiffarray = new double[nbin];
			for (int i = 0; i < nbin; i++){	hdiffarray[i] = 0;}
			for (int i = 0; i < FinalBin; i++){hdiffarray[i] = Fit1 -> Eval(i*BinWidth);}

			//turn that array into a histogram
			TH1D *HextendFit1 = new TH1D("HextendFit1","",nbin,0,cwidth);
			for (int i = 0; i < FinalBin; i++){	HextendFit1 -> SetBinContent(i,hdiffarray[i]);}

			//find the exponential fit from 10*Tw to 20*Tw, and expand the fake array and histogram to be based on that exponential
			TF1 *FextendFit1 = new TF1("FextendFit1","[0]*exp(-x*[1]/1e12)", 10*Tw, 20*Tw);
			FextendFit1->SetParameter(Scalar, 10000);
			TFitResultPtr r2 = HextendFit1 -> Fit(FextendFit1,"SQR0");
			for (int i = FinalBin; i < nbin; i++)
			{
				double Fextender = FextendFit1 -> Eval(i*BinWidth);
				HextendFit1 -> SetBinContent(i,Fextender);
				hdiffarray[i] = HextendFit1->GetBinContent(i);
			}

			auto shave1 = FirstBin;
			auto shave2 = 0;

			vector<thread> threads(thread_num);
			//perform the convolution between points for hdiff1 using the expanded fake array and create a replication for pomme's theory for second order "HORDER2"
			TH1D *HOrder2 = new TH1D("HOrder2","",nbin,0,cwidth);			
			for (int i = 0; i < thread_num; i++){threads[i] = thread(prob_2nd_sum,thread_num,FirstBin,UpTo9999,i,hdiffarray,HOrder2);}
			// Wait for all the threads to finish
			for (auto& th : threads) {th.join();}
	
			// Create a C++ array to store the histogram contents
			double *hdiff2array = new double[nbin];

			// Fill the array with bin contents from the histogram
			for (int hdiff2bin = 0; hdiff2bin < nbin; hdiff2bin++) {hdiff2array[hdiff2bin] = HOrder2->GetBinContent(hdiff2bin);}

			// create the next order distribution
			TH1D *HOrder3 = new TH1D("HOrder3","",nbin,0,cwidth);
			for (int i = 0; i < thread_num; i++){threads[i] = thread(prob_mix_sum,thread_num,FirstBin,UpTo9999,i,hdiffarray,hdiff2array,HOrder3);}
			// Wait for all the threads to finish
			for (auto& th : threads) {th.join();}

			//Create a histogram that shows the difference between theory and results
			TH1D *res1 = new TH1D("res1","",nbin,0,cwidth);
			res1 -> Add(hdiff1,1);
			res1 -> Add(HextendFit1,-1);
			res1 -> Scale(-1.0);
			res1 -> SetMinimum(0);
			//The HOrder2 distribution is just a probability distribution right now, a scalar multiplier (t1/t2) is included to make it comparable to the actual hdiff2 distribution
			double t1 = hdiff2   -> Integral(FirstBin*2,FinalBin);
			double t2 = HOrder2  -> Integral(FirstBin*2,FinalBin);
			HOrder2 -> Add(HOrder2,(t1/t2)-1);
			HOrder2 -> SetMinimum(0);
			HextendFit1 -> SetMinimum(0);

			//create a new histogram that is based on HORDER2 with changing sizes to match the residual differences
			TH1D *Fres1 = new TH1D("Fres1","",nbin,0,cwidth);
			t1 			= res1       -> GetMaximum();
			double tB 	= HOrder2    -> GetMaximum();
			t2 			= 0.5*(2*tB+3*TMath::Power(4*tB+9,0.5)+9);
			Fres1 -> Add(HOrder2,t1/t2);
			Fres1 -> SetMinimum(0);

			//Crazy stuff time, find the residual difference between res1 and HOrder2, make a res2 and use HOrder3
			TH1D *res2 = new TH1D("res2","",nbin,0,cwidth);
			res2 -> Add(Fres1,1);
			res2 -> Add(res1,-1);
			res2 -> SetMinimum(0);

			//create a new histogram that is based on HORDER2 with changing sizes to match the residual differences
			TH1D *Fres2 = new TH1D("Fres2","",nbin,0,cwidth);
			double res2_max = 0;
			double res2_maxbin = 0;
			for (int i = FirstBin*2;i < FinalBin; i++)
			{
				double nest = res2 -> GetBinContent(i);
				if ( nest > res2_max )
				{
					res2_max = nest;
					res2_maxbin = i;
				}
			}

			t1			= res2 -> Integral(2*FirstBin,-1);
			tB 			= HOrder3    -> GetMaximum();
			t2 			= 0.5*(2*tB+3*TMath::Power(4*tB+9,0.5)+9);
//			t2			= HOrder3 -> Integral(2*FirstBin,1);
			Fres2 -> Add(HOrder3,t1/t2);
			Fres2 -> SetMinimum(0);
			
			Fres1 -> Add(Fres2,-1);

			//adjust the original theory to the now HORDER2 matched to residual values histogram
			TH1D *ModHist1 = new TH1D("ModHist1","",nbin,0,cwidth);
			ModHist1 -> Add(HextendFit1);
			ModHist1 -> Add(Fres1,-1);
			
			ModHist1 -> SetMinimum(0);
			
//			HextendFit1		-> SetLineColor(kBlack);
//			HextendFit1 	-> Draw("same");
//			HOrder2			-> SetLineColor(kGreen);
//			HOrder2			-> Draw("same");
			res1			-> SetLineColor(kViolet);
			res1 			-> Draw("same");
//			res2			-> SetLineColor(kRed);
//			res2 			-> Draw("same");
			Fres1			-> SetLineColor(kGreen);
			Fres1 			-> Draw("same");
			Fres2			-> Draw("same");
			HOrder3			-> Draw("same");
			ModHist1        	-> SetLineColor(kBlue);
			ModHist1		-> Draw("same");
			HextendFit1 		-> Write();
			HOrder2 		-> Write();
			res1 			-> Write();
			res2			-> Write();
			Fres1 			-> Write();
			Fres2			-> Write();
			HOrder3			-> Write();
			ModHist1 		-> Write();
			cs -> Write();
			
			double HOrder2Con = Fres1 -> Integral(FirstBin,UpTo9999) / (HextendFit1->Integral(FirstBin,UpTo9999) + Fres1 -> Integral(FirstBin,UpTo9999) + Fres2 -> Integral(FirstBin,UpTo9999) );
			double HOrder3Con = Fres2 -> Integral(FirstBin,UpTo9999) / (HextendFit1->Integral(FirstBin,UpTo9999) + Fres1 -> Integral(FirstBin,UpTo9999) + Fres2 -> Integral(FirstBin,UpTo9999) );
			cout << "1st Order Integral				: " << HextendFit1->Integral(FirstBin,UpTo9999) << endl;
			cout << "2nd Order Integral				: " << Fres1->Integral(FirstBin,UpTo9999) << endl;
			cout << "3rd Order Integral				: " << Fres2->Integral(FirstBin,UpTo9999) << endl;
			cout << "Convolution Factor 2nd-Order 	: " << HOrder2Con*100 << "%" << endl;
			cout << "Convolution Factor 3rd-Order 	: " << HOrder3Con*100 << "%" << endl;
			double estimated_ICR = (r1->Parameter(1)/r1->Parameter(2))*(1e12/r1->Parameter(2)) * (1 + (1*HOrder2Con) + (2*HOrder3Con));
			cout << "Estimated ICR			:" << estimated_ICR << endl;
			cout << fixed << setprecision(1) << endl;
			cout << "Estimated Deadtime 		:" << (1 - MCR/estimated_ICR)*100 << "%" << endl;
			HistToText(hdiff1,ModHist1,fname);

			cout << fixed << setprecision(2) << endl;

			Chi_Test_HistHist(hdiff1,ModHist1,0.98);
		}
		else if (MODA == 3)
		{
			TFitResultPtr r1 = hdiff1 -> Fit(Fit1,"MSW0","L",Tw,Flat);

			//make a fake array that is based off of the Fit1 function		
			double *hdiffarray = new double[nbin];
			for (int i = 0; i < nbin; i++){	hdiffarray[i] = 0;	}
			for (int i = 0; i < FinalBin; i++){hdiffarray[i] = Fit1 -> Eval(i*BinWidth);}

			//turn that array into a histogram
			TH1D *HextendFit1 = new TH1D("HextendFit1","",nbin,0,cwidth);
			for (int i = 0; i < FinalBin; i++){	HextendFit1 -> SetBinContent(i,hdiffarray[i]);}

			//find the exponential fit from 10*Tw to 20*Tw, and expand the fake array and histogram to be based on that exponential
			TF1 *FextendFit1 = new TF1("FextendFit1","[0]*exp(-x*[1]/1e12)", 10*Tw, 20*Tw);
			FextendFit1->SetParameter(Scalar, 10000);
			TFitResultPtr r2 = HextendFit1 -> Fit(FextendFit1,"SQR0");
			for (int i = FinalBin; i < nbin; i++)
			{
				double Fextender = FextendFit1 -> Eval(i*BinWidth);
				HextendFit1 -> SetBinContent(i,Fextender);
				hdiffarray[i] = HextendFit1->GetBinContent(i);
			}

			auto shave1 = FirstBin;
			auto shave2 = 0;

			vector<thread> threads(thread_num);
			//perform the convolution between points for hdiff1 using the expanded fake array and create a replication for pomme's theory for second order "HORDER2"
			TH1D *HOrder2 = new TH1D("HOrder2","",nbin,0,cwidth);			
			for (int i = 0; i < thread_num; i++){threads[i] = thread(prob_2nd_sum,thread_num,FirstBin,UpTo9999,i,hdiffarray,HOrder2);}
			// Wait for all the threads to finish
			for (auto& th : threads) {th.join();}
	
			// Create a C++ array to store the histogram contents
			double *hdiff2array = new double[nbin];
			for (int hdiff2bin = 0; hdiff2bin < nbin; hdiff2bin++) {hdiff2array[hdiff2bin] = HOrder2->GetBinContent(hdiff2bin);}

			// create the next order distribution
			TH1D *HOrder3 = new TH1D("HOrder3","",nbin,0,cwidth);
			for (int i = 0; i < thread_num; i++){threads[i] = thread(prob_mix_sum,thread_num,FirstBin,UpTo9999,i,hdiffarray,hdiff2array,HOrder3);}
			// Wait for all the threads to finish
			for (auto& th : threads) {th.join();}

			// Fill the array with bin contents from the histogram
			double *hdiff3array = new double[nbin];
			for (int hdiff3bin = 0; hdiff3bin < nbin; hdiff3bin++) {hdiff3array[hdiff3bin] = HOrder3->GetBinContent(hdiff3bin);}

			TH1D *HOrder4 = new TH1D("HOrder4","",nbin,0,cwidth);
			for (int i = 0; i < thread_num; i++){threads[i] = thread(prob_mix_sum,thread_num,FirstBin,UpTo9999,i,hdiffarray,hdiff3array,HOrder4);}
			// Wait for all the threads to finish
			for (auto& th : threads) {th.join();}

			double *hdiff4array = new double[nbin];
			for (int hdiff4bin = 0; hdiff4bin < nbin; hdiff4bin++) {hdiff4array[hdiff4bin] = HOrder4->GetBinContent(hdiff4bin);}

			//Create a histogram that shows the difference between theory and results
			TH1D *res1 = new TH1D("res1","",nbin,0,cwidth);
			res1 -> Add(hdiff1,1);
			res1 -> Add(HextendFit1,-1);
			res1 -> Scale(-1.0);
			res1 -> SetMinimum(0);
			//The HOrder2 distribution is just a probability distribution right now, a scalar multiplier (t1/t2) is included to make it comparable to the actual hdiff2 distribution
			double t1 = hdiff2   -> Integral(FirstBin*2,FinalBin);
			double t2 = HOrder2  -> Integral(FirstBin*2,FinalBin);
			HOrder2 -> Add(HOrder2,(t1/t2)-1);
			HOrder2 -> SetMinimum(0);
			HextendFit1 -> SetMinimum(0);

			double step = 0;
			for (int i = 0; i < 1; i++)
			{
				double tc = 1.24e-6;
				double est_icr = 15300;
				double k0Integral = HextendFit1 -> Integral();
				double k1Integral = k0Integral * (pow(tc*est_icr,1) * exp(-tc*est_icr) / TMath::Factorial(1)) / exp(-tc*est_icr);
				double k2Integral = k0Integral * (pow(tc*est_icr,2) * exp(-tc*est_icr) / TMath::Factorial(2)) / exp(-tc*est_icr);
				double k3Integral = k0Integral * (pow(tc*est_icr,3) * exp(-tc*est_icr) / TMath::Factorial(3)) / exp(-tc*est_icr);
				HOrder2 -> Scale(k1Integral/HOrder2->Integral());
				HOrder3 -> Scale(k2Integral/HOrder3->Integral());
				HOrder4 -> Scale(k3Integral/HOrder4->Integral());
				HOrder2 -> SetMinimum(0);
				HOrder3 -> SetMinimum(0);
				HOrder4 -> SetMinimum(0);
				TH1D *ModHist1 = new TH1D("ModHist1","",nbin,0,cwidth);
				ModHist1 -> Add(HextendFit1);
				ModHist1 -> Add(HOrder2,-1);
				ModHist1 -> Add(HOrder3,1);
				ModHist1 -> Add(HOrder4,-1);
				
				hdiff1 -> Draw();
				HextendFit1 -> Draw("same");
				ModHist1 -> Draw("same");
				cs -> Write();
				
				cout << "1st Order Integral				: " << HextendFit1->Integral(FirstBin,UpTo9999) << endl;
				cout << "2nd Order Integral				: " << HOrder2->Integral(FirstBin,UpTo9999) << endl;
				cout << "3rd Order Integral				: " << HOrder3->Integral(FirstBin,UpTo9999) << endl;
				cout << "4rd Order Integral				: " << HOrder4->Integral(FirstBin,UpTo9999) << endl;
				double estimated_ICR;
				if (nano_time == 0){
					estimated_ICR = (r1->Parameter(1)/r1->Parameter(2))*(1e12/r1->Parameter(2)) / exp(-tc*est_icr);
				}
				else
				{
					estimated_ICR = (r1->Parameter(1)/r1->Parameter(2))*(1e9/r1->Parameter(2)) / exp(-tc*est_icr);
				}
				cout << "Estimated ICR			:" << estimated_ICR << endl;
				cout << fixed << setprecision(1) << endl;
				cout << "Estimated Deadtime 		:" << (1 - MCR/estimated_ICR)*100 << "%" << endl;
				HistToText(hdiff1,ModHist1,HextendFit1,HOrder2,HOrder3,HOrder4,fname);

				cout << fixed << setprecision(2) << endl;

				Chi_Test_HistHist(hdiff1,ModHist1,0.98);
			}
		}
		cout << "" << endl;
		cout << "" << endl;
	    }
	}
    }
    t.Stop();
    t.Print();
}
