double_t Func_Poisson(double_t *x,double_t *par)
{
    if (x[0]-par[0]*(par[1]+1) < 0) {return 0;}
    else 
    {
        double_t a = pow((x[0]-par[0]*(par[1]+1))/par[2],par[1])*exp(-(x[0]-par[0]*(par[1]+1))/par[2])/TMath::Factorial(par[1]);  
        return a;
    }
}

double_t Func_Poisson_PlusOne(double_t *x,double_t *par)
{
    if (x[0]-par[0]*(par[1]+1) < 0) {return 0;}
    else 
    {
        double_t new_k = par[1] + 1;
        double_t a = pow((x[0]-par[0]*(new_k))/par[2],new_k)*exp(-(x[0]-par[0]*(new_k))/par[2])/TMath::Factorial(new_k);  
        return a;
    }
}

double_t Func_Reverse(double_t *x,double_t *par)
{
    if (x[0]-par[0]*(par[1]+1) < 0) {return 0;}
    else 
    {
        double_t a = 1-pow((x[0]-par[0]*(par[1]+1))*par[3]/par[0],par[1])*exp(-(x[0]-par[0]*(par[1]+1))*par[3]/par[0])/TMath::Factorial(par[1]);
        return a;
    }
}

double_t Func_Product(double_t *x, double_t *par)
{
    double_t reverse_value = Func_Reverse(x, par);
    double_t exponential_value = Func_Poisson(x, par);
    double_t Final = reverse_value * exponential_value;
    return Final;
}


double_t Final_Conv(double_t *x, double_t *par)
{
    double_t Final = par[4]*Func_Product(x,par);
    Final += par[5]*Func_Poisson_PlusOne(x,par);
    return Final;
}

double Chi_Test_Func(TH1D* hist, TF1* func)
{
	// Perform the Chi-square test weighted
	int FirstBin = 3*hist-> FindFirstBinAbove(10,1,1-1);
    double binwidth = hist -> GetBinWidth(1);

    long long int Total = hist -> GetEntries();
	double_t UpTo99 = 1; 
	double Moving_Total = 0;
	for (int i = 1; i <= hist->GetNbinsX(); i++) 
	{
		Moving_Total += hist->GetBinContent(i);
		if (Moving_Total < 0.99 * Total) 
		{
			UpTo99 = i; // Update the last bin index
		}
	}

	double Chi2_W = 0;
	double Total_W = hist->Integral(FirstBin,UpTo99);
	double Avg_W = Total_W/(UpTo99-FirstBin);
	double Hold_A = 0;
	double Hold_B = 0;
	for (int i = FirstBin; i < UpTo99; i++)
	{
		Hold_A = hist->GetBinContent(i);
		Hold_B = func->Eval(i*binwidth);
        if (Hold_B != 0)
        {
            Chi2_W += (TMath::Power(Hold_A - Hold_B,2)/Hold_B)*(Hold_A/Avg_W);
        }
	}
	
	double Chi2_W_NDF = Chi2_W / (UpTo99 - FirstBin);
	cout << fixed << setprecision(2) <<"Non-Normalized Chi-square value :" << Chi2_W_NDF << endl;
	return Chi2_W_NDF;
}

double_t GiveMe99(TH1D* hist)
{
    long long int Total = hist -> GetEntries();
	double_t UpTo99 = 1; 
	double Moving_Total = 0;
	for (int i = 1; i <= hist->GetNbinsX(); i++) 
	{
		Moving_Total += hist->GetBinContent(i);
		if (Moving_Total < 0.999 * Total) 
		{
			UpTo99 = i; // Update the last bin index
		}
	}
    return UpTo99;
}

void SPoisson(const char *ext="Hist.root")
{
    TString dirname = gSystem->pwd();
    TStopwatch t;
    t.Start();
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
                TH1D* hdiff3 = (TH1D*)f1->Get("hdiff3");
                TH1D* hdiff4 = (TH1D*)f1->Get("hdiff4");				
                
                TAxis* axis = hdiff1 -> GetXaxis();
                double nbin = GiveMe99(hdiff1);
                double cwidth = axis->GetBinWidth(1)*GiveMe99(hdiff1);
                double BinWidth = axis->GetBinWidth(1);

                double secsum = 0;
                for (int i = 1; i <= nbin; i++) {
                    double freq = hdiff1->GetBinContent(i);
                    double width = hdiff1->GetBinWidth(i);
                    secsum += freq * width*i;
                }
                double hdiff1integral = hdiff1->Integral(1, nbin);
                double MCR = hdiff1integral / (secsum/1e12);

                double_t Tau= hdiff1->FindFirstBinAbove(10,1,1,-1)*BinWidth;
                double_t Mean = secsum / hdiff1->Integral(10,nbin);
                double_t Scalar = hdiff1 -> GetMaximum();
                double_t Scalar2 = Scalar*0.1;
                double_t Mean2 = Mean;
                double_t Z_value = 8;
                int Pcheck ((1e12/Mean)<(1e3));

                cout << "DataFile ID :" << fname << endl;
//                cout << "Histogram ID:" << "hdiff1" << endl;
//                TF1 *Poisson = new TF1("Poisson", "[0]*exp(-(x-[1])/[2])*(1-exp(-(x-[1])/(0.121*[1])))                         +           [3]*pow((x-[1])/[2],1)*exp(-(x-[1])/[2])/TMath::Factorial(1)", 0., 10.);
                TF1 *Poisson = new TF1("Poisson",Final_Conv,0.,10.,6);
                Poisson-> SetNpx(cwidth*0.1 / BinWidth);
                Poisson-> SetParameters(Tau,0,Mean,Z_value,Scalar,Scalar2);
                Poisson-> SetParNames("Tau","k-value","Mean","Z-Value","Scalar","Scalar_2");
                Poisson-> SetRange(Tau, cwidth);
                Poisson-> FixParameter(0,Tau);
                Poisson-> FixParameter(1,0);
                Poisson-> FixParameter(3,Z_value);
                Poisson-> SetParLimits(2,0.1*Mean,10*Mean);
                Poisson-> SetParLimits(4,0.1*Scalar,10*Scalar);
                Poisson-> SetParLimits(5,0,Scalar);
                if (Pcheck){
                    Poisson -> FixParameter(1,0);
                }
                //                Poisson->SetParLimits(1,0,Tau*2);
                //Poisson->SetParLimits(4,0,Mean2*2);
                TFitResultPtr r1 = hdiff1 -> Fit("Poisson","MSRWQ0");
//                cout << "" << endl;

//                cout << "Histogram ID:" << "hdiff2" << endl;

//                TF1 *Poisson2 = new TF1("Poisson2", "[0]*(((x-[1]*2))/[2])*exp(-((x-[1]*2))/[2])/TMath::Factorial(1)           +           [3]*pow((x-[1]*2)/[2],2)*exp(-((x-[1]*2))/[2])/TMath::Factorial(2)", 0., 10.);
                TF1 *Poisson2 = new TF1("Poisson2",Final_Conv,0.,10.,6);
                Poisson2-> SetNpx(cwidth*0.1 / BinWidth);
                Poisson2-> SetParameters(Tau,1,Mean,Z_value,Scalar,Scalar2);
                Poisson2-> SetParNames("Tau","k-value","Mean","Z-Value","Scalar","Scalar_2");
                Poisson2-> SetRange(2*Tau, cwidth);
                Poisson2-> FixParameter(0,Tau);
                Poisson2-> FixParameter(1,1);
                Poisson2-> FixParameter(3,Z_value);
                Poisson2-> SetParLimits(2,0.1*Mean,10*Mean);
                Poisson2-> SetParLimits(4,0.1*Scalar,10*Scalar);
                Poisson2-> SetParLimits(5,0,Scalar);
                if (Pcheck){
                    Poisson2 -> FixParameter(1,0);
                }
                TFitResultPtr r2 = hdiff2 -> Fit("Poisson2","MSRWQ0");
//                cout << "" << endl;

//                cout << "Histogram ID:" << "hdiff3" << endl;

//                TF1 *Poisson3 = new TF1("Poisson3", "[0]*pow((((x-[1]*3))/[2]),2)*exp(-((x-[1]*3))/[2])/TMath::Factorial(2)     +           [3]*pow((x-[1]*3)/[2],3)*exp(-((x-[1]*3))/[2])/TMath::Factorial(3)", 0., 10.);
                TF1 *Poisson3 = new TF1("Poisson3",Final_Conv,0.,10.,6);
                Poisson3-> SetNpx(cwidth*0.1 / BinWidth);
                Poisson3-> SetParameters(Tau,2,Mean,Z_value,Scalar,Scalar2);
                Poisson3-> SetParNames("Tau","k-value","Mean","Z-Value","Scalar","Scalar_2");
                Poisson3-> SetRange(3*Tau, cwidth);
                Poisson3-> FixParameter(0,Tau);
                Poisson3-> FixParameter(1,2);
                Poisson3-> FixParameter(3,Z_value);
                Poisson3-> SetParLimits(2,0.1*Mean,10*Mean);
                Poisson3-> SetParLimits(4,0.1*Scalar,10*Scalar);
                Poisson3-> SetParLimits(5,0,Scalar);
                if (Pcheck){
                    Poisson3 -> FixParameter(1,0);
                }
                TFitResultPtr r3 = hdiff3 -> Fit("Poisson3","MSRWQ0");
//                cout << "" << endl;

//                cout << "Histogram ID:" << "hdiff4" << endl;

//                TF1 *Poisson4 = new TF1("Poisson4", "[0]*pow((((x-[1]*4))/[2]),3)*exp(-((x-[1]*4))/[2])/TMath::Factorial(3)     +           [3]*pow((x-[1]*4)/[2],4)*exp(-((x-[1]*4))/[2])/TMath::Factorial(4)", 0., 10.);
                TF1 *Poisson4 = new TF1("Poisson4",Final_Conv,0.,10.,6);
                Poisson4-> SetNpx(cwidth*0.1 / BinWidth);
                Poisson4-> SetParameters(Tau,3,Mean,Z_value,Scalar,Scalar2);
                Poisson4-> SetParNames("Tau","k-value","Mean","Z-Value","Scalar","Scalar_2");
                Poisson4-> SetRange(4*Tau, cwidth);
                Poisson4-> FixParameter(0,Tau);
                Poisson4-> FixParameter(1,3);
                Poisson4-> FixParameter(3,Z_value);
                Poisson4-> SetParLimits(2,0.1*Mean,10*Mean);
                Poisson4-> SetParLimits(4,0.1*Scalar,10*Scalar);
                Poisson4-> SetParLimits(5,0,Scalar);
                if (Pcheck){
                    Poisson4 -> FixParameter(1,0);
                }
                TFitResultPtr r4 = hdiff4 -> Fit("Poisson4","MSRWQ0");
//                cout << "" << endl;

                // Unconvoluted models
                TF1 *Poisson1_Unconv = new TF1("Poisson1_Unconv", Final_Conv, 0., 10.,6);
                Poisson1_Unconv->SetNpx(cwidth*0.1 / BinWidth);
                Poisson1_Unconv->SetParameters(Tau,0,Mean,Z_value,Scalar,Scalar2);
                Poisson1_Unconv->SetParNames("Tau","k-value","Mean","Z-Value","Scalar","Scalar2");
                Poisson1_Unconv->SetRange(1*Tau,cwidth);
                Poisson1_Unconv->FixParameter(0,Tau);
                Poisson1_Unconv->FixParameter(1,0);
                Poisson1_Unconv->FixParameter(3,Z_value);
                Poisson1_Unconv->FixParameter(5,0);
                TFitResultPtr r1_U = hdiff1 -> Fit("Poisson1_Unconv","MSRWQ0");

                TF1 *Poisson2_Unconv = new TF1("Poisson2_Unconv", Final_Conv, 0., 10.,6);
                Poisson2_Unconv->SetNpx(cwidth*0.1 / BinWidth);
                Poisson2_Unconv->SetParameters(Tau,1,Mean,Z_value,Scalar,Scalar2);
                Poisson2_Unconv->SetParNames("Tau","k-value","Mean","Z-Value","Scalar","Scalar2");
                Poisson2_Unconv->SetRange(2*Tau,cwidth);
                Poisson2_Unconv->FixParameter(0,Tau);
                Poisson2_Unconv->FixParameter(1,1);
                Poisson2_Unconv->FixParameter(3,Z_value);
                Poisson2_Unconv->FixParameter(5,0);
                TFitResultPtr r2_U = hdiff2 -> Fit("Poisson2_Unconv","MSRWQ0");

                TF1 *Poisson3_Unconv = new TF1("Poisson3_Unconv", Final_Conv, 0., 10.,6);
                Poisson3_Unconv->SetNpx(cwidth*0.1 / BinWidth);
                Poisson3_Unconv->SetParameters(Tau,2,Mean,Z_value,Scalar,Scalar2);
                Poisson3_Unconv->SetParNames("Tau","k-value","Mean","Z-Value","Scalar","Scalar2");
                Poisson3_Unconv->SetRange(3*Tau,cwidth);
                Poisson3_Unconv->FixParameter(0,Tau);
                Poisson3_Unconv->FixParameter(1,2);
                Poisson3_Unconv->FixParameter(3,Z_value);
                Poisson3_Unconv->FixParameter(5,0);
                TFitResultPtr r3_U = hdiff3 -> Fit("Poisson3_Unconv","MSRWQ0");

                TF1 *Poisson4_Unconv = new TF1("Poisson4_Unconv", Final_Conv, 0., 10.,6);
                Poisson4_Unconv->SetNpx(cwidth*0.1 / BinWidth);
                Poisson4_Unconv->SetParameters(Tau,3,Mean,Z_value,Scalar,Scalar2);
                Poisson4_Unconv->SetParNames("Tau","k-value","Mean","Z-Value","Scalar","Scalar2");
                Poisson4_Unconv->SetRange(4*Tau,cwidth);
                Poisson4_Unconv->FixParameter(0,Tau);
                Poisson4_Unconv->FixParameter(1,3);
                Poisson4_Unconv->FixParameter(3,Z_value);
                Poisson4_Unconv->FixParameter(5,0);
                TFitResultPtr r4_U = hdiff4 -> Fit("Poisson4_Unconv","MSRWQ0");


                double t_mean = (r1->Parameter(2) + r2->Parameter(2) + r3->Parameter(2) + r4->Parameter(2)) / 4;
                double t_pileup = (r1->Parameter(5)/(r1->Parameter(4)+r1->Parameter(5))  +   r2->Parameter(5)/(r2->Parameter(4)+r2->Parameter(5))/2  +   r3->Parameter(5)/(r3->Parameter(4)+r3->Parameter(5))/3    +         r4->Parameter(5)/(r4->Parameter(4)+r4->Parameter(5))/4)/4;

                cout << "Entry Rate(cps)		:" << 1e12/Mean << endl;
                cout << "Poisson Rate(cps)	        :" << 1e12/t_mean << endl;
                cout << "Est. Deadtime		        :" << (1-t_mean/Mean)*100 << "%" << endl;
                cout << "Est Conv                       :" << t_pileup*100 << "%" << endl;
                cout << "Est. True Rate(cps)            :" << (1e12/t_mean)*(1+t_pileup) << endl;
                cout << "k = 0 uICR : " << (1e12/r1_U->Parameter(2)) << "   ± " << (r1_U -> ParError(2)/r1_U->Parameter(2)*(1e12/r1_U->Parameter(2))) << endl;
                cout << "k = 1 uICR : " << (1e12/r2_U->Parameter(2)) << "   ± " << (r2_U -> ParError(2)/r2_U->Parameter(2)*(1e12/r2_U->Parameter(2))) << endl;
                cout << "k = 2 uICR : " << (1e12/r3_U->Parameter(2)) << "   ± " << (r3_U -> ParError(2)/r3_U->Parameter(2)*(1e12/r3_U->Parameter(2))) << endl;
                cout << "k = 3 uICR : " << (1e12/r4_U->Parameter(2)) << "   ± " << (r4_U -> ParError(2)/r4_U->Parameter(2)*(1e12/r4_U->Parameter(2))) << endl;
                cout << "k = 0 Conv : " << (r1->Parameter(5)/(r1->Parameter(4)+r1->Parameter(5)))*100 << "%" << endl;
                cout << "k = 1 Conv : " << (r2->Parameter(5)/(r2->Parameter(4)+r2->Parameter(5)))*100 << "%" << endl;
                cout << "k = 2 Conv : " << (r3->Parameter(5)/(r3->Parameter(4)+r3->Parameter(5)))*100 << "%" << endl;
                cout << "k = 3 Conv : " << (r4->Parameter(5)/(r4->Parameter(4)+r4->Parameter(5)))*100 << "%" << endl;
                cout << "k = 0 ICR  : " << (1e12/r1->Parameter(2))*(1+(r1->Parameter(5)/(r1->Parameter(4)+r1->Parameter(5)))/1) << "   ± " << sqrt(pow(r1 -> ParError(2)/r1->Parameter(2),2) + pow(r1 -> ParError(4)/r1->Parameter(4),2) + pow(r1 -> ParError(5)/r1->Parameter(5),2))*(1e12/r1->Parameter(2)*(1+(r1->Parameter(5)/(r1->Parameter(4)+r1->Parameter(5)))/1)) << endl;
                cout << "k = 1 ICR  : " << (1e12/r2->Parameter(2))*(1+(r2->Parameter(5)/(r2->Parameter(4)+r2->Parameter(5)))/2) << "   ± " << sqrt(pow(r2 -> ParError(2)/r2->Parameter(2),2) + pow(r2 -> ParError(4)/r2->Parameter(4),2) + pow(r2 -> ParError(5)/r2->Parameter(5),2))*(1e12/r2->Parameter(2)*(1+(r2->Parameter(5)/(r2->Parameter(4)+r2->Parameter(5)))/2)) << endl;
                cout << "k = 2 ICR  : " << (1e12/r3->Parameter(2))*(1+(r3->Parameter(5)/(r3->Parameter(4)+r3->Parameter(5)))/3) << "   ± " << sqrt(pow(r3 -> ParError(2)/r3->Parameter(2),2) + pow(r3 -> ParError(4)/r3->Parameter(4),2) + pow(r3 -> ParError(5)/r3->Parameter(5),2))*(1e12/r3->Parameter(2)*(1+(r3->Parameter(5)/(r3->Parameter(4)+r3->Parameter(5)))/3)) << endl;
                cout << "k = 3 ICR  : " << (1e12/r4->Parameter(2))*(1+(r4->Parameter(5)/(r4->Parameter(4)+r4->Parameter(5)))/4) << "   ± " << sqrt(pow(r4 -> ParError(2)/r4->Parameter(2),2) + pow(r4 -> ParError(4)/r4->Parameter(4),2) + pow(r4 -> ParError(5)/r4->Parameter(5),2))*(1e12/r4->Parameter(2)*(1+(r4->Parameter(5)/(r4->Parameter(4)+r4->Parameter(5)))/4)) << endl;
                cout << ((1e12/r1_U->Parameter(2)) + (1e12/r1_U->Parameter(2)) + (1e12/r1_U->Parameter(2)) + (1e12/r1_U->Parameter(2)))/4 << endl;
                cout << "" << endl;

                //make a fake array that is based off of the Fit1/2/3/4 function		
                double *hdiffarray1 = new double[nbin];
                double *hdiffarray2 = new double[nbin];
                double *hdiffarray3 = new double[nbin];
                double *hdiffarray4 = new double[nbin];

                for (int i = 0; i < nbin; i++)
                {
                    hdiffarray1[i] = Poisson -> Eval(i*BinWidth);
                    hdiffarray2[i] = Poisson2-> Eval(i*BinWidth);
                    hdiffarray3[i] = Poisson3-> Eval(i*BinWidth);
                    hdiffarray4[i] = Poisson4-> Eval(i*BinWidth);
                }

                TH1F *res1 = new TH1F("res1","",nbin,0,cwidth);
                TH1F *res2 = new TH1F("res2","",nbin,0,cwidth);
                TH1F *res3 = new TH1F("res3","",nbin,0,cwidth);
                TH1F *res4 = new TH1F("res4","",nbin,0,cwidth);
                for (Int_t i=1;i<=nbin;i++) {
                    Double_t res1diff = Poisson->Eval(hdiff1->GetBinCenter(i))-hdiff1->GetBinContent(i);
                    res1->SetBinContent(i,res1diff);
                    Double_t res2diff = Poisson2->Eval(hdiff2->GetBinCenter(i))-hdiff2->GetBinContent(i);
                    res2->SetBinContent(i,res2diff);
                    Double_t res3diff = Poisson3->Eval(hdiff3->GetBinCenter(i))-hdiff3->GetBinContent(i);
                    res3->SetBinContent(i,res3diff);
                    Double_t res4diff = Poisson4->Eval(hdiff4->GetBinCenter(i))-hdiff4->GetBinContent(i);
                    res4->SetBinContent(i,res4diff);
                }
                res1 ->Draw();
                res2 ->Draw();
                res3 ->Draw();
                res4 ->Draw();

                int bin99 = GiveMe99(hdiff1);
                cout << "UnConvoluted Chi-Square" << endl;
                Chi_Test_Func(hdiff1,Poisson1_Unconv);
                Chi_Test_Func(hdiff2,Poisson2_Unconv);
                Chi_Test_Func(hdiff3,Poisson3_Unconv);
                Chi_Test_Func(hdiff4,Poisson4_Unconv);
                cout << "Convoluted Chi-Square" << endl;
                Chi_Test_Func(hdiff1,Poisson );
                Chi_Test_Func(hdiff2,Poisson2);
                Chi_Test_Func(hdiff3,Poisson3);
                Chi_Test_Func(hdiff4,Poisson4);

                hdiff1->Draw();
                hdiff2->SetLineColor(kRed);
                hdiff2->Draw("same");
                hdiff3->SetLineColor(kGreen);
                hdiff3->Draw("same");
                hdiff4->SetLineColor(kOrange);
                hdiff4->Draw("same");
                hdiff1->SetLineWidth(5);
                hdiff2->SetLineWidth(5);
                hdiff3->SetLineWidth(5);
                hdiff4->SetLineWidth(5);
                Poisson1_Unconv->SetLineColor(kCyan);
                Poisson1_Unconv->Draw("same");
                Poisson2_Unconv->SetLineColor(kCyan);
                Poisson2_Unconv->Draw("same");
                Poisson3_Unconv->SetLineColor(kCyan);
                Poisson3_Unconv->Draw("same");
                Poisson4_Unconv->SetLineColor(kCyan);
                Poisson4_Unconv->Draw("same");
                Poisson->SetLineColor(kBlack);
                Poisson->Draw("same");
                Poisson2->SetLineColor(kBlack);
                Poisson2->Draw("same");
                Poisson3->SetLineColor(kBlack);
                Poisson3->Draw("same");
                Poisson4->SetLineColor(kBlack);
                Poisson4->Draw("same");
                hdiff1 -> SetTitle("Exp. k = 0");
                hdiff2 -> SetTitle("Exp. k = 1");
                hdiff3 -> SetTitle("Exp. k = 2");
                hdiff4 -> SetTitle("Exp. k = 3");
                Poisson -> SetTitle("Poisson Conv. k = 0");
                Poisson2 -> SetTitle("Poisson Conv. k = 1");
                Poisson3 -> SetTitle("Poisson Conv. k = 2");
                Poisson4 -> SetTitle("Poisson Conv. k = 3");
                Poisson1_Unconv -> SetTitle("Poisson Unconv. k = 0");
                Poisson2_Unconv -> SetTitle("Poisson Unconv. k = 1");
                Poisson3_Unconv -> SetTitle("Poisson Unconv. k = 2");
                Poisson4_Unconv -> SetTitle("Poisson Unconv. k = 3");
                hdiff1 -> GetXaxis() -> SetRangeUser(0,GiveMe99(hdiff1)*hdiff1->GetBinWidth(1));

                fname.ReplaceAll("DataR_CH0@DT5781_1967_","");
                fname.ReplaceAll(".root","");
                fname.Append("_Results.root");
                TFile f(fname,"RECREATE");
                hdiff1 -> Write();
                hdiff2 -> Write();
                hdiff3 -> Write();
                hdiff4 -> Write();
                Poisson -> Write();
                Poisson2 -> Write();
                Poisson3 -> Write();
                Poisson4 -> Write();
                Poisson1_Unconv -> Write();
                Poisson2_Unconv -> Write();
                Poisson3_Unconv -> Write();
                Poisson4_Unconv -> Write();
                r1 -> Write();
                r2 -> Write();
                r3 -> Write();
                r4 -> Write();
                res1 -> Write();
                res2 -> Write();
                res3 -> Write();
                res4 -> Write();
                cs -> Write();
                cs -> BuildLegend(0.75,0.75,1,1);

//                ofstream output;
//                fname.ReplaceAll(".root","");
//
//                output.open(fname + "data.txt");
//                char strng[400];
//                for (int i = 0; i < hdiff1->GetNbinsX(); i++)
//                {
//                    sprintf(strng,"   %10.4f   %10.4f   %10.4f   %10.4f   %10.4f   \n",hdiff1->GetBinCenter(i), hdiff1->GetBinContent(i), hdiff2->GetBinContent(i), hdiff3->GetBinContent(i), hdiff4->GetBinContent(i));
//                    output << strng;
//                }            
//                output.close();    
//                output.open(fname + "fit.txt");
//                for (int i = 0; i < hdiff1->GetNbinsX(); i++)
//                {
//                    sprintf(strng,"   %10.4f   %10.4f   %10.4f   %10.4f   %10.4f   \n",hdiff1->GetBinCenter(i), hdiffarray1[i], hdiffarray2[i], hdiffarray3[i], hdiffarray4[i]);
//                    output << strng;
//                }
//                output.close();    
            }
        }
    }
    t.Stop();
    t.Print();
    cout << "" << endl;
}

