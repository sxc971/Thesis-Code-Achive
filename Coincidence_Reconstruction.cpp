void Coincidence_Reconstruction(TString fname_0,TString fname_1)
{
    // prelimary data gathering and library creation component
    TStopwatch t;
    t.Start();
//    freopen("Terminal.txt", "w", stdout);
//    TCanvas *cs = new TCanvas("cs","cs",1920,1080);
//  the actual fun stuff, we assume Ch1 is the HPGe(reference channel swap ch0File = fname; to ch1File = fname and vise versa otherwise)
    
    cout << fname_0 << endl;
    cout << fname_1 << endl;

    TString ch0File = "";
    TString ch1File = "";   
    //This is how we grab two file names at the same time, kinda jank
    ch0File = fname_0;
    ch1File = fname_1;
    TFile *g0 = new TFile(ch0File.Data());
    TFile *g1 = new TFile(ch1File.Data());
    ULong64_t rt_0{0};
    ULong64_t rt_1{0};  
    UShort_t Eg_0 {0};
    UShort_t Eg_1 {0};
    TTree *kListData_0 = (TTree*)g0->Get("Data_R");         
    TTree *kListData_1 = (TTree*)g1->Get("Data_R"); 
    TString fname = fname_0;
    fname.ReplaceAll("CH0","");
    fname.ReplaceAll("CH1","");
    fname.ReplaceAll("CH","");
    fname.ReplaceAll("@","");
    fname.ReplaceAll(".root","");
    fname.Append("Coincidence.root");
    TFile *f = new TFile(fname,"RECREATE");

    // this part isn't nessesary but I found it helps with processing speed
    kListData_0->SetBranchStatus("*",0);
    kListData_0->SetBranchStatus("Timestamp",1);
    kListData_0->SetBranchStatus("Energy",1);
    kListData_0->SetBranchAddress("Energy",&Eg_0); 
    kListData_0->SetBranchAddress("Timestamp",&rt_0);       
    kListData_1->SetBranchStatus("*",0);
    kListData_1->SetBranchStatus("Timestamp",1);
    kListData_1->SetBranchStatus("Energy",1);
    kListData_1->SetBranchAddress("Energy",&Eg_1); 
    kListData_1->SetBranchAddress("Timestamp",&rt_1);

    // end unnessesary part 
    const int Total_Count_0 = kListData_0->GetEntries();
    const int Total_Count_1 = kListData_1->GetEntries();
    kListData_0->GetEntry(Total_Count_0-1);
    kListData_1->GetEntry(Total_Count_1-1);
    double_t Total_Time_0 = rt_0;
    double_t Total_Time_1 = rt_1;   
    // Loop through timestamps of Ch1 data, do a window search from Ch0 data. Keep track of a minimum starting position for Ch0 to lower processing time
    // if Ch0 timestamp exceeds Ch1(stamp) + Window end early, the variables are in pico-seconds
    double Static_Adjust = 0;
    double window = 3e6;
    TH1D *Time_Coincidence_0        = new TH1D("Time_Coincidence_0","",4*window/10000,-2*window,2*window);
    TH1D *Time_AntiCoincidence_0    = new TH1D("Time_AntiCoincidence_0","",4*window/10000,-2*window,2*window);
    TH1D *Hist_CH0_Energy           = new TH1D("Hist_CH0_Energy","",pow(2,14)+1,0,pow(2,14)*0.2);
    TH1D *Hist_CH1_Energy           = new TH1D("Hist_CH1_Energy","",pow(2,14)+1,0,pow(2,14)*0.2);
    TH1D *Hist_Coin_Energy          = new TH1D("Hist_Coin_Energy","",pow(2,14)+1,0,pow(2,14)*0.2);
    TH1D *Hist_Anti_CH0_Energy      = new TH1D("Hist_Anti_CH0_Energy","",pow(2,14)+1,0,pow(2,14)*0.2);
    TH1D *Hist_Anti_CH1_Energy      = new TH1D("Hist_Anti_CH1_Energy","",pow(2,14)+1,0,pow(2,14)*0.2);
    TH2D *EE_Coincidence            = new TH2D("EE_Coincidence","",8192,0,3276.6,8192,0,3276.6);
    TTree *Energy_Coincidence       = new TTree("Energy_Coincidence","");

    Double_t CH0_Branch_Entry;
    Double_t CH1_Branch_Entry;
    double_t CH0_Time;
    double_t CH1_Time;
    TBranch *Energy_CH0 = Energy_Coincidence->Branch("Energy_CH0", &CH0_Branch_Entry,"Energy_CH0/D");
    TBranch *Energy_CH1 = Energy_Coincidence->Branch("Energy_CH1", &CH1_Branch_Entry,"Energy_CH1/D");
    TBranch *Time_CH0 = Energy_Coincidence->Branch("Time_CH0", &CH0_Time,"Time_CH0/D");
    TBranch *Time_CH1 = Energy_Coincidence->Branch("Time_CH1", &CH1_Time,"Time_CH1/D");

    double Search_Floor_0 = 0;
    double Search_Floor_1 = 0;
    double coincidence_count = 0;
    bool Coincidence_Success;
    double_t rt_0d;
    double_t rt_1d;
    int PhaseCounter = 1e7;

    //loop through the HPGe counts
    while (Search_Floor_1 < Total_Count_1)
    {
        if (static_cast<int>(Search_Floor_1) % 1000000 == 0)
        {
            float progress = ((float)coincidence_count / PhaseCounter) * 100.0f;
            std::cout << "Auto Phaseshifting: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
            std::cout << '\r';
        }
        kListData_1 -> GetEntry(Search_Floor_1);
        rt_1d = rt_1;
        rt_1d += Static_Adjust;
        // as the loop for HPGe counts goes, evaluate if there are any similiar values in the alt detector
        Coincidence_Success = false;
        while (Search_Floor_0 < Total_Count_0)
        {
            kListData_0 -> GetEntry(Search_Floor_0);
            rt_0d = rt_0;
            if (rt_0d - rt_1d > window || rt_0d - rt_1d < -1.*window)
            {
                // the count is outside the acceptable window
                // if alt detector time stamp is lower than the HPGe, move the starting search position 
                if (rt_0d - rt_1d < 0)  
                {
                    Search_Floor_0 += 1;
                }
                // if alt detector time stamp is higher than the HPGe, move to the next HPGe count
                else                                    
                {
                    break;
                }
            }
            else
            {
                // count is within the window'
                coincidence_count += 1;
                Time_Coincidence_0->Fill(rt_0d - rt_1d);
                Search_Floor_0 += 1;
                Coincidence_Success = true;
            };
        }
        // anti coincidence count
        if (Coincidence_Success == false)
        {
            Time_AntiCoincidence_0->Fill(rt_0d - rt_1d);
        }
        Search_Floor_1 += 1;
        if (rt_1d >= Total_Time_1 || coincidence_count > PhaseCounter)
        {
            break;
        }
    }

    cout << "" << endl;
    //cout << Time_Coincidence_0->GetMaximum() << endl;
    //cout << Time_Coincidence_0->GetMinimum(1) << endl;
    //cout << Time_Coincidence_0->GetBinCenter(Time_Coincidence_0->GetMaximumBin()) << endl;
    //cout << "" << endl;

    TF1 *gaussian = new TF1("gaussian", "gaus+[A]", -1*window+1, window-1);
    gaussian->SetParameters(Time_Coincidence_0->GetMaximum()-Time_Coincidence_0->GetMinimum(1), Time_Coincidence_0->GetBinCenter(Time_Coincidence_0->GetMaximumBin()), window/10,Time_Coincidence_0->GetMinimum(1));
	TFitResultPtr r0 = Time_Coincidence_0 -> Fit(gaussian,"SRQ0");

    // redo using the new parameter values
    Static_Adjust += r0->Parameter(1);
    window = TMath::Abs(r0->Parameter(2)*4);
    if (window < 4e6) {window = 4e6;}
    TH1D *Time_Coincidence_1        = new TH1D("Time_Coincidence_1","",4*window/10000,-2*window,2*window);
    TH1D *Time_AntiCoincidence_1    = new TH1D("Time_AntiCoincidence_1","",4*window/10000,-2*window,2*window);

    Search_Floor_0 = 0;
    Search_Floor_1 = 0;
    coincidence_count = 0;
    double anticoincidence_count = 0;
    //loop through the HPGe counts
    while (Search_Floor_1 < Total_Count_1)
    {
        if (static_cast<int>(Search_Floor_1) % 1000000 == 0)
        {
            float progress = ((float)Search_Floor_1 / Total_Count_1) * 100.0f;
            float Est_Rate = Energy_Coincidence -> GetEntries() / (rt_1d/1e12);
            std::cout << "Timing Progress _ Second Pass: " << std::fixed << std::setprecision(1) << progress << "%" << "        " << "Estimated Rate:" << std::fixed << std::setprecision(1) << Est_Rate << std::flush;
            std::cout << '\r';
        }
        kListData_1 -> GetEntry(Search_Floor_1);
        rt_1d = rt_1 + Static_Adjust;
        // as the loop for HPGe counts goes, evaluate if there are any similiar values in the alt detector
        Coincidence_Success = false;
        while (Search_Floor_0 < Total_Count_0)
        {
            kListData_0 -> GetEntry(Search_Floor_0);
            rt_0d = rt_0;
            if (rt_0d - rt_1d > window || rt_0d - rt_1d < -1.*window)
            {
                // the count is outside the acceptable window
                // if Alt detector time stamp is lower than the Main, move the starting search position 
                if (rt_0d - rt_1d < 0) 	{Search_Floor_0 += 1;}
                // if Alt detector time stamp is higher than the Main + Window, move to the next Main count
                else {break;}
            }
            else
            {
                // count is within the window, the inital value is the lowest rt_main - rt_alt + phase possible, but we want to find the closet match
                // we do this by making a two vector reference system each vector contains {Main Time - Alt Time + Phase, Main Count #, Alt Count #}
                // The first vector is used as a reference, and the second vector moves up the Alt detector entry up one, check to see if time diff has shrunk
                vector<double> TimeDiffVector_0;
                vector<double> TimeDiffVector_1;
                double Original_Search_0 = Search_Floor_0;
                double Original_Search_1 = Search_Floor_1;
                bool MoreThanOneMatch;
                while(Search_Floor_1 < Total_Count_1)
                {
                    TimeDiffVector_0 = {TMath::Abs(rt_0d - rt_1d),Search_Floor_0,Search_Floor_1};
                    kListData_1 -> GetEntry(Search_Floor_1 + 1);
                    rt_1d = rt_1 + Static_Adjust;
                    TimeDiffVector_1 = {TMath::Abs(rt_0d - rt_1d),Search_Floor_0,Search_Floor_1+1};
                    // Time difference was smaller, throw out the original match
                    if (TimeDiffVector_1[0] < TimeDiffVector_0[0])  
                    {
                        TimeDiffVector_0 = TimeDiffVector_1;
                        Search_Floor_1  += 1;
                        MoreThanOneMatch = true;
                    }
                    // Lowest absolute time difference has been found
                    else {break;}
                }
                // Now we just cross check to make sure that the next few counts in Main Det to see if they are a better match
                while(Search_Floor_0 < Total_Count_0)
                {
                    kListData_0 -> GetEntry(Search_Floor_0 + 1);
                    rt_0d = rt_0;
                    TimeDiffVector_1 = {TMath::Abs(rt_0d - rt_1d),Search_Floor_0 + 1,Search_Floor_1};
                    // Time difference was smaller, then that means there are actually two coincidences
                    // Original Main Count and Alt Det Main Count - 1 (that would've been the next closest), and the next closest that we are about to find
                    if (TimeDiffVector_1[0] < TimeDiffVector_0[0])  
                    {
                        if (Search_Floor_0 == Original_Search_0 && MoreThanOneMatch == true)
                        {
                            kListData_0 -> GetEntry(TimeDiffVector_0[1]);
                            kListData_1 -> GetEntry(TimeDiffVector_0[2]-1);
                            rt_0d = rt_0;
                            rt_1d = rt_1 + Static_Adjust;
                            EE_Coincidence      -> Fill(Eg_0*0.4,Eg_1*0.4);
                            Time_Coincidence_1  -> Fill(rt_0d - rt_1d);
                            Hist_Coin_Energy    -> Fill(Eg_0*0.2);

                            CH0_Branch_Entry    = Eg_0*0.2;
                            CH1_Branch_Entry    = Eg_1*0.2;
                            CH0_Time            = rt_0d;
                            CH1_Time            = rt_1d;
                            Energy_Coincidence  -> Fill();
                            MoreThanOneMatch    = false;
                        }
                        TimeDiffVector_0 = TimeDiffVector_1;
                        Search_Floor_0  += 1;
                    }
                    // Lowest absolute time difference has been found
                    else {break;}
                }
                // Best possible match has been made we'll use those values
                kListData_0 -> GetEntry(TimeDiffVector_0[1]);
                kListData_1 -> GetEntry(TimeDiffVector_0[2]);
                rt_0d = rt_0;
                rt_1d = rt_1 + Static_Adjust;
                EE_Coincidence      -> Fill(Eg_0*0.4,Eg_1*0.4);
                Time_Coincidence_1  -> Fill(rt_0d - rt_1d);
                Hist_Coin_Energy    -> Fill(Eg_0*0.2);

                CH0_Branch_Entry    = Eg_0*0.2;
                CH1_Branch_Entry    = Eg_1*0.2;
                CH0_Time            = rt_0d;
                CH1_Time            = rt_1d;
                Energy_Coincidence  -> Fill();
                coincidence_count 	+= 1;
                Search_Floor_0 		+= 1;
                Coincidence_Success = true;
            }
        }
        // non-coincidence count
        if (Coincidence_Success == false)
        {
            Time_AntiCoincidence_1  -> Fill(rt_0d - rt_1d);
            Hist_Anti_CH0_Energy    -> Fill(Eg_0*0.2);
            Hist_Anti_CH1_Energy    -> Fill(Eg_1*0.2);
            anticoincidence_count 	+= 1;
        }
        Search_Floor_1      += 1;
    }
    cout << "Coincidence Pairing Completed" << endl;
	
    double Integral_A = EE_Coincidence->Integral(1332*2.5,1340*2.5,1174*2.5,1180*2.5);
    double Integral_b = EE_Coincidence->Integral(1174*2.5,1180*2.5,1332*2.5,1340*2.5);
    cout << "" << endl;
    cout << Integral_A + Integral_b << endl;

    cout << "" << endl;

    TCanvas *canvas = new TCanvas("canvas", "TH2D Plot", 900, 900);
    EE_Coincidence 		-> GetXaxis() 	-> SetRangeUser(1170,1180);
    EE_Coincidence 		-> GetYaxis() 	-> SetRangeUser(1330,1340);
    EE_Coincidence      -> Draw("colz");
    canvas              -> Update();
    canvas              -> Draw();

    gaussian            -> Write();
    Hist_CH0_Energy     -> Write();
    Hist_CH1_Energy     -> Write();
    Hist_Anti_CH0_Energy-> Write();
    Hist_Anti_CH1_Energy-> Write();
    Hist_Coin_Energy    -> Write();
    Time_Coincidence_0  -> Write();
    Time_Coincidence_1  -> Write();
    Energy_Coincidence  -> Write();
    EE_Coincidence      -> Write();
    f->Close();


    t.Stop();
    t.Print();
}