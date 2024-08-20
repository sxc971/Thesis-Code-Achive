void SHist(TString fname)
{
    TFile *g0 = new TFile(fname.Data());
    ULong64_t rt{0};
    TTree *kListData = (TTree*)g0->Get("Data_R");
    kListData->SetBranchStatus("*",0);
    kListData->SetBranchStatus("Timestamp",1);
    kListData->SetBranchAddress("Timestamp",&rt);
    double Tcount = kListData->GetEntries();
    kListData->GetEntry(Tcount-1);
    ULong64_t Tlength = rt;

    ULong64_t width = rt/Tcount*10;
    double_t cwidth = TMath::Power(2,20)*40000;
    
    double bin = TMath::Power(2,20);

    TString Histname = fname;
    Histname.ReplaceAll("DataR_CH0@DT5781_1967_","");
    Histname.ReplaceAll(".root","");

    TH1D* hdiff1 = new TH1D("hdiff1","",bin,0,cwidth);
    TH1D* hdiff2 = new TH1D("hdiff2","",bin,0,cwidth);
    TH1D* hdiff3 = new TH1D("hdiff3","",bin,0,cwidth);
    TH1D* hdiff4 = new TH1D("hdiff4","",bin,0,cwidth);
    hdiff1->GetXaxis()->SetTitle("Time difference between counts (ps)");
    hdiff1->GetYaxis()->SetTitle("Frequency");

    for (Long64_t i = 0; i < kListData->GetEntries()-1; i++)
        //				for (Long64_t i = 0; i < 5e7+1; i++)
    {        
        if (i % 1000000 == 0)
        {
            float progress = (i / Tcount) * 100.0f;
            std::cout << "Timing Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
            std::cout << '\r';
        }
        kListData->GetEntry(i);
        double r1 = rt;
        double diff = 0;

        kListData->GetEntry(i+1);
        diff = rt - r1;
        hdiff1->Fill(diff);

        kListData->GetEntry(i+2);
        diff = rt - r1;
        hdiff2->Fill(diff);

        kListData->GetEntry(i+3);
        diff = rt - r1;
        hdiff3->Fill(diff);

        kListData->GetEntry(i+4);
        diff = rt - r1;
        hdiff4->Fill(diff);
    }

    double Integral = 0;
    for (int i = 0; i < bin; i++)
    {
        Integral += hdiff1 -> GetBinContent(i);
        if (Integral > 0.999 * Tcount)
        {
            TAxis* xAxis = hdiff1->GetXaxis();
            xAxis->SetRangeUser(0, i*hdiff1->GetXaxis()->GetBinWidth(1));
        }
    }

    cout << "" << endl;
    hdiff1->Draw();
    hdiff2->SetLineColor(kRed);
    hdiff2->Draw("same");
    hdiff3->SetLineColor(kGreen);
    hdiff3->Draw("same");
    hdiff4->SetLineColor(kMagenta);
    hdiff4->Draw("same");

    fname.ReplaceAll("DataR_CH0@DT5781_14469_","");
    fname.ReplaceAll(".root","");
    fname.Append("_Hist.root");
    TFile f(fname,"new");
    hdiff1 -> Write();
    hdiff2 -> Write();
    hdiff3 -> Write();
    hdiff4 -> Write();

    delete hdiff1;
    delete hdiff2;
    delete hdiff3;
    delete hdiff4;
    kListData->Delete("");
    f.Close();
}

void Custom_Deadtime(int Hist,const char *ext=".root")
{
    TStopwatch t;
    t.Start();
    TString dirname = gSystem -> WorkingDirectory();
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

                TFile *g0 = new TFile(fname.Data());
                TTree *kListData_0 = (TTree*)g0->Get("Data_R"); 

                ULong64_t rt{0};
                UShort_t Eg{0};
                TTree *kListData = (TTree*)g0->Get("Data_R");
                kListData->SetBranchAddress("Timestamp",&rt);  
                kListData->SetBranchAddress("Energy",&Eg);  

                ULong64_t Tcount = kListData->GetEntries();
                kListData->GetEntry(Tcount-1);
                ULong64_t Tlength = rt;

                double basePET = 7.0e6;
                double baseDeadtime = 8.5e6;
                std::vector<double> PET = {basePET,2*basePET,4*basePET,8*basePET,16*basePET,32*basePET,64*basePET,128*basePET};
                std::vector<double> Deadtime = {baseDeadtime,2*baseDeadtime,4*baseDeadtime,8*baseDeadtime,16*baseDeadtime,32*baseDeadtime,64*baseDeadtime,128*baseDeadtime};
                double est_icr_initial = Tlength / Tcount;
                double bin = TMath::Power(2,16);
                ULong64_t width = rt/Tcount*10;
                double_t cwidth = bin*10000;

                TTree* outputTree = kListData->CloneTree(0);
                std::vector<ULong64_t> CountRemover;
                int vector_size_cap = 1e6;

                for (int DeadtimeIter = 0; DeadtimeIter < PET.size(); DeadtimeIter++)
                {
                    TString newname = fname;
                    newname.ReplaceAll("DataR_CH0@DT5781_14469_","");
                    newname.ReplaceAll(".root", "");  // Remove the original extension
                    newname += Form("_PET_%i_Dead_%i.root", static_cast<int>(PET[DeadtimeIter]/1e6), static_cast<int>(Deadtime[DeadtimeIter]/1e6));
                    TFile newfile(newname, "recreate");

                    TTree* outputTree = kListData->CloneTree(0);
                    for (ULong64_t vector_number = 0; vector_number < Tcount/vector_size_cap + 1; vector_number++)
                    {
                        float progress = (static_cast<double>(vector_number*vector_size_cap) / static_cast<double>(Tcount)) * 100.0f;
                        std::cout << "Removal Flagging: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
                        std::cout << '\r';

                        CountRemover.clear();
                        ULong64_t Track = 0;

                        for (ULong64_t i = 0; i < vector_size_cap; i++)
                        {
                            kListData->GetEntry(i + vector_number*vector_size_cap);
                            double r1 = rt;
                            double diff = 0;

                            kListData->GetEntry(i + 1 + vector_number*vector_size_cap);
                            diff = rt - r1;

                            if (diff < PET[DeadtimeIter])
                            {
                                // remove both counts
                                CountRemover.push_back(i);
                                CountRemover.push_back(i + 1);
                            }
                            else if (diff < Deadtime[DeadtimeIter])
                            {
                                // remove last count
                                CountRemover.push_back(i + 1);
                            }
                            else if (diff > Deadtime[DeadtimeIter])
                            {
                                // keep both counts
                            }
                        }

                        // mark which entries should be ignored when rebuilding
                        // Create a set to store unique elements
                        std::set<ULong64_t> uniqueSet(CountRemover.begin(), CountRemover.end());

                        // Create a new vector with unique elements
                        std::vector<ULong64_t> uniqueVector(uniqueSet.begin(), uniqueSet.end());                
                    
                        for (Long64_t entry = 0; entry < vector_size_cap; entry++)
                        {
                            if (entry != uniqueVector[Track])
                            {
                                kListData -> GetEntry(vector_number*vector_size_cap + entry);
                                outputTree -> Fill();
                            }
                            else{Track += 1;}
                        }
                    }
                    outputTree->Write();
                    newfile.Close();

                    if (Hist != 0)
                    {
                        SHist(newname);
                    }
                }
            }
        }
    }
    t.Stop();
    t.Print();
}
