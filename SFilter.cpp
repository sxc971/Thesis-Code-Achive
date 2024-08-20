#include <TMath.h>

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

    ULong64_t width = rt/Tcount*10000;
    double_t cwidth = TMath::Power(2,20)*10000;
    
    double bin = TMath::Power(2,20);

    TString Histname = fname;
    Histname.ReplaceAll("DataR_CH0@DT5780_1967_","");
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

void SFilter(int Hist, const char *ext=".root")
{
    TStopwatch t;
    t.Start();
    TString dirname = gSystem -> WorkingDirectory();
    TCanvas *cs = new TCanvas("cs","cs",1920,1080);
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files) 
    {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) 
        {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext)) 
            {
                cout << fname.Data() << endl;
                TFile *g0 = new TFile(fname.Data());
                TTree *kListData_0 = (TTree*)g0->Get("Data_R");     

                fname.ReplaceAll(".root","_Filtered.root");
                TFile newfile(fname, "recreate");

                 // Create a new TTree with the same structure as the input tree
                TTree* outputTree = kListData_0->CloneTree(0);

                UInt_t Flag_Value {0};
                kListData_0 -> SetBranchAddress("Flags",&Flag_Value);
                long long int Total_Count = kListData_0 -> GetEntries();

                for (long long int i = 0; i < Total_Count; i++)
                {
                    if (i % 1000000 == 0)
                    {
                        float progress = ((float)i / Total_Count) * 100.0f;
                        std::cout << "Filter Progress: " << std::fixed << std::setprecision(1) << progress << "%" << std::flush;
                        std::cout << '\r';
                    }
                    kListData_0->GetEntry(i);
                    if (Flag_Value == 0)
                    {
                        outputTree -> Fill();
                    }
                }
                outputTree->Write();
                newfile.Close();
                if (Hist != 0)
                {
                    SHist(fname);
                }
            }
        }
    }
    t.Stop();
    t.Print();
}
