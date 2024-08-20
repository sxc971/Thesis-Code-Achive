void Merge_Fragments(TString MainFileName,int MaxFile, const char* ext = ".root")
{
    TChain Data_Set("Data_R");
    TString Og_Name = MainFileName;
    Og_Name.ReplaceAll(".root","_Merged.root");

    TString currentDir = gSystem->WorkingDirectory();
    
    for (int CurFile = 0; CurFile < MaxFile + 1; CurFile++)
    {
        Data_Set.Add(MainFileName);
        cout << MainFileName << endl;
        // Convert integer values to strings
        if (CurFile == 0)
        {
            MainFileName.ReplaceAll(".root","_0.root");
        }
        TString oldSubStr = Form("_%d.root", CurFile);
        TString newSubStr = Form("_%d.root", CurFile + 1);
        MainFileName.ReplaceAll(oldSubStr,newSubStr);
    }
    Data_Set.Merge(Og_Name);
}