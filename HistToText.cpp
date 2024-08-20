#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1D.h"

int HistToText(TString filename) {

  // Open the ROOT file
  TFile* file = TFile::Open(filename);

  // Get the histograms from the file
  TH1D* hdiff1 = (TH1D*)file->Get("Time_Coincidence_1");
  TH1D* Pomme0 = (TH1D*)file->Get("Time_Coincidence_1");

  // Open the output text file
  std::ofstream outfile(filename + ".txt");

  // Write the header line
  outfile << "bin center; hdiff1; Pomme0" << std::endl;

  // Loop over the bins of the histograms
  for (int i=1; i<=hdiff1->GetNbinsX(); i++) {

    // Get the bin center
    double bin_center = hdiff1->GetBinCenter(i);

    // Get the bin content for each histogram
    double hdiff1_content = hdiff1->GetBinContent(i);
    double hdiff2_content = Pomme0->GetBinContent(i);

    // Write the data to the output file
    outfile << bin_center << ";" << hdiff1_content << ";" << hdiff2_content << std::endl;
  }

  // Close the output file and the ROOT file
  outfile.close();
  file->Close();

  return 0;
}
