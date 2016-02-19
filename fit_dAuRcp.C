////////////////////////////////////////////////////////////////////////////////
//
// Fit the value of beta to the 0-20%/60-88% dAu Jet Rcp
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 18 Feb 2016
//
////////////////////////////////////////////////////////////////////////////////
//
// NOTE: Calculations must have been done by
//       calculate_RCP.C
//       prior to running this macro
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TBox.h>

#include "data.h"

#include <iostream>

using namespace std;

void fit_dAuRcp()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //=================================================//
  // SET RUNNING CONDITIONS
  //=================================================//

  const int NBETA = 22;
  double beta[] = {1.00, 1.30, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38,
                   1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48,
                   1.50, 1.80
                  };

  const char* fileNames[NBETA];
  for (int i = 0; i < NBETA; i++)
  {
    fileNames[i] = Form("Rcp_dAu_beta%03.0f.root", beta[i] * 100);
  }

  //=================================================//
  // DECLARE VARIABLES
  //=================================================//

  TGraph* grcp[NBETA];

  double chi2ndf[NBETA];

  TGraph* gchi2_beta;


  TFile* fin;

  //=================================================//
  // READ FROM FILES
  //=================================================//
  cout << endl;
  cout << "--> Reading Rcp results from files" << endl;

  for (int ibeta = 0; ibeta < NBETA; ibeta++)
  {
    cout << "----> Reading results from " << fileNames[ibeta] << endl;
    fin = TFile::Open(fileNames[ibeta]);
    if (!fin)
    {
      cout << "ERROR!! Unable to open " << fileNames[ibeta] << endl;
      return;
    }

    grcp[ibeta] = (TGraph*) fin->Get("gRcp_dAu_cent0");
    if (!grcp[ibeta])
    {
      cout << "ERROR!! Unable to find gRcp_dAu_cent0 in "
           << fileNames[ibeta] << endl;
      return;
    }
    // grcp[ibeta]->SetDirectory(0);
    grcp[ibeta]->SetName(Form("grcp_%i", ibeta));

    fin->Close();
    delete fin;

  }// ibeta


  //=================================================//
  // CALCULATE CHI2
  //=================================================//
  cout << endl;
  cout << "--> Calculating chi2 for each beta value" << endl;

  double pT_jet[NJET] = {0};
  double pTe_jet[NJET] = {0};
  double x_jet[NJET] = {0};
  double xl_jet[NJET] = {0};
  double xh_jet[NJET] = {0};
  double xe_jet[NJET] = {0};
  for (int i = 0; i < NJET; i++)
  {
    pT_jet[i] = 0.5 * (pTl_jet[i] + pTh_jet[i]);
    x_jet[i] = 2 * pT_jet[i] / 200;

    xl_jet[i] = 2 * pTl_jet[i] / 200;
    xh_jet[i] = 2 * pTh_jet[i] / 200;
  }

  double minbeta = 0;
  double minchi2 = 999;
  int minidx = 0;
  double NDF = NJET - 1; // data points - free parameters (beta)
  for (int ibeta = 0; ibeta < NBETA; ibeta++)
  {

    double chi2 = 0;
    for (int i = 0; i < NJET; i++)
    {
      double m = grcp[ibeta]->Eval(x_jet[i]);
      double sig = TMath::Power(Rcp_jet_A[0][i], 2);
      sig += TMath::Power(Rcp_jet_B[0][i], 2);

      double tmp = TMath::Power(Rcp_jet[0][i] - m, 2) / sig;

      chi2 += tmp;
    }

    chi2ndf[ibeta] = chi2 / NDF;
    cout << " beta:" << beta[ibeta]
         << " chi2/ndf:" << chi2ndf[ibeta]
         << endl;

    if (chi2ndf[ibeta] < minchi2)
    {
      minchi2 = chi2ndf[ibeta];
      minbeta = beta[ibeta];
      minidx = ibeta;
    }

  } // ibeta

  cout << "   found best        - "
       << " beta: " << beta[minidx]
       << " chi2/ndf: " << chi2ndf[minidx]
       << endl;


  gchi2_beta = new TGraph(NBETA, beta, chi2ndf);
  gchi2_beta->SetTitle(";#beta;#chi^2/NDF");
  gchi2_beta->SetLineColor(kBlue);


  // Search for the error bounds, where chi2 > minchi2 + 1
  // (before dividing by NDF)
  int minidxl = 0;
  int minidxh = 0;

  // First check to the right of the chosen value
  for (int i = minidx; i < NBETA; i++)
  {

    if (chi2ndf[i] > minchi2 + 1.0 / NDF)
    {

      // check that the pervious bin wasn't a better match
      double di = TMath::Abs(chi2ndf[i] * NDF - minchi2 * NDF);
      double dp = TMath::Abs(chi2ndf[i - 1] * NDF - minchi2 * NDF);
      if (TMath::Abs(dp - 1) < TMath::Abs(di - 1))
      {
        minidxh = i - 1;
      }
      else
      {
        minidxh = i;
      }

      cout << "   found upper bound - "
           << " beta: " << beta[minidxh]
           << " chi2: " << chi2ndf[minidxh] * NDF
           << " (vs " << minchi2 * NDF + 1
           << ")"
           << endl;
      break;
    }
  }
  // Now check to the left of the chosen value
  for (int i = minidx; i >= 0; i--)
  {
    if (chi2ndf[i] > minchi2 + 1.0 / NDF)
    {
      // check that the pervious bin wasn't a better match
      double di = TMath::Abs(chi2ndf[i] * NDF - minchi2 * NDF);
      double dp = TMath::Abs(chi2ndf[i + 1] * NDF - minchi2 * NDF);
      if (TMath::Abs(dp - 1) < TMath::Abs(di - 1))
      {
        minidxl = i + 1;
      }
      else
      {
        minidxl = i;
      }
      cout << "   found lower bound - "
           << " beta: " << beta[minidxl]
           << " chi2: " << chi2ndf[minidxl] * NDF
           << " (vs " << minchi2 * NDF + 1
           << ")"
           << endl;
      break;
    }
  }

  // cout << " minbeta: " << minbeta
  //      << " +" << minbetah << " -" << minbetal
  //      << " minchi2: " << minchi2 << endl;

  //=================================================//
  // PLOT OBJECTS
  //=================================================//
  cout << endl;

  TH1F *haxis = new TH1F("haxis", ";#beta;#chi^{2}/NDF", 100, 1.0, 1.8);
  haxis->SetMinimum(0);
  haxis->SetMaximum(5);


  TGraphErrors *gRCP_jet;
  TBox *bRcp_jet[NJET];
  gRCP_jet = new TGraphErrors(NJET,
                              x_jet, Rcp_jet[0],
                              xe_jet, Rcp_jet_A[0]);
  gRCP_jet->SetMarkerStyle(20);
  gRCP_jet->SetMarkerColor(kBlack);

  // systematic uncertainties
  for (int i = 0; i < NJET; i++)
  {
    double x1 = xl_jet[i];
    double x2 = xh_jet[i];

    double y1 = Rcp_jet[0][i] - Rcp_jet_B[0][i];
    double y2 = Rcp_jet[0][i] + Rcp_jet_B[0][i];

    bRcp_jet[i] = new TBox(x1, y1, x2, y2);
    bRcp_jet[i]->SetFillColorAlpha(kGray, 0.6);
  }

  TH1F* haxis_rcp = new TH1F("haxis_rcp",
                             ";x_{p};R_{CP}",
                             100, 0, 0.6);
  haxis_rcp->GetYaxis()->SetTitleOffset(1.2);
  haxis_rcp->GetYaxis()->CenterTitle();
  haxis_rcp->GetYaxis()->SetLabelSize(0.04);
  haxis_rcp->GetYaxis()->SetTitleSize(0.04);
  haxis_rcp->GetYaxis()->SetNdivisions(6, 3, 0);
  haxis_rcp->GetXaxis()->SetTitleOffset(0.9);
  haxis_rcp->GetXaxis()->SetLabelSize(0.04);
  haxis_rcp->GetXaxis()->SetTitleSize(0.04);
  haxis_rcp->GetXaxis()->SetNdivisions(5, 4, 0);
  haxis_rcp->SetMinimum(0.01);
  haxis_rcp->SetMaximum(1.19);

  TLine lmin;
  lmin.SetLineColor(kRed);
  lmin.SetLineStyle(2);

  TLatex tmin;
  tmin.SetTextColor(kRed);
  tmin.SetTextAlign(22);
  tmin.SetNDC();

  //=================================================//
  // PLOT
  //=================================================//

  TCanvas* cchi2 = new TCanvas("cchi2", "chi2", 800, 800);

  cchi2->cd(1);
  haxis->Draw();
  gchi2_beta->Draw("C");

  lmin.DrawLine(1.0, chi2ndf[minidx], 1.8, chi2ndf[minidx]);
  lmin.DrawLine(beta[minidx], 0, beta[minidx], 5);
  tmin.DrawLatex(0.7, 0.8,
                 Form("#beta = %.2f^{+%.2f}_{-%.2f}",
                      beta[minidx],
                      beta[minidxh] - beta[minidx],
                      beta[minidx] - beta[minidxl]));


  TCanvas* crcp = new TCanvas("crcp", "rcp", 800, 800);

  crcp->cd(1);
  haxis_rcp->Draw();

  for (int i = 0; i < NJET; i++)
    bRcp_jet[i]->Draw();
  gRCP_jet->Draw("P");


  grcp[minidxl]->SetLineColor(kBlue);
  grcp[minidxl]->Draw("C");

  grcp[minidxh]->SetLineColor(kBlue);
  grcp[minidxh]->Draw("C");

  grcp[minidx]->Draw("C");



}