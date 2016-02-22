////////////////////////////////////////////////////////////////////////////////
//
// Plot the resulting RpA and Rcp modifications
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 18 Feb 2016
//
////////////////////////////////////////////////////////////////////////////////
//
// File reads in results from
//  calculate_RCP.C
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TBox.h>
#include <TLine.h>
#include <TGaxis.h>

#include <data.h>

#include <iostream>
#include <iomanip>

using namespace std;


void plot_Results()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //=====================================================//
  // SET RUNNING CONDITIONS
  //=====================================================//

  bool printPlots = false;

  double beta = 1.38;

  const int NSYSTEMS = 3;
  const char *collSystem[] = {"pAu", "dAu", "3HeAu"};

  const char *inFile = Form("Rcp_systems_beta%03.0f.root", beta * 100);

  // Centrality bins
  const int NCENT = 4;
  double centl[NSYSTEMS][NCENT + 1] =
  {
    {0, 20, 40, 60, 84}, // pAu
    {0, 20, 40, 60, 88}, // dAu
    {0, 20, 40, 60, 88}, // He3Au
  };

  // line colors
  int colors[] = { kRed, kBlue, kGreen + 2};
  int lstyle[] = {4, 1, 7};
  int lineWidth = 3;

  // fill colors
  int fcolor_cent[NCENT] = {kBlue, kRed, kGreen + 2, kYellow + 2};

  //=====================================================//
  // DECLARE VARIABLES
  //=====================================================//

  TGraph *graa[NSYSTEMS][NCENT];
  TGraph *grcp[NSYSTEMS][NCENT];

  TGraph *graa_linMod[NSYSTEMS][NCENT];
  TGraph *grcp_linMod[NSYSTEMS][NCENT];

  TGraph *gQ[NSYSTEMS];

  // sigma modification function
  TF1* fsig = new TF1("fsig", "TMath::Power((1 + TMath::Exp(-1. *[0] * x)), 2) / 4.0", 0, 1);
  fsig->SetLineColor(kBlue);
  fsig->SetParameter(0, beta);

  char gname[500];

  //=====================================================//
  // READ RESULTS FROM FILE
  //=====================================================//
  cout << endl;
  cout << "--> Reading results from " << inFile << endl;

  TFile *fin = TFile::Open(inFile);
  if (!fin)
  {
    cout << "ERROR!! Unable to open " << inFile << endl;
    return;
  }

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    //mean Q vs x
    sprintf(gname, "gQ_%s", collSystem[isys]);
    gQ[isys] = (TGraph*) fin->Get(gname);
    if (!gQ[isys])
    {
      cout << "ERROR!! Unable to find graph " << gname
           << " in " << inFile << endl;
      return;
    }

    for (int icent = 0; icent < NCENT - 1; icent++)
    {
      // Rcp
      sprintf(gname, "gRcp_%s_cent%i", collSystem[isys], icent);
      grcp[isys][icent] = (TGraph*) fin->Get(gname);
      if (!grcp[isys][icent])
      {
        cout << "ERROR!! Unable to find graph " << gname
             << " in " << inFile << endl;
        return;
      }

      // Rcp w/ linear modification
      sprintf(gname, "gRcp_%s_cent%i_linMod", collSystem[isys], icent);
      grcp_linMod[isys][icent] = (TGraph*) fin->Get(gname);
      if (!grcp_linMod[isys][icent])
      {
        cout << "ERROR!! Unable to find graph " << gname
             << " in " << inFile << endl;
        return;
      }

    } // icent

    for (int icent = 0; icent < NCENT; icent++)
    {
      // raa
      sprintf(gname, "gRAA_%s_cent%i", collSystem[isys], icent);
      graa[isys][icent] = (TGraph*) fin->Get(gname);
      if (!graa[isys][icent])
      {
        cout << "ERROR!! Unable to find graph " << gname
             << " in " << inFile << endl;
        return;
      }

      // raa w/ linear modification
      sprintf(gname, "gRAA_%s_cent%i_linMod", collSystem[isys], icent);
      graa_linMod[isys][icent] = (TGraph*) fin->Get(gname);
      if (!graa_linMod[isys][icent])
      {
        cout << "ERROR!! Unable to find graph " << gname
             << " in " << inFile << endl;
        return;
      }

    }// icent
  }// isys



  //=====================================================//
  // Run 8 d+Au Jet Rcp
  //=====================================================//
  cout << endl;
  cout << "--> Run 8 d+Au Jet Rcp" << endl;;

  // Data arrays included from data.h

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

  TGraphErrors *gRCP_jet[NCENT - 1];
  TBox *bRcp_jet[NCENT - 1][NJET];
  for (int icent = 0; icent < NCENT - 1; icent++)
  {
    gRCP_jet[icent] = new TGraphErrors(NJET,
                                       x_jet, Rcp_jet[icent],
                                       xe_jet, Rcp_jet_A[icent]);
    gRCP_jet[icent]->SetMarkerStyle(20);
    gRCP_jet[icent]->SetMarkerColor(kBlack);

    // systematic uncertainties
    for (int i = 0; i < NJET; i++)
    {
      double x1 = xl_jet[i];
      double x2 = xh_jet[i];

      double y1 = Rcp_jet[icent][i] - Rcp_jet_B[icent][i];
      double y2 = Rcp_jet[icent][i] + Rcp_jet_B[icent][i];

      bRcp_jet[icent][i] = new TBox(x1, y1, x2, y2);
      bRcp_jet[icent][i]->SetFillColorAlpha(kGray, 0.6);
    }
  }

  TGraphErrors *gRdAu_jet[NCENT];
  TBox *bRdAu_jet[NCENT][NJET];
  for (int icent = 0; icent < NCENT; icent++)
  {
    gRdAu_jet[icent] = new TGraphErrors(NJET,
                                        x_jet, RdAu_jet[icent],
                                        xe_jet, RdAu_jet_A[icent]);
    gRdAu_jet[icent]->SetMarkerStyle(20);
    gRdAu_jet[icent]->SetMarkerColor(kBlack);

    // systematic uncertainties
    for (int i = 0; i < NJET; i++)
    {
      double x1 = xl_jet[i];
      double x2 = xh_jet[i];

      double y1 = RdAu_jet[icent][i] - RdAu_jet_B[icent][i];
      double y2 = RdAu_jet[icent][i] + RdAu_jet_B[icent][i];

      bRdAu_jet[icent][i] = new TBox(x1, y1, x2, y2);
      bRdAu_jet[icent][i]->SetFillColorAlpha(kGray, 0.6);
    }
  }



  //=====================================================//
  // DI-HADRON CORRELATIONS (PPG128)
  //=====================================================//
  cout << endl;
  cout << "--> Run 8 d+Au dihadron correlations" << endl;

  // data arrays included from data.h
  double xp_JdA[NPPG128] = {};
  double xe_JdA[NPPG128] = {0};


  double zfrag = 0.6;
  for (int j = 0; j < NPPG128; j++)
  {
    double xAufrag = (pT_trg[j] * TMath::Exp(-1.*eta_trg[j]) + pT_assoc[j] * TMath::Exp(-1.*eta_assoc)) / 200.;
    double xpfrag = (pT_trg[j] * TMath::Exp(eta_trg[j]) + pT_assoc[j] * TMath::Exp(eta_assoc)) / 200.;

    xp_JdA[j] = xpfrag / zfrag;

    cout << " " << j << " "  << xp_JdA[j] << " " << JdA[0][j] << endl;
  }

  // calcualte Jcp and ratio to model
  double Jcp[NPPG128] = {0};
  double Jcp_A[NPPG128] = {0};
  double Jcp_Bl[NPPG128] = {0};
  double Jcp_Bh[NPPG128] = {0};

  double JdA_rat[2][NPPG128];
  double JdA_rat_A[2][NPPG128];
  double JdA_rat_Bl[2][NPPG128];
  double JdA_rat_Bh[2][NPPG128];

  for (int j = 0; j < NPPG128; j++)
  {
    //Jcp
    Jcp[j] = JdA[0][j] / JdA[1][j];
    Jcp_A[j] = TMath::Sqrt(TMath::Power(JdA_A[0][j] / JdA[0][j], 2) +
                           TMath::Power(JdA_A[1][j] / JdA[1][j], 2));
    Jcp_Bl[j] = TMath::Sqrt(TMath::Power(JdA_Bl[0][j] / JdA[0][j], 2) +
                            TMath::Power(JdA_Bl[1][j] / JdA[1][j], 2));
    Jcp_Bh[j] = TMath::Sqrt(TMath::Power(JdA_Bh[0][j] / JdA[0][j], 2) +
                            TMath::Power(JdA_Bh[1][j] / JdA[1][j], 2));

    Jcp_A[j] = Jcp[j] * Jcp_A[j];
    Jcp_Bl[j] = Jcp[j] * Jcp_Bl[j];
    Jcp_Bh[j] = Jcp[j] * Jcp_Bh[j];

    // ratio to the model
    double mod_cent = graa[1][0]->Eval(xp_JdA[j]);
    double mod_periph = graa[1][NCENT - 1]->Eval(xp_JdA[j]);

    JdA_rat[0][j] = JdA[0][j] / mod_cent;
    JdA_rat_A[0][j] = JdA_A[0][j] / mod_cent;
    JdA_rat_Bl[0][j] = JdA_Bl[0][j] / mod_cent;
    JdA_rat_Bh[0][j] = JdA_Bh[0][j] / mod_cent;

    JdA_rat[1][j] = JdA[1][j] / mod_periph;
    JdA_rat_A[1][j] = JdA_A[1][j] / mod_periph;
    JdA_rat_Bl[1][j] = JdA_Bl[1][j] / mod_periph;
    JdA_rat_Bh[1][j] = JdA_Bh[1][j] / mod_periph;

  }



  const int NXBINS = 3;
  float xbins[] = {1e-4, 1e-3, 1e-2, 1e-1};
  float xmin_xbins[NXBINS] = {0};
  float xmax_xbins[NXBINS] = {0};
  int xbinColor[] = {kBlack, kMagenta + 2, kOrange + 2};
  TGraphErrors *gJdA[NXBINS][2];
  TGraphErrors *gJcp[NXBINS];
  TBox *bJdA[NXBINS][2][NPPG128];
  TBox *bJcp[NXBINS][NPPG128];
  int nxjda[NXBINS] = {0};
  for (int ix = 0; ix < NXBINS; ix++)
  {
    //-- JdA
    for (int k = 0; k < 2; k++)
    {
      nxjda[ix] = 0;

      gJdA[ix][k] = new TGraphErrors();
      gJdA[ix][k]->SetMarkerStyle(20 + ix);
      gJdA[ix][k]->SetMarkerColor(xbinColor[ix]);

      for (int j = 0; j < NPPG128; j++)
      {
        double x1 = xp_JdA[j] - 0.01;// * xp_JdA[j];
        double x2 = xp_JdA[j] + 0.01;// * xp_JdA[j];
        double y1 = JdA[k][j] - JdA_Bl[k][j];
        double y2 = JdA[k][j] + JdA_Bh[k][j];

        if (xAu_frag[j] > xbins[ix] && xAu_frag[j] < xbins[ix + 1])
        {
          gJdA[ix][k]->SetPoint(nxjda[ix], xp_JdA[j], JdA[k][j]);
          gJdA[ix][k]->SetPointError(nxjda[ix], 0, JdA_A[k][j]);

          bJdA[ix][k][nxjda[ix]] = new TBox(x1, y1, x2, y2);
          bJdA[ix][k][nxjda[ix]]->SetFillColorAlpha(kGray, 0.6);

          nxjda[ix]++;


        }
      } // j
    } // k

    //-- Jcp
    gJcp[ix] = new TGraphErrors();
    gJcp[ix]->SetMarkerStyle(20 + ix);
    gJcp[ix]->SetMarkerColor(xbinColor[ix]);

    nxjda[ix] = 0;

    for (int j = 0; j < NPPG128; j++)
    {
      double x1 = xp_JdA[j] - 0.01;// * xp_JdA[j];
      double x2 = xp_JdA[j] + 0.01;// * xp_JdA[j];
      double y1 = Jcp[j] - Jcp_Bl[j];
      double y2 = Jcp[j] + Jcp_Bh[j];

      if (xAu_frag[j] > xbins[ix] && xAu_frag[j] < xbins[ix + 1])
      {
        gJcp[ix]->SetPoint(nxjda[ix], xp_JdA[j], Jcp[j]);
        gJcp[ix]->SetPointError(nxjda[ix], 0, Jcp_A[j]);

        bJcp[ix][nxjda[ix]] = new TBox(x1, y1, x2, y2);
        // bJcp[ix][nxjda[ix]]->SetFillColorAlpha(kMagenta - 3 + ix, 0.6);
        bJcp[ix][nxjda[ix]]->SetFillColorAlpha(kGray, 0.6);

        nxjda[ix]++;
      }
    } // j

    // find the min/max xAu values
    for (int j = 0; j < NPPG128; j++)
    {
      if (j == 0) xmin_xbins[ix] = 1;

      if (xAu_frag[j] > xbins[ix] && xAu_frag[j] < xbins[ix + 1])
      {
        if (xAu_frag[j] < xmin_xbins[ix]) xmin_xbins[ix] = xAu_frag[j];
        if (xAu_frag[j] > xmax_xbins[ix]) xmax_xbins[ix] = xAu_frag[j];
      }
    }

  } // ix


  TGraphErrors *gJdA_rat[2];
  TBox *bJdA_rat[2][NPPG128];
  for (int k = 0; k < 2; k++)
  {

    gJdA_rat[k] = new TGraphErrors(NPPG128, xp_JdA, JdA_rat[k], xe_JdA, JdA_rat_A[k]);
    gJdA_rat[k]->SetMarkerStyle(20);
    gJdA_rat[k]->SetMarkerColor(kBlack);

    for (int j = 0; j < NPPG128; j++)
    {
      double x1 = xp_JdA[j] - 0.01;// * xp_JdA[j];
      double x2 = xp_JdA[j] + 0.01;// * xp_JdA[j];
      double y1 = JdA_rat[k][j] - JdA_rat_Bl[k][j];
      double y2 = JdA_rat[k][j] + JdA_rat_Bh[k][j];
      bJdA_rat[k][j] = new TBox(x1, y1, x2, y2);
      bJdA_rat[k][j]->SetFillColorAlpha(kGray, 0.6);
    }
  }



  // //=====================================================//
  // // PLOT OBJECTS
  // //=====================================================//
  // cout << endl;
  // cout << "--> Plotting ..." << endl;

  // TLatex label;
  // label.SetNDC();
  // label.SetTextAlign(22);

  // TH1F* haxis_rcp_pipT = new TH1F("haxis_rcp_pipT",
  //                                 ";#pi p_{T};R_{CP}",
  //                                 100, 0, 20);
  // haxis_rcp_pipT->GetYaxis()->SetTitleOffset(1.3);
  // haxis_rcp_pipT->GetXaxis()->SetTitleOffset(1.3);
  // haxis_rcp_pipT->GetYaxis()->CenterTitle();
  // haxis_rcp_pipT->SetMinimum(0.01);
  // haxis_rcp_pipT->SetMaximum(1.09);

  // TLine l1;
  // l1.SetLineStyle(2);

  // int xplot = 2;

  // TLegend *legBBC[NSYSTEMS];
  // TLegend *legNcollMB[NSYSTEMS];
  // TLegend *legNcoll[NSYSTEMS][NCENT];
  // for (int isys = 0; isys < NSYSTEMS; isys++)
  // {
  //   //-- BBCs charge
  //   hBBCscNcollA[isys][xplot]->SetLineColor(kBlack);
  //   hBBCscModNcollA[isys][xplot]->SetLineColor(kBlue);
  //   hBBCscNcollModA[isys][xplot]->SetLineColor(kGreen + 2);
  //   hBBCscModNcollModA[isys][xplot]->SetLineColor(kRed);

  //   hBBCscNcollA[isys][xplot]->SetLineWidth(lineWidth);
  //   hBBCscModNcollModA[isys][xplot]->SetLineWidth(lineWidth);

  //   legBBC[isys] = new TLegend(0.15, 0.15, 0.4, 0.4,
  //                              Form("%s @ 200 GeV, x = %.2f",
  //                                   collSystem[isys], x[xplot]));
  //   legBBC[isys]->SetFillStyle(0);
  //   legBBC[isys]->SetBorderSize(0);
  //   legBBC[isys]->SetTextSize(0.04);
  //   legBBC[isys]->AddEntry(hBBCscNcollA[isys][xplot], "Unmodified", "L");
  //   // legBBC[isys]->AddEntry(hBBCscModNcollA[isys][xplot], "Mod BBCs chg", "L");
  //   // legBBC[isys]->AddEntry(hBBCscNcollModA[isys][xplot], "Mod N_{coll}", "L");
  //   legBBC[isys]->AddEntry(hBBCscModNcollModA[isys][xplot], "Mod BBCs chg and N_{coll}", "L");

  //   //-- MB Ncoll
  //   hNcoll_MB[isys][xplot]->SetLineColor(kBlue);
  //   hNcollMod_MB[isys][xplot]->SetLineColor(kRed);

  //   hNcoll_MB[isys][xplot]->SetLineWidth(lineWidth);
  //   hNcollMod_MB[isys][xplot]->SetLineWidth(lineWidth);

  //   if (isys == 0)
  //   {
  //     legNcollMB[isys] = new TLegend(0.5, 0.6, 0.75, 0.85,
  //                                    Form("%s @ 200 GeV, x = %.2f",
  //                                         collSystem[isys], x[xplot]));
  //   }
  //   else
  //   {
  //     legNcollMB[isys] = new TLegend(0.15, 0.15, 0.4, 0.4,
  //                                    Form("%s @ 200 GeV, x = %.2f",
  //                                         collSystem[isys], x[xplot]));
  //   }
  //   legNcollMB[isys]->SetFillStyle(0);
  //   legNcollMB[isys]->SetBorderSize(0);
  //   legNcollMB[isys]->SetTextSize(0.04);
  //   legNcollMB[isys]->AddEntry((TObject*)0,
  //                              Form("<N_{coll}> = %.2f (PHENIX)",
  //                                   Ncoll_PHENIX[isys][NCENT]), "");
  //   legNcollMB[isys]->AddEntry(hNcoll_MB[isys][xplot],
  //                              Form("<N_{coll}> = %.2f",
  //                                   hNcoll_MB[isys][xplot]->GetMean())
  //                              , "L");
  //   legNcollMB[isys]->AddEntry(hNcollMod_MB[isys][xplot],
  //                              Form("<N_{coll}^{mod}> = %.2f",
  //                                   hNcollMod_MB[isys][xplot]->GetMean())
  //                              , "L");

  //   //-- centrality dep Ncoll
  //   for (int icent = 0; icent < NCENT; icent++)
  //   {
  //     hNcoll_cent[isys][xplot][icent]->SetLineColor(kBlue);
  //     hNcollMod_cent[isys][xplot][icent]->SetLineColor(kRed);

  //     hNcoll_cent[isys][xplot][icent]->SetLineWidth(lineWidth);
  //     hNcollMod_cent[isys][xplot][icent]->SetLineWidth(lineWidth);

  //     if (isys == 0 || icent >= NCENT - 2)
  //     {
  //       legNcoll[isys][icent] = new TLegend(0.6, 0.5, 0.9, 0.85,
  //                                           Form("%s, x = %.2f, %.0f-%.0f%%",
  //                                               collSystem[isys], x[xplot],
  //                                               centl[isys][icent], centl[isys][icent + 1]));
  //     }
  //     else
  //     {
  //       legNcoll[isys][icent] = new TLegend(0.15, 0.12, 0.4, 0.45,
  //                                           Form("%s @ 200 GeV, x = %.2f",
  //                                               collSystem[isys], x[xplot]));
  //     }
  //     legNcoll[isys][icent]->SetFillStyle(0);
  //     legNcoll[isys][icent]->SetBorderSize(0);
  //     legNcoll[isys][icent]->SetTextSize(0.06);
  //     legNcoll[isys][icent]->AddEntry((TObject*)0,
  //                                     Form("<N_{coll}> = %.2f (PHENIX)",
  //                                          Ncoll_PHENIX[isys][icent]), "");
  //     legNcoll[isys][icent]->AddEntry(hNcoll_cent[isys][xplot][icent],
  //                                     Form("<N_{coll}> = %.2f",
  //                                          hNcoll_cent[isys][xplot][icent]->GetMean())
  //                                     , "L");
  //     legNcoll[isys][icent]->AddEntry(hNcollMod_cent[isys][xplot][icent],
  //                                     Form("<N_{coll}^{mod}> = %.2f",
  //                                          hNcollMod_cent[isys][xplot][icent]->GetMean())
  //                                     , "L");


  //   }

  // }

  // TLegend *legx = new TLegend(0.8, 0.5, 0.98, 0.98);
  // legx->SetFillStyle(0);
  // legx->SetBorderSize(0);
  // legx->SetTextSize(0.04);
  // legx->AddEntry(hBBCscNcollA[0][0], "x_{p} = 0.00", "L");
  // for (int ix = 0; ix < NX; ix++)
  // {
  //   legx->AddEntry(hBBCscModNcollModA[0][ix],
  //                  Form("x_{p} = %.2f", x[ix]),
  //                  "L");
  // }

  // TLegend *legmod[NSYSTEMS];
  // for (int isys = 0; isys < NSYSTEMS; isys++)
  // {

  //   legmod[isys] = new TLegend(0.2, 0.7, 0.98, 0.98);
  //   legmod[isys]->SetFillStyle(0);
  //   legmod[isys]->SetBorderSize(0);
  //   legmod[isys]->SetTextSize(0.08);
  //   legmod[isys]->SetNColumns(3);
  //   legmod[isys]->AddEntry(graa_MB[isys],
  //                          "N_{coll}^{Mod} & BBCsc^{Mod}",
  //                          "L");
  //   legmod[isys]->AddEntry(graa_NcollModABBCsc_MB[isys],
  //                          "N_{coll}^{Mod} & BBCsc",
  //                          "L");
  //   legmod[isys]->AddEntry(graa_NcollABBCscMod_MB[isys],
  //                          "N_{coll} & BBCsc^{Mod}",
  //                          "L");
  // }


  // //=====================================================//
  // // PLOT
  // //=====================================================//

  // TCanvas *crcppipt = new TCanvas("crcppipt", "RCP", 1400, 400);
  // crcppipt->SetTopMargin(0.0);
  // crcppipt->SetRightMargin(0.0);
  // crcppipt->SetBottomMargin(0.0);
  // crcppipt->SetLeftMargin(0.0);
  // crcppipt->Divide(NCENT - 1, 1, 0, 0);

  // for (int icent = 0; icent < NCENT - 1; icent++)
  // {
  //   crcppipt->GetPad(icent + 1)->SetTopMargin(0.02);
  //   crcppipt->GetPad(icent + 1)->SetRightMargin(0.02);
  //   crcppipt->GetPad(icent + 1)->SetBottomMargin(0.10);
  //   crcppipt->GetPad(icent + 1)->SetLeftMargin(0.10);
  //   crcppipt->GetPad(icent + 1)->SetTicks(1, 1);

  //   crcppipt->cd(icent + 1);
  //   haxis_rcp_pipT->GetXaxis()->SetRangeUser(0, 20);
  //   haxis_rcp_pipT->Draw();

  //   for (int isys = 0; isys < NSYSTEMS; isys++)
  //     grcp_NcollModABBCscMod_pipT[isys][icent]->Draw("C");

  //   l1.DrawLine(0, 1, 20., 1);
  // }

  // TCanvas *cbbc = new TCanvas("cbbc", "bbc", 600, 1200);
  // cbbc->SetTopMargin(0.00);
  // cbbc->SetRightMargin(0.00);
  // cbbc->SetBottomMargin(0.00);
  // cbbc->SetLeftMargin(0.00);
  // cbbc->Divide(1, NSYSTEMS, 0, 0);

  // for (int isys = 0; isys < NSYSTEMS; isys++)
  // {
  //   cbbc->GetPad(isys + 1)->SetTopMargin(0.02);
  //   cbbc->GetPad(isys + 1)->SetRightMargin(0.02);
  //   cbbc->GetPad(isys + 1)->SetBottomMargin(0.10);
  //   cbbc->GetPad(isys + 1)->SetLeftMargin(0.10);
  //   cbbc->GetPad(isys + 1)->SetTicks(1, 1);

  //   cbbc->cd(isys + 1);
  //   gPad->SetLogy();
  //   hBBCsc[isys][0]->GetYaxis()->SetRangeUser(1e2, 2e4);
  //   hBBCsc[isys][0]->GetXaxis()->SetRangeUser(1, 80);
  //   hBBCsc[isys][0]->SetTitle(";Q_{BBC,Au}");
  //   hBBCsc[isys][0]->DrawCopy("hist");

  //   for (int i = 0; i < NCENT; i++)
  //     hBBCsc_cent[isys][i]->Draw("hist same");

  //   label.DrawLatex(0.7, 0.7, collSystem[isys]);
  // }

  // TCanvas *cyield = new TCanvas("cyield", "yield", 1400, 400);
  // cyield->Divide(NSYSTEMS, 1, 0, 0);
  // label.DrawLatex(0.5, 0.96, "BBCsc charge for events with high-p_{T} particle");
  // for (int isys = 0; isys < NSYSTEMS; isys++)
  // {
  //   cyield->GetPad(isys + 1)->SetTopMargin(0.02);
  //   cyield->GetPad(isys + 1)->SetRightMargin(0.02);
  //   cyield->GetPad(isys + 1)->SetBottomMargin(0.10);
  //   cyield->GetPad(isys + 1)->SetLeftMargin(0.10);
  //   cyield->GetPad(isys + 1)->SetTicks(1, 1);

  //   cyield->cd(isys + 1);
  //   gPad->SetLogy();
  //   hBBCscNcollA[isys][xplot]->GetYaxis()->SetRangeUser(1e2, 2e4);
  //   hBBCscNcollA[isys][xplot]->GetXaxis()->SetRangeUser(1, 80);
  //   hBBCscNcollA[isys][xplot]->SetTitle(";Q_{BBC,Au}");
  //   hBBCscNcollA[isys][xplot]->DrawCopy("hist");
  //   // hBBCscNcollModA[isys][xplot]->Draw("same");
  //   // hBBCscModNcollA[isys][xplot]->Draw("same");
  //   hBBCscModNcollModA[isys][xplot]->DrawCopy("hist same");
  //   legBBC[isys]->Draw("same");
  // }

  // TCanvas* cyieldpau = new TCanvas("cyieldpau", "yield pAu", 600, 400);
  // cyieldpau->SetTopMargin(0.02);
  // cyieldpau->SetRightMargin(0.02);
  // cyieldpau->SetBottomMargin(0.10);
  // cyieldpau->SetLeftMargin(0.10);
  // cyieldpau->SetTicks(1, 1);

  // cyieldpau->cd(1);
  // // gPad->SetLogy();
  // // hBBCscNcollA[0][0]->GetYaxis()->SetRangeUser(1e2, 2e4);
  // hBBCscNcollA[0][0]->GetXaxis()->SetRangeUser(1, 30);
  // hBBCscNcollA[0][0]->SetLineWidth(lineWidth);
  // hBBCscNcollA[0][0]->SetLineColor(kBlack);
  // hBBCscNcollA[0][0]->SetTitle(";Q_{BBC,Au}");
  // hBBCscNcollA[0][0]->DrawCopy("hist");
  // for (int ix = 0; ix < NX; ix++)
  // {
  //   hBBCscModNcollModA[0][ix]->SetLineWidth(1);
  //   hBBCscModNcollModA[0][ix]->SetLineColor(2 + ix);
  //   hBBCscModNcollModA[0][ix]->DrawCopy("hist same");
  // }
  // legx->Draw("same");

  // TCanvas *cncollmb = new TCanvas("cncollmb", "MB Ncoll", 1400, 400);
  // cncollmb->Divide(NSYSTEMS, 1, 0, 0);
  // for (int isys = 0; isys < NSYSTEMS; isys++)
  // {
  //   cncollmb->GetPad(isys + 1)->SetTopMargin(0.02);
  //   cncollmb->GetPad(isys + 1)->SetRightMargin(0.02);
  //   cncollmb->GetPad(isys + 1)->SetBottomMargin(0.10);
  //   cncollmb->GetPad(isys + 1)->SetLeftMargin(0.10);
  //   cncollmb->GetPad(isys + 1)->SetTicks(1, 1);

  //   cncollmb->cd(isys + 1);
  //   gPad->SetLogy();
  //   hNcoll_MB[isys][xplot]->GetXaxis()->SetRangeUser(0, 50);
  //   hNcoll_MB[isys][xplot]->Draw();
  //   hNcollMod_MB[isys][xplot]->Draw("same");
  //   legNcollMB[isys]->Draw("same");
  // }

  // TCanvas *cncoll = new TCanvas("cncoll", "Ncoll", 400 * NSYSTEMS, 200 * NCENT);
  // cncoll->Divide(NSYSTEMS, NCENT, 0, 0);
  // for (int isys = 0; isys < NSYSTEMS; isys++)
  // {
  //   for (int icent = 0; icent < NCENT; icent++)
  //   {
  //     cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetTopMargin(0.02);
  //     cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetRightMargin(0.02);
  //     cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetBottomMargin(0.10);
  //     cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetLeftMargin(0.10);
  //     cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetTicks(1, 1);

  //     cncoll->cd(icent * NSYSTEMS + isys + 1);
  //     gPad->SetLogy();
  //     hNcoll_cent[isys][xplot][icent]->GetXaxis()->SetRangeUser(0, 50);
  //     hNcoll_cent[isys][xplot][icent]->GetYaxis()->SetRangeUser(1, 1e5);
  //     hNcoll_cent[isys][xplot][icent]->Draw();
  //     hNcollMod_cent[isys][xplot][icent]->Draw("same");
  //     legNcoll[isys][icent]->Draw("same");
  //   }
  // }

  //=====================================================//
  // FIGURE OBJECTS FOR PAPER
  //=====================================================//

  const char *lsystem[] = {"p+Au", "d+Au", "^{3}He+Au"};
  TLegend *legrcp_paper[NCENT - 1];
  for (int icent = 0; icent < NCENT - 1; icent++)
  {
    double x1, y1;
    if (icent == NCENT - 2)
    {
      x1 = 0.25;
      y1 = 0.2;
    }
    else
    {
      x1 = 0.25;
      y1 = 0.05;
    }
    legrcp_paper[icent] = new TLegend(x1, y1, x1 + 0.25, y1 + 0.4,
                                      Form("(%.0f-%.0f%%) / (%.0f-%.0f%%)",
                                           centl[1][icent], centl[1][icent + 1],
                                           centl[1][NCENT - 1], centl[1][NCENT]));
    legrcp_paper[icent]->SetFillStyle(0);
    legrcp_paper[icent]->SetBorderSize(0);
    legrcp_paper[icent]->SetTextSize(0.06);
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      legrcp_paper[icent]->AddEntry(grcp[isys][icent],
                                    lsystem[isys], "L");
    }
    legrcp_paper[icent]->AddEntry(gRCP_jet[icent], "PHENIX d+Au", "P");
    legrcp_paper[icent]->AddEntry((TObject*)0, "arXiv:1509.04657", "");

  }

  TLegend *legraa_paper = new TLegend(0.2, 0.05, 0.45, 0.4);
  legraa_paper->SetFillStyle(0);
  legraa_paper->SetBorderSize(0);
  legraa_paper->SetTextSize(0.06);
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    legraa_paper->AddEntry(graa[isys][0],
                           lsystem[isys], "L");
  }
  legraa_paper->AddEntry(gRdAu_jet[0], "PHENIX d+Au", "P");
  legraa_paper->AddEntry((TObject*)0, " arXiv:1509.04657", "");


  TLegend *legjda_paper = new TLegend(0.4, 0.6, 0.95, 0.92);
  legjda_paper->SetFillStyle(0);
  legjda_paper->SetBorderSize(0);
  legjda_paper->SetTextSize(0.04);
  legjda_paper->AddEntry(graa[1][0], lsystem[1], "L");
  legjda_paper->AddEntry(graa_linMod[1][0],
                         Form("%s #times Linear CNM Effects", lsystem[1]),
                         "L");
  // legjda_paper->AddEntry((TObject*)0, "PHENIX d+Au", "P");
  legjda_paper->AddEntry((TObject*)0, "Phys.Rev.Lett. 107 (2011) 172301", "");
  for (int ix = 0; ix < NXBINS; ix++)
  {
    legjda_paper->AddEntry(gJdA[ix][0],
                           // Form("PHENIX d+Au %.1e < x^{frag}_{Au} < %.1e", xbins[ix], xbins[ix + 1]),
                           Form("PHENIX d+Au %.1e < x^{frag}_{Au} < %.1e", xmin_xbins[ix], xmax_xbins[ix]),
                           "P");
  }
  // legjda_paper->AddEntry(gJdA_hpt[0], "PHENIX d+Au, p_{T}^{trig}>1 & p_{T}^{assoc}>1", "P");


  TLegend *legQ_paper = new TLegend(0.20, 0.25, 0.4, 0.5);
  legQ_paper->SetFillStyle(0);
  legQ_paper->SetBorderSize(0);
  legQ_paper->SetTextSize(0.05);
  for (int isys = 0; isys < NSYSTEMS; isys++)
    legQ_paper->AddEntry(gQ[isys], lsystem[isys], "L");



  TH1F *haxis_sig = new TH1F("haxis_sig", ";x_{p};#sigma(x_{p}) / #sigma_{NN}", 100, 0, 1);
  haxis_sig->SetMinimum(0);
  haxis_sig->SetMaximum(1);
  haxis_sig->GetYaxis()->SetTitleOffset(1.15);
  haxis_sig->GetYaxis()->CenterTitle();
  haxis_sig->GetYaxis()->SetLabelSize(0.05);
  haxis_sig->GetYaxis()->SetTitleSize(0.06);
  haxis_sig->GetYaxis()->SetNdivisions(6, 3, 0);
  haxis_sig->GetXaxis()->SetTitleOffset(0.9);
  haxis_sig->GetXaxis()->SetLabelSize(0.05);
  haxis_sig->GetXaxis()->SetTitleSize(0.06);
  haxis_sig->GetXaxis()->SetNdivisions(5, 4, 0);


  TH1F* haxis_rcp_paper = new TH1F("haxis_rcp_paper",
                                   ";x_{p};R_{CP}",
                                   100, 0, 1);
  haxis_rcp_paper->GetYaxis()->SetTitleOffset(0.95);
  haxis_rcp_paper->GetYaxis()->CenterTitle();
  haxis_rcp_paper->GetYaxis()->SetLabelSize(0.06);
  haxis_rcp_paper->GetYaxis()->SetTitleSize(0.08);
  haxis_rcp_paper->GetYaxis()->SetNdivisions(6, 3, 0);
  haxis_rcp_paper->GetXaxis()->SetTitleOffset(0.9);
  haxis_rcp_paper->GetXaxis()->SetLabelSize(0.06);
  haxis_rcp_paper->GetXaxis()->SetTitleSize(0.08);
  haxis_rcp_paper->GetXaxis()->SetNdivisions(5, 4, 0);
  haxis_rcp_paper->SetMinimum(0.01);
  haxis_rcp_paper->SetMaximum(1.19);

  TH1F* haxis_raa_paper = new TH1F("haxis_raa_paper",
                                   ";x_{p};R_{AA}",
                                   100, 0, 1);
  haxis_raa_paper->GetYaxis()->SetTitleOffset(0.9);
  haxis_raa_paper->GetYaxis()->CenterTitle();
  haxis_raa_paper->GetYaxis()->SetLabelSize(0.06);
  haxis_raa_paper->GetYaxis()->SetTitleSize(0.08);
  haxis_raa_paper->GetYaxis()->SetNdivisions(6, 3, 0);
  haxis_raa_paper->GetXaxis()->SetTitleOffset(0.9);
  haxis_raa_paper->GetXaxis()->SetLabelSize(0.06);
  haxis_raa_paper->GetXaxis()->SetTitleSize(0.08);
  haxis_raa_paper->GetXaxis()->SetNdivisions(5, 4, 0);
  haxis_raa_paper->SetMinimum(0.01);
  haxis_raa_paper->SetMaximum(1.19);


  TH1F* haxis_jda_paper = new TH1F("haxis_jda_paper",
                                   ";x_{p};J_{dA}",
                                   100, 0, 1);
  haxis_jda_paper->GetYaxis()->SetTitleOffset(0.9);
  haxis_jda_paper->GetYaxis()->CenterTitle();
  haxis_jda_paper->GetYaxis()->SetLabelSize(0.06);
  haxis_jda_paper->GetYaxis()->SetTitleSize(0.08);
  haxis_jda_paper->GetYaxis()->SetNdivisions(6, 3, 0);
  haxis_jda_paper->GetXaxis()->SetTitleOffset(0.9);
  haxis_jda_paper->GetXaxis()->SetLabelSize(0.06);
  haxis_jda_paper->GetXaxis()->SetTitleSize(0.08);
  haxis_jda_paper->GetXaxis()->SetNdivisions(5, 4, 0);
  haxis_jda_paper->SetMinimum(0.01);
  haxis_jda_paper->SetMaximum(1.99);

  TH1F* haxis_jcp_paper = new TH1F("haxis_jcp_paper",
                                   ";x_{p};J_{cp}",
                                   100, 0, 1);
  haxis_jcp_paper->GetYaxis()->SetTitleOffset(0.9);
  haxis_jcp_paper->GetYaxis()->CenterTitle();
  haxis_jcp_paper->GetYaxis()->SetLabelSize(0.06);
  haxis_jcp_paper->GetYaxis()->SetTitleSize(0.08);
  haxis_jcp_paper->GetYaxis()->SetNdivisions(6, 3, 0);
  haxis_jcp_paper->GetXaxis()->SetTitleOffset(0.9);
  haxis_jcp_paper->GetXaxis()->SetLabelSize(0.06);
  haxis_jcp_paper->GetXaxis()->SetTitleSize(0.08);
  haxis_jcp_paper->GetXaxis()->SetNdivisions(5, 4, 0);
  haxis_jcp_paper->SetMinimum(0.01);
  haxis_jcp_paper->SetMaximum(1.99);

  TGaxis *gx_paper = new TGaxis(xp_JdA[0], 1.09, xp_JdA[8], 1.09,
                                xAu_frag[0], xAu_frag[8],
                                510, "-S");
  gx_paper->SetTickLength(0.15);
  gx_paper->SetTitle("x_{Au}^{frag}");
  gx_paper->SetLineColor(kRed);
  gx_paper->SetTitleColor(kRed);
  gx_paper->SetLabelColor(kRed);
  gx_paper->SetLabelSize(0.05);
  gx_paper->SetTitleSize(0.06);
  gx_paper->SetLabelOffset(-0.01);
  gx_paper->SetTitleOffset(0.8);
  gx_paper->CenterTitle();
  gx_paper->SetNoExponent(kTRUE);
  gx_paper->SetNdivisions(802);


  TH1F* haxis_Q_paper = new TH1F("haxis_Q_paper",
                                 ";x_{p};<dN^{hard}/dQ(x_{p})>/<dN^{hard}/dQ(x_{p}=0)>",
                                 100, 0, 1);
  haxis_Q_paper->GetYaxis()->SetTitleOffset(0.9);
  haxis_Q_paper->GetYaxis()->CenterTitle();
  haxis_Q_paper->GetYaxis()->SetLabelSize(0.04);
  haxis_Q_paper->GetYaxis()->SetTitleSize(0.05);
  haxis_Q_paper->GetYaxis()->SetNdivisions(10, 2, 0);
  haxis_Q_paper->GetXaxis()->SetTitleOffset(0.9);
  haxis_Q_paper->GetXaxis()->SetLabelSize(0.04);
  haxis_Q_paper->GetXaxis()->SetTitleSize(0.05);
  haxis_Q_paper->GetXaxis()->SetNdivisions(10, 2, 0);
  haxis_Q_paper->SetMinimum(0.41);
  haxis_Q_paper->SetMaximum(1.09);

  TLatex lpaper;
  lpaper.SetNDC();
  lpaper.SetTextAlign(22);
  lpaper.SetTextSize(0.08);

  double topMargin = 0.05;
  double bottomMargin = 0.15;
  double yl, yh, fracArea;

  const char *plabel[4] = {"(a)", "(b)", "(c)", "(d)"};


  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);

  TLine l1;
  l1.SetLineStyle(2);

  //=====================================================//
  // FIGURES FOR PAPER
  //=====================================================//

  TCanvas *crcp_paper = new TCanvas("crcp_paper", "RCP", 500, 1200);
  crcp_paper->SetTopMargin(0.0);
  crcp_paper->SetRightMargin(0.0);
  crcp_paper->SetBottomMargin(0.0);
  crcp_paper->SetLeftMargin(0.0);

  TPad *prcp_paper[NCENT - 1];
  fracArea = 1. / (float)((NCENT - 1) - 2 + 1. / (1 - topMargin) + 1. / (1 - bottomMargin));
  for (int icent = 0; icent < NCENT - 1; icent++)
  {
    if (icent == 0)
    {
      yl = 1.0 - fracArea * (1. / (1. - topMargin));
      yh = 1;
    }
    else if (icent == NCENT - 2)
    {
      yl = 0;
      yh = 1.0 - (1. / (1. - topMargin) + icent - 1) * fracArea;
    }
    else
    {
      yl = 1.0 - (1. / (1. - topMargin) + icent) * fracArea;
      yh = 1.0 - (1. / (1. - topMargin) + icent - 1) * fracArea;
    }
    cout << " " << icent << " yl:" << yl
         << " yh:" << yh
         << " h:" << yh - yl
         << endl;

    prcp_paper[icent] = new TPad(Form("prcp_paper_%i", icent), "", 0, yl, 1, yh);
    prcp_paper[icent]->SetLeftMargin(0.15);
    prcp_paper[icent]->SetRightMargin(0.05);
    if (icent == 0)
      prcp_paper[icent]->SetTopMargin(topMargin);
    else
      prcp_paper[icent]->SetTopMargin(0.0);
    if (icent == NCENT - 2)
      prcp_paper[icent]->SetBottomMargin(bottomMargin);
    else
      prcp_paper[icent]->SetBottomMargin(0.0);
    prcp_paper[icent]->SetTicks(1, 1);
    prcp_paper[icent]->Draw();

    prcp_paper[icent]->cd();
    haxis_rcp_paper->GetXaxis()->SetRangeUser(0, 0.6);
    haxis_rcp_paper->Draw();

    for (int i = 0; i < NJET; i++)
      bRcp_jet[icent][i]->Draw();
    gRCP_jet[icent]->Draw("P");


    for (int isys = 0; isys < NSYSTEMS; isys++)
      grcp[isys][icent]->Draw("C");


    l1.DrawLine(0, 1, 0.5, 1);
    lpaper.SetTextSize(0.06 * fracArea / (yh - yl));
    if (icent == NCENT - 2)
      legrcp_paper[icent]->Draw("same");
    else if (icent == 0)
    {
      lpaper.DrawLatex(0.38, 0.88, Form("(%.0f-%.0f%%) / (%.0f-%.0f%%)",
                                        centl[1][icent], centl[1][icent + 1],
                                        centl[1][NCENT - 1], centl[1][NCENT]));
    }
    else
    {
      lpaper.DrawLatex(0.38, 0.90, Form("(%.0f-%.0f%%) / (%.0f-%.0f%%)",
                                        centl[1][icent], centl[1][icent + 1],
                                        centl[1][NCENT - 1], centl[1][NCENT]));
    }
    if (icent == 0)
      lpaper.DrawLatex(0.9, 0.88, plabel[icent]);
    else
      lpaper.DrawLatex(0.9, 0.93, plabel[icent]);

    crcp_paper->cd();
  }


  TCanvas *craa_paper = new TCanvas("craa_paper", "RAA", 400, 1000);
  craa_paper->SetTopMargin(0.0);
  craa_paper->SetRightMargin(0.0);
  craa_paper->SetBottomMargin(0.0);
  craa_paper->SetLeftMargin(0.0);

  TPad *praa_paper[NCENT];
  fracArea = 1. / (float)(NCENT - 2 + 1. / (1 - topMargin) + 1. / (1 - bottomMargin));
  double raal[NCENT] = {0.01, 0.41, 0.61, 0.76};
  double raah[NCENT] = {1.19, 1.59, 1.59, 2.49};
  for (int icent = 0; icent < NCENT; icent++)
  {

    if (icent == 0)
    {
      yl = 1.0 - fracArea * (1. / (1. - topMargin));
      yh = 1;
    }
    else if (icent == NCENT - 1)
    {
      yl = 0;
      yh = 1.0 - (1. / (1. - topMargin) + icent - 1) * fracArea;
    }
    else
    {
      yl = 1.0 - (1. / (1. - topMargin) + icent) * fracArea;
      yh = 1.0 - (1. / (1. - topMargin) + icent - 1) * fracArea;
    }
    cout << " " << icent << " yl:" << yl
         << " yh:" << yh
         << " h:" << yh - yl
         << endl;

    praa_paper[icent] = new TPad(Form("praa_paper_%i", icent), "", 0, yl, 1, yh);
    praa_paper[icent]->SetLeftMargin(0.15);
    praa_paper[icent]->SetRightMargin(0.05);
    if (icent == 0)
      praa_paper[icent]->SetTopMargin(topMargin);
    else
      praa_paper[icent]->SetTopMargin(0.0);
    if (icent == NCENT - 1)
      praa_paper[icent]->SetBottomMargin(bottomMargin);
    else
      praa_paper[icent]->SetBottomMargin(0.0);
    praa_paper[icent]->SetTicks(1, 1);
    praa_paper[icent]->Draw();

    praa_paper[icent]->cd();
    haxis_raa_paper->GetYaxis()->SetRangeUser(raal[icent], raah[icent]);
    haxis_raa_paper->GetXaxis()->SetRangeUser(0, 0.6);
    haxis_raa_paper->DrawCopy();

    for (int i = 0; i < NJET; i++)
      bRdAu_jet[icent][i]->Draw();
    gRdAu_jet[icent]->Draw("P");

    for (int isys = 0; isys < NSYSTEMS; isys++)
      graa[isys][icent]->Draw("C");

    l1.DrawLine(0, 1, 1.0, 1);
    if (icent == 0)
      legraa_paper->Draw("same");
    lpaper.SetTextSize(0.06 * fracArea / (yh - yl));
    if (icent == NCENT - 1)
      lpaper.DrawLatex(0.9, 0.22, plabel[icent]);
    else
      lpaper.DrawLatex(0.9, 0.07, plabel[icent]);
    if (icent == 0)
      lpaper.DrawLatex(0.3, 0.85, Form("%.0f-%.0f%%",
                                       centl[1][icent], centl[1][icent + 1]));
    else
      lpaper.DrawLatex(0.3, 0.90, Form("%.0f-%.0f%%",
                                       centl[1][icent], centl[1][icent + 1]));
    craa_paper->cd();

  }

  TCanvas *csig_paper = new TCanvas("csig_paper", "sig mod", 500, 500);
  csig_paper->SetTopMargin(0.05);
  csig_paper->SetRightMargin(0.05);
  csig_paper->SetBottomMargin(0.15);
  csig_paper->SetLeftMargin(0.15);
  csig_paper->SetTicks(1, 1);

  csig_paper->cd(1);
  haxis_sig->Draw();
  fsig->Draw("same");


  TCanvas *cjda_paper = new TCanvas("cjda_paper", "jda", 500, 1000);
  cjda_paper->SetTopMargin(0);
  cjda_paper->SetRightMargin(0);
  cjda_paper->SetBottomMargin(0);
  cjda_paper->SetLeftMargin(0);

  TPad *pjda_paper[3];
  fracArea = 1. / (float)(3 - 2 + 1. / (1 - topMargin) + 1. / (1 - bottomMargin));

  //-- 0-20%
  cjda_paper->cd();
  yl = 1.0 - fracArea * (1. / (1. - topMargin));
  yh = 1;
  cout << " jdA 0 - yl:" << yl << " yh:" << yh <<  "fracArea:" << fracArea << endl;
  pjda_paper[0] = new TPad("pjda_paper_0", "", 0, yl, 1, yh);
  pjda_paper[0]->SetLeftMargin(0.15);
  pjda_paper[0]->SetRightMargin(0.05);
  pjda_paper[0]->SetTopMargin(topMargin);
  pjda_paper[0]->SetBottomMargin(0);
  pjda_paper[0]->SetTicks(1, 1);
  pjda_paper[0]->Draw();
  pjda_paper[0]->cd();
  haxis_jda_paper->GetYaxis()->SetRangeUser(0.01, 1.09);
  haxis_jda_paper->DrawCopy();

  for (int ix = 0; ix < NXBINS; ix++)
  {
    for (int j = 0; j < nxjda[ix]; j++)
      bJdA[ix][0][j]->Draw("same");
    gJdA[ix][0]->Draw("P");
  }

  graa[1][0]->Draw("C");
  graa_linMod[1][0]->Draw("C");

  l1.DrawLine(0, 1, 1.0, 1);
  label.SetTextSize(0.06 * fracArea / (yh - yl));
  label.DrawLatex(0.3, 0.1, "0-20%");
  label.DrawLatex(0.18, 0.95 - topMargin, "(a)");
  legjda_paper->Draw("same");

  //-- 60-88%
  cjda_paper->cd();
  yl = 1.0 - (1. / (1. - topMargin) + 1) * fracArea;
  yh = 1.0 - (1. / (1. - topMargin) + 1 - 1) * fracArea;
  cout << " jdA 1 - yl:" << yl << " yh:" << yh <<  "fracArea:" << fracArea << endl;
  pjda_paper[1] = new TPad("pjda_paper_1", "", 0, yl, 1, yh);
  pjda_paper[1]->SetLeftMargin(0.15);
  pjda_paper[1]->SetRightMargin(0.05);
  pjda_paper[1]->SetTopMargin(0);
  pjda_paper[1]->SetBottomMargin(0);
  pjda_paper[1]->SetTicks(1, 1);
  pjda_paper[1]->Draw();
  pjda_paper[1]->cd();
  haxis_jda_paper->GetYaxis()->SetRangeUser(0.01, 1.99);
  haxis_jda_paper->DrawCopy();

  for (int ix = 0; ix < NXBINS; ix++)
  {
    for (int j = 0; j < nxjda[ix]; j++)
      bJdA[ix][1][j]->Draw("same");
    gJdA[ix][1]->Draw("P");
  }

  graa[1][NCENT - 1]->Draw("C");
  graa_linMod[1][NCENT - 1]->Draw("C");

  l1.DrawLine(0, 1, 1.0, 1);
  label.SetTextSize(0.06 * fracArea / (yh - yl));
  label.DrawLatex(0.3, 0.1, "60-88%");
  label.DrawLatex(0.18, 0.95, "(b)");

  //-- Jcp
  cjda_paper->cd();
  yl = 0;
  yh = 1.0 - (1. / (1. - topMargin) + 2 - 1) * fracArea;
  cout << " jdA 2 - yl:" << yl << " yh:" << yh <<  "fracArea:" << fracArea << endl;
  pjda_paper[2] = new TPad("pjda_paper_2", "", 0, yl, 1, yh);
  pjda_paper[2]->SetLeftMargin(0.15);
  pjda_paper[2]->SetRightMargin(0.05);
  pjda_paper[2]->SetTopMargin(0);
  pjda_paper[2]->SetBottomMargin(bottomMargin);
  pjda_paper[2]->SetTicks(1, 1);
  pjda_paper[2]->Draw();
  pjda_paper[2]->cd();
  haxis_jcp_paper->GetYaxis()->SetRangeUser(0.01, 1.09);
  haxis_jcp_paper->DrawCopy();

  for (int ix = 0; ix < NXBINS; ix++)
  {
    for (int j = 0; j < nxjda[ix]; j++)
      bJcp[ix][j]->Draw("same");
    gJcp[ix]->Draw("P");
  }


  grcp[1][0]->Draw("same");
  grcp_linMod[1][0]->Draw("same");

  l1.DrawLine(0, 1, 1.0, 1);
  label.SetTextSize(0.06 * fracArea / (yh - yl));
  label.DrawLatex(0.3, 0.1 + bottomMargin, "(0-20%)/(60-88%)");
  label.DrawLatex(0.18, 0.95, "(c)");



  TCanvas *cq_paper = new TCanvas("cq_paper", "<Q>", 900, 600);
  cq_paper->SetTopMargin(0.02);
  cq_paper->SetRightMargin(0.02);
  cq_paper->SetBottomMargin(0.10);
  cq_paper->SetLeftMargin(0.10);
  cq_paper->SetTicks(1, 1);

  cq_paper->cd(1);
  haxis_Q_paper->Draw();

  for (int isys = 0; isys < NSYSTEMS; isys++)
    gQ[isys]->Draw("C");

  l1.DrawLine(0, 1, 1.0, 1);
  legQ_paper->Draw("same");


//=====================================================//
// SAVE
//=====================================================//
  if (printPlots)
  {
    cout << endl;
    cout << "--> Printing plots to pdf" << endl;

    crcp_paper->Print("Rcp_paper.pdf");
    craa_paper->Print("RAA_paper.pdf");
    csig_paper->Print("sigmod.pdf");
    cq_paper->Print("Q_paper.pdf");
    cjda_paper->Print("JdA_paper.pdf");
    // cjda_rat_paper->Print("JdA_rat_paper.pdf");
  }

}