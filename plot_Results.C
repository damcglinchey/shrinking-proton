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

  bool printPlots = true;

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


  //<rT> values for each centrality bin
  // (From Jamie)
  float mean_rT[] = {3.3302, 4.12443, 4.83264, 5.74247};

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

  const int NXBINS = 3;
  float xbins[] = {1e-4, 1e-3, 1e-2, 1e-1};
  float xmin_xbins[NXBINS] = {0};
  float xmax_xbins[NXBINS] = {0};
  int xbinColor[] = {kBlack, kMagenta + 2, kOrange + 2};
  TGraphErrors *gJdA[NXBINS][4];
  TBox *bJdA[NXBINS][4][NPPG128];
  int nxjda[NXBINS] = {0};
  for (int ix = 0; ix < NXBINS; ix++)
  {
    //-- JdA
    for (int k = 0; k < 4; k++)
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
  }// ix

  //find the min and max x bins
  for (int ib = 0; ib < NXBINS; ib++)
  {
    //get the max
    xmax_xbins[ib] = 0;
    for (int ix = 0; ix < NPPG128; ix++)
    {
      if (xAu_frag[ix] > xbins[ib] && xAu_frag[ix] < xbins[ib + 1])
      {
        if (xAu_frag[ix] > xmax_xbins[ib])
          xmax_xbins[ib] = xAu_frag[ix];
      }
    }

    //get the min
    xmin_xbins[ib] = 1;
    for (int ix = 0; ix < NPPG128; ix++)
    {
      if (xAu_frag[ix] > xbins[ib] && xAu_frag[ix] < xbins[ib + 1])
      {
        if (xAu_frag[ix] < xmin_xbins[ib])
          xmin_xbins[ib] = xAu_frag[ix];
      }
    }
  }

  //=====================================================//
  // FIGURE OBJECTS FOR PAPER
  //=====================================================//

  int fontSizeLabel = 16;
  int fontSizeTitle = 20;
  int fontSizeLegend = 16;

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
    legrcp_paper[icent] = new TLegend(x1, y1, x1 + 0.25, y1 + 0.4);
    legrcp_paper[icent]->SetFillStyle(0);
    legrcp_paper[icent]->SetBorderSize(0);
    legrcp_paper[icent]->SetTextFont(63);
    legrcp_paper[icent]->SetTextSize(fontSizeTitle);
    legrcp_paper[icent]->SetHeader(Form("(%.0f-%.0f%%) / (%.0f-%.0f%%)",
                                        centl[1][icent], centl[1][icent + 1],
                                        centl[1][NCENT - 1], centl[1][NCENT]));
    // legrcp_paper[icent]->SetTextSize(0.06);
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      legrcp_paper[icent]->AddEntry(grcp[isys][icent],
                                    lsystem[isys], "L");
    }
    legrcp_paper[icent]->AddEntry(gRCP_jet[icent], "PHENIX d+Au", "P");
    // legrcp_paper[icent]->AddEntry((TObject*)0, "arXiv:1509.04657", "");

  }

  TLegend *legraa_paper = new TLegend(0.2, 0.05, 0.45, 0.45);
  legraa_paper->SetFillStyle(0);
  legraa_paper->SetBorderSize(0);
  legraa_paper->SetTextFont(63);
  legraa_paper->SetTextSize(fontSizeLegend);
  // legraa_paper->SetTextSize(0.06);
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    legraa_paper->AddEntry(graa[isys][0],
                           lsystem[isys], "L");
  }
  legraa_paper->AddEntry(gRdAu_jet[0], "PHENIX d+Au", "P");
  // legraa_paper->AddEntry((TObject*)0, " arXiv:1509.04657", "");


  TLegend *legjda_model_paper = new TLegend(0.45, 0.55, 0.90, 0.95);
  legjda_model_paper->SetFillStyle(0);
  legjda_model_paper->SetBorderSize(0);
  legjda_model_paper->SetTextFont(63);
  legjda_model_paper->SetTextSize(fontSizeLegend);
  legjda_model_paper->SetHeader("d+Au #sqrt{s_{NN}}=200 GeV");
  legjda_model_paper->AddEntry(graa[1][0], "High-x_{p}", "L");
  legjda_model_paper->AddEntry(graa_linMod[1][0],
                               "High-x_{p} #times low-x_{Au}",
                               "L");

  TLegend *legjda_data_paper = new TLegend(0.32, 0.6, 0.96, 0.85);
  legjda_data_paper->SetFillStyle(0);
  legjda_data_paper->SetBorderSize(0);
  legjda_data_paper->SetTextFont(63);
  legjda_data_paper->SetTextSize(fontSizeLegend);
  // legjda_data_paper->SetTextSize(14);
  // legjda_data_paper->SetHeader("Phys.Rev.Lett. 107 (2011) 172301");
  // legjda_data_paper->SetHeader("PHENIX (PRL 107 (2011) 172301)");
  // legjda_data_paper->SetHeader("PHENIX");
  legjda_data_paper->SetNColumns(2);
  legjda_data_paper->SetEntrySeparation(0.02);
  // legjda_data_paper->AddEntry((TObject*)0, "PHENIX", "");
  // legjda_data_paper->AddEntry((TObject*)0, "Phys.Rev.Lett. 107 (2011) 172301", "");
  // for (int ix = 0; ix < NXBINS; ix++)
  // {
  //   legjda_data_paper->AddEntry(gJdA[ix][0],
  //                               // Form("%.1e < x^{frag}_{Au} < %.1e", xmin_xbins[ix], xmax_xbins[ix]),
  //                               Form("x^{frag}_{Au} ~ %.0e", 0.5*(xmin_xbins[ix] + xmax_xbins[ix])),
  //                               "P");
  // }
  legjda_data_paper->AddEntry(gJdA[0][0], "x_{Au}^{frag} ~ 5x10^{4}", "P");
  legjda_data_paper->AddEntry(gJdA[1][0], "x_{Au}^{frag} ~ 5x10^{3}", "P");
  legjda_data_paper->AddEntry(gJdA[2][0], "x_{Au}^{frag} ~ 2x10^{2}", "P");

  TLegend *legQ_paper = new TLegend(0.20, 0.25, 0.4, 0.5);
  legQ_paper->SetFillStyle(0);
  legQ_paper->SetBorderSize(0);
  legQ_paper->SetTextFont(63);
  legQ_paper->SetTextSize(fontSizeLegend);
  // legQ_paper->SetTextSize(0.05);
  for (int isys = 0; isys < NSYSTEMS; isys++)
    legQ_paper->AddEntry(gQ[isys], lsystem[isys], "L");



  TH1F *haxis_sig = new TH1F("haxis_sig", ";x_{p};#sigma(x_{p}) / #sigma_{NN}", 100, 0, 1);
  haxis_sig->SetMinimum(0);
  haxis_sig->SetMaximum(1);
  haxis_sig->GetYaxis()->SetTitleOffset(1.5);
  haxis_sig->GetYaxis()->CenterTitle();
  haxis_sig->GetYaxis()->SetTitleFont(63);
  haxis_sig->GetYaxis()->SetLabelFont(63);
  haxis_sig->GetYaxis()->SetLabelSize(fontSizeLabel);
  haxis_sig->GetYaxis()->SetTitleSize(fontSizeTitle);
  haxis_sig->GetYaxis()->SetNdivisions(6, 3, 0);
  haxis_sig->GetXaxis()->SetTitleOffset(1.5);
  haxis_sig->GetXaxis()->SetTitleFont(63);
  haxis_sig->GetXaxis()->SetLabelFont(63);
  haxis_sig->GetXaxis()->SetLabelSize(fontSizeLabel);
  haxis_sig->GetXaxis()->SetTitleSize(fontSizeTitle);
  haxis_sig->GetXaxis()->SetNdivisions(5, 4, 0);


  TH1F* haxis_rcp_paper = new TH1F("haxis_rcp_paper",
                                   ";x_{p};R_{CP}",
                                   100, 0, 1);
  haxis_rcp_paper->GetYaxis()->SetTitleOffset(2.5);
  haxis_rcp_paper->GetYaxis()->CenterTitle();
  haxis_rcp_paper->GetYaxis()->SetTitleFont(63);
  haxis_rcp_paper->GetYaxis()->SetLabelFont(63);
  haxis_rcp_paper->GetYaxis()->SetLabelSize(fontSizeLabel);
  haxis_rcp_paper->GetYaxis()->SetTitleSize(fontSizeTitle);
  haxis_rcp_paper->GetYaxis()->SetNdivisions(6, 3, 0);
  haxis_rcp_paper->GetXaxis()->SetTitleOffset(2.5);
  haxis_rcp_paper->GetXaxis()->SetTitleFont(63);
  haxis_rcp_paper->GetXaxis()->SetLabelFont(63);
  haxis_rcp_paper->GetXaxis()->SetLabelSize(fontSizeLabel);
  haxis_rcp_paper->GetXaxis()->SetTitleSize(fontSizeTitle);
  haxis_rcp_paper->GetXaxis()->SetNdivisions(5, 4, 0);
  haxis_rcp_paper->SetMinimum(0.01);
  haxis_rcp_paper->SetMaximum(1.19);

  TH1F* haxis_raa_paper = new TH1F("haxis_raa_paper",
                                   ";x_{p} = 2p_{T}^{jet}/#sqrt{s_{NN}};R_{p/d/^{3}He+Au}",
                                   100, 0, 1);
  haxis_raa_paper->GetYaxis()->SetTitleOffset(2.5);
  haxis_raa_paper->GetYaxis()->CenterTitle();
  haxis_raa_paper->GetYaxis()->SetTitleFont(63);
  haxis_raa_paper->GetYaxis()->SetLabelFont(63);
  haxis_raa_paper->GetYaxis()->SetLabelSize(fontSizeLabel);
  haxis_raa_paper->GetYaxis()->SetTitleSize(fontSizeTitle);
  haxis_raa_paper->GetYaxis()->SetNdivisions(6, 3, 0);
  // haxis_raa_paper->GetXaxis()->SetTitleOffset(2.5);
  haxis_raa_paper->GetXaxis()->SetTitleOffset(3.5);
  haxis_raa_paper->GetXaxis()->SetTitleFont(63);
  haxis_raa_paper->GetXaxis()->SetLabelFont(63);
  haxis_raa_paper->GetXaxis()->SetLabelSize(fontSizeLabel);
  haxis_raa_paper->GetXaxis()->SetTitleSize(fontSizeTitle);
  haxis_raa_paper->GetXaxis()->SetNdivisions(5, 4, 0);
  haxis_raa_paper->SetMinimum(0.01);
  haxis_raa_paper->SetMaximum(1.19);


  TH1F* haxis_jda_paper = new TH1F("haxis_jda_paper",
                                   ";x_{p};J_{d+Au}",
                                   100, 0, 1);
  haxis_jda_paper->GetYaxis()->CenterTitle();
  haxis_jda_paper->GetYaxis()->SetTitleOffset(2.5);
  haxis_jda_paper->GetYaxis()->SetTitleFont(63);
  haxis_jda_paper->GetYaxis()->SetLabelFont(63);
  haxis_jda_paper->GetYaxis()->SetLabelSize(fontSizeLabel);
  haxis_jda_paper->GetYaxis()->SetTitleSize(fontSizeTitle);
  haxis_jda_paper->GetYaxis()->SetNdivisions(4, 5, 0);
  haxis_jda_paper->GetXaxis()->SetTitleOffset(3.5);
  haxis_jda_paper->GetXaxis()->SetTitleFont(63);
  haxis_jda_paper->GetXaxis()->SetLabelFont(63);
  haxis_jda_paper->GetXaxis()->SetLabelSize(fontSizeLabel);
  haxis_jda_paper->GetXaxis()->SetTitleSize(fontSizeTitle);
  haxis_jda_paper->GetXaxis()->SetNdivisions(5, 4, 0);
  haxis_jda_paper->SetMinimum(0.01);
  haxis_jda_paper->SetMaximum(1.99);


  TH1F* haxis_Q_paper = new TH1F("haxis_Q_paper",
                                 ";x_{p};<dN^{hard}/dQ(x_{p})>/<dN^{hard}/dQ(x_{p}=0)>",
                                 100, 0, 1);
  haxis_Q_paper->GetYaxis()->CenterTitle();
  haxis_Q_paper->GetYaxis()->SetTitleOffset(1.5);
  haxis_Q_paper->GetYaxis()->SetTitleFont(63);
  haxis_Q_paper->GetYaxis()->SetLabelFont(63);
  haxis_Q_paper->GetYaxis()->SetLabelSize(fontSizeLabel);
  haxis_Q_paper->GetYaxis()->SetTitleSize(fontSizeTitle);
  haxis_Q_paper->GetYaxis()->SetNdivisions(10, 2, 0);
  haxis_Q_paper->GetXaxis()->SetTitleOffset(1.5);
  haxis_Q_paper->GetXaxis()->SetTitleFont(63);
  haxis_Q_paper->GetXaxis()->SetLabelFont(63);
  haxis_Q_paper->GetXaxis()->SetLabelSize(fontSizeLabel);
  haxis_Q_paper->GetXaxis()->SetTitleSize(fontSizeTitle);
  haxis_Q_paper->GetXaxis()->SetNdivisions(10, 2, 0);
  haxis_Q_paper->SetMinimum(0.41);
  haxis_Q_paper->SetMaximum(1.09);

  TLatex lpaper;
  lpaper.SetNDC();
  lpaper.SetTextAlign(22);
  lpaper.SetTextFont(63);
  lpaper.SetTextSize(fontSizeTitle);
  // lpaper.SetTextSize(0.08);

  double topMargin = 0.05;
  // double bottomMargin = 0.15;
  double bottomMargin = 0.2;
  double yl, yh, fracArea;

  const char *plabel[4] = {"(a)", "(b)", "(c)", "(d)"};


  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);
  label.SetTextFont(63);
  label.SetTextSize(fontSizeLegend);

  TLine l1;
  l1.SetLineStyle(2);

  //=====================================================//
  // FIGURES FOR PAPER
  //=====================================================//

  //-------------------------------------//
  //-- Rcp vs xp for 3 centrality bins --//
  //-------------------------------------//
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


    l1.DrawLine(0, 1, 0.6, 1);
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



  //-------------------------------------//
  //-- RAA vs xp for 4 centrality bins --//
  //-------------------------------------//
  TCanvas *craa_paper = new TCanvas("craa_paper", "RAA", 400, 1000);
  craa_paper->SetTopMargin(0.0);
  craa_paper->SetRightMargin(0.0);
  craa_paper->SetBottomMargin(0.0);
  craa_paper->SetLeftMargin(0.0);

  TPad *praa_paper[NCENT];
  fracArea = 1. / (float)(NCENT - 2 + 1. / (1 - topMargin) + 1. / (1 - bottomMargin));
  // double raal[NCENT] = {0.01, 0.41, 0.61, 0.76};
  // double raah[NCENT] = {1.19, 1.59, 1.59, 2.49};
  double raal[NCENT] = {0.31, 0.51, 0.71, 0.89};
  double raah[NCENT] = {1.19, 1.39, 1.59, 2.19};
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

    l1.DrawLine(0, 1, 0.6, 1);
    if (icent == 0)
      legraa_paper->Draw("same");
    if (icent == NCENT - 1)
      lpaper.DrawLatex(0.9, 0.26, plabel[icent]);
    else
      lpaper.DrawLatex(0.9, 0.07, plabel[icent]);
    if (icent == 0)
      lpaper.DrawLatex(0.3, 0.85, Form("%.0f-%.0f%%",
                                       centl[1][icent], centl[1][icent + 1]));
    else
      lpaper.DrawLatex(0.3, 0.90, Form("%.0f-%.0f%%",
                                       centl[1][icent], centl[1][icent + 1]));
    if (icent == 1)
      lpaper.DrawLatex(0.35, 0.07, "#sqrt{s_{NN}}=200 GeV");
    craa_paper->cd();

  }


  //------------------------------------//
  //-- JdA vs xp in 4 centrality bins --//
  //------------------------------------//
  TCanvas *cjda_paper = new TCanvas("cjda_paper", "JdAu", 400, 1000);
  cjda_paper->SetTopMargin(0.0);
  cjda_paper->SetRightMargin(0.0);
  cjda_paper->SetBottomMargin(0.0);
  cjda_paper->SetLeftMargin(0.0);

  TPad *pjda_paper[NCENT];
  fracArea = 1. / (float)(NCENT - 2 + 1. / (1 - topMargin) + 1. / (1 - bottomMargin));
  double jdal[NCENT] = {0.01, 0.01, 0.01, 0.01};
  double jdah[NCENT] = {1.99, 1.99, 1.99, 1.99};
  // double jdah[NCENT] = {1.19, 1.19, 1.99, 1.99};
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

    pjda_paper[icent] = new TPad(Form("pjda_paper_%i", icent), "", 0, yl, 1, yh);
    pjda_paper[icent]->SetLeftMargin(0.15);
    pjda_paper[icent]->SetRightMargin(0.05);
    if (icent == 0)
      pjda_paper[icent]->SetTopMargin(topMargin);
    else
      pjda_paper[icent]->SetTopMargin(0.0);
    if (icent == NCENT - 1)
      pjda_paper[icent]->SetBottomMargin(bottomMargin);
    else
      pjda_paper[icent]->SetBottomMargin(0.0);
    pjda_paper[icent]->SetTicks(1, 1);
    pjda_paper[icent]->Draw();

    pjda_paper[icent]->cd();
    haxis_jda_paper->GetYaxis()->SetRangeUser(jdal[icent], jdah[icent]);
    haxis_jda_paper->DrawCopy();

    for (int ix = 0; ix < NXBINS; ix++)
    {
      for (int j = 0; j < nxjda[ix]; j++)
        bJdA[ix][icent][j]->Draw("same");
      gJdA[ix][icent]->Draw("P");
    }

    graa[1][icent]->Draw("C");
    graa_linMod[1][icent]->Draw("C");

    l1.DrawLine(0, 1, 1.0, 1);
    if (icent == 0)
      legjda_model_paper->Draw("same");
    if (icent == 1)
    {
      label.DrawLatex(0.65, 0.89, "PHENIX");
      legjda_data_paper->Draw("same");
    }
    // if (icent == 1)
    //   lpaper.DrawLatex(0.7, 0.90, "d+Au #sqrt{s_{NN}}=200 GeV");
    if (icent == NCENT - 1)
      lpaper.DrawLatex(0.2, 0.26, plabel[icent]);
    else
      lpaper.DrawLatex(0.2, 0.07, plabel[icent]);
    if (icent == 0)
      lpaper.DrawLatex(0.25, 0.85, Form("%.0f-%.0f%%",
                                        centl[1][icent], centl[1][icent + 1]));
    else
      lpaper.DrawLatex(0.25, 0.90, Form("%.0f-%.0f%%",
                                        centl[1][icent], centl[1][icent + 1]));
    cjda_paper->cd();

  }





  // //------------------------------------//
  // //-- JdA vs xp in 2 centrality bins --//
  // //------------------------------------//
  // TCanvas *cjda_paper = new TCanvas("cjda_paper", "jda", 500, 600);
  // cjda_paper->SetTopMargin(0);
  // cjda_paper->SetRightMargin(0);
  // cjda_paper->SetBottomMargin(0);
  // cjda_paper->SetLeftMargin(0);

  // TPad *pjda_paper[3];
  // fracArea = 1. / (float)(2 - 2 + 1. / (1 - topMargin) + 1. / (1 - bottomMargin));

  // //-- 0-20%
  // cjda_paper->cd();
  // yl = 1.0 - fracArea * (1. / (1. - topMargin));
  // yh = 1;
  // cout << " jdA 0 - yl:" << yl << " yh:" << yh <<  "fracArea:" << fracArea << endl;
  // pjda_paper[0] = new TPad("pjda_paper_0", "", 0, yl, 1, yh);
  // pjda_paper[0]->SetLeftMargin(0.15);
  // pjda_paper[0]->SetRightMargin(0.05);
  // pjda_paper[0]->SetTopMargin(topMargin);
  // pjda_paper[0]->SetBottomMargin(0);
  // pjda_paper[0]->SetTicks(1, 1);
  // pjda_paper[0]->Draw();
  // pjda_paper[0]->cd();
  // haxis_jda_paper->GetYaxis()->SetRangeUser(0.01, 1.19);
  // haxis_jda_paper->DrawCopy();

  // for (int ix = 0; ix < NXBINS; ix++)
  // {
  //   for (int j = 0; j < nxjda[ix]; j++)
  //     bJdA[ix][0][j]->Draw("same");
  //   gJdA[ix][0]->Draw("P");
  // }

  // graa[1][0]->Draw("C");
  // graa_linMod[1][0]->Draw("C");

  // l1.DrawLine(0, 1, 1.0, 1);
  // lpaper.DrawLatex(0.23, 0.06, "0-20%");
  // lpaper.DrawLatex(0.18, 0.95 - topMargin, "(a)");
  // legjda_paper->Draw("same");

  // //-- 60-88%
  // cjda_paper->cd();
  // // yl = 1.0 - (1. / (1. - topMargin) + 1) * fracArea;
  // yl = 0.0;
  // yh = 1.0 - (1. / (1. - topMargin) + 1 - 1) * fracArea;
  // cout << " jdA 1 - yl:" << yl << " yh:" << yh <<  "fracArea:" << fracArea << endl;
  // pjda_paper[1] = new TPad("pjda_paper_1", "", 0, yl, 1, yh);
  // pjda_paper[1]->SetLeftMargin(0.15);
  // pjda_paper[1]->SetRightMargin(0.05);
  // pjda_paper[1]->SetTopMargin(0);
  // pjda_paper[1]->SetBottomMargin(bottomMargin);
  // pjda_paper[1]->SetTicks(1, 1);
  // pjda_paper[1]->Draw();
  // pjda_paper[1]->cd();
  // haxis_jda_paper->GetYaxis()->SetRangeUser(0.01, 1.99);
  // haxis_jda_paper->DrawCopy();

  // for (int ix = 0; ix < NXBINS; ix++)
  // {
  //   for (int j = 0; j < nxjda[ix]; j++)
  //     bJdA[ix][NCENT - 1][j]->Draw("same");
  //   gJdA[ix][NCENT - 1]->Draw("P");
  // }

  // graa[1][NCENT - 1]->Draw("C");
  // graa_linMod[1][NCENT - 1]->Draw("C");

  // l1.DrawLine(0, 1, 1.0, 1);
  // lpaper.DrawLatex(0.23, 0.06 + bottomMargin, "60-88%");
  // lpaper.DrawLatex(0.18, 0.95, "(b)");



  //-----------------//
  //-- sigma vs xp --//
  //-----------------//
  TCanvas *csig_paper = new TCanvas("csig_paper", "sig mod", 500, 500);
  csig_paper->SetTopMargin(0.05);
  csig_paper->SetRightMargin(0.05);
  csig_paper->SetBottomMargin(0.15);
  csig_paper->SetLeftMargin(0.15);
  csig_paper->SetTicks(1, 1);

  csig_paper->cd(1);
  haxis_sig->Draw();
  fsig->Draw("same");


  //---------------//
  //-- <Q> vs xp --//
  //---------------//
  TCanvas *cq_paper = new TCanvas("cq_paper", "<Q>", 500, 500);
  cq_paper->SetTopMargin(0.02);
  cq_paper->SetRightMargin(0.02);
  cq_paper->SetBottomMargin(0.15);
  cq_paper->SetLeftMargin(0.15);
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