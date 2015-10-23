////////////////////////////////////////////////////////////////////////////////
//
// Calculate RCP for small systems - (p, d, He3)+Au - from Glauber model
// output.
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// Oct 15 2015
//
////////////////////////////////////////////////////////////////////////////////
//
// Useful References:
//  -- pAu --
// https://www.phenix.bnl.gov/WWW/p/draft/nagle/PHENIX/nagle_run15_pau_update_05282015.pdf
//  -- dAu --
// http://www.phenix.bnl.gov/phenix/WWW/p/info/an/900/Run8_dAu_200GeV_Centrality_Categorization.pdf
// http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1087/Run8_dAu_200GeV_Centrality_Addendum-01.pdf
//  -- He3Au --
// http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1207/Run14_3HeAu_200GeV_Centrality_Categorization.pdf
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
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

#include <iostream>

using namespace std;

Double_t NBD(Double_t *x, Double_t *par)
{
  // Negative Binomial Distribution
  double P = 0;
  double n = x[0];
  double mu = par[0];
  double k = par[1];

  if (n + k > 100.0) {
    // log method for handling large numbers
    P  = TMath::LnGamma(n + k) - TMath::LnGamma(n + 1.) - TMath::LnGamma(k);
    P  += n * TMath::Log(mu / k) - (n + k) * TMath::Log(1.0 + mu / k);
    P = TMath::Exp(P);
  } else {
    P = TMath::Gamma(n + k);
    P /= TMath::Gamma(n + 1) * TMath::Gamma(k);
    P *= TMath::Exp(n * TMath::Log(mu / k) - (n + k) * TMath::Log(1 + mu / k));
  }
  return P;
}

void calculate_RCP()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //=====================================================//
  // SET RUNNING CONDITIONS
  //=====================================================//

  const int NX = 6;             // Number of x values
  double x[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};         // x values

  const int NSYSTEMS = 3;
  const char *collSystem[] = {"pAu", "dAu", "He3Au"};

  bool saveRcp = true;
  const char *outFile = "Rcp_systems.root";


  // Negative Binomial Distribution parameters for each system
  double NBD_mu[] =
  {
    3.14, // pAu
    3.04, // dAu
    2.91, // He3Au
  };
  double NBD_k[] =
  {
    0.47, // pAu
    0.46, // dAu
    0.55, // He3Au
  };

  // Files for each x value for each system
  const char *xFiles[NSYSTEMS][NX] =
  {
    { // pAu Files
      "glauber_pau_snn42_x0_ntuple_100k.root",
      "glauber_pau_snn42_x01_ntuple_100k.root",
      "glauber_pau_snn42_x02_ntuple_100k.root",
      "glauber_pau_snn42_x03_ntuple_100k.root",
      "glauber_pau_snn42_x04_ntuple_100k.root",
      "glauber_pau_snn42_x05_ntuple_100k.root",
    },
    { // dAu Files
      "glauber_dau_snn42_x0_ntuple_100k.root",
      "glauber_dau_snn42_x01_ntuple_100k.root",
      "glauber_dau_snn42_x02_ntuple_100k.root",
      "glauber_dau_snn42_x03_ntuple_100k.root",
      "glauber_dau_snn42_x04_ntuple_100k.root",
      "glauber_dau_snn42_x05_ntuple_100k.root",
    },
    { // He3Au Files
      "glauber_he3au_snn42_x0_ntuple_100k.root",
      "glauber_he3au_snn42_x01_ntuple_100k.root",
      "glauber_he3au_snn42_x02_ntuple_100k.root",
      "glauber_he3au_snn42_x03_ntuple_100k.root",
      "glauber_he3au_snn42_x04_ntuple_100k.root",
      "glauber_he3au_snn42_x05_ntuple_100k.root",
    },

  };

  // Ntuple name stored in files for each system
  const char *ntpName[NSYSTEMS] =
  {
    "nt_p_Au", // pAu
    "nt_dh_Au", // dAu
    "nt_He3_Au", // He3Au
  };

  // Centrality bins
  // double dcent_cent   = 0.05; // 0-5%
  // double dcent_periph = 0.30; // 70-100%
  double dcent_cent   = 0.20; // 0-20%
  double dcent_periph = 0.40; // 60-100%


  // const int NCENT = 4;
  // double centl[] = { 0, 20, 40,  60};
  // double centh[] = {20, 40, 60, 100};
  const int NCENT = 9;
  double centl[] = { 0,  5, 10, 20, 30, 40, 50, 60,  70};
  double centh[] = { 5, 10, 20, 30, 40, 50, 60, 70, 100};

  // line colors
  int colors[NSYSTEMS] = { kBlue, kRed, kGreen + 2};

  //=====================================================//
  // DECLARE VARIABLES
  //=====================================================//

  // collision system variables
  double BBC_central = 0;
  double BBC_periph = 0;

  TF1 *fNBD = new TF1("fNBD", NBD, 0, 200, 2);

  // For running over ntuples
  TFile *fin;
  TNtuple *ntp;
  Float_t Ncoll;

  TH1D *hNcoll[NSYSTEMS][NX];
  TH1D *hBBCs[NSYSTEMS][NX];
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    for (int i = 0; i < NX; i++)
    {
      hBBCs[isys][i] = new TH1D(Form("hBBCs_%i_%i", isys, i), ";BBCs charge", 800, 0, 200);
      hBBCs[isys][i]->SetLineColor(i + 1);
    }
  }

  TH1D* htmp;


  // For counting
  double Nevent_central[NSYSTEMS][NX] = {0};
  double Nevent_periph[NSYSTEMS][NX] = {0};

  // Calculating RCP
  TGraph *grcp[NSYSTEMS];


  double Nevent_cent[NSYSTEMS][NX][NCENT] = {0};
  double BBC_cent[NCENT] = {0};
  TGraph *grcp_cent[NSYSTEMS][NX];

  //=====================================================//
  // GET NTUPLE(S) FROM FILE AND CALCULATE YIELD
  //=====================================================//
  cout << endl;
  cout << "--> Reading Ntuples from files ..." << endl;

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cout << "------------------------------" << endl;
    cout << "          " << collSystem[isys] << endl;
    cout << "------------------------------" << endl;

    // print values
    cout << "  NBD : mu=" << NBD_mu[isys] << " k=" << NBD_k[isys] << endl;
    cout << "  Ntuple name: " << ntpName[isys] << endl;

    fNBD->SetParameters(NBD_mu[isys], NBD_k[isys]);
    // feff->SetParameters(trig_eff_params[0], trig_eff_params[1]);


    for (int ix = 0; ix < NX; ix++)
    {
      fin = TFile::Open(xFiles[isys][ix]);
      if (!fin)
      {
        cout << "ERROR!! Unable to open " << xFiles[isys][ix] << endl;
        return;
      }
      else
        cout << "----> Reading " << xFiles[isys][ix] << " for x=" << x[ix] << endl;

      ntp = (TNtuple*) fin->Get(ntpName[isys]);
      if (!ntp)
      {
        cout << "ERROR!! Unable to find " << ntpName[isys] << " in " << xFiles[isys][ix] << endl;
        return;
      }

      ntp->Draw("Ncoll>>htmp(50, 0.5, 50.5)", "", "goff");

      hNcoll[isys][ix] = (TH1D*) gDirectory->FindObject("htmp");
      hNcoll[isys][ix]->SetDirectory(0);
      hNcoll[isys][ix]->SetName(Form("hNcoll_%i_%i", isys, ix));
      hNcoll[isys][ix]->SetTitle(";N_{coll}");
      hNcoll[isys][ix]->SetLineColor(ix + 1);

      cout << "   <Ncoll>: " << hNcoll[isys][ix]->GetMean() << endl;

      //Convolve with NBD to get BBC charge
      for (int i = 1; i <= hNcoll[isys][ix]->GetNbinsX(); i++)
      {
        double ncoll = hNcoll[isys][ix]->GetBinCenter(i);
        double w = hNcoll[isys][ix]->GetBinContent(i);

        fNBD->SetParameters(ncoll * NBD_mu[isys], ncoll * NBD_k[isys]);

        for (int j = 1; j <= hBBCs[isys][ix]->GetNbinsX(); j++)
        {
          double c = hBBCs[isys][ix]->GetBinCenter(j);
          double bc = hBBCs[isys][ix]->GetBinContent(j);

          double v = w * fNBD->Eval(c);

          hBBCs[isys][ix]->SetBinContent(j, v + bc);

          // cout << " ncoll:" << ncoll << " BBCc:" << c
          //      << " w:" << w << " NBD:" << fNBD->Eval(c)
          //      << " v:" << v << " bc:" << bc
          //      << " new bc:" << v + bc
          //      << endl;
        }
      }

      // // apply the trigger efficiency
      // for (int j = 1; j <= hBBCs[ix]->GetNbinsX(); j++)
      // {
      //   double c = hBBCs[ix]->GetBinCenter(j);
      //   double bc = hBBCs[ix]->GetBinContent(j);
      //   double eff = feff->Eval(c);
      //   hBBCs[ix]->SetBinContent(c, bc * eff);
      // }

      // check the quantiles
      double xq[2];
      double yq[2] = {0};
      xq[0] = dcent_periph;
      xq[1] = 1. - dcent_cent;
      hBBCs[isys][ix]->GetQuantiles(2, yq, xq);
      cout << " quantiles: "
           << xq[0] << ": " << yq[0] << "  "
           << xq[1] << ": " << yq[1] << "  "
           << endl;
      if (ix == 0)
      {
        BBC_periph = yq[0];
        BBC_central = yq[1];
      }

      // Count the number of events
      Nevent_central[isys][ix] = hBBCs[isys][ix]->Integral(hBBCs[isys][ix]->FindBin(BBC_central),
                                 hBBCs[isys][ix]->GetNbinsX());
      Nevent_periph[isys][ix] = hBBCs[isys][ix]->Integral(1,
                                hBBCs[isys][ix]->FindBin(BBC_periph));


      cout << "   N central    = " << Nevent_central[isys][ix] << endl;
      cout << "   N peripheral = " << Nevent_periph[isys][ix] << endl;



      // get the number of events for each centrality bin
      if (ix == 0)
      {
        double xq_c[NCENT - 1];
        double yq_c[NCENT - 1] = {0};
        for (int icent = 0; icent < NCENT - 1; icent++)
          xq_c[icent] = (100. - centh[icent]) / 100.;

        hBBCs[isys][ix]->GetQuantiles(NCENT - 1, yq_c, xq_c);

        for (int icent = 0; icent < NCENT - 1; icent++)
        {
          BBC_cent[icent] = yq_c[icent];
          cout << "   " << icent << " "
               << xq_c[icent] << " "
               << yq_c[icent] << endl;
        }
      }

      for (int icent = 0; icent < NCENT; icent++)
      {
        if (icent == 0)
          Nevent_cent[isys][ix][icent] = hBBCs[isys][ix]->Integral(hBBCs[isys][ix]->FindBin(BBC_cent[icent]),
                                         hBBCs[isys][ix]->GetNbinsX());
        else if (icent == NCENT - 1)
          Nevent_cent[isys][ix][icent] = hBBCs[isys][ix]->Integral(1,
                                         hBBCs[isys][ix]->FindBin(BBC_cent[icent - 1]));
        else
          Nevent_cent[isys][ix][icent] = hBBCs[isys][ix]->Integral(hBBCs[isys][ix]->FindBin(BBC_cent[icent]),
                                         hBBCs[isys][ix]->FindBin(BBC_cent[icent - 1]));

        cout << "   " << isys << " " << ix << " " << icent
             << " " << Nevent_cent[isys][ix][icent] << endl;
      }

      // delete htmp;
      delete ntp;
      fin->Close();
      delete fin;

    }

  }
  //=====================================================//
  // CALCULATE RCP
  //=====================================================//
  cout << endl;
  cout << "--> Calculating RCP ..." << endl;

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    grcp[isys] = new TGraph();
    grcp[isys]->SetName(Form("grcp_%i", isys));
    grcp[isys]->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    grcp[isys]->SetLineStyle(1);
    grcp[isys]->SetLineColor(colors[isys]);
    // grcp[isys]->SetPoint(0, 0, 1);

    for (int ix = 0; ix < NX; ix++)
    {
      // x = pT / sqrt(s)
      double pT = 0.5 * x[ix] * 200;

      double rcp = Nevent_central[isys][ix] / dcent_cent;
      rcp /= Nevent_periph[isys][ix] / dcent_periph;

      cout << " " << x[ix] << " " << pT << " " << rcp << endl;

      grcp[isys]->SetPoint(ix + 1, x[ix], rcp);
    }
  }

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    for (int ix = 0; ix < NX; ix++)
    {
      grcp_cent[isys][ix] = new TGraph();
      grcp_cent[isys][ix]->SetName(Form("grcp_cent_%i_%i", isys, ix));
      grcp_cent[isys][ix]->SetTitle("; centrality ; R ");
      grcp_cent[isys][ix]->SetLineStyle(1);
      grcp_cent[isys][ix]->SetLineColor(1 + ix);
      grcp_cent[isys][ix]->SetLineWidth(2);

      for (int icent = 0; icent < NCENT; icent++)
      {
        double rcp = Nevent_cent[isys][ix][icent] / (centh[icent] - centl[icent]);
        rcp /= Nevent_cent[isys][ix][4] / (centh[4] - centl[4]);


        // cout << "  " << isys << " " << x[ix]
        //      << " " << icent
        //      << " " << Nevent_cent[isys][ix][icent] << " " << (centh[icent] - centl[icent])
        //      << " " << Nevent_cent[isys][ix][4] << " " << (centh[4] - centl[4])
        //      << " " << rcp
        //      << endl;

        grcp_cent[isys][ix]->SetPoint(icent, icent, rcp);
      }
    }
  }

  //=====================================================//
  // Run 8 d+Au Jet Rcp
  //=====================================================//
  cout << endl;
  cout << "--> Run 8 d+Au Jet Rcp";

  // Run 8 d+Au Jet RCP
  // http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1222/dvp-dAuJet-FinalResults-Note-03.pdf
  // 0-20% / 60-88%
  const int NJET = 8;
  double pTl_jet[] =   {12.1, 14.5, 17.3, 20.7, 24.7, 29.4, 35.1, 41.9};
  double pTh_jet[] =   {14.5, 17.3, 20.7, 24.7, 29.4, 35.1, 41.9, 50.0};
  double Rcp_jet[] =   {0.73, 0.71, 0.66, 0.61, 0.57, 0.54, 0.52, 0.50};
  double Rcp_jet_A[] = {0.01, 0.01, 0.02, 0.02, 0.02, 0.04, 0.05, 0.06};
  double Rcp_jet_B[] = {0.06, 0.03, 0.02, 0.02, 0.03, 0.05, 0.08, 0.09};

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

  TGraphErrors *gRCP_jet = new TGraphErrors(NJET,
      x_jet, Rcp_jet,
      xe_jet, Rcp_jet_A);
  gRCP_jet->SetMarkerStyle(20);
  gRCP_jet->SetMarkerColor(kBlack);


  // systematic uncertainties
  TBox *bRcp_jet[NJET];
  for (int i = 0; i < NJET; i++)
  {
    double x1 = xl_jet[i];
    double x2 = xh_jet[i];

    double y1 = Rcp_jet[i] - Rcp_jet_B[i];
    double y2 = Rcp_jet[i] + Rcp_jet_B[i];

    bRcp_jet[i] = new TBox(x1, y1, x2, y2);
    bRcp_jet[i]->SetFillColor(kBlack);
    bRcp_jet[i]->SetFillStyle(3003);
  }

  //=====================================================//
  // PLOT OBJECTS
  //=====================================================//
  cout << endl;
  cout << "--> Plotting ..." << endl;

  TLegend *legNcoll[NSYSTEMS];
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    if (isys == 0)
      legNcoll[isys] = new TLegend(0.5, 0.5, 0.9, 0.9);
    else
      legNcoll[isys] = new TLegend(0.15, 0.15, 0.5, 0.6);
    legNcoll[isys]->SetFillStyle(0);
    legNcoll[isys]->SetBorderSize(0);
    for (int ix = 0; ix < NX; ix++)
    {
      legNcoll[isys]->AddEntry(hNcoll[isys][ix],
                               Form("x=%.2f <N_{coll}>=%.2f", x[ix], hNcoll[isys][ix]->GetMean()),
                               "L");
    }
  }

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);

  TH1F* haxis_rcp = new TH1F("haxis_rcp",
                             Form(";x (x_{jet}=2 * p_{T}^{jet} / #sqrt{s_{NN}});R_{CP}(0-%.0f%% / %.0f-100%%)",
                                  dcent_cent * 100, (1 - dcent_periph) * 100),
                             100, 0, 1);
  haxis_rcp->SetMinimum(0);
  haxis_rcp->SetMaximum(1.1);

  TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.83);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  for (int isys = 0; isys < NSYSTEMS; isys++)
    leg->AddEntry(grcp[isys], collSystem[isys], "L");

  TLine l1;
  l1.SetLineStyle(2);

  TH1F* haxis_rcpcent = new TH1F("haxis_rcpcent",
                                 ";centrality; R",
                                 NCENT, 0, NCENT - 1);
  haxis_rcpcent->SetMinimum(0);
  haxis_rcpcent->SetMaximum(2.0);
  // haxis_rcpcent->GetXaxis()->SetNdivisions(NCENT);
  for (int ibin = 1; ibin <= haxis_rcpcent->GetNbinsX(); ibin++)
  {
    haxis_rcpcent->GetXaxis()->SetBinLabel(ibin,
                                           Form("%.0f - %.0f",
                                               centl[ibin - 1],
                                               centh[ibin - 1]));
  }
  haxis_rcpcent->LabelsOption("d");

  TLegend *legx = new TLegend(0.6, 0.15, 0.9, 0.5);
  legx->SetFillStyle(0);
  legx->SetBorderSize(0);
  for (int ix = 0; ix < NX; ix++)
  {
    legx->AddEntry(grcp_cent[0][ix], Form("x = %.2f", x[ix]), "L");
  }

  //=====================================================//
  // PLOT
  //=====================================================//

  TCanvas *crcp = new TCanvas("crcp", "RCP", 800, 800);
  crcp->cd(1);
  haxis_rcp->GetXaxis()->SetRangeUser(0, 0.5);
  haxis_rcp->Draw();
  for (int isys = 0; isys < NSYSTEMS; isys++)
    grcp[isys]->Draw("L");

  gRCP_jet->Draw("P");
  for (int i = 0; i < NJET; i++)
    bRcp_jet[i]->Draw();

  l1.DrawLine(0, 1, 0.5, 1);
  leg->Draw("same");


  TCanvas *ctest = new TCanvas("ctest", "test", 2000, 1000);
  ctest->Divide(2, NSYSTEMS);

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    ctest->cd(2 * isys + 1);
    gPad->SetLogy();
    hNcoll[isys][0]->Draw();
    for (int ix = 1; ix < NX; ix++)
      hNcoll[isys][ix]->Draw("same");
    legNcoll[isys]->Draw("same");
    label.DrawLatex(0.5, 0.95, collSystem[isys]);

    ctest->cd(2 * isys + 2);
    gPad->SetLogy();
    hBBCs[isys][0]->Draw();
    for (int ix = 1; ix < NX; ix++)
      hBBCs[isys][ix]->Draw("same");
    label.DrawLatex(0.5, 0.95, collSystem[isys]);
  }

  TCanvas *crcpcent = new TCanvas("crcpcent", "rcp cent", 1200, 600);
  crcpcent->Divide(NSYSTEMS, 1);


  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    crcpcent->GetPad(isys + 1)->SetTicks(1, 1);
    crcpcent->cd(isys + 1);
    haxis_rcpcent->Draw();
    for (int ix = 0; ix < NX; ix++)
      grcp_cent[isys][ix]->Draw("L");

    l1.DrawLine(0, 1, NCENT-1, 1);
    label.DrawLatex(0.5, 0.95, collSystem[isys]);
    if (isys == NSYSTEMS - 1) legx->Draw("same");
  }

  //=====================================================//
  // SAVE
  //=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving RCP to " << outFile << endl;

    ctest->Print("NcollBBCcharge.pdf");
    crcp->Print("Rcp_systems.pdf");
    crcpcent->Print("R_systems.pdf");

    // TFile *fout = new TFile(outFile, "UPDATE");

    // grcp->Write(Form("Rcp_%s", collSystem));

    // fout->Close();
  }
}