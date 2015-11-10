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

#include <iostream>
#include <iomanip>

using namespace std;

double evalNBD(double n, double mu, double k)
{
  // Negative Binomial Distribution
  double P = 0;

  if (mu < 0 || k < 0) {
    P = 0;
  } else if (n + k > 100.0) {
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
  double x[] = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5};

  double beta = 1.5;

  const int NSYSTEMS = 3;
  const char *collSystem[] = {"pAu", "dAu", "3HeAu"};

  bool saveRcp = true;
  const char *outFile = "Rcp_systems.root";

  // Files for each x value for each system
  const char *xFiles[NSYSTEMS][NX] =
  {
    { // pAu Files
      "rootfiles/glauber_pau_snn42_x001_ntuple_100k.root",
      "rootfiles/glauber_pau_snn42_x01_ntuple_100k.root",
      "rootfiles/glauber_pau_snn42_x02_ntuple_100k.root",
      "rootfiles/glauber_pau_snn42_x03_ntuple_100k.root",
      "rootfiles/glauber_pau_snn42_x04_ntuple_100k.root",
      "rootfiles/glauber_pau_snn42_x05_ntuple_100k.root",
    },
    { // dAu Files
      "rootfiles/glauber_dau_snn42_x001_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_x01_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_x02_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_x03_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_x04_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_x05_ntuple_100k.root",
    },
    { // He3Au Files
      "rootfiles/glauber_he3au_snn42_x001_ntuple_100k.root",
      "rootfiles/glauber_he3au_snn42_x01_ntuple_100k.root",
      "rootfiles/glauber_he3au_snn42_x02_ntuple_100k.root",
      "rootfiles/glauber_he3au_snn42_x03_ntuple_100k.root",
      "rootfiles/glauber_he3au_snn42_x04_ntuple_100k.root",
      "rootfiles/glauber_he3au_snn42_x05_ntuple_100k.root",
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
  const int NCENT = 4;
  double centl[NSYSTEMS][NCENT + 1] =
  {
    {0, 20, 40, 60, 84}, // pAu
    {0, 20, 40, 60, 88}, // dAu
    {0, 20, 40, 60, 88}, // He3Au
  };

  // PHENIX values (last entry is 0-100%)
  double Ncoll_PHENIX[NSYSTEMS][NCENT + 1] =
  {
    { 8.20, 6.06, 4.43, 2.62, 4.67}, // pAu
    {15.1, 10.2, 6.6, 3.2, 7.6}, // dAu
    {22.37, 14.71, 8.28, 3.38, 10.45}, //HeAu
  };

  // Negative Binomial Distribution parameters for each system
  double NBD_par[NSYSTEMS][2] =
  {
    {3.14, 0.47}, //pAu {mu, k}
    {3.04, 0.46}, //dAu {mu, k}
    {2.91, 0.55}, //3HeAu {mu, k}
  };

  // MB trigger efficiency function parameters
  double eff_par[NSYSTEMS][2] =
  {
    {1.07552e+00, 6.02328e-01}, // pAu
    {0.897, 0.612}, // dAu (from D. McGlinchey thesis)
    {1.22134e+00, 5.10114e-01}, // He3Au
  };

  // line colors
  int colors[NSYSTEMS] = { kRed, kBlue, kGreen + 2};
  // int lstyle[NSYSTEMS] = {    3,     5,          7};
  int lstyle[NSYSTEMS] = {1, 1, 1};

  // fill colors
  int fcolor_cent[NCENT] = {kBlue, kRed, kGreen + 2, kYellow + 2};

  //=====================================================//
  // DECLARE VARIABLES
  //=====================================================//

  // For running over ntuples
  TFile *fin;
  TNtuple *ntp;

  TH2D *hNcoll_NcollMod[NSYSTEMS][NX];
  TH2D *hNcoll_NcollA[NSYSTEMS][NX];
  TH2D *hNcollMod_NcollModA[NSYSTEMS][NX];

  TH2D *hNcoll_NcollModA[NSYSTEMS][NX];
  TH2D *hNcollMod_NcollA[NSYSTEMS][NX];



  // Calculated histograms
  TH1D *hBBCsc[NSYSTEMS][NX];

  TH2D *hNcoll_BBCsc[NSYSTEMS][NX];
  TH2D *hNcollMod_BBCscMod[NSYSTEMS][NX];

  TH1D* hBBCscNcollA[NSYSTEMS][NX];
  TH1D* hBBCscNcollModA[NSYSTEMS][NX];
  TH1D* hBBCscModNcollA[NSYSTEMS][NX];
  TH1D* hBBCscModNcollModA[NSYSTEMS][NX];

  for (int i = 0; i < NSYSTEMS; i++)
  {
    for (int j = 0; j < NX; j++)
    {
      // BBCs charge (unmodified)
      hBBCsc[i][j] = new TH1D(Form("hBBCsc_%i_%i", i, j),
                              ";BBCs charge",
                              1596, 1, 400);

      // Ncoll vs BBCs charge
      hNcoll_BBCsc[i][j] = new TH2D(Form("hNcoll_BBCsc_%i_%i", i, j),
                                    ";BBCs charge; N_{coll}",
                                    1596, 1, 400,
                                    101, -0.5, 100.5);

      hNcollMod_BBCscMod[i][j] = new TH2D(Form("hNcollMod_BBCscMod_%i_%i", i, j),
                                          ";BBCs charge Mod; N_{coll}^{mod}",
                                          1596, 1, 400,
                                          101, -0.5, 100.5);


      // Ncoll weighted BBC charge distributions
      hBBCscNcollA[i][j] = new TH1D(Form("hBBCscNcollA_%i_%i", i, j),
                                    ";N_{coll}(A) #time BBCs charge",
                                    1596, 1, 400);

      hBBCscModNcollModA[i][j] = new TH1D(Form("hBBCscModNcollModA_%i_%i", i, j),
                                          ";N_{coll}^{mod}(A) #time BBCs charge Mod",
                                          1596, 1, 400);

      hBBCscNcollModA[i][j] = new TH1D(Form("hBBCscNcollModA_%i_%i", i, j),
                                       ";N_{coll}^{mod}(A) #time BBCs charge",
                                       1596, 1, 400);

      hBBCscModNcollA[i][j] = new TH1D(Form("hBBCscModNcollA_%i_%i", i, j),
                                       ";N_{coll}(A) #time BBCs charge Mod",
                                       1596, 1, 400);


    } // j
  } // i


  TH1D *hBBCsc_cent[NSYSTEMS][NCENT]; // distributions for each centrality bin


  TH1D *hNcoll_MB[NSYSTEMS][NX];
  TH1D *hNcollMod_MB[NSYSTEMS][NX];

  TH1D *hNcoll_cent[NSYSTEMS][NX][NCENT];
  TH1D *hNcollMod_cent[NSYSTEMS][NX][NCENT];


  TH1D* htmp;

  // trigger efficiencies
  TF1 *ftrigeff = new TF1("ftrigeff", "1.0-TMath::Exp(-pow((x/[0]), [1]))", 0.0, 200.0);


  // for Rcp
  double bias_NcollModABBCsc[NSYSTEMS][NX][NCENT];
  double bias_NcollABBCscMod[NSYSTEMS][NX][NCENT];
  double bias_NcollModABBCscMod[NSYSTEMS][NX][NCENT];

  double rcp_NcollModABBCscMod[NSYSTEMS][NCENT - 1][NX];
  TGraph *grcp_NcollModABBCscMod[NSYSTEMS][NCENT - 1];
  TGraph *grcp_NcollModABBCscMod_pipT[NSYSTEMS][NCENT - 1];

  // for mean BBCs charge
  double mean_BBCs[NSYSTEMS][NX];

  TGraph *gmean_BBCs[NSYSTEMS];
  TGraph *gmean_BBCs_pipT[NSYSTEMS];

  // sigma modification function
  TF1* fsig = new TF1("fsig", "TMath::Power((1 + TMath::Exp(-1. *[0] * x)), 2) / 4.0", 0, 1);
  fsig->SetLineColor(kBlue);
  fsig->SetParameter(0, beta);

  double pi0_pT[NX] = {0};
  for (int ix = 0; ix < NX; ix++)
    pi0_pT[ix] = 0.6 * 0.5 * x[ix] * 200;

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
    cout << "  NBD : mu=" << NBD_par[isys][0] << " k=" << NBD_par[isys][1] << endl;
    cout << "  Ntuple name: " << ntpName[isys] << endl;

    ftrigeff->SetParameters(eff_par[isys][0], eff_par[isys][1]);

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

      ntp->Draw("Sum$(NcollA):Sum$(NcollModA) >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
                "", "goff");
      hNcoll_NcollMod[isys][ix] = (TH2D*) gDirectory->FindObject("htmp");
      hNcoll_NcollMod[isys][ix]->SetDirectory(0);
      hNcoll_NcollMod[isys][ix]->SetName(Form("hNcoll_NcollMod_%i_%i", isys, ix));
      hNcoll_NcollMod[isys][ix]->SetTitle(";N_{coll}^{mod};N_{coll}");

      ntp->Draw("Sum$(NcollA):NcollA[modA] >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
                "", "goff");
      hNcoll_NcollA[isys][ix] = (TH2D*) gDirectory->FindObject("htmp");
      hNcoll_NcollA[isys][ix]->SetDirectory(0);
      hNcoll_NcollA[isys][ix]->SetName(Form("hNcoll_NcollA_%i_%i", isys, ix));
      hNcoll_NcollA[isys][ix]->SetTitle(";N_{coll}(A);N_{coll}");

      ntp->Draw("Sum$(NcollA):NcollModA[modA] >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
                "", "goff");
      hNcoll_NcollModA[isys][ix] = (TH2D*) gDirectory->FindObject("htmp");
      hNcoll_NcollModA[isys][ix]->SetDirectory(0);
      hNcoll_NcollModA[isys][ix]->SetName(Form("hNcoll_NcollModA_%i_%i", isys, ix));
      hNcoll_NcollModA[isys][ix]->SetTitle(";N_{coll}^{mod}(A);N_{coll}");

      ntp->Draw("Sum$(NcollModA):NcollModA[modA] >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
                "", "goff");
      hNcollMod_NcollModA[isys][ix] = (TH2D*) gDirectory->FindObject("htmp");
      hNcollMod_NcollModA[isys][ix]->SetDirectory(0);
      hNcollMod_NcollModA[isys][ix]->SetName(Form("hNcollMod_NcollModA_%i_%i", isys, ix));
      hNcollMod_NcollModA[isys][ix]->SetTitle(";N_{coll}^{mod}(A);N_{coll}^{mod}");

      ntp->Draw("Sum$(NcollModA):NcollA[modA] >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
                "", "goff");
      hNcollMod_NcollA[isys][ix] = (TH2D*) gDirectory->FindObject("htmp");
      hNcollMod_NcollA[isys][ix]->SetDirectory(0);
      hNcollMod_NcollA[isys][ix]->SetName(Form("hNcollMod_NcollA_%i_%i", isys, ix));
      hNcollMod_NcollA[isys][ix]->SetTitle(";N_{coll}(A);N_{coll}^{mod}");

      // we're done with the file, might as well close it
      delete ntp;
      fin->Close();
      delete fin;

      // Get the MB NColl distributions
      hNcoll_MB[isys][ix] = (TH1D*) hNcoll_NcollMod[isys][ix]->ProjectionY(
                              Form("hNcoll_MB_%i_%i", isys, ix),
                              1,
                              hNcoll_NcollMod[isys][ix]->GetNbinsX());

      hNcollMod_MB[isys][ix] = (TH1D*) hNcoll_NcollMod[isys][ix]->ProjectionX(
                                 Form("hNcollMod_MB_%i_%i", isys, ix),
                                 1,
                                 hNcoll_NcollMod[isys][ix]->GetNbinsY());


      // convolve with NBD to get BBCs charge & apply trigger efficiency
      for (int ibx = 1; ibx <= hNcoll_NcollMod[isys][ix]->GetNbinsX(); ibx++)
      {
        for (int iby = 1; iby <= hNcoll_NcollMod[isys][ix]->GetNbinsY(); iby++)
        {

          float w = hNcoll_NcollMod[isys][ix]->GetBinContent(ibx, iby);
          float Ncoll = hNcoll_NcollMod[isys][ix]->GetYaxis()->GetBinCenter(iby);
          float NcollMod = hNcoll_NcollMod[isys][ix]->GetXaxis()->GetBinCenter(ibx);

          if (w <= 0) continue;

          for (int j = 1; j <= hNcoll_BBCsc[isys][ix]->GetNbinsX(); j++)
          {
            double c = hNcoll_BBCsc[isys][ix]->GetXaxis()->GetBinCenter(j);
            double NBD = evalNBD(c, Ncoll * NBD_par[isys][0],
                                 Ncoll * NBD_par[isys][1]);
            double NBDMod = evalNBD(c, NcollMod * NBD_par[isys][0],
                                    NcollMod * NBD_par[isys][1]);
            double eff = ftrigeff->Eval(c);

            // NBD=nan for mu/k = 0
            if (Ncoll > 0)
            {
              hNcoll_BBCsc[isys][ix]->Fill(c, Ncoll, w * NBD * eff);
              hBBCsc[isys][ix]->Fill(c, w * NBD * eff);
            }
            if (NcollMod > 0)
            {
              hNcollMod_BBCscMod[isys][ix]->Fill(c, NcollMod, w * NBDMod * eff);
            }
          } // j
        } // iby
      } // ibx

      // Calculate the BBC charge weighted by Ncoll of the modified nucleon
      for (int ibx = 1; ibx <= hNcoll_NcollA[isys][ix]->GetNbinsX(); ibx++)
      {
        for (int iby = 1; iby <= hNcoll_NcollA[isys][ix]->GetNbinsY(); iby++)
        {
          double w = hNcoll_NcollA[isys][ix]->GetBinContent(ibx, iby);
          double wMod = hNcollMod_NcollModA[isys][ix]->GetBinContent(ibx, iby);
          double wModA = hNcoll_NcollModA[isys][ix]->GetBinContent(ibx, iby);
          double wModBBCsc = hNcollMod_NcollA[isys][ix]->GetBinContent(ibx, iby);

          if (!(w > 0 || wMod > 0 || wModA > 0 || wModBBCsc > 0)) continue;

          double NcollA = hNcoll_NcollA[isys][ix]->GetXaxis()->GetBinCenter(ibx);
          double Ncoll = hNcoll_NcollA[isys][ix]->GetYaxis()->GetBinCenter(iby);

          for (int j = 1; j <= hBBCscNcollA[isys][ix]->GetNbinsX(); j++)
          {
            double c = hBBCscNcollA[isys][ix]->GetXaxis()->GetBinCenter(j);
            double NBD = evalNBD(c, Ncoll * NBD_par[isys][0], Ncoll * NBD_par[isys][1]);
            double eff = ftrigeff->Eval(c);

            if (Ncoll > 0)
            {
              hBBCsc[isys][ix]->Fill(c, w * NBD * eff);

              hBBCscNcollA[isys][ix]->Fill(c, NcollA * w * NBD * eff);
              hBBCscModNcollModA[isys][ix]->Fill(c, NcollA * wMod * NBD * eff);

              hBBCscNcollModA[isys][ix]->Fill(c, NcollA * wModA * NBD * eff);
              hBBCscModNcollA[isys][ix]->Fill(c, NcollA * wModBBCsc * NBD * eff);
            }
          } // j
        } // iby
      } // ibx

      //-- Deal with charge distributions --//
      mean_BBCs[isys][ix] = hBBCscModNcollModA[isys][ix]->GetMean();

      //-- Calculate RCP --//
      // get the BBCcs limits for each centrality bin
      double xq_c[NCENT];
      double yq_c[NCENT] = {0};
      double blim[NCENT][2];
      for (int icent = 0; icent < NCENT; icent++)
        xq_c[icent] = (centl[isys][NCENT] - centl[isys][icent + 1]) / centl[isys][NCENT];

      hBBCsc[isys][ix]->GetQuantiles(NCENT, yq_c, xq_c);

      for (int icent = 0; icent < NCENT; icent++)
      {
        blim[icent][0] = hBBCsc[isys][ix]->FindBin(yq_c[icent]);
        if (icent > 0)
          blim[icent][1] = hBBCsc[isys][ix]->FindBin(yq_c[icent - 1]) - 1;
        else
          blim[icent][1] = hBBCsc[isys][ix]->GetNbinsX();

        cout << "  xq[" << icent << "]: " << xq_c[icent]
             << " yq[" << icent << "]: " << yq_c[icent] << endl;
        cout << "    blim[" << icent << "]:"
             << " [" << blim[icent][0] << ", " << blim[icent][1] << "]"
             << " -> [" << hBBCsc[isys][ix]->GetBinLowEdge(blim[icent][0]) << ","
             << hBBCsc[isys][ix]->GetBinLowEdge(blim[icent][1] + 1) << "]"
             << endl;
        double frac = hBBCsc[isys][ix]->Integral(blim[icent][0], blim[icent][1]);
        frac /= hBBCsc[isys][ix]->Integral();
        cout << "    frac: " << frac << endl;

      }


      double yieldA[NCENT] = {0};
      double yieldA_NcollModABBCscMod[NCENT] = {0};
      double yieldA_NcollModABBCsc[NCENT] = {0};
      double yieldA_NcollABBCscMod[NCENT] = {0};

      // calculate all the bias factors
      cout << endl;
      for (int icent = 0; icent < NCENT; icent++)
      {
        // get Ncoll distributions for each centrality bin
        hNcoll_cent[isys][ix][icent] = (TH1D*) hNcoll_BBCsc[isys][ix]->ProjectionY(
                                         Form("hNcoll_cent_%i_%i_%i", isys, ix, icent),
                                         blim[icent][0], blim[icent][1]);

        hNcollMod_cent[isys][ix][icent] = (TH1D*) hNcollMod_BBCscMod[isys][ix]->ProjectionY(
                                            Form("hNcollMod_cent_%i_%i_%i", isys, ix, icent),
                                            blim[icent][0], blim[icent][1]);

        // get "jet" yields
        yieldA[icent] = hBBCscNcollA[isys][ix]->Integral(blim[icent][0], blim[icent][1]);
        yieldA_NcollModABBCscMod[icent] = hBBCscModNcollModA[isys][ix]->Integral(blim[icent][0], blim[icent][1]);
        yieldA_NcollModABBCsc[icent] = hBBCscNcollModA[isys][ix]->Integral(blim[icent][0], blim[icent][1]);
        yieldA_NcollABBCscMod[icent] = hBBCscModNcollA[isys][ix]->Integral(blim[icent][0], blim[icent][1]);

        // Calculate the bias factor for each centrality bin
        bias_NcollModABBCscMod[isys][ix][icent] = yieldA_NcollModABBCscMod[icent] / yieldA[icent];
        bias_NcollModABBCsc[isys][ix][icent] = yieldA_NcollModABBCsc[icent] / yieldA[icent];
        bias_NcollABBCscMod[isys][ix][icent] = yieldA_NcollABBCscMod[icent] / yieldA[icent];

        // print
        cout << " cent: " << centl[isys][icent] << " - " << centl[isys][icent + 1] << endl;
        cout << "    <Ncoll>   : " << hNcoll_cent[isys][ix][icent]->GetMean() << endl;
        cout << "    <NcollMod>: " << hNcollMod_cent[isys][ix][icent]->GetMean() << endl;
        cout << "    Bias (both): " << bias_NcollModABBCscMod[isys][ix][icent] << endl;
        cout << "    Bias (Mod A): " << bias_NcollModABBCsc[isys][ix][icent] << endl;
        cout << "    Bias (Mod B): " << bias_NcollABBCscMod[isys][ix][icent] << endl;

      } //icent

      // Calculate Rcp
      cout << endl;
      cout << " Rcp: " << endl;
      for (int icent = 0; icent < NCENT - 1; icent++)
      {
        rcp_NcollModABBCscMod[isys][icent][ix] = bias_NcollModABBCscMod[isys][ix][icent];
        rcp_NcollModABBCscMod[isys][icent][ix] /= bias_NcollModABBCscMod[isys][ix][NCENT - 1];

        cout << "   " << centl[isys][icent] << " - " << centl[isys][icent + 1] << ": "
             << rcp_NcollModABBCscMod[isys][icent][ix]
             << endl;
      }


      if (ix == 0)
      {
        for (int icent = 0; icent < NCENT; icent++)
        {
          hBBCsc_cent[isys][icent] = (TH1D*) hBBCsc[isys][ix]->Clone(
                                       Form("hBBCsc_cent_%i_%i", isys, icent));
          hBBCsc_cent[isys][icent]->SetFillColorAlpha(fcolor_cent[icent], 0.8);
          hBBCsc_cent[isys][icent]->Reset();

          for (int ib = blim[icent][0]; ib <= blim[icent][1]; ib++)
          {
            double bc = hBBCsc[isys][ix]->GetBinContent(ib);
            hBBCsc_cent[isys][icent]->SetBinContent(ib, bc);
          }
        }
      }

    } // ix

    for (int icent = 0; icent < NCENT - 1; icent++)
    {
      grcp_NcollModABBCscMod[isys][icent] = new TGraph(NX, x, rcp_NcollModABBCscMod[isys][icent]);
      grcp_NcollModABBCscMod[isys][icent]->SetName(Form("grcp_NcollModABBCscMod_%i_%i", isys, icent));
      grcp_NcollModABBCscMod[isys][icent]->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
      grcp_NcollModABBCscMod[isys][icent]->SetLineStyle(1);
      grcp_NcollModABBCscMod[isys][icent]->SetLineWidth(2);
      grcp_NcollModABBCscMod[isys][icent]->SetLineColor(colors[isys]);
      grcp_NcollModABBCscMod[isys][icent]->SetLineStyle(lstyle[isys]);

      // vs pion pT
      grcp_NcollModABBCscMod_pipT[isys][icent] = new TGraph(NX, pi0_pT, rcp_NcollModABBCscMod[isys][icent]);
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetName(Form("grcp_NcollModABBCscMod_pipT_%i_%i", isys, icent));
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetTitle(";#pi p_{T};R_{CP}");
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetLineStyle(1);
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetLineWidth(2);
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetLineColor(colors[isys]);
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetLineStyle(lstyle[isys]);

    }

    cout << "<Q>" << endl;
    double denom = mean_BBCs[isys][1];
    for (int ix = 0; ix < NX; ix++)
      mean_BBCs[isys][ix] = mean_BBCs[isys][ix] / denom;

    gmean_BBCs[isys] = new TGraph(NX, x, mean_BBCs[isys]);
    gmean_BBCs[isys]->SetName(Form("gmean_BBCs_%i", isys));
    gmean_BBCs[isys]->SetTitle(Form(";x;<Q_{BBC, Au}> / <Q_{BBC,Au}(x=%.2f)", x[1]));
    gmean_BBCs[isys]->SetLineWidth(2);
    gmean_BBCs[isys]->SetLineColor(colors[isys]);
    gmean_BBCs[isys]->SetLineStyle(lstyle[isys]);

    gmean_BBCs_pipT[isys] = new TGraph(NX, pi0_pT, mean_BBCs[isys]);
    gmean_BBCs_pipT[isys]->SetName(Form("gmean_BBCs_%i", isys));
    gmean_BBCs_pipT[isys]->SetTitle(Form(";#pi p_{T} = 0.6#frac{#sqrt{s}*x}{2};<Q_{BBC, Au}> / <Q_{BBC,Au}(x=%.2f)", x[1]));
    gmean_BBCs_pipT[isys]->SetLineWidth(2);
    gmean_BBCs_pipT[isys]->SetLineColor(colors[isys]);
    gmean_BBCs_pipT[isys]->SetLineStyle(lstyle[isys]);

  } // isys

//=====================================================//
// DO SOME PRINTING
//=====================================================//
  cout << endl;
  cout << "--> Printing ... " << endl;

  cout << setprecision(2) << fixed;
  cout << endl;
  cout << " ncoll MB " << endl;
  cout << "x";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " & " << collSystem[isys];
  cout << " \\\\" << endl;
  cout << "0.00";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " & " << hNcoll_MB[isys][0]->GetMean();
  cout << "\\\\" << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    for (int isys = 0; isys < NSYSTEMS; isys++)
      cout << " & " << hNcollMod_MB[isys][ix]->GetMean();
    cout << "\\\\" << endl;
  }

  cout << endl;
  cout << "-- ncoll ==" << endl;
  cout << "x sig(x)";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  cout << endl;
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cout << " & 0-100";
    for (int icent = 0; icent < NCENT; icent++)
      cout << " & " << centl[isys][icent] << "-" << centl[isys][icent + 1];
  }
  cout << " \\\\" << endl;
  cout << "PHENIX & 42.00";
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cout << " & " << Ncoll_PHENIX[isys][NCENT];
    for (int icent = 0; icent < NCENT; icent++)
      cout << " & " << Ncoll_PHENIX[isys][icent];
  }
  cout << "0.00 & 42.00";
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cout << " & " << hNcoll_MB[isys][0]->GetMean();
    for (int icent = 0; icent < NCENT; icent++)
      cout << " & " << hNcoll_cent[isys][0][icent]->GetMean();
  }
  cout << " \\\\" << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      cout << " & " << hNcollMod_MB[isys][ix]->GetMean();
      for (int icent = 0; icent < NCENT; icent++)
        cout << " & " << hNcollMod_cent[isys][ix][icent]->GetMean();
    }
    cout << " \\\\" << endl;
  }

  cout << endl;
  cout << " Bias factors (both mod)" << endl;
  cout << "x & sig(x)";
  for (int icent = 0; icent < NCENT; icent++)
  {
    cout << " & " << centl[1][icent] << "-" << centl[1][icent + 1];
  }
  cout << "\\\\" << endl;
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cout << " & & \\multicolumn{" << NCENT << "}{|c}{"
         << collSystem[isys] << "} \\\\" << endl;
    for (int ix = 0; ix < NX; ix++)
    {
      cout << x[ix];
      cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
      for (int icent = 0; icent < NCENT; icent++)
        cout << " & " << bias_NcollModABBCscMod[isys][ix][icent];
      cout << "\\\\" << endl;
    }
  }

  cout << endl;
  cout << " Bias factors (mod A)" << endl;
  cout << "x";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " & " << collSystem[isys];
  cout << "\\\\" << endl;
// for (int isys = 0; isys < NSYSTEMS; isys++)
// {
//   for (int icent = 0; icent < NCENT; icent++)
//     cout << " & " << centl[isys][icent] << "-" << centl[isys][icent +1];
// }
  cout << " \\\\" << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    // cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      // for (int icent = 0; icent < NCENT; icent++)
      cout << " & " << bias_NcollModABBCsc[isys][ix][0];
    }
    cout << " | ";
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      // for (int icent = 0; icent < NCENT; icent++)
      cout << " & " << bias_NcollModABBCsc[isys][ix][NCENT - 1];
    }
    cout << " \\\\" << endl;
  }

  cout << endl;
  cout << " Bias factors (mod BBCsc)" << endl;
  cout << "x";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " & " << collSystem[isys];
  cout << "\\\\" << endl;
// for (int isys = 0; isys < NSYSTEMS; isys++)
// {
//   for (int icent = 0; icent < NCENT; icent++)
//     cout << " & " << centl[isys][icent] << "-" << centl[isys][icent + 1];
// }
  cout << " \\\\" << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    // cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      // for (int icent = 0; icent < NCENT; icent++)
      cout << " & " << bias_NcollABBCscMod[isys][ix][0];
    }
    cout << " | ";
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      // for (int icent = 0; icent < NCENT; icent++)
      cout << " & " << bias_NcollABBCscMod[isys][ix][NCENT - 1];
    }
    cout << " \\\\" << endl;
  }


  cout << endl;
  cout << " Rcp cent " << endl;
  cout << "x ";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " & " << collSystem[isys];
  cout << "\\\\" << endl;
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    for (int icent = 0; icent < NCENT; icent++)
      cout << " & " << centl[isys][icent] << "-" << centl[isys][icent + 1];
  }
  cout << " \\\\" << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    // cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      for (int icent = 0; icent < NCENT; icent++)
        cout << " & " << rcp_NcollModABBCscMod[isys][icent][ix];
    }
    cout << " \\\\" << endl;
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
  double pTl_jet[NJET] = {12.1, 14.5, 17.3, 20.7, 24.7, 29.4, 35.1, 41.9};
  double pTh_jet[NJET] = {14.5, 17.3, 20.7, 24.7, 29.4, 35.1, 41.9, 50.0};
  double Rcp_jet[NCENT - 1][NJET] = {
    {0.73, 0.71, 0.66, 0.61, 0.57, 0.54, 0.52, 0.50}, // 00-20% / 60-88%
    {0.92, 0.88, 0.85, 0.82, 0.76, 0.71, 0.67, 0.63}, // 20-40% / 60-88%
    {0.98, 0.96, 0.95, 0.93, 0.91, 0.88, 0.86, 0.84}, // 40-60% / 60-88%
  };
  double Rcp_jet_A[NCENT - 1][NJET] = {
    {0.01, 0.01, 0.02, 0.02, 0.02, 0.04, 0.05, 0.06}, // 00-20% / 60-88%
    {0.01, 0.02, 0.02, 0.02, 0.04, 0.05, 0.07, 0.08}, // 20-40% / 60-88%
    {0.02, 0.02, 0.02, 0.03, 0.04, 0.07, 0.09, 0.11}, // 40-60% / 60-88%
  };
  double Rcp_jet_B[NCENT - 1][NJET] = {
    {0.06, 0.03, 0.02, 0.02, 0.03, 0.05, 0.08, 0.09}, // 00-20% / 60-88%
    {0.05, 0.04, 0.03, 0.04, 0.07, 0.11, 0.16, 0.19}, // 20-40% / 60-88%
    {0.04, 0.05, 0.05, 0.06, 0.05, 0.07, 0.10, 0.12}, // 40-60% / 60-88%
  };

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

//=====================================================//
// PLOT OBJECTS
//=====================================================//
  cout << endl;
  cout << "--> Plotting ..." << endl;

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);

  TH1F* haxis_rcp = new TH1F("haxis_rcp",
                             ";x (x_{jet}=2 * p_{T}^{jet} / #sqrt{s_{NN}});R_{CP}",
                             100, 0, 1);
  haxis_rcp->GetYaxis()->SetTitleOffset(1.3);
  haxis_rcp->GetXaxis()->SetTitleOffset(1.3);
  haxis_rcp->GetYaxis()->CenterTitle();
  haxis_rcp->SetMinimum(0.01);
  haxis_rcp->SetMaximum(1.09);

  TH1F* haxis_rcp_pipT = new TH1F("haxis_rcp_pipT",
                                  ";#pi p_{T};R_{CP}",
                                  100, 0, 20);
  haxis_rcp_pipT->GetYaxis()->SetTitleOffset(1.3);
  haxis_rcp_pipT->GetXaxis()->SetTitleOffset(1.3);
  haxis_rcp_pipT->GetYaxis()->CenterTitle();
  haxis_rcp_pipT->SetMinimum(0.01);
  haxis_rcp_pipT->SetMaximum(1.09);

  TLine l1;
  l1.SetLineStyle(2);

  TLegend *legrcp[NCENT - 1];
  for (int icent = 0; icent < NCENT - 1; icent++)
  {
    legrcp[icent] = new TLegend(0.15, 0.15, 0.4, 0.4,
                                Form("(%.0f-%.0f%%) / (%.0f-%.0f%%)",
                                     centl[1][icent], centl[1][icent + 1],
                                     centl[1][NCENT - 1], centl[1][NCENT]));
    legrcp[icent]->SetFillStyle(0);
    legrcp[icent]->SetBorderSize(0);
    legrcp[icent]->SetTextSize(0.04);
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      legrcp[icent]->AddEntry(grcp_NcollModABBCscMod[isys][icent],
                              collSystem[isys], "L");
    }
    legrcp[icent]->AddEntry(gRCP_jet[icent], "PHENIX, arXiv:1509.04657", "P");

  }

  int xplot = 2;

  TLegend *legBBC[NSYSTEMS];
  TLegend *legNcollMB[NSYSTEMS];
  TLegend *legNcoll[NSYSTEMS][NCENT];
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    //-- BBCs charge
    hBBCscNcollA[isys][xplot]->SetLineColor(kBlack);
    hBBCscModNcollA[isys][xplot]->SetLineColor(kBlue);
    hBBCscNcollModA[isys][xplot]->SetLineColor(kGreen + 2);
    hBBCscModNcollModA[isys][xplot]->SetLineColor(kRed);

    hBBCscNcollA[isys][xplot]->SetLineWidth(2);
    hBBCscModNcollModA[isys][xplot]->SetLineWidth(2);

    legBBC[isys] = new TLegend(0.15, 0.15, 0.4, 0.4,
                               Form("%s @ 200 GeV, x = %.2f",
                                    collSystem[isys], x[xplot]));
    legBBC[isys]->SetFillStyle(0);
    legBBC[isys]->SetBorderSize(0);
    legBBC[isys]->SetTextSize(0.04);
    legBBC[isys]->AddEntry(hBBCscNcollA[isys][xplot], "Unmodified", "L");
    // legBBC[isys]->AddEntry(hBBCscModNcollA[isys][xplot], "Mod BBCs chg", "L");
    // legBBC[isys]->AddEntry(hBBCscNcollModA[isys][xplot], "Mod N_{coll}", "L");
    legBBC[isys]->AddEntry(hBBCscModNcollModA[isys][xplot], "Mod BBCs chg and N_{coll}", "L");

    //-- MB Ncoll
    hNcoll_MB[isys][xplot]->SetLineColor(kBlue);
    hNcollMod_MB[isys][xplot]->SetLineColor(kRed);

    hNcoll_MB[isys][xplot]->SetLineWidth(2);
    hNcollMod_MB[isys][xplot]->SetLineWidth(2);

    if (isys == 0)
    {
      legNcollMB[isys] = new TLegend(0.5, 0.6, 0.75, 0.85,
                                     Form("%s @ 200 GeV, x = %.2f",
                                          collSystem[isys], x[xplot]));
    }
    else
    {
      legNcollMB[isys] = new TLegend(0.15, 0.15, 0.4, 0.4,
                                     Form("%s @ 200 GeV, x = %.2f",
                                          collSystem[isys], x[xplot]));
    }
    legNcollMB[isys]->SetFillStyle(0);
    legNcollMB[isys]->SetBorderSize(0);
    legNcollMB[isys]->SetTextSize(0.04);
    legNcollMB[isys]->AddEntry((TObject*)0,
                               Form("<N_{coll}> = %.2f (PHENIX)",
                                    Ncoll_PHENIX[isys][NCENT]), "");
    legNcollMB[isys]->AddEntry(hNcoll_MB[isys][xplot],
                               Form("<N_{coll}> = %.2f",
                                    hNcoll_MB[isys][xplot]->GetMean())
                               , "L");
    legNcollMB[isys]->AddEntry(hNcollMod_MB[isys][xplot],
                               Form("<N_{coll}^{mod}> = %.2f",
                                    hNcollMod_MB[isys][xplot]->GetMean())
                               , "L");

    //-- centrality dep Ncoll
    for (int icent = 0; icent < NCENT; icent++)
    {
      hNcoll_cent[isys][xplot][icent]->SetLineColor(kBlue);
      hNcollMod_cent[isys][xplot][icent]->SetLineColor(kRed);

      hNcoll_cent[isys][xplot][icent]->SetLineWidth(2);
      hNcollMod_cent[isys][xplot][icent]->SetLineWidth(2);

      if (isys == 0 || icent >= NCENT - 2)
      {
        legNcoll[isys][icent] = new TLegend(0.6, 0.5, 0.9, 0.85,
                                            Form("%s, x = %.2f, %.0f-%.0f%%",
                                                collSystem[isys], x[xplot],
                                                centl[isys][icent], centl[isys][icent + 1]));
      }
      else
      {
        legNcoll[isys][icent] = new TLegend(0.15, 0.12, 0.4, 0.45,
                                            Form("%s @ 200 GeV, x = %.2f",
                                                collSystem[isys], x[xplot]));
      }
      legNcoll[isys][icent]->SetFillStyle(0);
      legNcoll[isys][icent]->SetBorderSize(0);
      legNcoll[isys][icent]->SetTextSize(0.06);
      legNcoll[isys][icent]->AddEntry((TObject*)0,
                                      Form("<N_{coll}> = %.2f (PHENIX)",
                                           Ncoll_PHENIX[isys][icent]), "");
      legNcoll[isys][icent]->AddEntry(hNcoll_cent[isys][xplot][icent],
                                      Form("<N_{coll}> = %.2f",
                                           hNcoll_cent[isys][xplot][icent]->GetMean())
                                      , "L");
      legNcoll[isys][icent]->AddEntry(hNcollMod_cent[isys][xplot][icent],
                                      Form("<N_{coll}^{mod}> = %.2f",
                                           hNcollMod_cent[isys][xplot][icent]->GetMean())
                                      , "L");


    }

  }

  TLegend *legQ = new TLegend(0.15, 0.15, 0.5, 0.4);
  legQ->SetFillStyle(0);
  legQ->SetBorderSize(0);
  legQ->SetTextSize(0.04);
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    legQ->AddEntry(gmean_BBCs[isys], collSystem[isys], "L");
  }

  TLegend *legx = new TLegend(0.8, 0.5, 0.98, 0.98);
  legx->SetFillStyle(0);
  legx->SetBorderSize(0);
  legx->SetTextSize(0.04);
  legx->AddEntry(hBBCscNcollA[0][0], "x_{p} = 0.00", "L");
  for (int ix = 0; ix < NX; ix++)
  {
    legx->AddEntry(hBBCscModNcollModA[0][ix],
                   Form("x_{p} = %.2f", x[ix]),
                   "L");
  }

  //=====================================================//
  // PLOT
  //=====================================================//

  TCanvas *crcp = new TCanvas("crcp", "RCP", 1400, 400);
  crcp->SetTopMargin(0.0);
  crcp->SetRightMargin(0.0);
  crcp->SetBottomMargin(0.0);
  crcp->SetLeftMargin(0.0);
  crcp->Divide(NCENT - 1, 1, 0, 0);

  for (int icent = 0; icent < NCENT - 1; icent++)
  {
    crcp->GetPad(icent + 1)->SetTopMargin(0.02);
    crcp->GetPad(icent + 1)->SetRightMargin(0.02);
    crcp->GetPad(icent + 1)->SetBottomMargin(0.10);
    crcp->GetPad(icent + 1)->SetLeftMargin(0.10);
    crcp->GetPad(icent + 1)->SetTicks(1, 1);

    crcp->cd(icent + 1);
    haxis_rcp->GetXaxis()->SetRangeUser(0, 0.5);
    haxis_rcp->Draw();

    for (int i = 0; i < NJET; i++)
      bRcp_jet[icent][i]->Draw();
    gRCP_jet[icent]->Draw("P");


    for (int isys = 0; isys < NSYSTEMS; isys++)
      grcp_NcollModABBCscMod[isys][icent]->Draw("C");


    l1.DrawLine(0, 1, 0.5, 1);
    legrcp[icent]->Draw("same");
  }

  TCanvas *crcppipt = new TCanvas("crcppipt", "RCP", 1400, 400);
  crcppipt->SetTopMargin(0.0);
  crcppipt->SetRightMargin(0.0);
  crcppipt->SetBottomMargin(0.0);
  crcppipt->SetLeftMargin(0.0);
  crcppipt->Divide(NCENT - 1, 1, 0, 0);

  for (int icent = 0; icent < NCENT - 1; icent++)
  {
    crcppipt->GetPad(icent + 1)->SetTopMargin(0.02);
    crcppipt->GetPad(icent + 1)->SetRightMargin(0.02);
    crcppipt->GetPad(icent + 1)->SetBottomMargin(0.10);
    crcppipt->GetPad(icent + 1)->SetLeftMargin(0.10);
    crcppipt->GetPad(icent + 1)->SetTicks(1, 1);

    crcppipt->cd(icent + 1);
    haxis_rcp_pipT->GetXaxis()->SetRangeUser(0, 20);
    haxis_rcp_pipT->Draw();

    for (int isys = 0; isys < NSYSTEMS; isys++)
      grcp_NcollModABBCscMod_pipT[isys][icent]->Draw("C");

    l1.DrawLine(0, 1, 20., 1);
    legrcp[icent]->Draw("same");
  }

  TCanvas *cbbc = new TCanvas("cbbc", "bbc", 600, 1200);
  cbbc->SetTopMargin(0.00);
  cbbc->SetRightMargin(0.00);
  cbbc->SetBottomMargin(0.00);
  cbbc->SetLeftMargin(0.00);
  cbbc->Divide(1, NSYSTEMS, 0, 0);

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cbbc->GetPad(isys + 1)->SetTopMargin(0.02);
    cbbc->GetPad(isys + 1)->SetRightMargin(0.02);
    cbbc->GetPad(isys + 1)->SetBottomMargin(0.10);
    cbbc->GetPad(isys + 1)->SetLeftMargin(0.10);
    cbbc->GetPad(isys + 1)->SetTicks(1, 1);

    cbbc->cd(isys + 1);
    gPad->SetLogy();
    hBBCsc[isys][0]->GetYaxis()->SetRangeUser(1e2, 2e4);
    hBBCsc[isys][0]->GetXaxis()->SetRangeUser(1, 80);
    hBBCsc[isys][0]->SetTitle(";Q_{BBC,Au}");
    hBBCsc[isys][0]->DrawCopy("hist");

    for (int i = 0; i < NCENT; i++)
      hBBCsc_cent[isys][i]->Draw("hist same");

    label.DrawLatex(0.7, 0.7, collSystem[isys]);
  }

  TCanvas *cyield = new TCanvas("cyield", "yield", 1400, 400);
  cyield->Divide(NSYSTEMS, 1, 0, 0);
  label.DrawLatex(0.5, 0.96, "BBCsc charge for events with high-p_{T} particle");
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cyield->GetPad(isys + 1)->SetTopMargin(0.02);
    cyield->GetPad(isys + 1)->SetRightMargin(0.02);
    cyield->GetPad(isys + 1)->SetBottomMargin(0.10);
    cyield->GetPad(isys + 1)->SetLeftMargin(0.10);
    cyield->GetPad(isys + 1)->SetTicks(1, 1);

    cyield->cd(isys + 1);
    gPad->SetLogy();
    hBBCscNcollA[isys][xplot]->GetYaxis()->SetRangeUser(1e2, 2e4);
    hBBCscNcollA[isys][xplot]->GetXaxis()->SetRangeUser(1, 80);
    hBBCscNcollA[isys][xplot]->SetTitle(";Q_{BBC,Au}");
    hBBCscNcollA[isys][xplot]->DrawCopy("hist");
    // hBBCscNcollModA[isys][xplot]->Draw("same");
    // hBBCscModNcollA[isys][xplot]->Draw("same");
    hBBCscModNcollModA[isys][xplot]->DrawCopy("hist same");
    legBBC[isys]->Draw("same");
  }

  TCanvas* cyieldpau = new TCanvas("cyieldpau", "yield pAu", 600, 400);
  cyieldpau->SetTopMargin(0.02);
  cyieldpau->SetRightMargin(0.02);
  cyieldpau->SetBottomMargin(0.10);
  cyieldpau->SetLeftMargin(0.10);
  cyieldpau->SetTicks(1, 1);

  cyieldpau->cd(1);
  // gPad->SetLogy();
  // hBBCscNcollA[0][0]->GetYaxis()->SetRangeUser(1e2, 2e4);
  hBBCscNcollA[0][0]->GetXaxis()->SetRangeUser(1, 30);
  hBBCscNcollA[0][0]->SetLineWidth(2);
  hBBCscNcollA[0][0]->SetLineColor(kBlack);
  hBBCscNcollA[0][0]->SetTitle(";Q_{BBC,Au}");
  hBBCscNcollA[0][0]->DrawCopy("hist");
  for (int ix = 0; ix < NX; ix++)
  {
    hBBCscModNcollModA[0][ix]->SetLineWidth(1);
    hBBCscModNcollModA[0][ix]->SetLineColor(2 + ix);
    hBBCscModNcollModA[0][ix]->DrawCopy("hist same");
  }
  legx->Draw("same");

  TCanvas *cncollmb = new TCanvas("cncollmb", "MB Ncoll", 1400, 400);
  cncollmb->Divide(NSYSTEMS, 1, 0, 0);
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cncollmb->GetPad(isys + 1)->SetTopMargin(0.02);
    cncollmb->GetPad(isys + 1)->SetRightMargin(0.02);
    cncollmb->GetPad(isys + 1)->SetBottomMargin(0.10);
    cncollmb->GetPad(isys + 1)->SetLeftMargin(0.10);
    cncollmb->GetPad(isys + 1)->SetTicks(1, 1);

    cncollmb->cd(isys + 1);
    gPad->SetLogy();
    hNcoll_MB[isys][xplot]->GetXaxis()->SetRangeUser(0, 50);
    hNcoll_MB[isys][xplot]->Draw();
    hNcollMod_MB[isys][xplot]->Draw("same");
    legNcollMB[isys]->Draw("same");
  }

  TCanvas *cncoll = new TCanvas("cncoll", "Ncoll", 400 * NSYSTEMS, 200 * NCENT);
  cncoll->Divide(NSYSTEMS, NCENT, 0, 0);
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    for (int icent = 0; icent < NCENT; icent++)
    {
      cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetTopMargin(0.02);
      cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetRightMargin(0.02);
      cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetBottomMargin(0.10);
      cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetLeftMargin(0.10);
      cncoll->GetPad(icent * NSYSTEMS + isys + 1)->SetTicks(1, 1);

      cncoll->cd(icent * NSYSTEMS + isys + 1);
      gPad->SetLogy();
      hNcoll_cent[isys][xplot][icent]->GetXaxis()->SetRangeUser(0, 50);
      hNcoll_cent[isys][xplot][icent]->GetYaxis()->SetRangeUser(1, 1e5);
      hNcoll_cent[isys][xplot][icent]->Draw();
      hNcollMod_cent[isys][xplot][icent]->Draw("same");
      legNcoll[isys][icent]->Draw("same");
    }
  }

  TCanvas *cq = new TCanvas("cq", "<Q>", 900, 600);
  cq->SetTopMargin(0.02);
  cq->SetRightMargin(0.02);
  cq->SetBottomMargin(0.10);
  cq->SetLeftMargin(0.10);
  cq->SetTicks(1, 1);

  cq->cd(1);
  gmean_BBCs_pipT[0]->GetYaxis()->SetTitleOffset(1.30);
  gmean_BBCs_pipT[0]->GetYaxis()->SetRangeUser(0.7, 1.2);
  gmean_BBCs_pipT[0]->GetXaxis()->SetRangeUser(0, 20);
  gmean_BBCs_pipT[0]->Draw("AC");
  for (int isys = 1; isys < NSYSTEMS; isys++)
    gmean_BBCs_pipT[isys]->Draw("C");
  l1.DrawLine(0, 1, 0.5, 1);
  legQ->Draw("same");



//=====================================================//
// FIGURES FOR PAPER
//=====================================================//

  const char *lsystem[] = {"p+Au", "d+Au", "^{3}He+Au"};
  TLegend *legrcp_paper[NCENT - 1];
  for (int icent = 0; icent < NCENT - 1; icent++)
  {
    double x1, y1;
    if (icent == NCENT - 2)
    {
      x1 = 0.2;
      y1 = 0.2;
    }
    else
    {
      x1 = 0.2;
      y1 = 0.05;
    }
    legrcp_paper[icent] = new TLegend(x1, y1, x1 + 0.25, y1 + 0.35,
                                      Form("(%.0f-%.0f%%) / (%.0f-%.0f%%)",
                                           centl[1][icent], centl[1][icent + 1],
                                           centl[1][NCENT - 1], centl[1][NCENT]));
    legrcp_paper[icent]->SetFillStyle(0);
    legrcp_paper[icent]->SetBorderSize(0);
    legrcp_paper[icent]->SetTextSize(0.06);
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      legrcp_paper[icent]->AddEntry(grcp_NcollModABBCscMod[isys][icent],
                                    lsystem[isys], "L");
    }
    legrcp_paper[icent]->AddEntry(gRCP_jet[icent], "PHENIX, arXiv:1509.04657", "P");

  }

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
  haxis_rcp_paper->GetYaxis()->SetTitleOffset(0.9);
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

  TLatex lpaper;
  lpaper.SetNDC();
  lpaper.SetTextAlign(22);
  lpaper.SetTextSize(0.08);


  TCanvas *crcp_paper = new TCanvas("crcp_paper", "RCP", 400, 1200);
  crcp_paper->SetTopMargin(0.0);
  crcp_paper->SetRightMargin(0.0);
  crcp_paper->SetBottomMargin(0.0);
  crcp_paper->SetLeftMargin(0.0);
  crcp_paper->Divide(1, NCENT - 1, 0, 0);

  const char *plabel[4] = {"(a)", "(b)", "(c)", "(d)"};
  for (int icent = 0; icent < NCENT - 1; icent++)
  {
    if (icent == 0)
      crcp_paper->GetPad(icent + 1)->SetTopMargin(0.05);
    else
      crcp_paper->GetPad(icent + 1)->SetTopMargin(0.00);
    if (icent == NCENT - 2)
      crcp_paper->GetPad(icent + 1)->SetBottomMargin(0.15);
    else
      crcp_paper->GetPad(icent + 1)->SetBottomMargin(0.00);

    crcp_paper->GetPad(icent + 1)->SetRightMargin(0.05);
    crcp_paper->GetPad(icent + 1)->SetLeftMargin(0.15);
    crcp_paper->GetPad(icent + 1)->SetTicks(1, 1);

    crcp_paper->cd(icent + 1);
    haxis_rcp_paper->GetXaxis()->SetRangeUser(0, 0.5);
    haxis_rcp_paper->Draw();

    for (int i = 0; i < NJET; i++)
      bRcp_jet[icent][i]->Draw();
    gRCP_jet[icent]->Draw("P");


    for (int isys = 0; isys < NSYSTEMS; isys++)
      grcp_NcollModABBCscMod[isys][icent]->Draw("C");


    l1.DrawLine(0, 1, 0.5, 1);
    legrcp_paper[icent]->Draw("same");
    if (icent == 0)
      lpaper.DrawLatex(0.9, 0.88, plabel[icent]);
    else
      lpaper.DrawLatex(0.9, 0.93, plabel[icent]);
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


//=====================================================//
// SAVE
//=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving RCP to " << outFile << endl;

    crcp->Print("Rcp_systems.pdf");
    cyield->Print("yield_systems.pdf");
    cncoll->Print("Ncoll_systems.pdf");
    cncollmb->Print("MBNcoll_systems.pdf");
    cq->Print("Q_systems.pdf");
    cyieldpau->Print("yield_pAu.pdf");
    crcppipt->Print("Rcp_pipT_systems.pdf");
    cbbc->Print("BBCQ_systems.pdf");

    crcp_paper->Print("Rcp_paper.pdf");
    csig_paper->Print("sigmod.pdf");
  }
}