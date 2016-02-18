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
#include <TGaxis.h>

#include <data.h>

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

  // const int NX = 21;             // Number of x values
  // double x[] = {0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
  //               0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99
  //              };
  const int NX = 7;
  double x[NX] = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

  double beta = 1.5;

  const int NSYSTEMS = 1;
  const char *collSystem[] = {"dAu", "dAu", "3HeAu"};
  // const char *collSystem[] = {"pAu", "dAu", "3HeAu"};

  bool saveRcp = true;
  const char *outFile = "Rcp_systems.root";
  bool printPlots = false;

  // Files for each x value for each system
  const char *xFiles[NSYSTEMS][NX] =
  {
    {
      "rootfiles/glauber_dau_snn42_beta15_x001_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_beta15_x010_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_beta15_x020_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_beta15_x030_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_beta15_x040_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_beta15_x050_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_beta15_x060_ntuple_100k.root",
    },
  };
  // for (int isys = 0; isys < NSYSTEMS; isys++)
  // {
  //   for (int ix = 0; ix < NX; ix++)
  //   {
  //     xFiles[isys][ix] = Form("rootfiles/glauber_dau_snn42_beta%.0f_x%03.0f_ntuple_100k.root", beta * 10, x[ix] * 100);
  //     cout << xFiles[isys][ix] << endl;
  //   }
  // }
  // const char *xFiles[NSYSTEMS][NX] =
  // {
  //   { // pAu Files
  //     "rootfiles/glauber_pau_snn42_x001_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x005_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x01_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x015_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x02_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x025_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x03_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x035_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x04_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x045_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x05_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x055_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x06_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x065_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x07_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x075_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x08_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x085_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x09_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x095_ntuple_100k.root",
  //     "rootfiles/glauber_pau_snn42_x099_ntuple_100k.root",
  //   },
  //   { // dAu Files
  //     "rootfiles/glauber_dau_snn42_x001_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x005_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x01_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x015_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x02_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x025_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x03_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x035_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x04_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x045_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x05_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x055_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x06_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x065_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x07_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x075_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x08_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x085_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x09_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x095_ntuple_100k.root",
  //     "rootfiles/glauber_dau_snn42_x099_ntuple_100k.root",
  //   },
  //   { // He3Au Files
  //     "rootfiles/glauber_he3au_snn42_x001_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x005_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x01_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x015_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x02_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x025_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x03_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x035_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x04_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x045_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x05_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x055_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x06_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x065_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x07_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x075_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x08_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x085_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x09_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x095_ntuple_100k.root",
  //     "rootfiles/glauber_he3au_snn42_x099_ntuple_100k.root",
  //   },

  // };

  // Ntuple name stored in files for each system
  const char *ntpName[NSYSTEMS] =
  {
    // "nt_p_Au", // pAu
    "nt_dh_Au", // dAu
    // "nt_He3_Au", // He3Au
  };

  // Centrality bins
  const int NCENT = 4;
  double centl[NSYSTEMS][NCENT + 1] =
  {
    // {0, 20, 40, 60, 84}, // pAu
    {0, 20, 40, 60, 88}, // dAu
    // {0, 20, 40, 60, 88}, // He3Au
  };

  // PHENIX values (last entry is 0-100%)
  double Ncoll_PHENIX[NSYSTEMS][NCENT + 1] =
  {
    // { 8.20, 6.06, 4.43, 2.62, 4.67}, // pAu
    {15.1, 10.2, 6.6, 3.2, 7.6}, // dAu
    // {22.37, 14.71, 8.28, 3.38, 10.45}, //HeAu
  };

  // Negative Binomial Distribution parameters for each system
  double NBD_par[NSYSTEMS][2] =
  {
    // {3.14, 0.47}, //pAu {mu, k}
    {3.04, 0.46}, //dAu {mu, k}
    // {2.91, 0.55}, //3HeAu {mu, k}
  };

  // MB trigger efficiency function parameters
  double eff_par[NSYSTEMS][2] =
  {
    // {1.07552e+00, 6.02328e-01}, // pAu
    {0.897, 0.612}, // dAu (from D. McGlinchey thesis)
    // {1.22134e+00, 5.10114e-01}, // He3Au
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

  TH1D* hBBCscNcollA_noTrig[NSYSTEMS][NX];
  TH1D* hBBCscModNcollModA_noTrig[NSYSTEMS][NX];
  TH1D* hBBCscNcollModA_noTrig[NSYSTEMS][NX];
  TH1D* hBBCscModNcollA_noTrig[NSYSTEMS][NX];

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


      hBBCscNcollA_noTrig[i][j] = new TH1D(Form("hBBCscNcollA_noTrig_%i_%i", i, j),
                                           ";N_{coll}(A) #time BBCs charge",
                                           1596, 1, 400);

      hBBCscModNcollModA_noTrig[i][j] = new TH1D(Form("hBBCscModNcollModA_noTrig_%i_%i", i, j),
          ";N_{coll}^{mod}(A) #time BBCs charge Mod",
          1596, 1, 400);

      hBBCscNcollModA_noTrig[i][j] = new TH1D(Form("hBBCscNcollModA_noTrig_%i_%i", i, j),
                                              ";N_{coll}^{mod}(A) #time BBCs charge",
                                              1596, 1, 400);

      hBBCscModNcollA_noTrig[i][j] = new TH1D(Form("hBBCscModNcollA_noTrig_%i_%i", i, j),
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


  // bias factors
  double bias_NcollModABBCscMod_MB[NSYSTEMS][NX];
  double bias_NcollModABBCsc_MB[NSYSTEMS][NX];
  double bias_NcollABBCscMod_MB[NSYSTEMS][NX];

  double bias_NcollModABBCsc[NSYSTEMS][NX][NCENT];
  double bias_NcollABBCscMod[NSYSTEMS][NX][NCENT];
  double bias_NcollModABBCscMod[NSYSTEMS][NX][NCENT];

  // for RAA
  double raa_NcollModABBCscMod[NSYSTEMS][NCENT][NX];
  double raa_NcollModABBCsc[NSYSTEMS][NCENT][NX];
  double raa_NcollABBCscMod[NSYSTEMS][NCENT][NX];

  TGraph *graa_NcollModABBCscMod[NSYSTEMS][NCENT];
  TGraph *graa_NcollModABBCsc[NSYSTEMS][NCENT];
  TGraph *graa_NcollABBCscMod[NSYSTEMS][NCENT];

  TGraph *graa_NcollModABBCscMod_MB[NSYSTEMS];
  TGraph *graa_NcollModABBCsc_MB[NSYSTEMS];
  TGraph *graa_NcollABBCscMod_MB[NSYSTEMS];




  // for Rcp
  double rcp_NcollModABBCscMod[NSYSTEMS][NCENT - 1][NX];
  double rcp_NcollModABBCsc[NSYSTEMS][NCENT - 1][NX];
  double rcp_NcollABBCscMod[NSYSTEMS][NCENT - 1][NX];
  TGraph *grcp_NcollModABBCscMod[NSYSTEMS][NCENT - 1];
  TGraph *grcp_NcollModABBCsc[NSYSTEMS][NCENT - 1];
  TGraph *grcp_NcollABBCscMod[NSYSTEMS][NCENT - 1];
  TGraph *grcp_NcollModABBCscMod_pipT[NSYSTEMS][NCENT - 1];

  // for mean BBCs charge
  double mean_BBCs_NcollModABBCscMod[NSYSTEMS][NX];
  double mean_BBCs_NcollModABBCsc[NSYSTEMS][NX];
  double mean_BBCs_NcollABBCscMod[NSYSTEMS][NX];

  TGraph *gmean_BBCs_NcollModABBCscMod[NSYSTEMS];
  TGraph *gmean_BBCs_NcollModABBCsc[NSYSTEMS];
  TGraph *gmean_BBCs_NcollABBCscMod[NSYSTEMS];
  TGraph *gmean_BBCs_pipT[NSYSTEMS];

  // sigma modification function
  TF1* fsig = new TF1("fsig", "TMath::Power((1 + TMath::Exp(-1. *[0] * x)), 2) / 4.0", 0, 1);
  fsig->SetLineColor(kBlue);
  fsig->SetParameter(0, beta);

  double pi0_pT[NX] = {0};
  for (int ix = 0; ix < NX; ix++)
    pi0_pT[ix] = 0.6 * 0.5 * x[ix] * 200;


  // for testing linear modification due to CNM effects
  double linCentMod = 0.20;
  double linPeriphMod = 0.60;
  TGraph *graa_NcollABBCscMod_linMod[NSYSTEMS][NCENT];
  TGraph *grcp_NcollABBCscMod_linMod[NSYSTEMS][NCENT - 1];


  //=====================================================//
  // GET NTUPLE(S) FROM FILE AND CALCULATE YIELD
  //=====================================================//
  cout << endl;
  cout << "--> Reading Ntuples from files ..." << endl;

  // loop over each collision system
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cout << "------------------------------" << endl;
    cout << "          " << collSystem[isys] << endl;
    cout << "------------------------------" << endl;

    // print values
    cout << "  NBD : mu=" << NBD_par[isys][0] << " k=" << NBD_par[isys][1] << endl;
    cout << "  Ntuple name: " << ntpName[isys] << endl;

    ftrigeff->SetParameters(eff_par[isys][0], eff_par[isys][1]);

    // loop over the glauber output for each x in a given collision system
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
      // These represent the "yield" as a function of BBC south charge
      for (int ibx = 1; ibx <= hNcoll_NcollA[isys][ix]->GetNbinsX(); ibx++)
      {
        for (int iby = 1; iby <= hNcoll_NcollA[isys][ix]->GetNbinsY(); iby++)
        {
          // In all cases, the x and y axes have the same binning and range
          // We can therefore do all permutations of modifying the
          // nucleon Ncoll and/or event Ncoll at the same time
          //   ibx: Nucleon Ncoll (containing high-x quark)
          //   iby: Event Ncoll

          // Unmodified case:
          double w = hNcoll_NcollA[isys][ix]->GetBinContent(ibx, iby);
          // Modified case:
          double wMod = hNcollMod_NcollModA[isys][ix]->GetBinContent(ibx, iby);
          // Modified nucleon Ncoll:
          double wModA = hNcoll_NcollModA[isys][ix]->GetBinContent(ibx, iby);
          // Modified collision Ncoll:
          double wModBBCsc = hNcollMod_NcollA[isys][ix]->GetBinContent(ibx, iby);

          // make sure that we have entries in at least one of the cases
          if (!(w > 0 || wMod > 0 || wModA > 0 || wModBBCsc > 0)) continue;

          double NcollA = hNcoll_NcollA[isys][ix]->GetXaxis()->GetBinCenter(ibx);
          double Ncoll = hNcoll_NcollA[isys][ix]->GetYaxis()->GetBinCenter(iby);

          // Loop over BBC south charge (our centrality estimator)
          for (int j = 1; j <= hBBCscNcollA[isys][ix]->GetNbinsX(); j++)
          {
            double c = hBBCscNcollA[isys][ix]->GetXaxis()->GetBinCenter(j);
            double NBD = evalNBD(c,
                                 Ncoll * NBD_par[isys][0],
                                 Ncoll * NBD_par[isys][1]);
            double eff = ftrigeff->Eval(c);

            if (Ncoll > 0)
            {
              hBBCsc[isys][ix]->Fill(c, w * NBD * eff);

              hBBCscNcollA[isys][ix]->Fill(c, NcollA * w * NBD * eff);
              hBBCscModNcollModA[isys][ix]->Fill(c, NcollA * wMod * NBD * eff);
              hBBCscNcollModA[isys][ix]->Fill(c, NcollA * wModA * NBD * eff);
              hBBCscModNcollA[isys][ix]->Fill(c, NcollA * wModBBCsc * NBD * eff);

              hBBCscNcollA_noTrig[isys][ix]->Fill(c, NcollA * w * NBD);
              hBBCscModNcollModA_noTrig[isys][ix]->Fill(c, NcollA * wMod * NBD);
              hBBCscNcollModA_noTrig[isys][ix]->Fill(c, NcollA * wModA * NBD);
              hBBCscModNcollA_noTrig[isys][ix]->Fill(c, NcollA * wModBBCsc * NBD);

            }
          } // j
        } // iby
      } // ibx

      //-----------------------------------------------------------------
      // From here on out, we'll consider three cases:
      // 1) NcollModABBCscMod:
      //        Use the modified projectile nucleon Ncoll for weighting
      //        Use the modified BBC south charge
      // 1) NcollModABBCsc:
      //        Use the modified projectile nucleon Ncoll for weighting
      //        Use the unmodified BBC south charge
      // 1) NcollABBCscMod:
      //        Use the unmodified projectile nucleon Ncoll for weighting
      //        Use the modified BBC south charge
      //-----------------------------------------------------------------


      //-- Deal with charge distributions --//
      mean_BBCs_NcollModABBCscMod[isys][ix] =
        hBBCscModNcollModA[isys][ix]->GetMean();
      mean_BBCs_NcollModABBCsc[isys][ix] =
        hBBCscNcollModA[isys][ix]->GetMean();
      mean_BBCs_NcollABBCscMod[isys][ix] =
        hBBCscModNcollA[isys][ix]->GetMean();

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
        hNcoll_cent[isys][ix][icent] =
          (TH1D*) hNcoll_BBCsc[isys][ix]->ProjectionY(
            Form("hNcoll_cent_%i_%i_%i", isys, ix, icent),
            blim[icent][0], blim[icent][1]);

        hNcollMod_cent[isys][ix][icent] =
          (TH1D*) hNcollMod_BBCscMod[isys][ix]->ProjectionY(
            Form("hNcollMod_cent_%i_%i_%i", isys, ix, icent),
            blim[icent][0], blim[icent][1]);

        // get "jet" yields
        yieldA[icent] =
          hBBCscNcollA[isys][ix]->Integral(blim[icent][0], blim[icent][1]);
        yieldA_NcollModABBCscMod[icent] =
          hBBCscModNcollModA[isys][ix]->Integral(blim[icent][0], blim[icent][1]);
        yieldA_NcollModABBCsc[icent] =
          hBBCscNcollModA[isys][ix]->Integral(blim[icent][0], blim[icent][1]);
        yieldA_NcollABBCscMod[icent] =
          hBBCscModNcollA[isys][ix]->Integral(blim[icent][0], blim[icent][1]);

        // Calculate the bias factor for each centrality bin
        bias_NcollModABBCscMod[isys][ix][icent] =
          yieldA_NcollModABBCscMod[icent] / yieldA[icent];
        bias_NcollModABBCsc[isys][ix][icent] =
          yieldA_NcollModABBCsc[icent] / yieldA[icent];
        bias_NcollABBCscMod[isys][ix][icent] =
          yieldA_NcollABBCscMod[icent] / yieldA[icent];

        // print
        // cout << " cent: " << centl[isys][icent] << " - " << centl[isys][icent + 1] << endl;
        // cout << "    <Ncoll>   : " << hNcoll_cent[isys][ix][icent]->GetMean() << endl;
        // cout << "    <NcollMod>: " << hNcollMod_cent[isys][ix][icent]->GetMean() << endl;
        // cout << "    Bias (both): " << bias_NcollModABBCscMod[isys][ix][icent] << endl;
        // cout << "    Bias (Mod A): " << bias_NcollModABBCsc[isys][ix][icent] << endl;
        // cout << "    Bias (Mod B): " << bias_NcollABBCscMod[isys][ix][icent] << endl;

      } //icent

      // Calculate Rcp
      for (int icent = 0; icent < NCENT - 1; icent++)
      {
        rcp_NcollModABBCscMod[isys][icent][ix] = bias_NcollModABBCscMod[isys][ix][icent];
        rcp_NcollModABBCscMod[isys][icent][ix] /= bias_NcollModABBCscMod[isys][ix][NCENT - 1];

        rcp_NcollModABBCsc[isys][icent][ix] = bias_NcollModABBCsc[isys][ix][icent];
        rcp_NcollModABBCsc[isys][icent][ix] /= bias_NcollModABBCsc[isys][ix][NCENT - 1];

        rcp_NcollABBCscMod[isys][icent][ix] = bias_NcollABBCscMod[isys][ix][icent];
        rcp_NcollABBCscMod[isys][icent][ix] /= bias_NcollABBCscMod[isys][ix][NCENT - 1];
      }

      //-- Calculate the 0-100% modification and RAA --//
      double yieldA_MB = hBBCscNcollA_noTrig[isys][ix]->Integral();
      double yieldA_NcollModABBCscMod_MB = hBBCscModNcollModA_noTrig[isys][ix]->Integral();
      double yieldA_NcollModABBCsc_MB = hBBCscNcollModA_noTrig[isys][ix]->Integral();
      double yieldA_NcollABBCscMod_MB = hBBCscModNcollA_noTrig[isys][ix]->Integral();

      bias_NcollModABBCscMod_MB[isys][ix] = yieldA_NcollModABBCscMod_MB / yieldA_MB;
      bias_NcollModABBCsc_MB[isys][ix] = yieldA_NcollModABBCsc_MB / yieldA_MB;
      bias_NcollABBCscMod_MB[isys][ix] = yieldA_NcollABBCscMod_MB / yieldA_MB;

      // Scale by the MB modification
      //  argument is that this modification will be present in pp, and
      //  therefore should cancel in MB RAA. This scaling enforces that
      //  assumption
      for (int icent = 0; icent < NCENT; icent++)
      {
        raa_NcollModABBCscMod[isys][icent][ix] = bias_NcollModABBCscMod[isys][ix][icent];
        raa_NcollModABBCscMod[isys][icent][ix] /= bias_NcollModABBCscMod_MB[isys][ix];

        raa_NcollModABBCsc[isys][icent][ix] = bias_NcollModABBCsc[isys][ix][icent];
        raa_NcollModABBCsc[isys][icent][ix] /= bias_NcollModABBCsc_MB[isys][ix];

        raa_NcollABBCscMod[isys][icent][ix] = bias_NcollABBCscMod[isys][ix][icent];
        raa_NcollABBCscMod[isys][icent][ix] /= bias_NcollABBCscMod_MB[isys][ix];
      }


      //-- Make the centrality selected charge distribution for later plotting
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


    //-- Define MB RAA graphs
    graa_NcollModABBCscMod_MB[isys] = new TGraph(NX,
        x,
        bias_NcollModABBCscMod_MB[isys]);
    graa_NcollModABBCscMod_MB[isys]->SetName(
      Form("graa_NcollModABBCscMod_MB_%i", isys));
    graa_NcollModABBCscMod_MB[isys]->SetTitle(";x_{p};R_{AA}");
    graa_NcollModABBCscMod_MB[isys]->SetLineStyle(1);
    graa_NcollModABBCscMod_MB[isys]->SetLineWidth(lineWidth);
    graa_NcollModABBCscMod_MB[isys]->SetLineColor(colors[isys]);
    graa_NcollModABBCscMod_MB[isys]->SetLineStyle(lstyle[isys]);

    graa_NcollModABBCsc_MB[isys] = new TGraph(NX,
        x,
        bias_NcollModABBCsc_MB[isys]);
    graa_NcollModABBCsc_MB[isys]->SetName(
      Form("graa_NcollModABBCsc_MB_%i", isys));
    graa_NcollModABBCsc_MB[isys]->SetTitle(";x_{p};R_{AA}");
    graa_NcollModABBCsc_MB[isys]->SetLineStyle(1);
    graa_NcollModABBCsc_MB[isys]->SetLineWidth(lineWidth);
    graa_NcollModABBCsc_MB[isys]->SetLineColor(colors[isys]);
    graa_NcollModABBCsc_MB[isys]->SetLineStyle(4);

    graa_NcollABBCscMod_MB[isys] = new TGraph(NX,
        x,
        bias_NcollABBCscMod_MB[isys]);
    graa_NcollABBCscMod_MB[isys]->SetName(
      Form("graa_NcollABBCscMod_MB_%i", isys));
    graa_NcollABBCscMod_MB[isys]->SetTitle(";x_{p};R_{AA}");
    graa_NcollABBCscMod_MB[isys]->SetLineStyle(1);
    graa_NcollABBCscMod_MB[isys]->SetLineWidth(lineWidth);
    graa_NcollABBCscMod_MB[isys]->SetLineColor(colors[isys]);
    graa_NcollABBCscMod_MB[isys]->SetLineStyle(7);




    //-- Define RAA centrality graphs
    for (int icent = 0; icent < NCENT; icent++)
    {
      graa_NcollModABBCscMod[isys][icent] = new TGraph(NX,
          x,
          raa_NcollModABBCscMod[isys][icent]);
      graa_NcollModABBCscMod[isys][icent]->SetName(
        Form("graa_NcollModABBCscMod_%i_%i", isys, icent));
      graa_NcollModABBCscMod[isys][icent]->SetTitle(";x_{p};R_{AA}");
      graa_NcollModABBCscMod[isys][icent]->SetLineStyle(1);
      graa_NcollModABBCscMod[isys][icent]->SetLineWidth(lineWidth);
      graa_NcollModABBCscMod[isys][icent]->SetLineColor(colors[isys]);
      graa_NcollModABBCscMod[isys][icent]->SetLineStyle(lstyle[isys]);

      graa_NcollModABBCsc[isys][icent] = new TGraph(NX,
          x,
          raa_NcollModABBCsc[isys][icent]);
      graa_NcollModABBCsc[isys][icent]->SetName(
        Form("graa_NcollModABBCsc_%i_%i", isys, icent));
      graa_NcollModABBCsc[isys][icent]->SetTitle(";x_{p};R_{AA}");
      graa_NcollModABBCsc[isys][icent]->SetLineStyle(1);
      graa_NcollModABBCsc[isys][icent]->SetLineWidth(lineWidth);
      graa_NcollModABBCsc[isys][icent]->SetLineColor(colors[isys]);
      graa_NcollModABBCsc[isys][icent]->SetLineStyle(lstyle[isys]);

      graa_NcollABBCscMod[isys][icent] = new TGraph(NX,
          x,
          raa_NcollABBCscMod[isys][icent]);
      graa_NcollABBCscMod[isys][icent]->SetName(
        Form("graa_NcollABBCscMod_%i_%i", isys, icent));
      graa_NcollABBCscMod[isys][icent]->SetTitle(";x_{p};R_{AA}");
      graa_NcollABBCscMod[isys][icent]->SetLineStyle(1);
      graa_NcollABBCscMod[isys][icent]->SetLineWidth(lineWidth);
      graa_NcollABBCscMod[isys][icent]->SetLineColor(colors[isys]);
      graa_NcollABBCscMod[isys][icent]->SetLineStyle(lstyle[isys]);


    }

    //-- Define Rcp graphs
    for (int icent = 0; icent < NCENT - 1; icent++)
    {
      grcp_NcollModABBCscMod[isys][icent] = new TGraph(NX,
          x,
          rcp_NcollModABBCscMod[isys][icent]);
      grcp_NcollModABBCscMod[isys][icent]->SetName(
        Form("grcp_NcollModABBCscMod_%i_%i", isys, icent));
      grcp_NcollModABBCscMod[isys][icent]->SetTitle(
        ";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
      grcp_NcollModABBCscMod[isys][icent]->SetLineStyle(1);
      grcp_NcollModABBCscMod[isys][icent]->SetLineWidth(lineWidth);
      grcp_NcollModABBCscMod[isys][icent]->SetLineColor(colors[isys]);
      grcp_NcollModABBCscMod[isys][icent]->SetLineStyle(lstyle[isys]);

      grcp_NcollModABBCsc[isys][icent] = new TGraph(NX,
          x,
          rcp_NcollModABBCsc[isys][icent]);
      grcp_NcollModABBCsc[isys][icent]->SetName(
        Form("grcp_NcollModABBCsc_%i_%i", isys, icent));
      grcp_NcollModABBCsc[isys][icent]->SetTitle(
        ";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
      grcp_NcollModABBCsc[isys][icent]->SetLineStyle(1);
      grcp_NcollModABBCsc[isys][icent]->SetLineWidth(lineWidth);
      grcp_NcollModABBCsc[isys][icent]->SetLineColor(colors[isys]);
      grcp_NcollModABBCsc[isys][icent]->SetLineStyle(lstyle[isys]);

      grcp_NcollABBCscMod[isys][icent] = new TGraph(NX,
          x,
          rcp_NcollABBCscMod[isys][icent]);
      grcp_NcollABBCscMod[isys][icent]->SetName(
        Form("grcp_NcollABBCscMod_%i_%i", isys, icent));
      grcp_NcollABBCscMod[isys][icent]->SetTitle(
        ";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
      grcp_NcollABBCscMod[isys][icent]->SetLineStyle(1);
      grcp_NcollABBCscMod[isys][icent]->SetLineWidth(lineWidth);
      grcp_NcollABBCscMod[isys][icent]->SetLineColor(colors[isys]);
      grcp_NcollABBCscMod[isys][icent]->SetLineStyle(lstyle[isys]);

      // vs pion pT
      grcp_NcollModABBCscMod_pipT[isys][icent] = new TGraph(NX,
          pi0_pT,
          rcp_NcollModABBCscMod[isys][icent]);
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetName(
        Form("grcp_NcollModABBCscMod_pipT_%i_%i", isys, icent));
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetTitle(";#pi p_{T};R_{CP}");
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetLineStyle(1);
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetLineWidth(lineWidth);
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetLineColor(colors[isys]);
      grcp_NcollModABBCscMod_pipT[isys][icent]->SetLineStyle(lstyle[isys]);

    }

    //-- scale the <BBCsc> by the unmodified value (x=0)
    double denom = mean_BBCs_NcollModABBCscMod[isys][0];
    for (int ix = NX - 1; ix >= 0; ix--)
    {
      mean_BBCs_NcollModABBCscMod[isys][ix] /= mean_BBCs_NcollModABBCscMod[isys][0];
      mean_BBCs_NcollModABBCsc[isys][ix] /= mean_BBCs_NcollModABBCsc[isys][0];
      mean_BBCs_NcollABBCscMod[isys][ix] /= mean_BBCs_NcollABBCscMod[isys][0];
    }

    gmean_BBCs_NcollModABBCscMod[isys] = new TGraph(NX, x, mean_BBCs_NcollModABBCscMod[isys]);
    gmean_BBCs_NcollModABBCscMod[isys]->SetName(Form("gmean_BBCs_NcollModABBCscMod_%i", isys));
    gmean_BBCs_NcollModABBCscMod[isys]->SetTitle(Form(";x;<Q_{BBC, Au}> / <Q_{BBC,Au}(x=%.2f)", x[1]));
    gmean_BBCs_NcollModABBCscMod[isys]->SetLineWidth(lineWidth);
    gmean_BBCs_NcollModABBCscMod[isys]->SetLineColor(colors[isys]);
    gmean_BBCs_NcollModABBCscMod[isys]->SetLineStyle(lstyle[isys]);

    gmean_BBCs_NcollModABBCsc[isys] = new TGraph(NX, x, mean_BBCs_NcollModABBCsc[isys]);
    gmean_BBCs_NcollModABBCsc[isys]->SetName(Form("gmean_BBCs_NcollModABBCsc_%i", isys));
    gmean_BBCs_NcollModABBCsc[isys]->SetTitle(Form(";x;<Q_{BBC, Au}> / <Q_{BBC,Au}(x=%.2f)", x[1]));
    gmean_BBCs_NcollModABBCsc[isys]->SetLineWidth(lineWidth);
    gmean_BBCs_NcollModABBCsc[isys]->SetLineColor(colors[isys]);
    gmean_BBCs_NcollModABBCsc[isys]->SetLineStyle(lstyle[isys]);

    gmean_BBCs_NcollABBCscMod[isys] = new TGraph(NX, x, mean_BBCs_NcollABBCscMod[isys]);
    gmean_BBCs_NcollABBCscMod[isys]->SetName(Form("gmean_BBCs_NcollABBCscMod_%i", isys));
    gmean_BBCs_NcollABBCscMod[isys]->SetTitle(Form(";x;<Q_{BBC, Au}> / <Q_{BBC,Au}(x=%.2f)", x[1]));
    gmean_BBCs_NcollABBCscMod[isys]->SetLineWidth(lineWidth);
    gmean_BBCs_NcollABBCscMod[isys]->SetLineColor(colors[isys]);
    gmean_BBCs_NcollABBCscMod[isys]->SetLineStyle(lstyle[isys]);

    gmean_BBCs_pipT[isys] = new TGraph(NX,
                                       pi0_pT,
                                       mean_BBCs_NcollModABBCscMod[isys]);
    gmean_BBCs_pipT[isys]->SetName(Form("gmean_BBCs_%i", isys));
    gmean_BBCs_pipT[isys]->SetTitle(
      Form(";#pi p_{T} = 0.6#frac{#sqrt{s}*x}{2};<Q_{BBC, Au}> / <Q_{BBC,Au}(x=%.2f)", x[1]));
    gmean_BBCs_pipT[isys]->SetLineWidth(lineWidth);
    gmean_BBCs_pipT[isys]->SetLineColor(colors[isys]);
    gmean_BBCs_pipT[isys]->SetLineStyle(lstyle[isys]);


    //-- test linear modification due to CNM effects
    graa_NcollABBCscMod_linMod[isys][0] = new TGraph();
    graa_NcollABBCscMod_linMod[isys][0]->SetName(
      Form("graa_NcollABBCscMod_linMod_%i_%i", isys, 0));
    graa_NcollABBCscMod_linMod[isys][0]->SetTitle(
      ";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    graa_NcollABBCscMod_linMod[isys][0]->SetLineStyle(1);
    graa_NcollABBCscMod_linMod[isys][0]->SetLineWidth(lineWidth);
    graa_NcollABBCscMod_linMod[isys][0]->SetLineColor(colors[isys]);
    graa_NcollABBCscMod_linMod[isys][0]->SetLineStyle(7);

    graa_NcollABBCscMod_linMod[isys][NCENT - 1] = new TGraph();
    graa_NcollABBCscMod_linMod[isys][NCENT - 1]->SetName(
      Form("graa_NcollABBCscMod_linMod_%i_%i", isys, NCENT - 1));
    graa_NcollABBCscMod_linMod[isys][NCENT - 1]->SetTitle(
      ";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    graa_NcollABBCscMod_linMod[isys][NCENT - 1]->SetLineStyle(1);
    graa_NcollABBCscMod_linMod[isys][NCENT - 1]->SetLineWidth(lineWidth);
    graa_NcollABBCscMod_linMod[isys][NCENT - 1]->SetLineColor(colors[isys]);
    graa_NcollABBCscMod_linMod[isys][NCENT - 1]->SetLineStyle(7);

    grcp_NcollABBCscMod_linMod[isys][0] = new TGraph();
    grcp_NcollABBCscMod_linMod[isys][0]->SetName(
      Form("grcp_NcollABBCscMod_linMod_%i_%i", isys, 0));
    grcp_NcollABBCscMod_linMod[isys][0]->SetTitle(
      ";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    grcp_NcollABBCscMod_linMod[isys][0]->SetLineStyle(1);
    grcp_NcollABBCscMod_linMod[isys][0]->SetLineWidth(lineWidth);
    grcp_NcollABBCscMod_linMod[isys][0]->SetLineColor(colors[isys]);
    grcp_NcollABBCscMod_linMod[isys][0]->SetLineStyle(7);

    int itr = 0;
    for (int ix = 0; ix < NX; ix++)
    {
      if (x[ix] > 0.4)
      {
        graa_NcollABBCscMod_linMod[isys][0]->SetPoint(itr,
            x[ix],
            linCentMod * raa_NcollABBCscMod[isys][0][ix]);

        graa_NcollABBCscMod_linMod[isys][NCENT - 1]->SetPoint(itr,
            x[ix],
            linPeriphMod * raa_NcollABBCscMod[isys][NCENT - 1][ix]);

        grcp_NcollABBCscMod_linMod[isys][0]->SetPoint(itr,
            x[ix],
            linCentMod / linPeriphMod * rcp_NcollABBCscMod[isys][0][ix]);

        itr++;
      }
    }

  } // isys

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


  // calculate the chi2 between the 0-20% RCP and the calculation
  double chi2 = 0;
  double NDF = NJET - 1; // data points - free parameters (beta)
  for (int i = 0; i < NJET; i++)
  {
    double m = grcp_NcollABBCscMod[0][0]->Eval(x_jet[i]);
    double tmp = TMath::Power(Rcp_jet[0][i] - m, 2);
    tmp /= TMath::Power(Rcp_jet_A[0][i], 2);

    chi2 += tmp;
  }
  cout << "-------------" << endl;
  cout << " chi2/ndf = " << chi2 / NDF << endl;
  cout << "-------------" << endl;


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
    double mod_cent = graa_NcollModABBCscMod[1][0]->Eval(xp_JdA[j]);
    double mod_periph = graa_NcollModABBCscMod[1][NCENT - 1]->Eval(xp_JdA[j]);

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



  //=====================================================//
  // PLOT OBJECTS
  //=====================================================//
  cout << endl;
  cout << "--> Plotting ..." << endl;

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);

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

    hBBCscNcollA[isys][xplot]->SetLineWidth(lineWidth);
    hBBCscModNcollModA[isys][xplot]->SetLineWidth(lineWidth);

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

    hNcoll_MB[isys][xplot]->SetLineWidth(lineWidth);
    hNcollMod_MB[isys][xplot]->SetLineWidth(lineWidth);

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

      hNcoll_cent[isys][xplot][icent]->SetLineWidth(lineWidth);
      hNcollMod_cent[isys][xplot][icent]->SetLineWidth(lineWidth);

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

  TLegend *legmod[NSYSTEMS];
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {

    legmod[isys] = new TLegend(0.2, 0.7, 0.98, 0.98);
    legmod[isys]->SetFillStyle(0);
    legmod[isys]->SetBorderSize(0);
    legmod[isys]->SetTextSize(0.08);
    legmod[isys]->SetNColumns(3);
    legmod[isys]->AddEntry(graa_NcollModABBCscMod_MB[isys],
                           "N_{coll}^{Mod} & BBCsc^{Mod}",
                           "L");
    legmod[isys]->AddEntry(graa_NcollModABBCsc_MB[isys],
                           "N_{coll}^{Mod} & BBCsc",
                           "L");
    legmod[isys]->AddEntry(graa_NcollABBCscMod_MB[isys],
                           "N_{coll} & BBCsc^{Mod}",
                           "L");
  }


  //=====================================================//
  // PLOT
  //=====================================================//

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
  hBBCscNcollA[0][0]->SetLineWidth(lineWidth);
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
      legrcp_paper[icent]->AddEntry(grcp_NcollModABBCscMod[isys][icent],
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
    legraa_paper->AddEntry(graa_NcollModABBCscMod[isys][0],
                           lsystem[isys], "L");
  }
  legraa_paper->AddEntry(gRdAu_jet[0], "PHENIX d+Au", "P");
  legraa_paper->AddEntry((TObject*)0, " arXiv:1509.04657", "");


  TLegend *legjda_paper = new TLegend(0.4, 0.6, 0.95, 0.92);
  legjda_paper->SetFillStyle(0);
  legjda_paper->SetBorderSize(0);
  legjda_paper->SetTextSize(0.04);
  legjda_paper->AddEntry(graa_NcollABBCscMod[1][0], lsystem[1], "L");
  legjda_paper->AddEntry(graa_NcollABBCscMod_linMod[1][0],
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
    legQ_paper->AddEntry(gmean_BBCs_NcollABBCscMod[isys], lsystem[isys], "L");



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
      grcp_NcollABBCscMod[isys][icent]->Draw("C");
    // grcp_NcollModABBCscMod[isys][icent]->Draw("C");


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
      graa_NcollABBCscMod[isys][icent]->Draw("C");
    // graa_NcollModABBCscMod[isys][icent]->Draw("C");

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

  graa_NcollABBCscMod[1][0]->Draw("C");
  graa_NcollABBCscMod_linMod[1][0]->Draw("C");
  // graa_NcollModABBCscMod[1][0]->Draw("C");

  l1.DrawLine(0, 1, 1.0, 1);
  label.SetTextSize(0.06 * fracArea / (yh - yl));
  label.DrawLatex(0.3, 0.1, "0-20%");
  // label.DrawLatex(0.9, 0.95 - topMargin, "(a)");
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

  // graa_NcollModABBCscMod[1][NCENT - 1]->Draw("C");
  graa_NcollABBCscMod[1][NCENT - 1]->Draw("C");
  graa_NcollABBCscMod_linMod[1][NCENT - 1]->Draw("C");

  l1.DrawLine(0, 1, 1.0, 1);
  label.SetTextSize(0.06 * fracArea / (yh - yl));
  label.DrawLatex(0.3, 0.1, "60-88%");
  // label.DrawLatex(0.9, 0.95, "(b)");
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


  // grcp_NcollModABBCscMod[1][0]->Draw("same");
  grcp_NcollABBCscMod[1][0]->Draw("same");
  grcp_NcollABBCscMod_linMod[1][0]->Draw("same");

  l1.DrawLine(0, 1, 1.0, 1);
  label.SetTextSize(0.06 * fracArea / (yh - yl));
  label.DrawLatex(0.3, 0.1 + bottomMargin, "(0-20%)/(60-88%)");
  // label.DrawLatex(0.9, 0.95, "(c)");
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
    gmean_BBCs_NcollABBCscMod[isys]->Draw("C");

  l1.DrawLine(0, 1, 1.0, 1);
  legQ_paper->Draw("same");


//=====================================================//
// SAVE
//=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving RCP to " << outFile << endl;
  }
  if (printPlots)
  {
    cout << endl;
    cout << "--> Printing plots to pdf" << endl;

    cyield->Print("yield_systems.pdf");
    cncoll->Print("Ncoll_systems.pdf");
    cncollmb->Print("MBNcoll_systems.pdf");
    cyieldpau->Print("yield_pAu.pdf");
    crcppipt->Print("Rcp_pipT_systems.pdf");
    cbbc->Print("BBCQ_systems.pdf");

    crcp_paper->Print("Rcp_paper.pdf");
    craa_paper->Print("RAA_paper.pdf");
    csig_paper->Print("sigmod.pdf");
    cq_paper->Print("Q_paper.pdf");
    cjda_paper->Print("JdA_paper.pdf");
    // cjda_rat_paper->Print("JdA_rat_paper.pdf");
  }

//   //=====================================================//
//   // DO SOME PRINTING
//   //=====================================================//
//   cout << endl;
//   cout << "--> Printing ... " << endl;

//   cout << setprecision(2) << fixed;
//   cout << endl;
//   cout << " ncoll MB " << endl;
//   cout << "x";
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//     cout << " & " << collSystem[isys];
//   cout << " \\\\" << endl;
//   cout << "0.00";
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//     cout << " & " << hNcoll_MB[isys][0]->GetMean();
//   cout << "\\\\" << endl;
//   for (int ix = 0; ix < NX; ix++)
//   {
//     cout << x[ix];
//     for (int isys = 0; isys < NSYSTEMS; isys++)
//       cout << " & " << hNcollMod_MB[isys][ix]->GetMean();
//     cout << "\\\\" << endl;
//   }

//   cout << endl;
//   cout << "-- ncoll ==" << endl;
//   cout << "x sig(x)";
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//     cout << " " << collSystem[isys];
//   cout << endl;
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//   {
//     cout << " & 0-100";
//     for (int icent = 0; icent < NCENT; icent++)
//       cout << " & " << centl[isys][icent] << "-" << centl[isys][icent + 1];
//   }
//   cout << " \\\\" << endl;
//   cout << "PHENIX & 42.00";
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//   {
//     cout << " & " << Ncoll_PHENIX[isys][NCENT];
//     for (int icent = 0; icent < NCENT; icent++)
//       cout << " & " << Ncoll_PHENIX[isys][icent];
//   }
//   cout << "0.00 & 42.00";
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//   {
//     cout << " & " << hNcoll_MB[isys][0]->GetMean();
//     for (int icent = 0; icent < NCENT; icent++)
//       cout << " & " << hNcoll_cent[isys][0][icent]->GetMean();
//   }
//   cout << " \\\\" << endl;
//   for (int ix = 0; ix < NX; ix++)
//   {
//     cout << x[ix];
//     cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
//     for (int isys = 0; isys < NSYSTEMS; isys++)
//     {
//       cout << " & " << hNcollMod_MB[isys][ix]->GetMean();
//       for (int icent = 0; icent < NCENT; icent++)
//         cout << " & " << hNcollMod_cent[isys][ix][icent]->GetMean();
//     }
//     cout << " \\\\" << endl;
//   }

//   cout << endl;
//   cout << " Bias factors (both mod)" << endl;
//   cout << "x & sig(x)";
//   for (int icent = 0; icent < NCENT; icent++)
//   {
//     cout << " & " << centl[1][icent] << "-" << centl[1][icent + 1];
//   }
//   cout << " & 0-100";
//   cout << "\\\\" << endl;
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//   {
//     cout << " & & \\multicolumn{" << NCENT + 1 << "}{|c}{"
//          << collSystem[isys] << "} \\\\" << endl;
//     for (int ix = 0; ix < NX; ix++)
//     {
//       cout << x[ix];
//       cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
//       for (int icent = 0; icent < NCENT; icent++)
//         cout << " & " << bias_NcollModABBCscMod[isys][ix][icent];

//       cout << " & " << bias_NcollModABBCscMod_MB[isys][ix];
//       cout << "\\\\" << endl;
//     }
//   }

//   cout << endl;
//   cout << " Bias factors (mod A)" << endl;
//   cout << "x";
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//     cout << " & " << collSystem[isys];
//   cout << "\\\\" << endl;
// // for (int isys = 0; isys < NSYSTEMS; isys++)
// // {
// //   for (int icent = 0; icent < NCENT; icent++)
// //     cout << " & " << centl[isys][icent] << "-" << centl[isys][icent +1];
// // }
//   cout << " \\\\" << endl;
//   for (int ix = 0; ix < NX; ix++)
//   {
//     cout << x[ix];
//     // cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
//     for (int isys = 0; isys < NSYSTEMS; isys++)
//     {
//       // for (int icent = 0; icent < NCENT; icent++)
//       cout << " & " << bias_NcollModABBCsc[isys][ix][0];
//     }
//     cout << " | ";
//     for (int isys = 0; isys < NSYSTEMS; isys++)
//     {
//       // for (int icent = 0; icent < NCENT; icent++)
//       cout << " & " << bias_NcollModABBCsc[isys][ix][NCENT - 1];
//     }
//     cout << " \\\\" << endl;
//   }

//   cout << endl;
//   cout << " Bias factors (mod BBCsc)" << endl;
//   cout << "x";
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//     cout << " & " << collSystem[isys];
//   cout << "\\\\" << endl;
// // for (int isys = 0; isys < NSYSTEMS; isys++)
// // {
// //   for (int icent = 0; icent < NCENT; icent++)
// //     cout << " & " << centl[isys][icent] << "-" << centl[isys][icent + 1];
// // }
//   cout << " \\\\" << endl;
//   for (int ix = 0; ix < NX; ix++)
//   {
//     cout << x[ix];
//     // cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
//     for (int isys = 0; isys < NSYSTEMS; isys++)
//     {
//       // for (int icent = 0; icent < NCENT; icent++)
//       cout << " & " << bias_NcollABBCscMod[isys][ix][0];
//     }
//     cout << " | ";
//     for (int isys = 0; isys < NSYSTEMS; isys++)
//     {
//       // for (int icent = 0; icent < NCENT; icent++)
//       cout << " & " << bias_NcollABBCscMod[isys][ix][NCENT - 1];
//     }
//     cout << " \\\\" << endl;
//   }


//   cout << endl;
//   cout << " Rcp cent " << endl;
//   cout << "x ";
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//     cout << " & " << collSystem[isys];
//   cout << "\\\\" << endl;
//   for (int isys = 0; isys < NSYSTEMS; isys++)
//   {
//     for (int icent = 0; icent < NCENT - 1; icent++)
//       cout << " & " << centl[isys][icent] << "-" << centl[isys][icent + 1];
//   }
//   cout << " \\\\" << endl;
//   for (int ix = 0; ix < NX; ix++)
//   {
//     cout << x[ix];
//     // cout << " & " << 42 * 0.25 * TMath::Power(1 + TMath::Exp(-0.75 * x[ix]), 2);
//     for (int isys = 0; isys < NSYSTEMS; isys++)
//     {
//       for (int icent = 0; icent < NCENT - 1; icent++)
//         cout << " & " << rcp_NcollModABBCscMod[isys][icent][ix];
//     }
//     cout << " \\\\" << endl;
//   }









}