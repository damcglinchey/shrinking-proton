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
  // NOTE only central and peripheral values real, others set ot 1
  // to make things easier
  // double linCentMod = 0.20;
  // double linPeriphMod = 0.60;
  double linMod[NCENT] = {0.20, 1.0, 1.0, 0.60};
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
    for (int icent = 0; icent < NCENT; icent++)
    {
      graa_NcollABBCscMod_linMod[isys][icent] = new TGraph();
      graa_NcollABBCscMod_linMod[isys][icent]->SetName(
        Form("graa_NcollABBCscMod_linMod_%i_%i", isys, icent));
      graa_NcollABBCscMod_linMod[isys][icent]->SetTitle(
        ";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
      graa_NcollABBCscMod_linMod[isys][icent]->SetLineStyle(1);
      graa_NcollABBCscMod_linMod[isys][icent]->SetLineWidth(lineWidth);
      graa_NcollABBCscMod_linMod[isys][icent]->SetLineColor(colors[isys]);
      graa_NcollABBCscMod_linMod[isys][icent]->SetLineStyle(7);

      if (icent < NCENT - 1)
      {
        grcp_NcollABBCscMod_linMod[isys][icent] = new TGraph();
        grcp_NcollABBCscMod_linMod[isys][icent]->SetName(
          Form("grcp_NcollABBCscMod_linMod_%i_%i", isys, icent));
        grcp_NcollABBCscMod_linMod[isys][icent]->SetTitle(
          ";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
        grcp_NcollABBCscMod_linMod[isys][icent]->SetLineStyle(1);
        grcp_NcollABBCscMod_linMod[isys][icent]->SetLineWidth(lineWidth);
        grcp_NcollABBCscMod_linMod[isys][icent]->SetLineColor(colors[isys]);
        grcp_NcollABBCscMod_linMod[isys][icent]->SetLineStyle(7);
      }

      int itr = 0;
      for (int ix = 0; ix < NX; ix++)
      {
        if (x[ix] > 0.4)
        {
          graa_NcollABBCscMod_linMod[isys][icent]->SetPoint(itr,
              x[ix],
              linMod[icent] * raa_NcollABBCscMod[isys][icent][ix]);

          if (icent < NCENT - 1)
          {
            grcp_NcollABBCscMod_linMod[isys][icent]->SetPoint(itr,
                x[ix],
                linMod[icent] / linMod[NCENT - 1] * rcp_NcollABBCscMod[isys][icent][ix]);
          }
          itr++;
        }
      }
    }// icent

  } // isys


  //=====================================================//
  // SAVE RESULTS TO FILE
  //=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving results to " << outFile << endl;

    TFile *fout = new TFile(outFile, "UPDATE");

    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      // write Rcp
      for (int icent = 0; icent < NCENT - 1; icent++)
      {
        grcp_NcollABBCscMod[isys][icent]->Write(
          Form("gRcp_%s_cent%i",
               collSystem[isys],
               icent));

        grcp_NcollABBCscMod_linMod[isys][icent]->Write(
          Form("gRcp_%s_cent%i_linMod",
               collSystem[isys],
               icent));
      }

      // write RpA
      for (int icent = 0; icent < NCENT; icent++)
      {
        graa_NcollABBCscMod[isys][icent]->Write(
          Form("gRAA_%s_cent%i",
               collSystem[isys],
               icent));

        graa_NcollABBCscMod_linMod[isys][icent]->Write(
          Form("gRAA_%s_cent%i_linMod",
               collSystem[isys],
               icent));
      }

      // BBC charge dist
      gmean_BBCs_NcollABBCscMod[isys]->Write(
        Form("gQ_%s", collSystem[isys]));


      // BBC
      for (int ix = 0; ix < NX; ix++)
      {
        hBBCscModNcollA[isys][ix]->Write(
          Form("yield_vs_BBCQ_%s_x%03.0f",
               collSystem[isys],
               x[ix] * 1000));

        hBBCsc[isys][ix]->Write(
          Form("BBCQ_MB_%s_x%03.0f",
               collSystem[isys],
               x[ix] * 1000));
      }

    }// isys
  } // saveRcp

  //=====================================================//
  // DO SOME PRINTING
  //=====================================================//
  //   cout << endl;
  //   cout << "-- > Printing ... " << endl;


}