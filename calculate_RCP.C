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

  const int NX = 4;             // Number of x values
  double x[] = {0.01, 0.10, 0.30, 0.50};

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
      "rootfiles/glauber_pau_snn42_x03_ntuple_100k.root",
      "rootfiles/glauber_pau_snn42_x05_ntuple_100k.root",
    },
    { // dAu Files
      "rootfiles/glauber_dau_snn42_x001_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_x01_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_x03_ntuple_100k.root",
      "rootfiles/glauber_dau_snn42_x05_ntuple_100k.root",
    },
    { // He3Au Files
      "rootfiles/glauber_he3au_snn42_x001_ntuple_100k.root",
      "rootfiles/glauber_he3au_snn42_x01_ntuple_100k.root",
      "rootfiles/glauber_he3au_snn42_x03_ntuple_100k.root",
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
  double dcent_cent[NSYSTEMS] = {0.20, 0.20, 0.20}; // 0-20%
  double dcent_periph[NSYSTEMS] = {0.24, 0.28, 0.28}; // 60-84% (pAu), 60-88% (dAu, 3HeAu)
  double cent_tot[NSYSTEMS] = {0.84, 0.88, 0.88};
  // double dcent_cent[NSYSTEMS] = {0.05, 0.05, 0.05}; // 0-20%
  // double dcent_periph[NSYSTEMS] = {0.14, 0.18, 0.18}; // 60-84% (pAu), 60-88% (dAu, 3HeAu)
  // double cent_tot[NSYSTEMS] = {0.84, 0.88, 0.88};

  const int NCENT = 9;
  double centl[NSYSTEMS][NCENT + 1] =
  {
    {0, 5, 10, 20, 30, 40, 50, 60, 70, 84}, // pAu
    {0, 5, 10, 20, 30, 40, 50, 60, 70, 88}, // dAu
    {0, 5, 10, 20, 30, 40, 50, 60, 70, 88}, // He3Au
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

  //=====================================================//
  // DECLARE VARIABLES
  //=====================================================//

  // For running over ntuples
  TFile *fin;
  TNtuple *ntp;
  TH2D *hNcoll_NcollMod[NSYSTEMS][NX];

  // Calculated histograms
  TH1D *hBBCs[NSYSTEMS][NX];
  TH1D *hBBCsMod[NSYSTEMS][NX];

  TH1D *hNcoll_MB[NSYSTEMS][NX];
  TH1D *hNcollMod_MB[NSYSTEMS][NX];

  TH2D *hNcollMod_BBCsc[NSYSTEMS][NX];
  TH2D *hNcoll_BBCsc[NSYSTEMS][NX];
  TH2D *hNcollMod_BBCscMod[NSYSTEMS][NX];
  TH2D *hNcoll_BBCscMod[NSYSTEMS][NX];

  TH1D *hBBCscNcoll[NSYSTEMS][NX];
  TH1D* hBBCscModNcoll[NSYSTEMS][NX];
  TH1D *hBBCscNcollMod[NSYSTEMS][NX];
  TH1D *hBBCscModNcollMod[NSYSTEMS][NX];
  for (int i = 0; i < NSYSTEMS; i++)
  {
    for (int j = 0; j < NX; j++)
    {

      hNcollMod_BBCsc[i][j] = new TH2D(Form("hNcollMod_BBCsc_%i_%i", i, j),
                                       ";BBCs charge; N_{coll}^{mod}",
                                       1596, 1, 400,
                                       101, -0.5, 100.5);

      hNcoll_BBCsc[i][j] = new TH2D(Form("hNcoll_BBCsc_%i_%i", i, j),
                                    ";BBCs charge; N_{coll}",
                                    1596, 1, 400,
                                    101, -0.5, 100.5);

      hNcollMod_BBCscMod[i][j] = new TH2D(Form("hNcollMod_BBCscMod_%i_%i", i, j),
                                          ";BBCs charge Mod; N_{coll}^{mod}",
                                          1596, 1, 400,
                                          101, -0.5, 100.5);

      hNcoll_BBCscMod[i][j] = new TH2D(Form("hNcoll_BBCscMod_%i_%i", i, j),
                                       ";BBCs charge Mod; N_{coll}",
                                       1596, 1, 400,
                                       101, -0.5, 100.5);

      // Ncoll weighted BBC charge distributions
      hBBCscNcoll[i][j] = new TH1D(Form("hBBCscNcoll_%i_%i", i, j),
                                   ";N_{coll} #time BBCs charge",
                                   1596, 1, 400);

      hBBCscModNcoll[i][j] = new TH1D(Form("hBBCscModNcoll_%i_%i", i, j),
                                      ";N_{coll} #time BBCs charge Mod",
                                      1596, 1, 400);

      hBBCscNcollMod[i][j] = new TH1D(Form("hBBCscNcollMod_%i_%i", i, j),
                                      ";N_{coll}^{mod} #time BBCs charge",
                                      1596, 1, 400);

      hBBCscModNcollMod[i][j] = new TH1D(Form("hBBCscModNcollMod_%i_%i", i, j),
                                         ";N_{coll}^{mod} #time BBCs charge Mod",
                                         1596, 1, 400);

    }
  }


  TH1D* htmp;

  // NBD
  // TF1 *fNBD = new TF1("fNBD", NBD, 0, 200, 2);

  // trigger efficiencies
  TF1 *ftrigeff = new TF1("ftrigeff", "1.0-TMath::Exp(-pow((x/[0]), [1]))", 0.0, 200.0);


  // for Rcp
  TH1D *hNcoll[NSYSTEMS][NX][2];
  TH1D *hNcollMod[NSYSTEMS][NX][2];
  TH1D *hNcollModBBCscMod[NSYSTEMS][NX][2];
  TH1D *hNcollBBCscMod[NSYSTEMS][NX][2];

  double bias_BBCcsMod[NSYSTEMS][NX][2];
  double bias_NcollMod[NSYSTEMS][NX][2];
  double bias_NcollModBBCcsMod[NSYSTEMS][NX][2];

  double rcp_BBCscMod[NSYSTEMS][NX];
  double rcp_NcollMod[NSYSTEMS][NX];
  double rcp_NcollModBBCscMod[NSYSTEMS][NX];

  TGraph *grcp_NcollMod[NSYSTEMS];
  TGraph *grcp_BBCscMod[NSYSTEMS];
  TGraph *grcp_NcollModBBCscMod[NSYSTEMS];

  // R (x/mid-central)
  double bias_BBCcsMod_cent[NSYSTEMS][NX][NCENT];
  double r_BBCscMod_cent[NSYSTEMS][NX][NCENT];
  TGraph *gr_BBBCscMod_cent[NSYSTEMS][NX];


  // alt
  double altyield_BBCscNcoll[NSYSTEMS][NX][2];
  double altyield_BBCscModNcoll[NSYSTEMS][NX][2];
  double altyield_BBCscModNcollMod[NSYSTEMS][NX][2];

  double altrcp_BBCscNcoll[NSYSTEMS][NX];
  double altrcp_BBCscModNcoll[NSYSTEMS][NX];
  double altrcp_BBCscModNcollMod[NSYSTEMS][NX];

  TGraph *galtrcp_BBCscModNcoll[NSYSTEMS];
  TGraph *galtrcp_BBCscModNcollMod[NSYSTEMS];

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

      // ntp->Draw("Sum$(NcollA):Sum$(NcollModA) >> htmp(100, 0.5, 100.5, 100, 0.5, 100.5)",
      //           "", "goff");
      ntp->Draw("Sum$(NcollA):Sum$(NcollModA) >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
                "", "goff");

      hNcoll_NcollMod[isys][ix] = (TH2D*) gDirectory->FindObject("htmp");
      hNcoll_NcollMod[isys][ix]->SetDirectory(0);
      hNcoll_NcollMod[isys][ix]->SetName(Form("hNcoll_NcollMod_%i_%i", isys, ix));
      hNcoll_NcollMod[isys][ix]->SetTitle(";N_{coll}^{mod};N_{coll}");

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
            // double eff = ftrigeff->Eval(c);
            double eff = 1;

            // NBD=nan for mu/k = 0
            if (Ncoll > 0)
            {
              hNcoll_BBCsc[isys][ix]->Fill(c, Ncoll, w * NBD * eff);
              hNcollMod_BBCsc[isys][ix]->Fill(c, NcollMod, w * NBD * eff);

              hBBCscNcoll[isys][ix]->Fill(c, Ncoll * w * NBD * eff);
              hBBCscNcollMod[isys][ix]->Fill(c, NcollMod * w * NBD * eff);
            }
            if (NcollMod > 0)
            {
              hNcoll_BBCscMod[isys][ix]->Fill(c, Ncoll, w * NBDMod * eff);
              hNcollMod_BBCscMod[isys][ix]->Fill(c, NcollMod, w * NBDMod * eff);

              hBBCscModNcoll[isys][ix]->Fill(c, Ncoll * w * NBDMod * eff);
              hBBCscModNcollMod[isys][ix]->Fill(c, NcollMod * w * NBDMod * eff);
            }
          } // j
        } // iby
      } // ibx


      // Project over unmodified/modified BBCcs (unmodified Ncoll)
      hBBCs[isys][ix] = (TH1D*)
                        hNcoll_BBCsc[isys][ix]->ProjectionX(
                          Form("hBBCs_%i_%i", isys, ix),
                          1,
                          hNcoll_BBCsc[isys][ix]->GetNbinsY());

      hBBCsMod[isys][ix] = (TH1D*)
                           hNcoll_BBCscMod[isys][ix]->ProjectionX(
                             Form("hBBCsMod_%i_%i", isys, ix),
                             1,
                             hNcoll_BBCscMod[isys][ix]->GetNbinsY());




      //-- Calculate RCP --//
      int nq = 2;
      double xq[2];
      double yq[2] = {0};
      xq[0] = dcent_periph[isys] / cent_tot[isys];
      xq[1] = (cent_tot[isys] - dcent_cent[isys]) / cent_tot[isys];
      hBBCs[isys][ix]->GetQuantiles(nq, yq, xq);
      cout << " quantiles: " << yq[0] << " " << yq[1] << endl;

      // get bin limits for each quantile
      int blim[2][2];
      // central
      blim[0][0] = hBBCs[isys][ix]->FindBin(yq[1]);
      blim[0][1] = hBBCs[isys][ix]->GetNbinsX();
      // peripheral
      blim[1][0] = 1;
      blim[1][1] = hBBCs[isys][ix]->FindBin(yq[0]);

      cout << "  blim[0]: [" << blim[0][0] << ", " << blim[0][1] << "]"
           << " -> [" << hBBCs[isys][ix]->GetBinLowEdge(blim[0][0]) << ","
           << hBBCs[isys][ix]->GetBinLowEdge(blim[0][1] + 1) << "]"
           << endl;
      cout << "  blim[1]: [" << blim[1][0] << ", " << blim[1][1] << "]"
           << " -> [" << hBBCs[isys][ix]->GetBinLowEdge(blim[1][0]) << ","
           << hBBCs[isys][ix]->GetBinLowEdge(blim[1][1] + 1) << "]"
           << endl;


      double Nevent[2];
      double NeventMod[2];
      double yield[2] = {0};
      double yield_NcollMod[2] = {0};
      double yield_BBCscMod[2] = {0};
      double yield_NcollModBBCscMod[2] = {0};


      for (int iq = 0; iq < 2; iq++)
      {
        // Calculate the unmodified number of events
        Nevent[iq] = hBBCs[isys][ix]->Integral(blim[iq][0], blim[iq][1]);

        // Calculate the number of modified events
        NeventMod[iq] = hBBCsMod[isys][ix]->Integral(blim[iq][0], blim[iq][1]);

        // Project the unmodified Ncoll (unmodified BBCs charge)
        hNcoll[isys][ix][iq] = (TH1D*) hNcoll_BBCsc[isys][ix]->ProjectionY(
                                 Form("hNcoll_%i_%i_%i", isys, ix, iq),
                                 blim[iq][0], blim[iq][1]);

        // Project the modified Ncoll (unmodified BBCs charge)
        hNcollMod[isys][ix][iq] = (TH1D*) hNcollMod_BBCsc[isys][ix]->ProjectionY(
                                    Form("hNcollMod_%i_%i_%i", isys, ix, iq),
                                    blim[iq][0], blim[iq][1]);

        // Project the unmodified Ncoll (modified BBCs charge)
        hNcollBBCscMod[isys][ix][iq] = (TH1D*) hNcoll_BBCscMod[isys][ix]->ProjectionY(
                                         Form("hNcollBBCscMod_%i_%i_%i", isys, ix, iq),
                                         blim[iq][0], blim[iq][1]);

        // Project the modified Ncoll for modified BBCs charge)
        hNcollModBBCscMod[isys][ix][iq] = (TH1D*) hNcollMod_BBCscMod[isys][ix]->ProjectionY(
                                            Form("hNcollModBBCscMod_%i_%i_%i", isys, ix, iq),
                                            blim[iq][0], blim[iq][1]);

        // Calculate yields (Ncoll weighted)
        for (int ib = 1; ib <= hNcoll[isys][ix][iq]->GetNbinsX(); ib++)
        {
          double ncoll = hNcoll[isys][ix][iq]->GetBinCenter(ib);

          yield[iq]                  += ncoll * hNcoll[isys][ix][iq]->GetBinContent(ib);
          yield_NcollMod[iq]         += ncoll * hNcollMod[isys][ix][iq]->GetBinContent(ib);
          yield_BBCscMod[iq]         += ncoll * hNcollBBCscMod[isys][ix][iq]->GetBinContent(ib);
          yield_NcollModBBCscMod[iq] += ncoll * hNcollModBBCscMod[isys][ix][iq]->GetBinContent(ib);

        }

        // Calculate the bias factor for each centrality bin
        bias_BBCcsMod[isys][ix][iq] = yield_BBCscMod[iq] / NeventMod[iq];
        bias_BBCcsMod[isys][ix][iq] /= yield[iq] / Nevent[iq];

        bias_NcollMod[isys][ix][iq] = yield_NcollMod[iq] / Nevent[iq];
        bias_NcollMod[isys][ix][iq] /= yield[iq] / Nevent[iq];

        bias_NcollModBBCcsMod[isys][ix][iq] = yield_NcollModBBCscMod[iq] / NeventMod[iq];
        bias_NcollModBBCcsMod[isys][ix][iq] /= yield[iq] / Nevent[iq];


        // Alternate calculation of the bias
        altyield_BBCscNcoll[isys][ix][iq] = hBBCscNcoll[isys][ix]->Integral(blim[iq][0], blim[iq][1]);
        altyield_BBCscModNcoll[isys][ix][iq] = hBBCscModNcoll[isys][ix]->Integral(blim[iq][0], blim[iq][1]);
        altyield_BBCscModNcollMod[isys][ix][iq] = hBBCscModNcollMod[isys][ix]->Integral(blim[iq][0], blim[iq][1]);

        cout << " alt[" << isys << "][" << ix << "][" << iq << "]: "
             << hBBCscNcoll[isys][ix]->Integral(blim[iq][0], blim[iq][1]) << " "
             << hBBCscModNcoll[isys][ix]->Integral(blim[iq][0], blim[iq][1]) << " "
             << hBBCscModNcollMod[isys][ix]->Integral(blim[iq][0], blim[iq][1]) << " "
             << endl;


      } // iq

      // alt rcp
      altrcp_BBCscModNcoll[isys][ix] = altyield_BBCscModNcoll[isys][ix][0] / altyield_BBCscNcoll[isys][ix][0];
      altrcp_BBCscModNcoll[isys][ix] /= altyield_BBCscModNcoll[isys][ix][1] / altyield_BBCscNcoll[isys][ix][1];

      altrcp_BBCscModNcollMod[isys][ix] = altyield_BBCscModNcollMod[isys][ix][0] / altyield_BBCscNcoll[isys][ix][0];
      altrcp_BBCscModNcollMod[isys][ix] /= altyield_BBCscModNcollMod[isys][ix][1] / altyield_BBCscNcoll[isys][ix][1];


      // Calculate Rcp
      rcp_BBCscMod[isys][ix] = bias_BBCcsMod[isys][ix][0] / bias_BBCcsMod[isys][ix][1];

      rcp_NcollMod[isys][ix] = bias_NcollMod[isys][ix][0] / bias_NcollMod[isys][ix][1];

      rcp_NcollModBBCscMod[isys][ix] = bias_NcollModBBCcsMod[isys][ix][0] / bias_NcollModBBCcsMod[isys][ix][1];



      cout << " x: " << x[ix] << endl;
      // cout << "   <Ncoll(0-20%)>       : " << hNcoll[isys][ix][0]->GetMean() << endl;
      // cout << "   <Ncoll^{mod}(0-20%)> : " << hNcollMod[isys][ix][0]->GetMean() << endl;
      // cout << "   Nevent(0-20%)        : " << Nevent[0] << endl;
      // cout << "   Nevent^{mod}(0-20%)  : " << NeventMod[0] << endl;
      // cout << "   yield(0-20%)         : " << yield[0] << endl;
      // cout << "   yield^{Nmod}(0-20%)  : " << yield_NcollMod[0] << endl;
      // cout << "   yield^{Bmod}(0-20%)  : " << yield_BBCscMod[0] << endl;
      // cout << "   yield^{NBmod}(0-20%) : " << yield_NcollModBBCscMod[0] << endl;
      // cout << endl;
      // cout << "   <Ncoll(60-88%)>      : " << hNcoll[isys][ix][1]->GetMean() << endl;
      // cout << "   <Ncoll^{mod}(60-88%)>: " << hNcollMod[isys][ix][1]->GetMean() << endl;
      // cout << "   Nevent(60-88%)       : " << Nevent[1] << endl;
      // cout << "   Nevent^{mod}(60-88%) : " << NeventMod[1] << endl;
      // cout << "   yield(60-88%)        : " << yield[1] << endl;
      // cout << "   yield^{Nmod}(60-88%) : " << yield_NcollMod[1] << endl;
      // cout << "   yield^{Bmod}(60-88%) : " << yield_BBCscMod[1] << endl;
      // cout << "   yield^{NBmod}(60-88%): " << yield_NcollModBBCscMod[1] << endl;
      // cout << endl;
      cout << "   Bias(0-20%) (BBCcs Mod) : " << bias_BBCcsMod[isys][ix][0] << endl;
      cout << "   Bias(0-20%) (Ncoll Mod) : " << bias_NcollMod[isys][ix][0] << endl;
      cout << "   Bias(0-20%) (both Mod)  : " << bias_NcollModBBCcsMod[isys][ix][0] << endl;
      cout << "   Bias(60-88%) (BBCcs Mod): " << bias_BBCcsMod[isys][ix][1] << endl;
      cout << "   Bias(60-88%) (Ncoll Mod): " << bias_NcollMod[isys][ix][1] << endl;
      cout << "   Bias(60-88%) (both Mod) : " << bias_NcollModBBCcsMod[isys][ix][1] << endl;
      cout << endl;
      cout << "   Rcp (BBCcs Mod) : " << rcp_BBCscMod[isys][ix] << endl;
      cout << "   Rcp (Ncoll Mod) : " << rcp_NcollMod[isys][ix] << endl;
      cout << "   Rcp (both Mod)  : " << rcp_NcollModBBCscMod[isys][ix] << endl;










      //-- Calculate R for each centrality bin --//

      // // get the BBCcs limits for each centrality bin
      // double xq_c[NCENT - 1];
      // double yq_c[NCENT - 1] = {0};
      // for (int icent = 0; icent < NCENT - 1; icent++)
      //   xq_c[icent] = (centl[isys][NCENT] - centl[isys][icent + 1]) / centl[isys][NCENT];

      // hBBCs[isys][ix]->GetQuantiles(NCENT - 1, yq_c, xq_c);

      // for (int icent = 0; icent < NCENT - 1; icent++)
      // {
      //   cout << "   " << icent << " "
      //        << xq_c[icent] << " "
      //        << yq_c[icent] << endl;
      // }














      //----- old ------//


      // // get the number of events for each centrality bin
      // if (ix == 0)
      // {
      //   double xq_c[NCENT - 1];
      //   double yq_c[NCENT - 1] = {0};
      //   for (int icent = 0; icent < NCENT - 1; icent++)
      //     xq_c[icent] = (centl[isys][NCENT] - centl[isys][icent + 1]) / centl[isys][NCENT];
      //   // xq_c[icent] = (100. - centh[icent]) / 100.;

      //   hBBCs[isys][ix]->GetQuantiles(NCENT - 1, yq_c, xq_c);

      //   for (int icent = 0; icent < NCENT - 1; icent++)
      //   {
      //     BBC_cent[icent] = yq_c[icent];
      //     cout << "   " << icent << " "
      //          << xq_c[icent] << " "
      //          << yq_c[icent] << endl;
      //   }
      // }

      // for (int icent = 0; icent < NCENT; icent++)
      // {
      //   if (icent == 0)
      //     Nevent_cent[isys][ix][icent] = hBBCs[isys][ix]->Integral(hBBCs[isys][ix]->FindBin(BBC_cent[icent]),
      //                                    hBBCs[isys][ix]->GetNbinsX());
      //   else if (icent == NCENT - 1)
      //     Nevent_cent[isys][ix][icent] = hBBCs[isys][ix]->Integral(1,
      //                                    hBBCs[isys][ix]->FindBin(BBC_cent[icent - 1]));
      //   else
      //     Nevent_cent[isys][ix][icent] = hBBCs[isys][ix]->Integral(hBBCs[isys][ix]->FindBin(BBC_cent[icent]),
      //                                    hBBCs[isys][ix]->FindBin(BBC_cent[icent - 1]));

      //   cout << "   " << isys << " " << ix << " " << icent
      //        << " " << Nevent_cent[isys][ix][icent] << endl;
      // }

      // // delete htmp;
      // delete ntp;
      // fin->Close();
      // delete fin;

    } // ix

    grcp_BBCscMod[isys] = new TGraph(NX, x, rcp_BBCscMod[isys]);
    grcp_BBCscMod[isys]->SetName(Form("grcp_BBCscMod_%i", isys));
    grcp_BBCscMod[isys]->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    grcp_BBCscMod[isys]->SetLineStyle(1);
    grcp_BBCscMod[isys]->SetLineWidth(2);
    grcp_BBCscMod[isys]->SetLineColor(colors[isys]);

    // grcp_NcollMod = new TGraph(NX, x, rcp_NcollMod);
    // grcp_NcollMod->SetName("grcp_NcollMod");
    // grcp_NcollMod->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    // grcp_NcollMod->SetLineStyle(3);
    // grcp_NcollMod->SetLineWidth(2);
    // grcp_NcollMod->SetLineColor(kGreen+2);

    // grcp_NcollModBBCscMod = new TGraph(NX, x, rcp_NcollModBBCscMod);
    // grcp_NcollModBBCscMod->SetName("grcp_NcollModBBCscMod");
    // grcp_NcollModBBCscMod->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    // grcp_NcollModBBCscMod->SetLineStyle(2);
    // grcp_NcollModBBCscMod->SetLineWidth(2);
    // grcp_NcollModBBCscMod->SetLineColor(kRed);

    galtrcp_BBCscModNcoll[isys] = new TGraph(NX, x, altrcp_BBCscModNcoll[isys]);
    galtrcp_BBCscModNcoll[isys]->SetName(Form("galtrcp_BBCscModNcoll_%i", isys));
    galtrcp_BBCscModNcoll[isys]->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    galtrcp_BBCscModNcoll[isys]->SetLineStyle(1);
    galtrcp_BBCscModNcoll[isys]->SetLineWidth(2);
    galtrcp_BBCscModNcoll[isys]->SetLineColor(colors[isys]);

    galtrcp_BBCscModNcollMod[isys] = new TGraph(NX, x, altrcp_BBCscModNcollMod[isys]);
    galtrcp_BBCscModNcollMod[isys]->SetName(Form("galtrcp_BBCscModNcollMod_%i", isys));
    galtrcp_BBCscModNcollMod[isys]->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
    galtrcp_BBCscModNcollMod[isys]->SetLineStyle(1);
    galtrcp_BBCscModNcollMod[isys]->SetLineWidth(2);
    galtrcp_BBCscModNcollMod[isys]->SetLineColor(colors[isys]);



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
    cout << " " << collSystem[isys];
  cout << endl;
  cout << "0.00";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << hNcoll_MB[isys][0]->GetMean();
  cout << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    for (int isys = 0; isys < NSYSTEMS; isys++)
      cout << " " << hNcollMod_MB[isys][ix]->GetMean();
    cout << endl;
  }

  cout << endl;
  cout << " ncoll cent ncoll periph " << endl;
  cout << "x sig(x)";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  cout << endl;
  cout << "0.00 42.0";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << hNcoll[isys][0][0]->GetMean();
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << hNcoll[isys][0][1]->GetMean();
  cout << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    cout << " " << 42 * TMath::Exp(-8.0 * x[ix]);
    for (int isys = 0; isys < NSYSTEMS; isys++)
      cout << " " << hNcollBBCscMod[isys][ix][0]->GetMean();
    for (int isys = 0; isys < NSYSTEMS; isys++)
      cout << " " << hNcollBBCscMod[isys][ix][1]->GetMean();
    cout << endl;
  }


  cout << endl;
  cout << " Bias(0-20%) & Bias(60-88%)" << endl;
  cout << "x sig(x)";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  cout << endl;
  cout << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    cout << " " << 42 * TMath::Exp(-8.0 * x[ix]);
    for (int isys = 0; isys < NSYSTEMS; isys++)
      cout << " " << bias_BBCcsMod[isys][ix][0];
    for (int isys = 0; isys < NSYSTEMS; isys++)
      cout << " " << bias_BBCcsMod[isys][ix][1];
    cout << endl;
  }

  cout << endl;
  cout << "-- Ncoll * BBCscMod --" << endl;
  cout << " Bias(0-20%) & Bias(60-88%)" << endl;
  cout << "x sig(x)";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  cout << endl;
  cout << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    cout << setw(7) << 42 * TMath::Exp(-8.0 * x[ix]) << " | ";
    for (int iq = 0; iq < 2; iq++)
    {
      for (int isys = 0; isys < NSYSTEMS; isys++)
        cout << setw(7) << altyield_BBCscModNcoll[isys][ix][iq] / altyield_BBCscNcoll[isys][ix][iq];
      cout << " | ";
    }
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      double rcp = altyield_BBCscModNcoll[isys][ix][0] / altyield_BBCscNcoll[isys][ix][0];
      rcp /= altyield_BBCscModNcoll[isys][ix][1] / altyield_BBCscNcoll[isys][ix][1];
      cout << setw(7) << rcp;
    }
    cout << endl;
  }

  cout << endl;
  cout << "-- NcollMod * BBCscMod --" << endl;
  cout << " Bias(0-20%) & Bias(60-88%)" << endl;
  cout << "x sig(x)";
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  for (int isys = 0; isys < NSYSTEMS; isys++)
    cout << " " << collSystem[isys];
  cout << endl;
  cout << endl;
  for (int ix = 0; ix < NX; ix++)
  {
    cout << x[ix];
    cout << setw(7) << 42 * TMath::Exp(-8.0 * x[ix]) << " | ";
    for (int iq = 0; iq < 2; iq++)
    {
      for (int isys = 0; isys < NSYSTEMS; isys++)
        cout << setw(7) << altyield_BBCscModNcollMod[isys][ix][iq] / altyield_BBCscNcoll[isys][ix][iq];
      cout << " | ";
    }
    for (int isys = 0; isys < NSYSTEMS; isys++)
    {
      double rcp = altyield_BBCscModNcollMod[isys][ix][0] / altyield_BBCscNcoll[isys][ix][0];
      rcp /= altyield_BBCscModNcollMod[isys][ix][1] / altyield_BBCscNcoll[isys][ix][1];
      cout << setw(7) << rcp;
    }
    cout << endl;
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

  TLegend *legncollmb = new TLegend(0.6, 0.2, 0.95, 0.95);
  legncollmb->SetFillStyle(0);
  legncollmb->SetBorderSize(0);
  legncollmb->SetTextSize(0.06);
  legncollmb->AddEntry(hNcoll_MB[0][0], "x = 0.00", "L");
  for (int i = 0; i < NX; i += 2)
  {
    legncollmb->AddEntry(hNcollMod_MB[0][i],
                         Form("x = %.2f", x[i]),
                         "L");
  }


  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);

  TH1F* haxis_rcp = new TH1F("haxis_rcp",
                             Form(";x (x_{jet}=2 * p_{T}^{jet} / #sqrt{s_{NN}});R_{CP}(0-%.0f%% / %.0f-88%%)",
                                  dcent_cent[1] * 100, (0.88 - dcent_periph[1]) * 100),
                             100, 0, 1);
  haxis_rcp->GetYaxis()->SetTitleOffset(1.3);
  haxis_rcp->GetXaxis()->SetTitleOffset(1.3);
  haxis_rcp->SetMinimum(0);
  haxis_rcp->SetMaximum(1.1);

  TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.83);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  for (int isys = 0; isys < NSYSTEMS; isys++)
    leg->AddEntry(grcp_BBCscMod[isys], collSystem[isys], "L");

  TLine l1;
  l1.SetLineStyle(2);

  // TH1F* haxis_rcpcent = new TH1F("haxis_rcpcent",
  //                                ";centrality; R",
  //                                NCENT, 0, NCENT - 1);
  // haxis_rcpcent->SetMinimum(0);
  // haxis_rcpcent->SetMaximum(2.0);
  // // haxis_rcpcent->GetXaxis()->SetNdivisions(NCENT);
  // for (int ibin = 1; ibin <= haxis_rcpcent->GetNbinsX(); ibin++)
  // {
  //   haxis_rcpcent->GetXaxis()->SetBinLabel(ibin,
  //                                          Form("%.0f - %.0f",
  //                                              centl[0][ibin - 1],
  //                                              centl[0][ibin]));
  // }
  // haxis_rcpcent->LabelsOption("d");

  // TLegend *legx = new TLegend(0.6, 0.15, 0.9, 0.5);
  // legx->SetFillStyle(0);
  // legx->SetBorderSize(0);
  // for (int ix = 0; ix < NX; ix++)
  // {
  //   legx->AddEntry(grcp_cent[0][ix], Form("x = %.2f", x[ix]), "L");
  // }

  int xplot = 2;

  TLegend *legNcoll[NSYSTEMS][2];
  const char* centlabel[2] = {"0-20", "60-88"};
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    for (int i = 0; i < 2; i++)
    {
      legNcoll[isys][i] = new TLegend(0.3, 0.7, 0.96, 0.96,
                                      Form("%s @ 200 GeV, x = %.2f, %s%%",
                                           collSystem[isys], x[xplot], centlabel[i]));
      legNcoll[isys][i]->SetFillStyle(0);
      legNcoll[isys][i]->SetBorderSize(0);
      legNcoll[isys][i]->SetTextSize(0.04);
      legNcoll[isys][i]->AddEntry(hNcoll[isys][xplot][i],
                                  Form("Unmodified <N_{coll}> = %.2f",
                                       hNcoll[isys][xplot][i]->GetMean()),
                                  "L"
                                 );
      legNcoll[isys][i]->AddEntry(hNcollBBCscMod[isys][xplot][i],
                                  Form("Mod BBCs chrg <N_{coll}> = %.2f",
                                       hNcollBBCscMod[isys][xplot][i]->GetMean()),
                                  "L"
                                 );
    }
  }

  TLegend *legBBC[NSYSTEMS];
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    legBBC[isys] = new TLegend(0.5, 0.7, 0.96, 0.96,
                               Form("%s @ 200 GeV, x = %.2f",
                                    collSystem[isys], x[xplot]));
    legBBC[isys]->SetFillStyle(0);
    legBBC[isys]->SetBorderSize(0);
    legBBC[isys]->SetTextSize(0.04);
    legBBC[isys]->AddEntry(hBBCs[isys][xplot], "Unmodified", "L");
    legBBC[isys]->AddEntry(hBBCsMod[isys][xplot], "Modified", "L");

  }
  //=====================================================//
  // PLOT
  //=====================================================//

  TCanvas *crcp = new TCanvas("crcp", "RCP", 800, 800);
  crcp->Divide(1, 1);
  crcp->GetPad(1)->SetTopMargin(0.02);
  crcp->GetPad(1)->SetRightMargin(0.02);
  crcp->GetPad(1)->SetBottomMargin(0.10);
  crcp->GetPad(1)->SetLeftMargin(0.10);
  crcp->GetPad(1)->SetTicks(1, 1);

  crcp->cd(1);
  haxis_rcp->GetXaxis()->SetRangeUser(0, 0.5);
  haxis_rcp->Draw();
  for (int isys = 0; isys < NSYSTEMS; isys++)
    grcp_BBCscMod[isys]->Draw("C");

  gRCP_jet->Draw("P");
  for (int i = 0; i < NJET; i++)
    bRcp_jet[i]->Draw();

  l1.DrawLine(0, 1, 0.5, 1);
  leg->Draw("same");


  TCanvas *cncollmb = new TCanvas("cncollmb", "ncoll MB", 1200, 400);
  cncollmb->Divide(NSYSTEMS, 1, 0, 0);
  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cncollmb->GetPad(isys + 1)->SetTopMargin(0.02);
    cncollmb->GetPad(isys + 1)->SetRightMargin(0.02);
    cncollmb->GetPad(isys + 1)->SetBottomMargin(0.10);
    cncollmb->GetPad(isys + 1)->SetLeftMargin(0.10);
    cncollmb->GetPad(isys + 1)->SetTicks(1, 1);

    hNcoll_MB[isys][0]->GetXaxis()->SetRangeUser(0, 50);
    hNcoll_MB[isys][0]->SetLineColor(kBlack);
    hNcoll_MB[isys][0]->SetLineWidth(2);

    cncollmb->cd(isys + 1);
    gPad->SetLogy();
    hNcoll_MB[isys][0]->Draw();
    int ic = 2;
    for (int ix = 0; ix < NX; ix += 2)
    {
      hNcollMod_MB[isys][ix]->SetLineColor(ic);
      hNcollMod_MB[isys][ix]->Draw("same");
      ic++;
    }
    if (isys == 0)
      legncollmb->Draw("same");

  }


  TCanvas *cncoll = new TCanvas("cncoll", "Ncoll", 1200, 800);
  cncoll->Divide(NSYSTEMS, 2, 0, 0);

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {

    cncoll->GetPad(isys + 1)->SetTopMargin(0.02);
    cncoll->GetPad(isys + 1)->SetRightMargin(0.02);
    cncoll->GetPad(isys + 1)->SetBottomMargin(0.10);
    cncoll->GetPad(isys + 1)->SetLeftMargin(0.10);
    cncoll->GetPad(isys + 1)->SetTicks(1, 1);

    cncoll->GetPad(isys + NSYSTEMS + 1)->SetTopMargin(0.02);
    cncoll->GetPad(isys + NSYSTEMS + 1)->SetRightMargin(0.02);
    cncoll->GetPad(isys + NSYSTEMS + 1)->SetBottomMargin(0.10);
    cncoll->GetPad(isys + NSYSTEMS + 1)->SetLeftMargin(0.10);
    cncoll->GetPad(isys + NSYSTEMS + 1)->SetTicks(1, 1);

    for (int i = 0; i < 2; i++)
    {
      cncoll->cd(isys + NSYSTEMS * i + 1);
      hNcoll[isys][xplot][i]->SetLineWidth(2);

      hNcoll[isys][xplot][i]->SetLineColor(kBlue);
      hNcollBBCscMod[isys][xplot][i]->SetLineColor(kRed);

      hNcoll[isys][xplot][i]->Scale(1. / hNcoll[isys][xplot][i]->Integral());
      hNcollBBCscMod[isys][xplot][i]->Scale(1. / hNcollBBCscMod[isys][xplot][i]->Integral());

      hNcoll[isys][xplot][i]->GetXaxis()->SetRangeUser(0, 60);
      hNcoll[isys][xplot][i]->GetYaxis()->SetRangeUser(1e-3, 1);

      gPad->SetLogy();
      hNcoll[isys][xplot][i]->Draw();
      hNcollBBCscMod[isys][xplot][i]->Draw("same");

      legNcoll[isys][i]->Draw("same");
    }

  }


  TCanvas *cbbc = new TCanvas("cbbc", "bbc", 1200, 400);
  cbbc->Divide(NSYSTEMS, 1, 0, 0);

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cbbc->GetPad(isys + 1)->SetTopMargin(0.02);
    cbbc->GetPad(isys + 1)->SetRightMargin(0.02);
    cbbc->GetPad(isys + 1)->SetBottomMargin(0.10);
    cbbc->GetPad(isys + 1)->SetLeftMargin(0.10);
    cbbc->GetPad(isys + 1)->SetTicks(1, 1);


    hBBCs[isys][xplot]->SetTitle(";BBCs charge");
    hBBCs[isys][xplot]->SetLineWidth(2);
    hBBCs[isys][xplot]->SetLineColor(kBlue);
    hBBCsMod[isys][xplot]->SetLineColor(kRed);

    hBBCs[isys][xplot]->GetXaxis()->SetRangeUser(0, 150);
    hBBCs[isys][xplot]->GetYaxis()->SetRangeUser(1, 5e3);

    cbbc->cd(isys + 1);
    gPad->SetLogy();
    hBBCs[isys][xplot]->Draw();
    hBBCsMod[isys][xplot]->Draw("same");
    legBBC[isys]->Draw("same");
  }

  TCanvas *cyield = new TCanvas("cyield", "yield", 1200, 400);
  cyield->Divide(NSYSTEMS, 1, 0, 0);

  for (int isys = 0; isys < NSYSTEMS; isys++)
  {
    cyield->GetPad(isys + 1)->SetTopMargin(0.02);
    cyield->GetPad(isys + 1)->SetRightMargin(0.02);
    cyield->GetPad(isys + 1)->SetBottomMargin(0.10);
    cyield->GetPad(isys + 1)->SetLeftMargin(0.10);
    cyield->GetPad(isys + 1)->SetTicks(1, 1);


    hBBCscNcoll[isys][xplot]->SetLineWidth(2);
    hBBCscNcoll[isys][xplot]->SetLineColor(kBlue);
    hBBCscModNcoll[isys][xplot]->SetLineColor(kRed);
    hBBCscModNcollMod[isys][xplot]->SetLineColor(kGreen + 2);

    hBBCscNcoll[isys][xplot]->GetXaxis()->SetRangeUser(0, 150);
    hBBCscNcoll[isys][xplot]->GetYaxis()->SetRangeUser(1, 5e5);

    cyield->cd(isys + 1);
    gPad->SetLogy();
    hBBCscNcoll[isys][xplot]->Draw("H");
    hBBCscModNcoll[isys][xplot]->Draw("H same");
    hBBCscModNcollMod[isys][xplot]->Draw("H same");
    legBBC[isys]->Draw("same");
  }


  TCanvas *caltrcp = new TCanvas("caltrcp", "RCP", 800, 800);
  caltrcp->Divide(2, 1);

  caltrcp->GetPad(1)->SetTopMargin(0.02);
  caltrcp->GetPad(1)->SetRightMargin(0.02);
  caltrcp->GetPad(1)->SetBottomMargin(0.10);
  caltrcp->GetPad(1)->SetLeftMargin(0.10);
  caltrcp->GetPad(1)->SetTicks(1, 1);

  caltrcp->GetPad(2)->SetTopMargin(0.02);
  caltrcp->GetPad(2)->SetRightMargin(0.02);
  caltrcp->GetPad(2)->SetBottomMargin(0.10);
  caltrcp->GetPad(2)->SetLeftMargin(0.10);
  caltrcp->GetPad(2)->SetTicks(1, 1);

  caltrcp->cd(1);
  haxis_rcp->GetXaxis()->SetRangeUser(0, 0.5);
  haxis_rcp->Draw();
  for (int isys = 0; isys < NSYSTEMS; isys++)
    galtrcp_BBCscModNcoll[isys]->Draw("C");

  gRCP_jet->Draw("P");
  for (int i = 0; i < NJET; i++)
    bRcp_jet[i]->Draw();

  l1.DrawLine(0, 1, 0.5, 1);
  leg->Draw("same");

  caltrcp->cd(2);
  haxis_rcp->GetXaxis()->SetRangeUser(0, 0.5);
  haxis_rcp->Draw();
  for (int isys = 0; isys < NSYSTEMS; isys++)
    galtrcp_BBCscModNcollMod[isys]->Draw("C");

  gRCP_jet->Draw("P");
  for (int i = 0; i < NJET; i++)
    bRcp_jet[i]->Draw();

  l1.DrawLine(0, 1, 0.5, 1);
  leg->Draw("same");



  // TCanvas *crcpcent = new TCanvas("crcpcent", "rcp cent", 1200, 600);
  // crcpcent->Divide(NSYSTEMS, 1);


  // for (int isys = 0; isys < NSYSTEMS; isys++)
  // {
  //   crcpcent->GetPad(isys + 1)->SetTicks(1, 1);
  //   crcpcent->cd(isys + 1);
  //   haxis_rcpcent->Draw();
  //   for (int ix = 0; ix < NX; ix++)
  //     grcp_cent[isys][ix]->Draw("L");

  //   l1.DrawLine(0, 1, NCENT - 1, 1);
  //   label.DrawLatex(0.5, 0.95, collSystem[isys]);
  //   if (isys == NSYSTEMS - 1) legx->Draw("same");
  // }

  // //=====================================================//
  // // SAVE
  // //=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving RCP to " << outFile << endl;

    crcp->Print("Rcp_systems.pdf");
    cncoll->Print("Ncoll_systems.pdf");
    cbbc->Print("BBCsc_systems.pdf");
    cncollmb->Print("MBNcoll_systems.pdf");

    // crcpcent->Print("R_systems.pdf");

    // TFile *fout = new TFile(outFile, "UPDATE");

    // grcp->Write(Form("Rcp_%s", collSystem));

    // fout->Close();
  }
}