////////////////////////////////////////////////////////////////////////////////
//
// Calculate RCP for small systems - (p, d, He3)+Au - from Glauber model
// output.
//
// Calculated explicitly for d+Au for comparison with the PHENIX Jet Rcp
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// Oct 15 2015
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

void calculate_RCP_dAuJet()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //=====================================================//
  // SET RUNNING CONDITIONS
  //=====================================================//

  const int NX = 11;             // Number of x values
  double x[] = {0.01, 0.05, 0.10, 0.15, 0.20,
                0.25, 0.30, 0.35, 0.40, 0.45,
                0.50
               };         // x values
  const char *xFiles[] = {    // Files for each x value
    "rootfiles/glauber_dau_snn42_x001_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x005_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x01_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x015_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x02_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x025_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x03_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x035_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x04_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x045_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x05_ntuple_100k.root",
  };

  // const int NX = 2;             // Number of x values
  // double x[] = {0.2, 0.5};         // x values
  // const char *xFiles[] = {    // Files for each x value
  //   "glauber_dau_snn42_x02_ntuple_100k.root",
  //   "glauber_dau_snn42_x05_ntuple_100k.root",
  // };

  // Set the collision system. Current options are:
  // "dAu"   - d+Au 200 GeV
  const char *collSystem = "dAu";

  bool saveRcp = true;

  //=====================================================//
  // DECLARE VARIABLES
  //=====================================================//

  // collision system variables
  double NBD_mu = 0;
  double NBD_k = 0;
  double trig_eff_params[2] = {0};

  TF1 *fNBD = new TF1("fNBD", NBD, 0, 200, 2);

  TF1* feff = new TF1("feff", "1.0 - TMath::Exp(-1.*TMath::Power(x/[0], [1]))", 0, 200);

  // For running over ntuples
  TFile *fin;
  TTree *ntp;
  const char *ntpName = "";
  Float_t Ncoll;


  TH2D *hNcoll_NcollMod[NX];
  TH2D *hNcollA_Ncoll[NX];
  TH2D *hNcollA_NcollMod[NX];
  TH2D *hNcollModA_Ncoll[NX];
  TH2D *hNcollModA_NcollMod[NX];

  TH2D *hNcollMod_BBCsc[NX];
  TH2D *hNcoll_BBCsc[NX];
  TH2D *hNcollMod_BBCscMod[NX];
  TH2D *hNcoll_BBCscMod[NX];
  TH2D *hNcollA_BBCsc[NX];
  TH2D *hNcollA_BBCscMod[NX];
  TH2D *hNcollModA_BBCsc[NX];
  TH2D *hNcollModA_BBCscMod[NX];
  for (int i = 0; i < NX; i++)
  {

    hNcollMod_BBCsc[i] = new TH2D(Form("hNcollMod_BBCsc_%i", i),
                                  ";BBCs charge; N_{coll}^{mod}",
                                  600, 0, 150,
                                  101, -0.5, 100.5);

    hNcoll_BBCsc[i] = new TH2D(Form("hNcoll_BBCsc_%i", i),
                               ";BBCs charge; N_{coll}",
                               600, 0, 150,
                               101, -0.5, 100.5);

    hNcollMod_BBCscMod[i] = new TH2D(Form("hNcollMod_BBCscMod_%i", i),
                                     ";BBCs charge Mod; N_{coll}^{mod}",
                                     600, 0, 150,
                                     101, -0.5, 100.5);

    hNcoll_BBCscMod[i] = new TH2D(Form("hNcoll_BBCscMod_%i", i),
                                  ";BBCs charge Mod; N_{coll}",
                                  600, 0, 150,
                                  101, -0.5, 100.5);


    hNcollA_BBCsc[i] = new TH2D(Form("hNcollA_BBCsc_%i", i),
                                ";BBCs charge; N_{coll}(A)",
                                600, 0, 150,
                                101, -0.5, 100.5);

    hNcollA_BBCscMod[i] = new TH2D(Form("hNcollA_BBCscMod_%i", i),
                                   ";BBCs charge Mod; N_{coll}(A)",
                                   600, 0, 150,
                                   101, -0.5, 100.5);

    hNcollModA_BBCsc[i] = new TH2D(Form("hNcollModA_BBCsc_%i", i),
                                   ";BBCs charge; N_{coll}^{mod}(A)",
                                   600, 0, 150,
                                   101, -0.5, 100.5);

    hNcollModA_BBCscMod[i] = new TH2D(Form("hNcollModA_BBCscMod_%i", i),
                                      ";BBCs charge Mod; N_{coll}^{mod}(A)",
                                      600, 0, 150,
                                      101, -0.5, 100.5);

  }


  TH1D *hNcoll[NX][2];
  TH1D *hNcollMod[NX][2];
  TH1D *hNcollModBBCscMod[NX][2];
  TH1D *hNcollBBCscMod[NX][2];
  TH1D *hBBCs[NX];
  TH1D *hBBCsMod[NX];

  TH1D *hNcollA[NX][2];
  TH1D *hNcollABBCscMod[NX][2];
  TH1D *hNcollModA[NX][2];
  TH1D *hNcollModABBCscMod[NX][2];

  TH1D* htmp;


  // For counting
  double Nevent[NX][2];
  double NeventMod[NX][2];

  double bias_BBCcsMod[NX][2];
  double bias_NcollMod[NX][2];
  double bias_NcollModBBCcsMod[NX][2];

  double rcp_BBCscMod[NX];
  double rcp_NcollMod[NX];
  double rcp_NcollModBBCscMod[NX];

  // Calculating RCP
  TGraph *grcp_NcollMod;
  TGraph *grcp_BBCscMod;
  TGraph *grcp_NcollModBBCscMod;

  //=====================================================//
  // SET UP COLLISION SYSTEM VALUES
  //=====================================================//
  cout << endl;
  cout << "--> Setting up the collision system ..." << endl;


  cout << "  System is d+Au" << endl;
  //http://www.phenix.bnl.gov/phenix/WWW/p/info/an/900/Run8_dAu_200GeV_Centrality_Categorization.pdf
  //http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1087/Run8_dAu_200GeV_Centrality_Addendum-01.pdf

  NBD_mu = 3.038;
  NBD_k  = 0.464;

  // currently from my thesis
  trig_eff_params[0] = 0.897;
  trig_eff_params[1] = 0.612;

  ntpName = "nt_dh_Au";

  // print values
  cout << "  NBD : mu=" << NBD_mu << " k=" << NBD_k << endl;
  cout << "  Trig Eff Pars: " << trig_eff_params[0] << ", " << trig_eff_params[1] << endl;
  cout << "  Ntuple name: " << ntpName << endl;

  fNBD->SetParameters(NBD_mu, NBD_k);
  feff->SetParameters(trig_eff_params[0], trig_eff_params[1]);


  //=====================================================//
  // GET NTUPLE(S) FROM FILE AND CALCULATE YIELD
  //=====================================================//
  cout << endl;
  cout << "--> Reading Ntuples from files ..." << endl;


  for (int ix = 0; ix < NX; ix++)
  {
    fin = TFile::Open(xFiles[ix]);
    if (!fin)
    {
      cout << "ERROR!! Unable to open " << xFiles[ix] << endl;
      return;
    }
    else
      cout << "----> Reading " << xFiles[ix] << " for x=" << x[ix] << endl;

    ntp = (TTree*) fin->Get(ntpName);
    if (!ntp)
    {
      cout << "ERROR!! Unable to find " << ntpName << " in " << xFiles[ix] << endl;
      return;
    }

    ntp->Draw("Sum$(NcollA):Sum$(NcollModA) >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
              "", "goff");
    hNcoll_NcollMod[ix] = (TH2D*) gDirectory->FindObject("htmp");
    hNcoll_NcollMod[ix]->SetDirectory(0);
    hNcoll_NcollMod[ix]->SetName(Form("hNcoll_NcollMod_%i", ix));
    hNcoll_NcollMod[ix]->SetTitle(";N_{coll}^{mod};N_{coll}");

    // for assuming hard process yield scales with Ncoll in only the "modified"
    // nucleon
    ntp->Draw("NcollA[modA]:Sum$(NcollA) >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
              "", "goff");
    hNcollA_Ncoll[ix] = (TH2D*) gDirectory->FindObject("htmp");
    hNcollA_Ncoll[ix]->SetDirectory(0);
    hNcollA_Ncoll[ix]->SetName(Form("hNcollModA_Ncoll_%i", ix));
    hNcollA_Ncoll[ix]->SetTitle(";N_{coll};N_{coll}^{mod}(A)");

    ntp->Draw("NcollA[modA]:Sum$(NcollModA) >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
              "", "goff");
    hNcollA_NcollMod[ix] = (TH2D*) gDirectory->FindObject("htmp");
    hNcollA_NcollMod[ix]->SetDirectory(0);
    hNcollA_NcollMod[ix]->SetName(Form("hNcollModA_Ncoll_%i", ix));
    hNcollA_NcollMod[ix]->SetTitle(";N_{coll}^{mod};N_{coll}^{mod}(A)");

    ntp->Draw("NcollModA[modA]:Sum$(NcollA) >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
              "", "goff");
    hNcollModA_Ncoll[ix] = (TH2D*) gDirectory->FindObject("htmp");
    hNcollModA_Ncoll[ix]->SetDirectory(0);
    hNcollModA_Ncoll[ix]->SetName(Form("hNcollModA_Ncoll_%i", ix));
    hNcollModA_Ncoll[ix]->SetTitle(";N_{coll};N_{coll}^{mod}(A)");

    ntp->Draw("NcollModA[modA]:Sum$(NcollModA) >> htmp(101, -0.5, 100.5, 101, -0.5, 100.5)",
              "", "goff");
    hNcollModA_NcollMod[ix] = (TH2D*) gDirectory->FindObject("htmp");
    hNcollModA_NcollMod[ix]->SetDirectory(0);
    hNcollModA_NcollMod[ix]->SetName(Form("hNcollA_NcollMod_%i", ix));
    hNcollModA_NcollMod[ix]->SetTitle(";N_{coll}^{mod};N_{coll}^{mod}(A)");

    // we're done with the file, might as well close it
    delete ntp;
    fin->Close();
    delete fin;



    // convolve with NBD to get BBCs charge & apply trigger efficiency
    for (int ibx = 1; ibx <= hNcoll_NcollMod[ix]->GetNbinsX(); ibx++)
    {
      for (int iby = 1; iby <= hNcoll_NcollMod[ix]->GetNbinsY(); iby++)
      {

        float w = hNcoll_NcollMod[ix]->GetBinContent(ibx, iby);
        if (w == 0) continue;

        float Ncoll = hNcoll_NcollMod[ix]->GetYaxis()->GetBinCenter(iby);
        float NcollMod = hNcoll_NcollMod[ix]->GetXaxis()->GetBinCenter(ibx);

        // Unmodified BBCs charge
        if (Ncoll > 0)
        {
          fNBD->SetParameters(Ncoll * NBD_mu, Ncoll * NBD_k);
          for (int j = 1; j <= hNcoll_BBCsc[ix]->GetNbinsX(); j++)
          {
            double c = hNcoll_BBCsc[ix]->GetXaxis()->GetBinCenter(j);
            double bc = hNcoll_BBCsc[ix]->GetBinContent(j, iby);
            double bcMod = hNcollMod_BBCsc[ix]->GetBinContent(j, ibx);

            double P = fNBD->Eval(c) * feff->Eval(c);

            hNcoll_BBCsc[ix]->SetBinContent(j, iby, w * P + bc);
            hNcollMod_BBCsc[ix]->SetBinContent(j, ibx, w * P + bcMod);
          }
        }
        // Modified BBCs charge
        if (NcollMod > 0)
        {
          fNBD->SetParameters(NcollMod * NBD_mu, NcollMod * NBD_k);
          for (int j = 1; j <= hNcoll_BBCscMod[ix]->GetNbinsX(); j++)
          {
            double c = hNcoll_BBCscMod[ix]->GetXaxis()->GetBinCenter(j);
            double bc = hNcoll_BBCscMod[ix]->GetBinContent(j, iby);
            double bcMod = hNcollMod_BBCscMod[ix]->GetBinContent(j, ibx);

            double P = fNBD->Eval(c) * feff->Eval(c);

            hNcoll_BBCscMod[ix]->SetBinContent(j, iby, w * P + bc);
            hNcollMod_BBCscMod[ix]->SetBinContent(j, ibx, w * P + bcMod);
          }
        }

      } // iby
    } // ibx


    // Project over unmodified/modified BBCcs (unmodified Ncoll)
    hBBCs[ix] = (TH1D*) hNcoll_BBCsc[ix]->ProjectionX(Form("hBBCs_%i", ix),
                1,
                hNcoll_BBCsc[ix]->GetNbinsY());

    hBBCsMod[ix] = (TH1D*) hNcoll_BBCscMod[ix]->ProjectionX(Form("hBBCsMod_%i", ix),
                   1,
                   hNcoll_BBCscMod[ix]->GetNbinsY());

    hBBCs[ix]->SetDirectory(0);
    hBBCsMod[ix]->SetDirectory(0);

    // Calculate the desired charge limits for each centrality bin
    int nq = 2;
    double xq[] = {(88. - 60.) / 88., (88. - 20.) / 88.}; // 60-88%, 0-20%
    double yq[2] = {0};
    hBBCs[ix]->GetQuantiles(nq, yq, xq);
    cout << "  quantiles: " << yq[0] << " " << yq[1] << endl;

    // get bin limits for each centrality bin
    double blim[2][2];
    //0-20%
    blim[0][0] = hBBCs[ix]->FindBin(yq[1]);
    blim[0][1] = hNcollA_BBCsc[ix]->GetNbinsX();
    //60-88%
    blim[1][0] = 1;
    blim[1][1] = hBBCs[ix]->FindBin(yq[0]);

    cout << "  blim[0]: [" << blim[0][0] << ", " << blim[0][1] << "]"
         << " -> [" << hBBCs[ix]->GetBinLowEdge(blim[0][0]) << ","
         << hBBCs[ix]->GetBinLowEdge(blim[0][1] + 1) << "]"
         << endl;
    cout << "  blim[1]: [" << blim[1][0] << ", " << blim[1][1] << "]"
         << " -> [" << hBBCs[ix]->GetBinLowEdge(blim[1][0]) << ","
         << hBBCs[ix]->GetBinLowEdge(blim[1][1] + 1) << "]"
         << endl;

    double yield[2] = {0};
    double yield_NcollMod[2] = {0};
    double yield_BBCscMod[2] = {0};
    double yield_NcollModBBCscMod[2] = {0};


    for (int iq = 0; iq < 2; iq++)
    {
      // Calculate the unmodified number of events
      Nevent[ix][iq] = hBBCs[ix]->Integral(blim[iq][0], blim[iq][1]);

      // Calculate the number of modified events
      NeventMod[ix][iq] = hBBCsMod[ix]->Integral(blim[iq][0], blim[iq][1]);

      // Project the unmodified Ncoll (unmodified BBCs charge)
      hNcoll[ix][iq] = (TH1D*) hNcoll_BBCsc[ix]->ProjectionY(
                        Form("hNcoll_%i_%i", ix, iq),
                        blim[iq][0], blim[iq][1]);

      // Project the modified Ncoll (unmodified BBCs charge)
      hNcollMod[ix][iq] = (TH1D*) hNcollMod_BBCsc[ix]->ProjectionY(
                           Form("hNcollMod_%i_%i", ix, iq),
                           blim[iq][0], blim[iq][1]);

      // Project the unmodified Ncoll (modified BBCs charge)
      hNcollBBCscMod[ix][iq] = (TH1D*) hNcoll_BBCscMod[ix]->ProjectionY(
                                Form("hNcollBBCscMod_%i_%i", ix, iq),
                                blim[iq][0], blim[iq][1]);

      // Project the modified Ncoll for modified BBCs charge)
      hNcollModBBCscMod[ix][iq] = (TH1D*) hNcollMod_BBCscMod[ix]->ProjectionY(
                                   Form("hNcollModBBCscMod_%i_%i", ix, iq),
                                   blim[iq][0], blim[iq][1]);



      // Calculate yields (Ncoll weighted)
      for (int ib = 1; ib <= hNcoll[ix][iq]->GetNbinsX(); ib++)
      {
        double ncoll = hNcoll[ix][iq]->GetBinCenter(ib);

        yield[iq]                  += ncoll * hNcoll[ix][iq]->GetBinContent(ib);
        yield_NcollMod[iq]         += ncoll * hNcollMod[ix][iq]->GetBinContent(ib);
        yield_BBCscMod[iq]         += ncoll * hNcollBBCscMod[ix][iq]->GetBinContent(ib);
        yield_NcollModBBCscMod[iq] += ncoll * hNcollModBBCscMod[ix][iq]->GetBinContent(ib);

      }


      // Think in terms of bias factors! (Drop <Ncoll>) and divide by unmodified
      // Calculate the bias factor for each centrality bin
      bias_BBCcsMod[ix][iq] = yield_BBCscMod[iq] / NeventMod[ix][iq];
      bias_BBCcsMod[ix][iq] /= yield[iq] / Nevent[ix][iq];

      bias_NcollMod[ix][iq] = yield_NcollMod[iq] / Nevent[ix][iq];
      bias_NcollMod[ix][iq] /= yield[iq] / Nevent[ix][iq];

      bias_NcollModBBCcsMod[ix][iq] = yield_NcollModBBCscMod[iq] / NeventMod[ix][iq];
      bias_NcollModBBCcsMod[ix][iq] /= yield[iq] / Nevent[ix][iq];



    } //iq

    cout << " x: " << x[ix] << endl;
    cout << "   <Ncoll(0-20%)>       : " << hNcoll[ix][0]->GetMean() << endl;
    cout << "   <Ncoll^{mod}(0-20%)> : " << hNcollMod[ix][0]->GetMean() << endl;
    cout << "   Nevent(0-20%)        : " << Nevent[ix][0] << endl;
    cout << "   Nevent^{mod}(0-20%)  : " << NeventMod[ix][0] << endl;
    cout << "   yield(0-20%)         : " << yield[0] << endl;
    cout << "   yield^{Nmod}(0-20%)  : " << yield_NcollMod[0] << endl;
    cout << "   yield^{Bmod}(0-20%)  : " << yield_BBCscMod[0] << endl;
    cout << "   yield^{NBmod}(0-20%) : " << yield_NcollModBBCscMod[0] << endl;
    cout << endl;
    cout << "   <Ncoll(60-88%)>      : " << hNcoll[ix][1]->GetMean() << endl;
    cout << "   <Ncoll^{mod}(60-88%)>: " << hNcollMod[ix][1]->GetMean() << endl;
    cout << "   Nevent(60-88%)       : " << Nevent[ix][1] << endl;
    cout << "   Nevent^{mod}(60-88%) : " << NeventMod[ix][1] << endl;
    cout << "   yield(60-88%)        : " << yield[1] << endl;
    cout << "   yield^{Nmod}(60-88%) : " << yield_NcollMod[1] << endl;
    cout << "   yield^{Bmod}(60-88%) : " << yield_BBCscMod[1] << endl;
    cout << "   yield^{NBmod}(60-88%): " << yield_NcollModBBCscMod[1] << endl;


    // Calculate Rcp
    rcp_BBCscMod[ix] = bias_BBCcsMod[ix][0] / bias_BBCcsMod[ix][1];

    rcp_NcollMod[ix] = bias_NcollMod[ix][0] / bias_NcollMod[ix][1];

    rcp_NcollModBBCscMod[ix] = bias_NcollModBBCcsMod[ix][0] / bias_NcollModBBCcsMod[ix][1];

    cout << endl;
    cout << "   Bias(0-20%) (BBCcs Mod) : " << bias_BBCcsMod[ix][0] << endl;
    cout << "   Bias(0-20%) (Ncoll Mod) : " << bias_NcollMod[ix][0] << endl;
    cout << "   Bias(0-20%) (both Mod)  : " << bias_NcollModBBCcsMod[ix][0] << endl;
    cout << "   Bias(60-88%) (BBCcs Mod): " << bias_BBCcsMod[ix][1] << endl;
    cout << "   Bias(60-88%) (Ncoll Mod): " << bias_NcollMod[ix][1] << endl;
    cout << "   Bias(60-88%) (both Mod) : " << bias_NcollModBBCcsMod[ix][1] << endl;


    cout << endl;
    cout << "   Rcp (BBCcs Mod) : " << rcp_BBCscMod[ix] << endl;
    cout << "   Rcp (Ncoll Mod) : " << rcp_NcollMod[ix] << endl;
    cout << "   Rcp (both Mod)  : " << rcp_NcollModBBCscMod[ix] << endl;




    //--- Considering only the "modified" nucleon ---//

    // Remember, the binning on all histograms is the same
    for (int iby = 1; iby <= hNcollA_Ncoll[ix]->GetNbinsY(); iby++)
    {
      for (int ibx = 1; ibx <= hNcollA_Ncoll[ix]->GetNbinsX(); ibx++)
      {

        double Ncoll = hNcollA_Ncoll[ix]->GetXaxis()->GetBinCenter(ibx);
        if (Ncoll == 0) continue;
        double NcollA = hNcollA_Ncoll[ix]->GetYaxis()->GetBinCenter(iby);
        fNBD->SetParameters(Ncoll * NBD_mu, Ncoll * NBD_k);

        double G_NcollA_Ncoll = hNcollA_Ncoll[ix]->GetBinContent(ibx, iby);
        double G_NcollModA_Ncoll = hNcollModA_Ncoll[ix]->GetBinContent(ibx, iby);
        double G_NcollA_NcollMod = hNcollA_NcollMod[ix]->GetBinContent(ibx, iby);
        double G_NcollModA_NcollMod = hNcollModA_NcollMod[ix]->GetBinContent(ibx, iby);

        //save some time if there are no entries
        if (! (G_NcollA_Ncoll > 0 || G_NcollModA_Ncoll > 0 ||
               G_NcollModA_Ncoll > 0 || G_NcollModA_NcollMod > 0))
          continue;

        for (int ic = 1; ic <= hNcollA_BBCsc[ix]->GetNbinsX(); ic++)
        {
          double c = hNcollA_BBCsc[ix]->GetXaxis()->GetBinCenter(ic);
          // double NBD = fNBD->Eval(c);
          double NBD = evalNBD(c, Ncoll * NBD_mu, Ncoll * NBD_k);
          double eff = feff->Eval(c);
          // double eff = 1;

          hNcollA_BBCsc[ix]->Fill(c, NcollA, G_NcollA_Ncoll * NBD * eff);
          hNcollModA_BBCsc[ix]->Fill(c, NcollA, G_NcollModA_Ncoll * NBD * eff);
          hNcollA_BBCscMod[ix]->Fill(c, NcollA, G_NcollA_NcollMod * NBD * eff);
          hNcollModA_BBCscMod[ix]->Fill(c, NcollA, G_NcollModA_NcollMod * NBD * eff);
        } // ic
      } // ibx
    } // iby

    cout << endl;
    cout << "Considering only the 'modified' nucleon" << endl;

    cout << "   <NcollA>   : "
         << hNcollA_Ncoll[ix]->ProjectionY("t", 1, hNcollA_Ncoll[ix]->GetNbinsX())->GetMean() << " "
         << hNcollA_NcollMod[ix]->ProjectionY("t", 1, hNcollA_NcollMod[ix]->GetNbinsX())->GetMean() << " "
         << hNcollA_BBCsc[ix]->ProjectionY("t", 1, hNcollA_BBCsc[ix]->GetNbinsX())->GetMean() << " "
         << hNcollA_BBCscMod[ix]->ProjectionY("t", 1, hNcollA_BBCscMod[ix]->GetNbinsX())->GetMean()
         << endl;
    cout << "   <NcollModA>: "
         << hNcollModA_Ncoll[ix]->ProjectionY("t", 1, hNcollModA_Ncoll[ix]->GetNbinsX())->GetMean() << " "
         << hNcollModA_NcollMod[ix]->ProjectionY("t", 1, hNcollModA_NcollMod[ix]->GetNbinsX())->GetMean() << " "
         << hNcollModA_BBCsc[ix]->ProjectionY("t", 1, hNcollModA_BBCsc[ix]->GetNbinsX())->GetMean() << " "
         << hNcollModA_BBCscMod[ix]->ProjectionY("t", 1, hNcollModA_BBCscMod[ix]->GetNbinsX())->GetMean()
         << endl;

    // Calculate yields (Ncoll weighted)
    double yield_NcollA_BBCsc[2] = {0};
    double yield_NcollA_BBCscMod[2] = {0};
    double yield_NcollModA_BBCsc[2] = {0};
    double yield_NcollModA_BBCscMod[2] = {0};


    for (int iq = 0; iq < nq; iq++)
    {
      for (int ib = 1; ib <= hNcollA_BBCsc[ix]->GetNbinsY(); ib++)
      {
        double NcollA = hNcollA_BBCsc[ix]->GetYaxis()->GetBinCenter(ib);

        yield_NcollA_BBCsc[iq]       +=
          NcollA * hNcollA_BBCsc[ix]->Integral(blim[iq][0], blim[iq][1], ib, ib);
        yield_NcollA_BBCscMod[iq]    +=
          NcollA * hNcollA_BBCscMod[ix]->Integral(blim[iq][0], blim[iq][1], ib, ib);
        yield_NcollModA_BBCsc[iq]    +=
          NcollA * hNcollModA_BBCsc[ix]->Integral(blim[iq][0], blim[iq][1], ib, ib);
        yield_NcollModA_BBCscMod[iq] +=
          NcollA * hNcollModA_BBCscMod[ix]->Integral(blim[iq][0], blim[iq][1], ib, ib);

      }

      cout << endl;
      cout << " yield[" << iq << "](NcollA, BBCsc)      : " << yield_NcollA_BBCsc[iq] << endl;
      cout << " yield[" << iq << "](NcollA, BBCscMod)   : " << yield_NcollA_BBCscMod[iq] << endl;
      cout << " yield[" << iq << "](NcollModA, BBCsc)   : " << yield_NcollModA_BBCsc[iq] << endl;
      cout << " yield[" << iq << "](NcollModA, BBCscMod): " << yield_NcollModA_BBCscMod[iq] << endl;


      // turn these into bias factors
      yield_NcollA_BBCsc[iq]       = yield_NcollA_BBCsc[iq] / Nevent[ix][iq];
      yield_NcollA_BBCscMod[iq]    = yield_NcollA_BBCscMod[iq] / NeventMod[ix][iq];
      yield_NcollModA_BBCsc[iq]    = yield_NcollModA_BBCsc[iq] / Nevent[ix][iq];
      yield_NcollModA_BBCscMod[iq] = yield_NcollModA_BBCscMod[iq] / NeventMod[ix][iq];

      cout << endl;
      cout << " yield[" << iq << "](NcollA, BBCsc)      : " << yield_NcollA_BBCsc[iq] << endl;
      cout << " yield[" << iq << "](NcollA, BBCscMod)   : " << yield_NcollA_BBCscMod[iq] << endl;
      cout << " yield[" << iq << "](NcollModA, BBCsc)   : " << yield_NcollModA_BBCsc[iq] << endl;
      cout << " yield[" << iq << "](NcollModA, BBCscMod): " << yield_NcollModA_BBCscMod[iq] << endl;

      yield_NcollA_BBCscMod[iq]    /= yield_NcollA_BBCsc[iq];
      yield_NcollModA_BBCsc[iq]    /= yield_NcollA_BBCsc[iq];
      yield_NcollModA_BBCscMod[iq] /= yield_NcollA_BBCsc[iq];
      yield_NcollA_BBCsc[iq]       /= yield_NcollA_BBCsc[iq];

      cout << endl;
      cout << " Bias[" << iq << "](NcollA, BBCsc)      : " << yield_NcollA_BBCsc[iq] << endl;
      cout << " Bias[" << iq << "](NcollA, BBCscMod)   : " << yield_NcollA_BBCscMod[iq] << endl;
      cout << " Bias[" << iq << "](NcollModA, BBCsc)   : " << yield_NcollModA_BBCsc[iq] << endl;
      cout << " Bias[" << iq << "](NcollModA, BBCscMod): " << yield_NcollModA_BBCscMod[iq] << endl;
    }


    cout << endl;
    cout << "   Rcp(NcollA, BBCsc)      : "
         << yield_NcollA_BBCsc[0] / yield_NcollA_BBCsc[1]
         << " (" << yield_NcollA_BBCsc[0] << " / " << yield_NcollA_BBCsc[1] << ")"
         << endl;
    cout << "   Rcp(NcollA, BBCscMod)   : "
         << yield_NcollA_BBCscMod[0] / yield_NcollA_BBCscMod[1]
         << " (" << yield_NcollA_BBCscMod[0] << " / " << yield_NcollA_BBCscMod[1] << ")"
         << endl;
    cout << "   Rcp(NcollModA, BBCsc)   : "
         << yield_NcollModA_BBCsc[0] / yield_NcollModA_BBCsc[1]
         << " (" << yield_NcollModA_BBCsc[0] << " / " << yield_NcollModA_BBCsc[1] << ")"
         << endl;
    cout << "   Rcp(NcollModA, BBCscMod): "
         << yield_NcollModA_BBCscMod[0] / yield_NcollModA_BBCscMod[1]
         << " (" << yield_NcollModA_BBCscMod[0] << " / " << yield_NcollModA_BBCscMod[1] << ")"
         << endl;

  }


  grcp_BBCscMod = new TGraph(NX, x, rcp_BBCscMod);
  grcp_BBCscMod->SetName("grcp_BBCscMod");
  grcp_BBCscMod->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
  grcp_BBCscMod->SetLineStyle(1);
  grcp_BBCscMod->SetLineWidth(2);
  grcp_BBCscMod->SetLineColor(kBlue);

  grcp_NcollMod = new TGraph(NX, x, rcp_NcollMod);
  grcp_NcollMod->SetName("grcp_NcollMod");
  grcp_NcollMod->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
  grcp_NcollMod->SetLineStyle(3);
  grcp_NcollMod->SetLineWidth(2);
  grcp_NcollMod->SetLineColor(kGreen + 2);

  grcp_NcollModBBCscMod = new TGraph(NX, x, rcp_NcollModBBCscMod);
  grcp_NcollModBBCscMod->SetName("grcp_NcollModBBCscMod");
  grcp_NcollModBBCscMod->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
  grcp_NcollModBBCscMod->SetLineStyle(2);
  grcp_NcollModBBCscMod->SetLineWidth(2);
  grcp_NcollModBBCscMod->SetLineColor(kRed);




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

  // TLegend *legNcoll;
  // if (collSystem == "pAu")
  //   legNcoll = new TLegend(0.5, 0.5, 0.9, 0.9);
  // else
  //   legNcoll = new TLegend(0.15, 0.15, 0.5, 0.6);
  // legNcoll->SetFillStyle(0);
  // legNcoll->SetBorderSize(0);
  // for (int ix = 0; ix < NX; ix++)
  // {
  //   legNcoll->AddEntry(hNcoll[ix],
  //                      Form("x=%.2f <N_{coll}>=%.2f", x[ix], hNcoll[ix]->GetMean()),
  //                      "L");
  // }

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);

  TLine l1;
  l1.SetLineColor(kBlack);
  l1.SetLineStyle(2);

  TH1F* haxis_rcp = new TH1F("haxis_rcp",
                             ";x (x_{jet}=2 * p_{T}^{jet} / #sqrt{s_{NN}});R_{CP} (0-20%/60-88%)",
                             100, 0, 1);
  haxis_rcp->SetMinimum(0);
  haxis_rcp->SetMaximum(1.1);

  TLegend *legRCP = new TLegend(0.15, 0.15, 0.5, 0.4);
  legRCP->SetFillStyle(0);
  legRCP->SetBorderSize(0);
  legRCP->AddEntry(gRCP_jet, "Run 8 Jet", "P");
  legRCP->AddEntry(grcp_BBCscMod, "Mod BBCs charge only", "L");
  // legRCP->AddEntry(grcp_NcollMod, "Mod N_{coll} only", "L");
  // legRCP->AddEntry(grcp_NcollModBBCscMod, "Mod both N_{coll} and BBCs chrg", "L");



  int xplot = 5;

  TLegend *legNcoll[2];
  const char* centl[2] = {"0-20", "60-88"};
  for (int i = 0; i < 2; i++)
  {
    legNcoll[i] = new TLegend(0.55, 0.7, 0.96, 0.96,
                              Form("d+Au @ 200 GeV, x = %.2f, %s%%",
                                   x[xplot], centl[i]));
    legNcoll[i]->SetFillStyle(0);
    legNcoll[i]->SetBorderSize(0);
    legNcoll[i]->AddEntry(hNcoll[xplot][i],
                          Form("Unmodified <N_{coll}> = %.2f",
                               hNcoll[xplot][i]->GetMean()),
                          "L"
                         );
    legNcoll[i]->AddEntry(hNcollBBCscMod[xplot][i],
                          Form("Mod BBCs chrg <N_{coll}> = %.2f",
                               hNcollBBCscMod[xplot][i]->GetMean()),
                          "L"
                         );
  }

  TLegend *legBBC = new TLegend(0.5, 0.7, 0.96, 0.96,
                                Form("d+Au @ 200 GeV, x = %.2f",
                                     x[xplot]));
  legBBC->SetFillStyle(0);
  legBBC->SetBorderSize(0);
  legBBC->AddEntry(hBBCs[xplot], "Unmodified", "L");
  legBBC->AddEntry(hBBCsMod[xplot], "Modified", "L");

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
  haxis_rcp->GetXaxis()->SetRangeUser(0, 0.6);
  haxis_rcp->Draw();
  grcp_BBCscMod->Draw("C");
  // grcp_NcollMod->Draw("C");
  // grcp_NcollModBBCscMod->Draw("C");
  legRCP->Draw("same");

  gRCP_jet->Draw("P");
  for (int i = 0; i < NJET; i++)
    bRcp_jet[i]->Draw();

  l1.DrawLine(0, 1, 0.6, 1);

  TCanvas *cncoll = new TCanvas("cncoll", "Ncoll", 1000, 500);
  cncoll->Divide(2, 1);

  cncoll->GetPad(1)->SetTopMargin(0.02);
  cncoll->GetPad(1)->SetRightMargin(0.02);
  cncoll->GetPad(1)->SetBottomMargin(0.10);
  cncoll->GetPad(1)->SetLeftMargin(0.10);
  cncoll->GetPad(1)->SetTicks(1, 1);

  cncoll->GetPad(2)->SetTopMargin(0.02);
  cncoll->GetPad(2)->SetRightMargin(0.02);
  cncoll->GetPad(2)->SetBottomMargin(0.10);
  cncoll->GetPad(2)->SetLeftMargin(0.10);
  cncoll->GetPad(2)->SetTicks(1, 1);

  for (int i = 0; i < 2; i++)
  {
    cncoll->cd(i + 1);
    hNcoll[xplot][i]->SetLineWidth(2);

    hNcoll[xplot][i]->SetLineColor(kBlue);
    hNcollBBCscMod[xplot][i]->SetLineColor(kRed);

    gPad->SetLogy();
    hNcoll[xplot][i]->Draw();
    hNcollBBCscMod[xplot][i]->Draw("same");

    legNcoll[i]->Draw("same");
  }


  TCanvas *cbbc = new TCanvas("cbbc", "bbc", 800, 600);
  cbbc->Divide(1, 1);
  cbbc->GetPad(1)->SetTopMargin(0.02);
  cbbc->GetPad(1)->SetRightMargin(0.02);
  cbbc->GetPad(1)->SetBottomMargin(0.10);
  cbbc->GetPad(1)->SetLeftMargin(0.10);
  cbbc->GetPad(1)->SetTicks(1, 1);

  cbbc->cd(1);
  hBBCs[xplot]->SetTitle(";BBCs charge");
  hBBCs[xplot]->SetLineWidth(2);
  hBBCs[xplot]->SetLineColor(kBlue);
  hBBCsMod[xplot]->SetLineColor(kRed);

  gPad->SetLogy();
  hBBCs[xplot]->Draw();
  hBBCsMod[xplot]->Draw("same");
  legBBC->Draw("same");


  //=====================================================//
  // SAVE
  //=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving plots" << endl;

    // ctest->Print("NcollBBCcharge_dAuJet.pdf");
    crcp->Print("Rcp_dAuJet.pdf");
    cncoll->Print("Ncoll_dAuJet.pdf");
    cbbc->Print("BBCs_dAuJet.pdf");
  }
}