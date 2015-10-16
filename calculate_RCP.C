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

#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>

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

  //=====================================================//
  // SET RUNNING CONDITIONS
  //=====================================================//

  const int NX = 3;             // Number of x values
  double x[] = {0.0, 0.05, 0.1};         // x values
  // const char *xFiles[NX] = {    // Files for each x value
  //   "glauber_pau_ntuple_100k.root",
  //   "glauber_pau_snn42_x005_ntuple_100k.root",
  //   "glauber_pau_snn42_x01_ntuple_100k.root",
  // };
  const char *xFiles[] = {    // Files for each x value
    "glauber_dau_snn42_x0_ntuple_100k.root",
    "glauber_dau_snn42_x005_ntuple_100k.root",
    "glauber_dau_snn42_x01_ntuple_100k.root",
  };
  // const char *xFiles[NX] = {    // Files for each x value
  //   "glauber_he3au_snn42_x0_ntuple_100k.root",
  //   "glauber_he3au_snn42_x005_ntuple_100k.root",
  //   "glauber_he3au_snn42_x01_ntuple_100k.root",
  // };

  // Set the collision system. Current options are:
  // "pAu"   - p+Au 200 GeV
  // "dAu"   - d+Au 200 GeV
  // "He3Au" - He3+Au 200 GeV
  const char *collSystem = "dAu";

  bool saveRcp = false;
  const char *outFile = "Rcp_systems.root";

  //=====================================================//
  // DECLARE VARIABLES
  //=====================================================//

  // collision system variables
  double NBD_mu = 0;
  double NBD_k = 0;
  double BBC_cent = 0;
  double BBC_periph = 0;
  double Ncoll_cent = 0;
  double Ncoll_periph = 0;
  double dcent_cent = 0;
  double dcent_periph = 0;
  double trig_eff_params[2] = {0};

  TF1 *fNBD = new TF1("fNBD", NBD, 0, 200, 2);

  TF1* feff = new TF1("feff", "1. - TMath::Exp(-1.*(x/[0])^{[1]})", 0, 200);

  // For running over ntuples
  TFile *fin;
  TNtuple *ntp;
  const char *ntpName = "";
  Float_t Ncoll;

  TH1D *hNcoll[NX];
  TH1D *hBBCs[NX];
  for (int i = 0; i < NX; i++)
  {
    hBBCs[i] = new TH1D(Form("hBBCs_%i", i), ";BBCs charge", 201, -0.5, 200.5);
    hBBCs[i]->SetLineColor(i + 1);
  }

  TH1D* htmp;


  // For counting
  double Nevent_cent[NX] = {0};
  double Nevent_periph[NX] = {0};

  // Calculating RCP
  TGraph *grcp;

  //=====================================================//
  // SET UP COLLISION SYSTEM VALUES
  //=====================================================//
  cout << endl;
  cout << "--> Setting up the collision system ..." << endl;

  if (collSystem == "pAu")
  {
    cout << "  System is p+Au" << endl;
    // https://www.phenix.bnl.gov/WWW/p/draft/nagle/PHENIX/nagle_run15_pau_update_05282015.pdf

    NBD_mu = 3.14;
    NBD_k  = 0.47;

    BBC_cent   = 48.816666; // 0-5% zvertex=0
    BBC_periph =  4.216667; // 70-84% zvertex=0

    Ncoll_cent   = 9.693; // 0-5%
    Ncoll_periph = 2.188;  // 70-84%

    dcent_cent   = 0.05;
    dcent_periph = 0.14;

    //currently for d+Au, need to fix
    trig_eff_params[0] = 0.897;
    trig_eff_params[1] = 0.612;

    ntpName = "nt_p_Au";
  }
  if (collSystem == "dAu")
  {
    cout << "  System is d+Au" << endl;
    //http://www.phenix.bnl.gov/phenix/WWW/p/info/an/900/Run8_dAu_200GeV_Centrality_Categorization.pdf
    //http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1087/Run8_dAu_200GeV_Centrality_Addendum-01.pdf

    NBD_mu = 3.038;
    NBD_k  = 0.464;

    BBC_cent   = 71.165; // 0-5% zvertex=0
    BBC_periph =  7.235; // 70-88% zvertex=0

    Ncoll_cent   = 18.115; // 0-5%
    Ncoll_periph = 2.638;  // 70-88%

    dcent_cent   = 0.05;
    dcent_periph = 0.18;

    // currently from my thesis
    trig_eff_params[0] = 0.897;
    trig_eff_params[1] = 0.612;

    ntpName = "nt_dh_Au";

  }
  if (collSystem == "He3Au")
  {
    cout << "  System is He3+Au" << endl;
    // http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1207/Run14_3HeAu_200GeV_Centrality_Categorization.pdf

    // WARNING!! THESE ARE TEMPORARILY COPIED FROM THE D+AU
    //           NEED TO GET THE REAL VALUES FROM JAMIE (NOT IN THE NOTE?)
    NBD_mu = 3.038;
    NBD_k  = 0.464;

    BBC_cent   = 87.4905; // 0-5% zvertex=0
    BBC_periph =  7.9065; // 70-83% zvertex=0

    Ncoll_cent   = 26.082; // 0-5%
    Ncoll_periph =  2.585; // 70-88%

    dcent_cent   = 0.05;
    dcent_periph = 0.18;

    // currently for d+Au, need to fix
    trig_eff_params[0] = 0.897;
    trig_eff_params[1] = 0.612;

    ntpName = "nt_He3_Au";

  }

  // print values
  cout << "  NBD : mu=" << NBD_mu << " k=" << NBD_k << endl;
  cout << "  BBCs: central>" << BBC_cent << " periph<" << BBC_periph << endl;
  cout << "  NColl: central=" << Ncoll_cent << " periph=" << Ncoll_periph << endl;
  cout << "  dcent: central=" << dcent_cent << " periph=" << dcent_periph << endl;
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

    ntp = (TNtuple*) fin->Get(ntpName);
    if (!ntp)
    {
      cout << "ERROR!! Unable to find " << ntpName << " in " << xFiles[ix] << endl;
      return;
    }

    ntp->SetBranchAddress("Ncoll", &Ncoll);

    // htmp = new TH1D("htmp", ";Ncoll", 51, -0.5, 50.5);
    ntp->Draw("Ncoll>>htmp(50, 0.5, 50.5)", "", "goff");

    hNcoll[ix] = (TH1D*) gDirectory->FindObject("htmp");
    hNcoll[ix]->SetDirectory(0);
    hNcoll[ix]->SetName(Form("hNcoll_%i", ix));
    hNcoll[ix]->SetLineColor(ix + 1);

    cout << "   <Ncoll>: " << hNcoll[ix]->GetMean() << endl;

    //Convolve with NBD to get BBC charge
    for (int i = 1; i <= hNcoll[ix]->GetNbinsX(); i++)
    {
      double ncoll = hNcoll[ix]->GetBinCenter(i);
      double w = hNcoll[ix]->GetBinContent(i);

      fNBD->SetParameters(ncoll * NBD_mu, ncoll * NBD_k);

      for (int j = 1; j <= hBBCs[ix]->GetNbinsX(); j++)
      {
        double c = hBBCs[ix]->GetBinCenter(j);
        double bc = hBBCs[ix]->GetBinContent(j);

        double v = w * fNBD->Eval(c);

        hBBCs[ix]->SetBinContent(j, v + bc);

        // cout << " ncoll:" << ncoll << " BBCc:" << c
        //      << " w:" << w << " NBD:" << fNBD->Eval(c)
        //      << " v:" << v << " bc:" << bc
        //      << " new bc:" << v + bc
        //      << endl;
      }
    }


    // apply the trigger efficiency
    for (int j = 1; j <= hBBCs[ix]->GetNbinsX(); j++)
    {
      double c = hBBCs[ix]->GetBinCenter(j);
      double bc = hBBCs[ix]->GetBinContent(j);
      double eff = feff->Eval(c);
      hBBCs[ix]->SetBinContent(c, bc * eff);
    }

    // check the quantiles
    int nq = 2;
    double xq[] = {18./88., 83./88.};
    double yq[2] = {0};
    hBBCs[ix]->GetQuantiles(nq, yq, xq);
    cout << " quantiles: " << yq[0] << " " << yq[1] << endl;
    if (ix == 0)
    {
      BBC_periph = yq[0];
      BBC_cent = yq[1];
    }

    // Count the number of events
    Nevent_cent[ix] = hBBCs[ix]->Integral(hBBCs[ix]->FindBin(BBC_cent),
                                          hBBCs[ix]->GetNbinsX());
    Nevent_periph[ix] = hBBCs[ix]->Integral(1,
                                            hBBCs[ix]->FindBin(BBC_periph));





    // // unsigned int Nentries = ntp->GetEntries();
    // unsigned int Nentries = 10000;

    // for (unsigned int ientry = 0; ientry < Nentries; ientry++)
    // {

    //   ntp->GetEntry(ientry);

    //   double BBCs = 0;
    //   fNBD->SetParameters(Ncoll * NBD_mu, Ncoll * NBD_k);
    //   BBCs = fNBD->GetRandom();
    //   // for (int i = 0; i < (int)Ncoll; i++)
    //   //   BBCs += fNBD->GetRandom();

    //   hBBCs[ix]->Fill(BBCs);

    //   if (BBCs > BBC_cent)
    //     Nevent_cent[ix]++;
    //   if (BBCs < BBC_periph)
    //     Nevent_periph[ix]++;
    // }

    cout << "   N central    = " << Nevent_cent[ix] << endl;
    cout << "   N peripheral = " << Nevent_periph[ix] << endl;

    // delete htmp;
    delete ntp;
    fin->Close();
    delete fin;

  }

  //=====================================================//
  // CALCULATE RCP
  //=====================================================//
  cout << endl;
  cout << "--> Calculating RCP ..." << endl;

  grcp = new TGraph();
  grcp->SetTitle(";p_{T} [GeV/c];R_{CP}");
  grcp->SetLineStyle(1);
  grcp->SetLineColor(kBlue);
  // grcp->SetPoint(0, 0, 1);

  for (int ix = 0; ix < NX; ix++)
  {
    // x = pT / sqrt(s)
    double pT = x[ix] * 200;

    double rcp = Nevent_cent[ix] / dcent_cent;
    rcp /= Nevent_periph[ix] / dcent_periph;
    // rcp *= Ncoll_periph / Ncoll_cent;

    cout << " " << pT << " " << rcp << endl;

    grcp->SetPoint(ix + 1, pT, rcp);
  }

  //=====================================================//
  // PLOT OBJECTS
  //=====================================================//
  cout << endl;
  cout << "--> Plotting ..." << endl;

  //=====================================================//
  // PLOT
  //=====================================================//

  TCanvas *crcp = new TCanvas("crcp", "RCP", 800, 800);
  crcp->cd(1);
  grcp->Draw("AL");


  TCanvas *ctest = new TCanvas("ctest", "test", 1500, 500);
  ctest->Divide(2, 1);

  ctest->cd(1);
  gPad->SetLogy();
  hNcoll[0]->Draw();
  for (int ix = 1; ix < NX; ix++)
    hNcoll[ix]->Draw("same");

  ctest->cd(2);
  gPad->SetLogy();
  hBBCs[0]->Draw();
  for (int ix = 1; ix < NX; ix++)
    hBBCs[ix]->Draw("same");

  //=====================================================//
  // SAVE
  //=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving RCP to " << outFile << endl;

    TFile *fout = new TFile(outFile, "UPDATE");

    grcp->Write(Form("Rcp_%s", collSystem));

    fout->Close();
  }
}