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
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TBox.h>

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

void calculate_RCP_dAuJet()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //=====================================================//
  // SET RUNNING CONDITIONS
  //=====================================================//

  const int NX = 6;             // Number of x values
  double x[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};         // x values
  const char *xFiles[] = {    // Files for each x value
    "glauber_dau_snn42_x0_ntuple_100k.root",
    // "glauber_dau_snn42_x005_ntuple_100k.root",
    "glauber_dau_snn42_x01_ntuple_100k.root",
    // "glauber_dau_snn42_x0125_ntuple_100k.root",
    "glauber_dau_snn42_x02_ntuple_100k.root",
    "glauber_dau_snn42_x03_ntuple_100k.root",
    "glauber_dau_snn42_x04_ntuple_100k.root",
    "glauber_dau_snn42_x05_ntuple_100k.root",
  };

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
  double BBC_cent = 0;
  double BBC_periph = 0;
  double Ncoll_cent = 0;
  double Ncoll_periph = 0;
  double dcent_cent = 0;
  double dcent_periph = 0;
  double trig_eff_params[2] = {0};

  TF1 *fNBD = new TF1("fNBD", NBD, 0, 200, 2);

  TF1* feff = new TF1("feff", "1.0 - TMath::Exp(-1.*TMath::Power(x/[0], [1]))", 0, 200);

  // For running over ntuples
  TFile *fin;
  TNtuple *ntp;
  const char *ntpName = "";
  Float_t Ncoll;

  TH1D *hNcoll[NX];
  TH1D *hBBCs[NX];
  for (int i = 0; i < NX; i++)
  {
    // hBBCs[i] = new TH1D(Form("hBBCs_%i", i), ";BBCs charge", 201, -0.5, 200.5);
    hBBCs[i] = new TH1D(Form("hBBCs_%i", i), ";BBCs charge", 800, 0, 200);
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

  // dcent_cent   = 0.05; // 0-5%
  // dcent_periph = 0.30; // 70-100%
  dcent_cent   = 0.2; // 0-20%
  dcent_periph = 0.28; // 60-88%


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
    hNcoll[ix]->SetTitle(";N_{coll}");
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

    // find the bin limits (for x=0 only)
    if (ix == 0)
    {
      // Find the central count
      int nbins = hBBCs[ix]->GetNbinsX();
      double tot = hBBCs[ix]->Integral(1, nbins);
      for (int ibin = hBBCs[ix]->GetNbinsX(); ibin > 0; ibin--)
      {
        double integral = hBBCs[ix]->Integral(ibin, nbins);

        //0 - 5%
        if (integral / tot > 0.20)
        {
          Nevent_cent[ix] = integral;
          cout << "     cent=" << hBBCs[ix]->GetBinCenter(ibin)
               << " frac=" << integral / tot
               << endl;
          break;
        }
      }

      // Find the peripheral count
      for (int ibin = 1; ibin <= nbins; ibin++)
      {
        double integral = hBBCs[ix]->Integral(1, ibin);

        // 70-100%
        if (integral / tot > 0.40)
        {
          Nevent_periph[ix] = integral;
          cout << "     periph=" << hBBCs[ix]->GetBinCenter(ibin)
               << " frac=" << integral / tot
               << endl;
          break;
        }
      }

      cout << "     Rcp= (" << Nevent_cent[ix] << " / " << Nevent_periph[ix] << ")"
           << " " << (Nevent_cent[ix] / 0.20) / (Nevent_periph[ix] / 0.40)
           << endl;

    }

    // apply the trigger efficiency
    for (int j = 1; j <= hBBCs[ix]->GetNbinsX(); j++)
    {
      double c = hBBCs[ix]->GetBinCenter(j);
      double bc = hBBCs[ix]->GetBinContent(j);
      double eff = feff->Eval(c);
      hBBCs[ix]->SetBinContent(j, bc * eff);

      // if (ix == 0)
      // {
      //   cout << "   " << j << " " << c << " " << bc << " " << eff << " " << bc*eff << endl;
      // }
    }

    // check the quantiles
    int nq = 2;
    // double xq[] = {18. / 88., 83. / 88.};
    // double xq[] = {0.30, 0.95};
    // double xq[] = {0.40, 0.80}; // 60-100%, 0-20%
    double xq[] = {(88.-60.)/88., (88.-20.)/88.}; // 60-88%, 0-20%
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
  grcp->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
  grcp->SetLineStyle(1);
  grcp->SetLineColor(kBlue);
  // grcp->SetPoint(0, 0, 1);

  for (int ix = 0; ix < NX; ix++)
  {
    // x = pT / sqrt(s)
    double pT = 0.5 * x[ix] * 200;

    double rcp = Nevent_cent[ix] / dcent_cent;
    rcp /= Nevent_periph[ix] / dcent_periph;
    // rcp *= Ncoll_periph / Ncoll_cent;

    cout << " " << x[ix] << " " << pT << " " << rcp << endl;

    grcp->SetPoint(ix + 1, x[ix], rcp);
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

  TLegend *legNcoll;
  if (collSystem == "pAu")
    legNcoll = new TLegend(0.5, 0.5, 0.9, 0.9);
  else
    legNcoll = new TLegend(0.15, 0.15, 0.5, 0.6);
  legNcoll->SetFillStyle(0);
  legNcoll->SetBorderSize(0);
  for (int ix = 0; ix < NX; ix++)
  {
    legNcoll->AddEntry(hNcoll[ix],
                       Form("x=%.2f <N_{coll}>=%.2f", x[ix], hNcoll[ix]->GetMean()),
                       "L");
  }

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);

  TH1F* haxis_rcp = new TH1F("haxis_rcp",
                             ";x (x_{jet}=2 * p_{T}^{jet} / #sqrt{s_{NN}});R_{CP} (0-20%/60-88%)",
                             100, 0, 1);
  haxis_rcp->SetMinimum(0);
  haxis_rcp->SetMaximum(1.1);

  //=====================================================//
  // PLOT
  //=====================================================//

  TCanvas *crcp = new TCanvas("crcp", "RCP", 800, 800);
  crcp->cd(1);
  haxis_rcp->GetXaxis()->SetRangeUser(0, 0.6);
  haxis_rcp->Draw();
  grcp->Draw("LS");

  gRCP_jet->Draw("P");
  for (int i = 0; i < NJET; i++)
    bRcp_jet[i]->Draw();


  TCanvas *ctest = new TCanvas("ctest", "test", 1500, 500);
  ctest->Divide(2, 1);

  ctest->cd(1);
  gPad->SetLogy();
  hNcoll[0]->Draw();
  for (int ix = 1; ix < NX; ix++)
    hNcoll[ix]->Draw("same");
  legNcoll->Draw("same");
  label.DrawLatex(0.5, 0.95, collSystem);

  ctest->cd(2);
  gPad->SetLogy();
  hBBCs[0]->Draw();
  for (int ix = 1; ix < NX; ix++)
    hBBCs[ix]->Draw("same");
  label.DrawLatex(0.5, 0.95, collSystem);

  //=====================================================//
  // SAVE
  //=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving plots" << endl;

    ctest->Print("NcollBBCcharge_dAuJet.pdf");
    crcp->Print("Rcp_dAuJet.pdf");
  }
}