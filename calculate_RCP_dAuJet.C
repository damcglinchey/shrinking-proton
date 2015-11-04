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

  // const int NX = 11;             // Number of x values
  // double x[] = {0.01, 0.05, 0.10, 0.15, 0.20,
  //               0.25, 0.30, 0.35, 0.40, 0.45,
  //               0.50
  //              };         // x values
  const int NX = 4;             // Number of x values
  double x[] = {0.01, 0.10, 0.30, 0.50};
  const char *xFiles[] = {    // Files for each x value
    "rootfiles/glauber_dau_snn42_x001_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x01_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x03_ntuple_100k.root",
    "rootfiles/glauber_dau_snn42_x05_ntuple_100k.root",
  };

  // Set the collision system. Current options are:
  // "dAu"   - d+Au 200 GeV
  const char *collSystem = "dAu";

  bool saveRcp = true;

  //=====================================================//
  // DECLARE VARIABLES
  //=====================================================//

  // collision system variables
  //http://www.phenix.bnl.gov/phenix/WWW/p/info/an/900/Run8_dAu_200GeV_Centrality_Categorization.pdf
  //http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1087/Run8_dAu_200GeV_Centrality_Addendum-01.pdf
  double NBD_mu = 3.038;
  double NBD_k = 0.464;
  double trig_eff_params[2] = {0.897, 0.612};

  TF1* feff = new TF1("feff", "1.0 - TMath::Exp(-1.*TMath::Power(x/[0], [1]))", 0, 200);
  feff->SetParameters(trig_eff_params[0], trig_eff_params[1]);

  // For running over ntuples
  TFile *fin;
  TTree *ntp;
  const char *ntpName = "nt_dh_Au";
  Float_t Ncoll;

  TH2D *hNcoll_NcollMod[NX];

  // calculated histograms
  TH2D *hNcoll_BBCsc[NX];
  TH2D *hNcoll_BBCscMod[NX];
  TH2D *hNcollMod_BBCscMod[NX];

  TH1D *hBBCsc[NX]; // unmodified BBC south charge for calculating cent limits

  TH1D* hBBCscNcoll[NX];
  TH1D* hBBCscModNcoll[NX];
  TH1D* hBBCscModNcollMod[NX];

  for (int i = 0; i < NX; i++)
  {

    hNcoll_BBCsc[i] = new TH2D(Form("hNcoll_BBCsc_%i", i),
                               ";BBCs charge; N_{coll}",
                               1596, 1, 400,
                               101, -0.5, 100.5);

    hNcoll_BBCscMod[i] = new TH2D(Form("hNcoll_BBCscMod_%i", i),
                                  ";BBCs charge Mod; N_{coll}",
                                  1596, 1, 400,
                                  101, -0.5, 100.5);

    hNcollMod_BBCscMod[i] = new TH2D(Form("hNcollMod_BBCscMod_%i", i),
                                     ";BBCs charge Mod; N_{coll}^{mod}",
                                     1596, 1, 400,
                                     101, -0.5, 100.5);

    // BBCs charge (unmodified)
    hBBCsc[i] = new TH1D(Form("hBBCsc_%i", i),
                         ";BBCs charge",
                         1596, 1, 400);


    // Ncoll weighted BBC charge distributions
    hBBCscNcoll[i] = new TH1D(Form("hBBCscNcoll_%i", i),
                              ";N_{coll} #time BBCs charge",
                              1596, 1, 400);

    hBBCscModNcoll[i] = new TH1D(Form("hBBCscModNcoll_%i", i),
                                 ";N_{coll} #time BBCs charge Mod",
                                 1596, 1, 400);

    hBBCscModNcollMod[i] = new TH1D(Form("hBBCscModNcollMod_%i", i),
                                    ";N_{coll}^{mod} #time BBCs charge Mod",
                                    1596, 1, 400);


  }

  TH1D* hNcoll_MB[NX];
  TH1D* hNcollMod_MB[NX];

  TH1D* hNcoll[NX][2];
  TH1D* hNcollMod[NX][2];

  double bias_BBCcsMod[NX][2];
  double bias_NcollModBBCcsMod[NX][2];

  double rcp_BBCscMod[NX];
  double rcp_NcollModBBCscMod[NX];

  TGraph *grcp_BBCscMod;
  TGraph *grcp_NcollModBBCscMod;

  //=====================================================//
  // PRINT RUNNING CONDITIONS
  //=====================================================//
  cout << endl;
  cout << "--> Running conditions ..." << endl;

  // print values
  cout << "  NBD : mu=" << NBD_mu << " k=" << NBD_k << endl;
  cout << "  Trig Eff Pars: " << trig_eff_params[0] << ", " << trig_eff_params[1] << endl;
  cout << "  Ntuple name: " << ntpName << endl;


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

    // we're done with the file, might as well close it
    delete ntp;
    fin->Close();
    delete fin;

    // project to get the MB Ncoll and NcollMod distributions
    hNcoll_MB[ix] = (TH1D*) hNcoll_NcollMod[ix]->ProjectionY(
                      Form("hNcoll_MB_%i", ix),
                      1, hNcoll_NcollMod[ix]->GetNbinsX());

    hNcollMod_MB[ix] = (TH1D*) hNcoll_NcollMod[ix]->ProjectionX(
                         Form("hNcollMod_MB_%i", ix),
                         1, hNcoll_NcollMod[ix]->GetNbinsY());

    cout << " <Ncoll>   : " << hNcoll_MB[ix]->GetMean() << endl;
    cout << " <NcollMod>: " << hNcollMod_MB[ix]->GetMean() << endl;

    // convolve with NBD to get BBCs charge & apply trigger efficiency
    for (int ibx = 1; ibx <= hNcoll_NcollMod[ix]->GetNbinsX(); ibx++)
    {
      for (int iby = 1; iby <= hNcoll_NcollMod[ix]->GetNbinsY(); iby++)
      {

        float w = hNcoll_NcollMod[ix]->GetBinContent(ibx, iby);
        if (w == 0) continue;

        float Ncoll = hNcoll_NcollMod[ix]->GetYaxis()->GetBinCenter(iby);
        float NcollMod = hNcoll_NcollMod[ix]->GetXaxis()->GetBinCenter(ibx);

        for (int j = 1; j <= hNcoll_BBCsc[ix]->GetNbinsX(); j++)
        {
          double c = hNcoll_BBCsc[ix]->GetXaxis()->GetBinCenter(j);
          double NBD = evalNBD(c, Ncoll * NBD_mu, Ncoll * NBD_k);
          double NBDMod = evalNBD(c, NcollMod * NBD_mu, NcollMod * NBD_k);
          double eff = feff->Eval(c);

          // Unmodified BBCs charge
          if (Ncoll > 0)
          {
            hNcoll_BBCsc[ix]->Fill(c, Ncoll, w * NBD * eff);

            hBBCsc[ix]->Fill(c, w * NBD * eff);

            hBBCscNcoll[ix]->Fill(c, Ncoll * w * NBD * eff);
          }
          // Modified BBCs charge
          if (NcollMod > 0)
          {
            hNcoll_BBCscMod[ix]->Fill(c, Ncoll, w * NBDMod * eff);
            hNcollMod_BBCscMod[ix]->Fill(c, NcollMod, w * NBDMod * eff);

            hBBCscModNcoll[ix]->Fill(c, Ncoll * w * NBDMod * eff);
            hBBCscModNcollMod[ix]->Fill(c, NcollMod * w * NBDMod * eff);
          }

        }// j
      } // iby
    } // ibx

    // Calculate the desired charge limits for each centrality bin
    int nq = 2;
    double xq[] = {(88. - 60.) / 88., (88. - 20.) / 88.}; // 60-88%, 0-20%
    double yq[2] = {0};
    hBBCsc[ix]->GetQuantiles(nq, yq, xq);
    cout << "  quantiles: " << yq[0] << " " << yq[1] << endl;

    // get bin limits for each centrality bin
    double blim[2][2];
    //0-20%
    blim[0][0] = hBBCsc[ix]->FindBin(yq[1]);
    blim[0][1] = hBBCsc[ix]->GetNbinsX();
    //60-88%
    blim[1][0] = 1;
    blim[1][1] = hBBCsc[ix]->FindBin(yq[0]);

    cout << "  blim[0]: [" << blim[0][0] << ", " << blim[0][1] << "]"
         << " -> [" << hBBCsc[ix]->GetBinLowEdge(blim[0][0]) << ","
         << hBBCsc[ix]->GetBinLowEdge(blim[0][1] + 1) << "]"
         << endl;
    cout << "  blim[1]: [" << blim[1][0] << ", " << blim[1][1] << "]"
         << " -> [" << hBBCsc[ix]->GetBinLowEdge(blim[1][0]) << ","
         << hBBCsc[ix]->GetBinLowEdge(blim[1][1] + 1) << "]"
         << endl;

    double yield[2] = {0};
    double yield_NcollMod[2] = {0};
    double yield_BBCscMod[2] = {0};
    double yield_NcollModBBCscMod[2] = {0};


    for (int iq = 0; iq < 2; iq++)
    {
      // get Ncoll distributions for each centrality bin
      hNcoll[ix][iq] = (TH1D*) hNcoll_BBCsc[ix]->ProjectionY(
                         Form("hNcoll_%i_%i", ix, iq),
                         blim[iq][0], blim[iq][1]);

      hNcollMod[ix][iq] = (TH1D*) hNcollMod_BBCscMod[ix]->ProjectionY(
                            Form("hNcollMod_%i_%i", ix, iq),
                            blim[iq][0], blim[iq][1]);

      // get "jet" yields
      yield[iq] = hBBCscNcoll[ix]->Integral(blim[iq][0], blim[iq][1]);
      yield_BBCscMod[iq] = hBBCscModNcoll[ix]->Integral(blim[iq][0], blim[iq][1]);
      yield_NcollModBBCscMod[iq] = hBBCscModNcollMod[ix]->Integral(blim[iq][0], blim[iq][1]);


      // Calculate the bias factor for each centrality bin
      bias_BBCcsMod[ix][iq] = yield_BBCscMod[iq] / yield[iq];

      bias_NcollModBBCcsMod[ix][iq] = yield_NcollModBBCscMod[iq] / yield[iq];

    } //iq

    // Calculate Rcp
    rcp_BBCscMod[ix] = bias_BBCcsMod[ix][0] / bias_BBCcsMod[ix][1];

    rcp_NcollModBBCscMod[ix] = bias_NcollModBBCcsMod[ix][0] / bias_NcollModBBCcsMod[ix][1];

    cout << endl;
    cout << "   <Ncoll(00-20%)>: " << hNcoll[ix][0]->GetMean() << " (vs 15.1)" << endl;
    cout << "   <Ncoll(60-88%)>: " << hNcoll[ix][1]->GetMean() << " (vs 3.2)" << endl;
    cout << "   Bias(00-20%) (BBCcs Mod): " << bias_BBCcsMod[ix][0] << endl;
    cout << "   Bias(00-20%) (both Mod) : " << bias_NcollModBBCcsMod[ix][0] << endl;
    cout << "   Bias(60-88%) (BBCcs Mod): " << bias_BBCcsMod[ix][1] << endl;
    cout << "   Bias(60-88%) (both Mod) : " << bias_NcollModBBCcsMod[ix][1] << endl;
    cout << "   Rcp (BBCcs Mod) : " << rcp_BBCscMod[ix] << endl;
    cout << "   Rcp (both Mod)  : " << rcp_NcollModBBCscMod[ix] << endl;


  }


  grcp_BBCscMod = new TGraph(NX, x, rcp_BBCscMod);
  grcp_BBCscMod->SetName("grcp_BBCscMod");
  grcp_BBCscMod->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
  grcp_BBCscMod->SetLineStyle(1);
  grcp_BBCscMod->SetLineWidth(2);
  grcp_BBCscMod->SetLineColor(kRed);

  grcp_NcollModBBCscMod = new TGraph(NX, x, rcp_NcollModBBCscMod);
  grcp_NcollModBBCscMod->SetName("grcp_NcollModBBCscMod");
  grcp_NcollModBBCscMod->SetTitle(";x (x=2 * p_{T} / #sqrt{s_{NN}});R_{CP}");
  grcp_NcollModBBCscMod->SetLineStyle(2);
  grcp_NcollModBBCscMod->SetLineWidth(2);
  grcp_NcollModBBCscMod->SetLineColor(kGreen + 2);




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

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(22);

  TLine l1;
  l1.SetLineColor(kBlack);
  l1.SetLineStyle(2);

  // set colors
  grcp_BBCscMod->SetLineColor(kRed);
  grcp_NcollModBBCscMod->SetLineColor(kGreen + 2);

  for (int ix = 0; ix < NX; ix++)
  {
    hBBCscNcoll[ix]->SetLineColor(kBlue);
    hBBCscModNcoll[ix]->SetLineColor(kRed);
    hBBCscModNcollMod[ix]->SetLineColor(kGreen+2);

    hBBCscNcoll[ix]->SetLineWidth(2);
    hBBCscModNcoll[ix]->SetLineWidth(2);
    hBBCscModNcollMod[ix]->SetLineWidth(2);
  }

  TH1F* haxis_rcp = new TH1F("haxis_rcp",
                             ";x (x_{jet}=2 * p_{T}^{jet} / #sqrt{s_{NN}});R_{CP} (0-20%/60-88%)",
                             100, 0, 1);
  haxis_rcp->SetMinimum(0);
  haxis_rcp->SetMaximum(1.1);
  haxis_rcp->GetYaxis()->SetTitleOffset(1.3);
  haxis_rcp->GetXaxis()->SetTitleOffset(1.3);

  TLegend *legRCP = new TLegend(0.15, 0.15, 0.5, 0.4);
  legRCP->SetFillStyle(0);
  legRCP->SetBorderSize(0);
  legRCP->AddEntry(gRCP_jet, "Run 8 Jet", "P");
  legRCP->AddEntry(grcp_BBCscMod, "Mod BBCs charge only", "L");
  legRCP->AddEntry(grcp_NcollModBBCscMod, "Mod BBCs chrg and N_{coll}", "L");



  int xplot = 2;

  TLegend *legBBC = new TLegend(0.5, 0.7, 0.96, 0.92,
                                Form("d+Au @ 200 GeV, x = %.2f",
                                     x[xplot]));
  legBBC->SetFillStyle(0);
  legBBC->SetBorderSize(0);
  legBBC->AddEntry(hBBCscNcoll[xplot], "Unmodified", "L");
  legBBC->AddEntry(hBBCscModNcoll[xplot], "Mod BBCs chg", "L");
  legBBC->AddEntry(hBBCscModNcollMod[xplot], "Mod BBCs chg and N_{coll}", "L");

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
  grcp_NcollModBBCscMod->Draw("C");
  legRCP->Draw("same");

  gRCP_jet->Draw("P");
  for (int i = 0; i < NJET; i++)
    bRcp_jet[i]->Draw();

  l1.DrawLine(0, 1, 0.6, 1);


  TCanvas *cbbc = new TCanvas("cbbc", "bbc", 800, 600);
  cbbc->Divide(1, 1);
  cbbc->GetPad(1)->SetTopMargin(0.06);
  cbbc->GetPad(1)->SetRightMargin(0.02);
  cbbc->GetPad(1)->SetBottomMargin(0.10);
  cbbc->GetPad(1)->SetLeftMargin(0.10);
  cbbc->GetPad(1)->SetTicks(1, 1);

  cbbc->cd(1);
  hBBCscNcoll[xplot]->SetTitle(";BBCs charge");

  gPad->SetLogy();
  hBBCscNcoll[xplot]->GetYaxis()->SetRangeUser(1, 2e5);
  hBBCscNcoll[xplot]->GetXaxis()->SetRangeUser(1, 200);
  hBBCscNcoll[xplot]->Draw();
  hBBCscModNcoll[xplot]->Draw("same");
  hBBCscModNcollMod[xplot]->Draw("same");
  legBBC->Draw("same");
  label.DrawLatex(0.5, 0.97, "BBCsc charge for events with high-p_{T} particle");


//=====================================================//
// SAVE
//=====================================================//
  if (saveRcp)
  {
    cout << endl;
    cout << "--> Saving plots" << endl;

    crcp->Print("Rcp_dAuJet.pdf");
    cbbc->Print("BBCs_dAuJet.pdf");
  }
}