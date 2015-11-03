PHOBOS MC Glauber NOTES - Started Oct 14 2015
=============================================
*Darren McGlinchey*

## Setup ##

Downloaded the source from
[hepforge](https://www.hepforge.org/downloads/tglaubermc)

See the reference documentation at
[arXiv:1408.2549](http://arxiv.org/pdf/1408.2549v2.pdf)

Note that the code does not compile with ROOT 6 out of the box. To avoid dealing with the compile issues (not clear to me), I also installed ROOT v5.34.34. Can access it by running the command `root5`. All subsequent ROOT sessions will also be ROOT 5.



Centrality Categorization NOTES - Started Oct 15 2015
=====================================================

## p+Au 200 GeV ##

Jamie has an initial centrality callibration explained here:
`https://www.phenix.bnl.gov/WWW/p/draft/nagle/PHENIX/nagle_run15_pau_update_05282015.pdf`

Best NBD parameters:

    mu = 3.14
    k  = 0.47

Trigger efficiency: 84% of inelastic c.s.

Efficiency parameters (`1.0-TMath::Exp(-pow((x/[0]), [1]))`) from Jamie:

    1.07552e+00, 6.02328e-01

When running with the Glauber, run with "p", "Au"

## d+Au 200 GeV ##

References: [AN980](http://www.phenix.bnl.gov/phenix/WWW/p/info/an/900/Run8_dAu_200GeV_Centrality_Categorization.pdf), [AN1087](http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1087/Run8_dAu_200GeV_Centrality_Addendum-01.pdf), [PPG160](http://journals.aps.org/prc/pdf/10.1103/PhysRevC.90.034902)

Best NBD parameters:

    mu = 3.04
    k  = 0.46

Trigger efficiency: 88% of inelastic c.s.

Efficiency parameters from my [thesis](http://www.phenix.bnl.gov/phenix/WWW/talk/archive/theses/2012/McGlinchey_Darren-McGlinchey_C_Dissertation_2012.pdf):

    0.897, 0.612

When running with the Glauber, run with "dh", "Au"

## 3He+Au 200 GeV ##

References: [AN1207](http://www.phenix.bnl.gov/phenix/WWW/p/info/an/1207/Run14_3HeAu_200GeV_Centrality_Categorization.pdf)

Best NBD parameters:

    mu = 2.91
    k  = 0.55

Trigger efficiency: 88% of inelastic c.s.

Efficiency parameters from Jamie:

    1.22134e+00, 5.10114e-01

When running with the Glauber, run with "He3", "Au"


Proton size fluctuation NOTES - Started Oct 15 2015
===================================================


## Modifying the nucleon-nucleon cross section ##

If you have a collision with a high-x parton in the projectile (p, d, He3), producing a high-pT particle at midrapidity, it's possible that the size of the nucleon is smaller. More of the nucleon's energy is contained in the high-x parton, thereby reducing the number of other sea quarks/gluons. This would cause a decrease in the nucleon-nucleon cross section (sigma_nn). 

Try modeling this in the glauber model as:

**Temporary**
    sigma_nn(x) = sigma_nn * e^{-8*x}

## Modifying the PHOBOS glauber to include this possibility ##


