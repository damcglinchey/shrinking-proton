PHOBOS MC Glauber NOTES - Started Oct 14 2015
=============================================
*Darren McGlinchey*

## Setup ##

Downloaded the source from
`https://www.hepforge.org/downloads/tglaubermc`

See the reference documentation at
`http://arxiv.org/pdf/1408.2549v2.pdf`

Note that the code does not compile with ROOT 6 out of the box. To avoid dealing with the compile issues (not clear to me), I also installed ROOT v5.34.34. Can access it by running the command `root5`. All subsequent ROOT sessions will also be ROOT 5.



Centrality Categorization NOTES - Started Oct 15 2015
=====================================================

## p+Au 200 GeV ##

Jamie has an initial centrality callibration explained here:
`https://www.phenix.bnl.gov/WWW/p/draft/nagle/PHENIX/nagle_run15_pau_update_05282015.pdf`

Best NBD parameters:

    mu = 3.14
    k = 0.47

Trigger efficiency: 84% of inelastic c.s.

Looking at the centrality categorization we get approximately

    [0-5%]   BBCs charge > 50
    [0-10%]  BBCs charge > 38
    [0-20%]  BBCs charge > 27
    [60-84%] BBCs charge < 7.5
    [70-84%] BBCs charge < 4.5



Proton size fluctuation NOTES - Started Oct 15 2015
===================================================

## Relation between pT and x ##

    x = sqrt(M^2 + pT^2) / sqrt(s) * exp(-y)

for pions assume M=0, and y=0 so

    x = pT / sqrt(s) = pT / 200 


## Modifying the nucleon-nucleon cross section ##

If you have a collision with a high-x parton in the projectile (p, d, He3), producing a high-pT particle at midrapidity, it's possible that the size of the nucleon is smaller. More of the nucleon's energy is contained in the high-x parton, thereby reducing the number of other sea quarks/gluons. This would cause a decrease in the nucleon-nucleon cross section (sigma_nn). 

Try modeling this in the glauber model as:

**Temporary**
    sigma_nn(x) = sigma_nn * (1-x)

## Modifying the PHOBOS glauber to include this possibility ##


## Run a grid of pion pT values ##

    pion pT    x     sigma_nn(x)/sigma_nn
    5        0.025         0.975
    10       0.05          0.950
    15       0.075         0.925
    20       0.1           0.900
