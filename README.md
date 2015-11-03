README
======

This repository is for investigating x-dependent proton fluctuations using small systems with a phenomenological model.

## Requirements ##

To run the PHOBOS Glauber MC, a root installation is required with the MathMore package included. Out of the box, the PHOBOS MC doesn't compile with ROOT 6.02/12, but does compile fine with ROOT 5.34/34. 

Other macros should be compatable with either ROOT 5 or ROOT 6.


## PHOBOS Glauber MC ##

The PHOBOS Glauber MC was downloaded from [hepforge](https://tglaubermc.hepforge.org/). The supporting documentation can be found at [arXiv:1408.2549](http://arxiv.org/abs/1408.2549)

The `runAndSaveNtuple()` function has been modified to take an input `x` value. This modifies the interaction cross section as implemented in `TGlauberMC::CalcEvent()`

The output of the `runAndSaveNtuple()` function has also been modified and simplified. It now outputs a `TTree`, the branches are defined in the `TGlauberMC::Run()` function.

To run the Glauber MC for 1k p+Au collisions with x=0.1:

    $ root
    root [0] gSystem->Load("libMathMore")
    (int)0
    root [1] .L runglauber_v2.2_C+
    root [2] runAndSaveNtuple(1000, "p", "Au", 42, -1, 0.4, 0.1, "testglauber.root")
    Setting up nucleus p
    Setting up nucleus Au
    Event # 950 x-sect = 1.79035 +- 0.0580562 b        
    Done!
    root [3] 

## calculate_RCP_dAuJet.C ##