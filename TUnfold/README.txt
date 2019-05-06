Unfolding example for charged particle multiplicities
=====================================================

24 April 2019   This is software is still under development and may
  change or may not compile.

Compiling:
==========
  copy the following files:
     src/*
     TUnfold_V17.7/*
     Makefile
     setup.sh

  setup working root version and c compiler
     source setup.sh

  complile
     make fillHist unfoldHist drawResult


Running:
========
  copy the files:
     minitree.dtd
     minitree.xml
     binning.dtd
     binning_etay.xml
  edit the XML files
     minitree.xml -> change location of data and MC, luminosity
     binning_etay.xml  -> change binning scheme and data luminosity

  set up environment if not yet done:
    source setup.sh  

  fill input histograms (loops over minitrees):
    ./fillHist -b binning_etay.xml -m minitree.xml

  unfold histograms:
    ./unfoldHist binning_etay.xml

  extract multiplicity histograms per (eta,Q2,y) bins
  extract covariance matrix per (eta,Q2,y) bin
    ./drawResults binning_etay.xml 

TODO:
=====
unfoldHist: use data uncertainties to unfold one MC from another
drawResults: extract covariance matrices between mixed kinematic bins
     (eta,Q2,y) and (eta',Q2',y')
use better command line for unfoldHist and drawResults
implement systematic uncertainties

Root files:
===========
fillHist
   writes a file  "unfoldingInput_yetalab.root"
   see tag <UnfoldingInput> in the XML file

unfoldHist
   reads  "unfoldingInput_yetalab.root"
   writes "unfoldingOutput_yetalab.root"
   see tags <UnfoldingInput> and <UnfoldingOutput> in the XML file

drawResult
   reads "unfoldingOutput_yetalab.root"
   writes new file "results_"+"unfoldingOutput_yetalab.root"
   see tag <UnfoldingOutput> in the XML file, the string "results" is
   fixed

Contents of "unfoldingInput_yetalab.root"
   binning schemes
      recQ2  (Q2rec,yrec,ntrack) binning
      genQ2  (Q2gen,ygen,ncharged) binning
      covBinning (eta,ntrack) binning
   TDirectories:
      data, djangoh, rapgap
   Content of the TDirectory:
      hist_rec_*  one histo per (eta).
      hist_recCovar_*  one histo per (Q2rec,yrec).
    only for for MC:
      hist_gen_*  one histo per (eta). 
      hist_genRec_*  one histo per (eta). Migration matrix
      hist_fake_*    one histo per (eta). Background from fakes 

Contents of "unfoldingOutput_yetalab.root"
   binning schemes
      copied from  "unfoldingInput_yetalab.root"
   TDirectories: each combination of MC with (data|MC)
      data_from_djangoh   -> unfold data with matrix taken from django
      rapgap_from_djangoh
      data_from_rapgap
      djangoh_from_rapgap
    NOTE: rapgap_from_djangoh and djangoh_from_rapgap are useless,
      because the MC statistics and error matrix is different from
      data
   Content of each TDirectory:
      hist_recCovar_* : copied from "unfoldingInput_yetalab.root"
      hist_rec_* : copied from "unfoldingInput_yetalab.root" DATA
      hist_genGec_* : copied from "unfoldingInput_yetalab.root" MC
      hist_fakes_* : unfolding background (c.f. hist_fake_* MC)
      hist_bias_* : unfolding bias (c.f. hist_gen_* MC)
      hist_unfolded_* : unfolding result for each eta bin
      hist_ematrixTotal_* : unfolding result covariance for each eta bin
      hist_rhoij_* : correlation coefficients
      hist_dxdy_* : matrix for error propagation
      Lcurve_* : L-curve
      logTauX_* : L-curve scan of X wrt log(tau)
      logTauY_* : L-curve scan of Y wrt log(tau)
      best*_* : points used for final unfolding
      LogTauCurvature : curvature of L-curve; the maximum is selected.

Contents of "results_unfoldingOutput_yetalab.root":
   binning schemes
      copied from  "unfoldingInput_yetalab.root"
   TDirectories: each combination of MC with (data|MC)
      data_from_djangoh   -> unfold data with matrix taken from django
      rapgap_from_djangoh
      data_from_rapgap
      djangoh_from_rapgap
    NOTE: rapgap_from_djangoh and djangoh_from_rapgap are useless,
      because the MC statistics and error matrix is different from
      data
   Content of each TDirectory:
      hist_unfolded_* : unfolding result in proper binning
      hist_ematrix_* : unfolding error matrix in proper binning
      hist_rhoij_* : unfolding correlations in proper binning
      hist_bias_* : truth of the MC which was used for unfolding
