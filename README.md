# MASFit_2
MASFit (Milano Antineutrino Spectrum Fitter)

MASFit is a code intended to simulate the antineutrino spectrum that JUNO will receive and evaluate the sensitivity of the experiment with a CHI squared test. The code will produce a simulated spectrum and fit it with both normal and inverted models, computing the difference between the two CHI squared.

In this README is explained how to run the code and how to use the input file. If you want further explanation on how the code works read the file *"INFO.pdf"*.

You can use the code *"MASFit_2.py"* with Python running the shell *"runMASFit.sh"*.
In these shells are conteined all the input file in the correct order.

The file *"inputFlux.txt"* contains the un-obscillated flux of anti-neutrinos from reactor. It has to have the same number of elements of the Nbin in the main code (you can set it from the input file).  
The files *"bkg_[bkg name].txt"* contains the energy spectrum of the varius background impemented: geo neutrinos, accidental backgrounds, Li9, (alpha,n) and fast neutrons. They must have the same number of elements of the Nbin in the main code (you can set it from the input file), and the same energies.   
The file *"inputCores.txt"* contains datas on the distances and powers of the 9 cores used in the simulation. You can change these values but you can't add new cores (because you would need to add also their systematics). If you want to run the code using the ideal baseline (only one core located at 52.5 km) you can choose it from the input file.  
The file *"MASFit_func_2.py"* is called in the main code and contains some of the functions used. 

The code will produce two different outputs: if you use an Asimov data-set (Fluctuations=0) the output will be a plot of the simulated data with the two fits (called "MASFit_plot.png") and a file with all the free parameters of the fit and their reconsructed values (called "MASFit_parameters.txt"). If the option Fluctuations is True (1) it will produce an istogram with the distribution of the Delta Chi Squared, doing M fits (called "MASFit_CHI.png").

This is an explanaion of the input file, in which you can control everything of the simulation (names and values as to be separated form "tab"):

**Fluctuations?** if the answer is 0 it will produce an Asimov data-set, if the answer i 1 it will introduce stathistical fluctuation (Poisson) on the Asimov data-set and it will do multiple fit and produce a plot of the distribution of the chi squared.

**M** Number of fits you want to do with statystical fluctuations (1000 fits takes nearly 1 hour)  
**Nbin** is the number of bin in which you divide the istogram (suggested 410, to have bins of 20 kev)  
**Emin(MeV)** is the lower energy for the simulation (don't go under 1.8)  
**Emax(MeV)** is the maximun energy for the simulation  
**Ncont** is the number of events you want to simulate (6y=100000, 20y=330000 etc..)

The following entries are the phisical parameters for neutrino oscillation

**Sin2Theta12**  
**Sin2Theta13_NO**  
**Sin2Theta13_IO**  
**DeltaM21**  
**DeltaM31_NO**  
**DeltaM32_IO**

The following lines are parameters for the energy resolution and the sistematic uncertainties

**a(%)** First term of the energy resolution (a/sqrt(E))  
**b(%)** Second term of the energy resolution (the constant one)  
**c(%)** Third therm in the energy resolution (c/E)  
**sigma_a(%)** Uncertainties on the terms of the energy resolution  
**sigma_b(%)** ""  
**sigma_c(%)** ""  
**sigma_alphaC(%)** Correlated reactor uncertainty  
**sigma_alphaD(%)** Detector uncertainty  
**sigma_b2b(%)** Bin to bin uncertainty  
**sigma_alphaR(%)** Reactor uncorrelated uncertainty  
**sigma_B_geo(%)** Systematic on background rate for geo-neutrini  
**sigma_sh_B_geo(%)** Systematic on background shape for geo-neutrini  
**sigma_B_acc(%)**	Systematic on background rate for accidentals  
**sigma_sh_B_acc(%)**	Systematic on background shape for accidentals  
**sigma_B_alpha(%)**	Systematic on background rate for alpha,n  
**sigma_sh_B_alpha(%)**	Systematic on background shape for alpha,n  
**sigma_B_fn(%)**	Systematic on background rate for fast neutrons  
**sigma_sh_B_fn(%)**	Systematic on background shape for fast neutrons  
**sigma_B_Li9(%)**	Systematic on background rate for Li9  
**sigma_sh_B_Li9(%)**	Systematic on background shape for Li9
**sigma_B_WR(%)**	Systematic on background rate for World reactors  
**sigma_sh_B_WR(%)**	Systematic on background shape for world reactors    

**Scan?** If the answer is 1 it will produce a plot with the chi squared scan in function of Delta m^2(3l), otherwise nothing will be produced.  
**Corr?** If the answer is 1 it will produce a correlation matrix of the free parameters of the fit, otherwise nothing will be produced.  

Here you have the possibility to chose which parameter of the fit are free anch which are fixed. If you put 0 the parameter will be free, if you put 1 it will be fixed

**Fix_M21**  
**Fix_M3l**  
**Fix_Theta13**  
**Fix_Theta12**  
**Fix_N**  
**Fix_a**  
**Fix_b**  
**Fix_c**

**IBD_rate**  Expected rate for inverse beta decay events  
**Geo_nu_rate** Expected rate for geo-neutrinos  
**Acc_bkg_rate** Expected rate for accidentals background  
**Li9_bkg_rate**  Expected rate for Li9 background  
**fn_bkg_rate** Expected rate for fast neutrons background  
**alpha_n_bkg_rate**  Expected rate for alpha,n bacground  
**WR_rate** Expected rate for world reactors background  

**Ideal_bsl[km]**   If the answer is 0 the real baseline will be used (taken from inputCores.txt). If the answer is different from 0 the code will ignore the file with the cores and simulate only one core located at the distance you have setted here

**Plot?** if the answer is 1 it will pop up and block the terminal until you have closed it, if it is 0 it won't appear, but will be saved anyways.
