# Wiser-fit
To install
`git clone https://github.com/JeffersonLab/Wiser-fit.git`

`cd Wiser-fit` 
`make`

> ratios

`E_beam,THETA_e,E_e,Q2,X,W2,PI-2E,K+2pi+,K-2pi-,DSIGX,SIGPIP,SIGPIM,SIGKAP,SIGKAM (last4 x-sections)



`!------------------------------------------------------------------------------`

` Calculate pi,K,p  cross section for electron beam on a proton target`

`! IntegrateQs over function WISER_FIT using integration routine QUADMO`

`! E0         is electron beam energy, OR max of Brem spectra`

`! P,E       is scattered particle  momentum,energy`

`! THETA_DEG  is kaon angle in degrees`

`! RAD_LEN (%)is the radiation length of target, including internal`

`!                (typically 5% for SLAC, 2% for JLab)`

`!               = .5 *(target radiation length in %) + 5. (+2 for JLab)`

`!       ***  =100. IF BREMSTRULUNG PHOTON BEAM OF 1 EQUIVIVENT QUANTA ***`

`! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar`

`! SIGMA      is output cross section in nanobars/GeV-str`

`!------------------------------------------------------------------------------`
 
 The code provided by Peter Bosted. The reference should be to Wiser's PhD thesis, as the data
were never published.  The actual fit and coding were done by Steve Rock and George
Chang in about 1982 or so. It was intended to be used by anyone who wanted it.

