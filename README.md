# Wiser-fit
To install
`git clone https://github.com/JeffersonLab/Wiser-fit.git`

`cd Wiser-fit` 
`make`



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

