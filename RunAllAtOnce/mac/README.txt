Description of my current setup to hopefully make your usage of full run script painless:


0) Get larlite, go to branch pi0reco.  Setup + build :
  > source config/setup_reco.sh
  > make -jN

1) Checkout David's hitremoval package, and move to tagged version v072017 :
   > git clone https://github.com/davidc1/HitRemoval
   > git checkout tags/v072017
   > make -jN

2) Checkout LArOpenCV (branch develop) and build
   > source $LARLITE_USERDEVDIR/LArOpenCV/setup_laropencv.sh
   > make -jN

3) Build the Filters package in this repository

4) Run script:
  > python both_cos_nu_removal_ahack.py name_of_sample /path/to/files.root
