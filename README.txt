Installation:
============
QMEAN needs the following external programs to be installed:

- dssp (mandatory)		either the "dssp"-executable needs to be in the PATH
  variable or the environment variable DSSP_EXECUTABLE must be defined
- R (for Z-score plots)		either the "R"-executable needs to be in the PATH
  variable or the environment variable R_EXECUTABLE must be defined 
  (Remark: in order to generate plots w/o X11, install the Cairo-package (in R:
  install.packages("Cairo"))
- PSIPRED and SSpro (see warning below); for SSpro use the following script: $SSPR_DIR/bin/predict_ss_sa.sh

Warning:
The QMEAN score is a linear combination of 6 terms: 4 statistical potentials and
2 agreement terms  (predicted vs calculated secondary structure by PSIPRED and
solvent accessibility by SSpro/ACCpro)
It is strongly recommended to provide a PSIPRED prediction file (horiz-format) 
and a SSpro prediction file (sspro_acc-format) of the model sequence for
reaching the full performance.


Example:
=======
The qmean executable can be found in the bin-subdirectory. 
The $QMEAN_ROOT/../data/examples have to be run from the qmean directory:

  0) to show options:
  ./bin/qmean -h

  1) to calculate QMEAN including Z-score plots (in R: without X11 install Cairo package in R and use option -Q 2)
  -> compare to directory "$QMEAN_ROOT/../data/examples/out_for_comparison/EXAMPLE1":
  $QMEAN_ROOT/bin/qmean -p $QMEAN_ROOT/../data/examples/example_data/BAKER-ROBETTA_CASP8_T0485_TS1 -S $QMEAN_ROOT/../data/examples/example_data/T0485.horiz -A $QMEAN_ROOT/../data/examples/example_data/T0485.sspro_acc -o $QMEAN_ROOT/../data/examples/out

  2) also calculate local scores and generate "energy profile":
  -> compare to directory "$QMEAN_ROOT/../data/examples/out_for_comparison/EXAMPLE2":
  $QMEAN_ROOT/bin/qmean -p $QMEAN_ROOT/../data/examples/example_data/BAKER-ROBETTA_CASP8_T0485_TS1 -S $QMEAN_ROOT/../data/examples/example_data/T0485.horiz -A $QMEAN_ROOT/../data/examples/example_data/T0485.sspro_acc -o $QMEAN_ROOT/../data/examples/out -L -P

  3) w/o Z-score plots:
  -> compare to directory "$QMEAN_ROOT/../data/examples/out_for_comparison/EXAMPLE3":
  $QMEAN_ROOT/bin/qmean -p $QMEAN_ROOT/../data/examples/example_data/BAKER-ROBETTA_CASP8_T0485_TS1 -S $QMEAN_ROOT/../data/examples/example_data/T0485.horiz -A $QMEAN_ROOT/../data/examples/example_data/T0485.sspro_acc -o $QMEAN_ROOT/../data/examples/out -n

  4) calculate scores for entire directory of structures
  -> compare to directory "$QMEAN_ROOT/../data/examples/out_for_comparison/EXAMPLE4":
  $QMEAN_ROOT/bin/qmean -d $QMEAN_ROOT/../data/examples/example_data/multiple_structures -S $QMEAN_ROOT/../data/examples/example_data/T0485.horiz -A $QMEAN_ROOT/../data/examples/example_data/T0485.sspro_acc -o $QMEAN_ROOT/../data/examples/out_multi -I my_example

  5) in analogy to 4) but merged output files:
  -> compare to directory "$QMEAN_ROOT/../data/examples/out_for_comparison/EXAMPLE5":
  $QMEAN_ROOT/bin/qmean -d $QMEAN_ROOT/../data/examples/example_data/multiple_structures -S $QMEAN_ROOT/../data/examples/example_data/T0485.horiz -A $QMEAN_ROOT/../data/examples/example_data/T0485.sspro_acc -o $QMEAN_ROOT/../data/examples/out_multi -I my_example -m -L


References:
==========
(1) QMEAN
Benkert, P., Tosatto, S.C.E. and Schomburg, D. (2008). 
"QMEAN: A comprehensive scoring function for model quality assessment." 
Proteins: Structure, Function, and Bioinformatics, 71(1):261-277. 

(2) OpenStruture
www.openstructure.org 


Documentation:
==============
For further documentation visit the QMEAN server help page: 
http://swissmodel.expasy.org/qmean/cgi/index.cgi?page=help


Contact: 
========
Dr. Pascal Benkert
Swiss Institute of Bioinformatics
Biozentrum, University of Basel
Klingelbergstrasse 50/70
CH-4056 Basel / Switzerland
pascal.benkert@unibas.ch
Tel. +41 61 267 15 80



