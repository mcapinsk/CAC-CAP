# CAC-CAP
The proof of chaos in the Standard Family is based on the CAPD library:

CAPD main webpage: http://capd.ii.uj.edu.pl

Full documentation: http://capd.ii.uj.edu.pl/html/

CAPD requirements: http://capd.ii.uj.edu.pl/html/capd_requirements.html

The proof also requires the OpenMP:

https://www.openmp.org

Please make sure you have it installed.

## Quick guide on how to build the CAPD library required for the proof

Clone the repository:

    git clone https://github.com/mcapinsk/CAC-CAP.git
    
Enter the repository, create the build folder, configure the library and then build:

    cd CAC-CAP
    mkdir build
    cd build
    cmake ..
    make

The above commands will build only the CAPD library.

For a detailed description on how to build the library see

http://capd.ii.uj.edu.pl/html/capd_compilation.html

## Building and running the proofs

To compile execute:

    cd ../CAC-CAP
    make

### To run the proof of Theorem 1.2 call

    ./CAC-CAP 0

### To run the proof of Theorem 1.3 using an N x N mesh call

    ./CAC-CAP 1 N

with N chosen as an integer. To produce the plot from the paper we have called

    nohup ./CAC-CAP 1 1000 > thm_1_3.txt

This computation took over four and a half hours running on 48 parallel threads.

### To run the proof of Theorem 1.4 choose k from 0 to 2 and call

    ./CAC-CAP 2 k

For the result from the paper we have called

    nohup ./CAC-CAP 2 2 > thm_1_4.txt

For different k we can obtain different accuracy:

k=0 gives a short computation and should validate 0.443 of the area.

k=1 should result in under an hour long computation and should validate 0.918 of the area. 

k=2 should take an hour or two on a desktop computer and should validate 0.9975 of the area. (Our computation took under 12 minutes running on 48 threads.)


## Results

After downloading the code files, the results from the computation which we have performed on our clusted will be appear in the files

CAC-CAP/CAC-CAP/thm_1_2.txt 

CAC-CAP/CAC-CAP/thm_1_3.txt 

CAC-CAP/CAC-CAP/thm_1_4.txt

The data files for the plot for the NTSF are included in the folder 

CAC-CAP/CAC-CAP/thm_1_3  

The plot can be made by calling 

    load 'plotNTSF.txt'

from gnuplot.

The data files for the plot for the dissipative case are included in the folder 

CAC-CAP/CAC-CAP/thm_1_4

The plot can be made by calling 

    load 'plot.txt'

from gnuplot.

When you execute the program the results for Theorems 1.2, 1.3 and 1.4 will be written out in terminal. The gnuplot files in the folders:

CAC-CAP/CAC-CAP/thm_1_3 

CAC-CAP/CAC-CAP/thm_1_4

will also be updated to the newes executed result. 

## Suggested order of reading the files

- validatedTrajectory.h/cpp
These files contain a method for validating that a trajectory goes above or below a prescribed level. The validation is performed by means of the Krawczyk parallel shooting method.

- conservativeMap.h/cpp
These files contain a computer assisted proof of chaos for conservative maps. The routine validateChaosInConservativeMaps() provides a computer assisted proof of Theorem 1.2 from the paper.

- nonTwistSF.h/cpp
These files contain a computer assisted proof of Theorem 1.3 from the paper. 

- dissipativeMap.h/cpp
These files contain the proof of Theorem 1.4.

## Authorship

Only the contents of the folder CAC-CAP/CAC-CAP constitute the proof from the paper, and the files included there have been written by Maciej J. Capinski. The remaining files are part of the CAPD library. These have been created by the CAPD Group.


