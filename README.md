# CAC-CAP
The proof of chaos in the Standarc Family. The proof is based on the CAPD library:

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

For detailed description on how to build the library see

http://capd.ii.uj.edu.pl/html/capd_compilation.html

## Building and running the Arnold diffusion proof

To compile and run the full proof of Arnold diffusion execute:

    cd ../CAC-CAP
    make
    ./CAC-CAP


## Results

During the computation the results are stored in the folder:

CAC-CAP/CAC-CAP/results
CAC-CAP/CAC-CAP/plots

## Suggested order of reading the files

- 
-

## Authorship

Only the contents of the folder CAC-CAP/CAC-CAP constitute the proof of Arnold diffusion in the 3bp, and the files included there have been written by Maciej J. Capinski. The remaining files are part of the CAPD library. These have been created by the CAPD Group.


