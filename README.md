# GOTHiC++

GOTHiC++ is a C++ implementation of the Bioconductor R package [GOTHiC](https://doi.org/doi:10.18129/B9.bioc.GOTHiC "GOTHiC on Bioconductor")

This is Hi-C analysis software which uses a cumulative binomial test to detect interactions between distal genomic loci that have significantly more reads than expected by chance in Hi-C experiments. It takes mapped paired NGS reads as input and gives back the list of significant interactions for a given bin size in the genome.

**GOTHiC Authors**: Borbala Mifsud and Robert Sugar

**GOTHiC Citation**: Mifsud B, Sugar R (2020). GOTHiC: Binomial test for Hi-C data analysis. R package version 1.24.0.

**GITHiC++ Authors**: Richard Thompson and Borbala Mifsud

**Contact Email**: ithompson[at]hbku.edu.qa

###Requirements
+ Cmake 3.4
+ C99 and C++17-compatible compiler
+ googletest

###Building and Installing

In order to compile GOTHiC++, the authors advise creating a build folder within the root:

```bash
mkdir build && cd build
```

Compiling GOTHiC++ follows a typical cmake pipeline:

```bash
cmake ..
make
```

it is possible that cmake will not use the correct compilers if you have multiple installes, or if it is not installed in a default location. The correct compiler can be passed to cmake using:

```bash
cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc -DCMAKE_CXX_COMPILER=/usr/local/bin/g++ ..
```

###Testing

Once compiled, the compilation can be tested using:

```bash
make test
```

###Running Single Sample GOTHiC++

```bash
<path/to>/gothic <path/to>gothic.conf
```

###Running Comparative GOTHiC++ 

In order to carry a comparative analysis, Samples must be individually run as described above for single samples, but with the "**Analysis: comparative**" option. This mode causes gothic to carry out the fragment identification, binning and frequency counting before outputing the interactions into a binary file. These files are then used as input for gothicomp.

```bash
<path/to>/gothicomp <path/to>gothicomp.conf
```

###Notes on config files


+ GOTHiC and GOTHiCOMP can carry out analysis using only cis (defined as on the same chromosome) or trans (defined as between chromosomes) or all. The config default is to analyse all interactions, but can be changed by altering the "CisTrans" option

+ By default, GOTHiC and GOTHiComp will remove diagnals. if this is not required, please uncomment the "#RemoveDiagonals: false" line

###Installing googletest

GOTHiC++ uses the googletest framework, this can be acquired using

```bash
git clone https://github.com/google/googletest.git
```
then compiled by... 

```bash
mkdir googletest/build && cd googletest/build
cmake ..
make
[sudo] make install
```

Google test uses c++11 standard; if make fails to compile googletest, try adding one of these commands before the target in googletest/CMakeLists.txt and re-running the cmake.

```
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
```

 