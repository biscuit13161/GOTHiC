# GOTHiC++

GOTHiC++ is a C++ implementation of the Bioconductor R package [GOTHiC](https://doi.org/doi:10.18129/B9.bioc.GOTHiC "GOTHiC on Bioconductor")

This is Hi-C analysis software which uses a cumulative binomial test to detect interactions between distal genomic loci that have significantly more reads than expected by chance in Hi-C experiments \[Mifsud *et al*,2017\]. It takes mapped paired NGS reads as input and gives back the list of significant interactions for a given bin size in the genome.


**GOTHiC Authors**: Borbala Mifsud and Robert Sugar

**GOTHiC Citation**: Mifsud B, Sugar R (2020). GOTHiC: Binomial test for Hi-C data analysis. R package version 1.24.0.

**GOTHiC++ Authors**: Richard Thompson, Elodie Darbo and Borbala Mifsud

**Contact Email**: ithompson[at]hbku.edu.qa

### Singularity Image

Users who don't have the access or experience to compile on the commandline,GOTHiC++ is available as a pre-built Singularity image which can be run locally or on a HPC

```bash
singularity pull library://gothic2020/default/gothic
```

#### Using Singularity to run GOTHiC++

GOTHiC++ runs using the same options described below when being run from the Singularity image. 

```bash
singularity exec <filename.sif> gothic <options>
singularity exec <filename.sif> gothicomp <optons>
```
Running the above commands without any arguments issues will print the usage, detailing the commandline options. 

### Running Single Sample GOTHiC++

GOTHiC++ does not, currently, readin HiCUP bam/sam files directly. Bam/sam files from HiCUP can be converted to txt files using the hicupToTable.sh script found in the src subdirectory where gothic is built.

```bash
src/hicupToTable.sh <input.bam>
```
Specifying bam/sam files in the config file, will result in gothic automatically converting `<input.bam>` to `<input.bam.txt>` if it does not exist. This file is then loaded by gothic.

```bash
<path/to>/gothic <path/to/gothic.conf>
```

or

```bash
gothic -i <filename> -s <name> -d <filename> [-t \#] [-r \#] [-o <dir>] [-c (all|trans|cis)]
      [-A (single|comparative)] [--verbose] [--no_rem_diag]

options:
    -i <filename>         Input filename
      --input <filename>
    -s <name>             Sample name
      --sample <name>
    -d <filename>         Digest of Restriction Enzyme, as used by HiCUP
      --digest <filename>
    -t #                  Num of threads to run, defaults to 1
      --threads #
    -r #                  Resolution in bases for bining interactions, defaults to 10000
      --res #
    -o <dir>              Output directory, defaults to './'
      --output <dir>
    -c (all|trans|cis)    Filter for Cis or Trans interactions,defaults to 'all'
      --cistrans (all|trans|cis)
    -A <option>           Analysis type, either 'single' or 'comparative'
      --analysis <option>
    --verbose             Print verbose output during run
    --no_rem_diag         Flag to stop removal of diagonals during analysis
```



### Running Comparative GOTHiC++

In order to carry a comparative analysis, Samples must be individually run as described above for single samples, but with the "**Analysis: comparative**" option. This mode causes gothic to carry out the fragment identification, binning and frequency counting before outputing the interactions into a binary file (`<SampleName>.inter.bin`). These files are then used as input for gothicomp.

```bash
<path/to>/gothicomp <path/to/gothicomp.conf>
```
 or 
 
 ```bash
gothicomp -c <filename> -s <filename> -n <name> -d <filename> [-b <filename>] [-t \#] [-r \#]
        [-o <dir>] [-a \#] [-A (bh|ihw)] [-C (all|trans|cis)] [--norandom] [--verbose|--debug]
        
options:
    -c <filename>         Control input filename
      --control <filename>
    -s <filename>         Sample input filename
      --sample <filename>
    -n <name>             Sample name
      --sample <name>
    -d <filename>         Digest of Restriction Enzyme, as used by HiCUP
      --digest <filename>
    -b <filename>         Baits file for analysis, as used by HiCUP
      --baits <filename>
    -t #                  Num of threads to run, defaults to 1
      --threads #
    -r #                  Resolution in bases for binning interactions, defaults to 10000
      --res #               - Only effective for binning the baits file, Sample data is binned using gothic.
    -o <dir>              Output directory, defaults to './'
      --output <dir>
    -C (all|trans|cis)    Filter for Cis or Trans interactions,defaults to 'all'
      --cistrans (all|trans|cis)
    -A (bh|ihw)           Algorithm for p-value correction, either 'bh' or 'ihw'
      --analysis (bh|ihw)
    -a #                  Alpha cutoff for p-value correction, defaults to '0.1'
      --alpha #
    --norandom            Turn off Random subsampling
    --verbose             Print verbose output during run
    --debug               Print very verbose output during run
 ```

- When there is a greater than 10% difference in read count between the samples, GOTHiComp will randomly subset the larger sample to match the size of the smaller sample before calculating the P values. This behaviour can be turned off with `--norandom` on the CLI, or setting `RandomSubset: false` if using a config file. 
- The resolution option of gothicomp only affects the bin size of the baits file, changes to the sample bin size require re-running gothic.

### Notes on config files

<ul>
<li>GOTHiC and GOTHiCOMP can carry out analysis using only 
	<ul>
	<li> cis (defined as both on the same chromosome)</li> 
  <li>trans (defined as between chromosomes)</li>
  <li>all</li>
  </ul> 
The config default is to analyse all interactions, but can be changed by altering the "CisTrans" option
</li>

<li>By default, GOTHiC and GOTHiComp will remove diagnals; if this is not required, please uncomment the "#RemoveDiagonals: false" line</li>
</ul>

### Requirements for building from Source
+ Cmake 3.4
+ C99 and C++17-compatible compiler
+ Boost C++ libraries (for unit testing)
+ tbb
+ SamTools (for Bam to txt conversion)
+ Rscript (for independent hypothesis weighting)
+ IHW BioConductor package (for independent hypothesis weighting)

### Building and Installing

In order to compile GOTHiC++, the authors advise creating a build folder within the root:

```bash
mkdir build && cd $_
```

Compiling GOTHiC++ follows a typical cmake pipeline:

```bash
cmake ..
make
```

it is possible that cmake will not use the correct compilers if you have multiple installed, or if it is not installed in a default location. The correct compiler can be passed to cmake using:

```bash
cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc -DCMAKE_CXX_COMPILER=/usr/local/bin/g++ ..
```

### Testing compiled code

Once compiled, the compilation can be tested using:

```bash
make test
```


### Installing Boost

If the Boost C++ libraries are not installed, please find the details at [boost.org][https://www.boost.org/]

### Sourcing Intel TBB

```bash
git clone https://github.com/oneapi-src/oneTBB.git
cd oneTBB
gmake
```

In order to locate TBB, Cmake utilises a `FindTBB.cmake` from [https://github.com/justusc/FindTBB](https://github.com/justusc/FindTBB).
Please note, it may be necessary to add the `-DTBB_DIR=<path/to>/oneTBB` option to the cmake command.

### Independent Hypothesis Weighting for Comparative Analyses
 
GOTHiC++ carries out Pvalue correction using either Benjamini Hochberg (BH) or Independent Hypothesis Weighting (IHW) \[Ignatiadis *et al*, 2016\], which is specified in the config file. Independent hypothesis weighting utilises the BioConductor IHW package, described from [https://bioconductor.org/packages/release/bioc/html/IHW.html](https://bioconductor.org/packages/release/bioc/html/IHW.html).


**N.B.** By default GOTHiC++ carries out Q value correction via the IHW R Bioconductor package. If IHW identifies only 1 bin, IHW attempts to use BH correction, however IHW does not have access to all the information it needs to accurately calculate BH correction. In this situation GOTHiC++ will revert to it's own BH correction.

### References

+ Ignatiadis, N., Klaus, B., Zaugg, J. *et al* (2016). Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. Nat Methods 13, 577–580. [https://doi.org/10.1038/nmeth.3885](https://doi.org/10.1038/nmeth.3885)
+ Mifsud B, Martincorena I, Darbo E, Sugar R, Schoenfelder S, Fraser P, *et al* (2017). GOTHiC, a probabilistic model to resolve complex biases and to identify real interactions in Hi-C data. PLoS ONE 12(4): e0174744 [https://doi.org/10.1371/journal.pone.0174744](https://doi.org/10.1371/journal.pone.0174744)
