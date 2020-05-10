# bayesLC
Code for calculating uncertainties on satellite based land cover, as described in the following papers:

[Cripps, E., O’Hagan, A., & Quaife, T. (2013). Quantifying uncertainty in remotely sensed land cover maps. Stochastic Environmental Research and Risk Assessment, 27(5), 1239-1251.](https://link.springer.com/article/10.1007/s00477-012-0660-3)

[Quaife, T., & Cripps, E. (2016). Bayesian analysis of uncertainty in the GlobCover 2009 land cover product at climate model grid scale. Remote Sensing, 8(4), 314.](https://www.mdpi.com/2072-4292/8/4/314)

## Usage

usage: bayeslc \[options\] control_file.txt

The control file must be the last word on the command line.

Options are:

-d %f        set the correlation factor to %f<br>
-l %d        set the length for the conditional auto regressive model to %d<br>
-n %d        set the number of simulations to %d (the higher the better)<br>
-u           display this message<br>

The format of the control file is a series of ascii text lines each specifying user
defined input and output filenames. Anything after a # character is treated as a comment.
Blank lines, or lines consisting only of a comment are ignored. Excluding these
the following information is expected, in order:


Line 1: a list of N files containing the input LC proportions<br>
Line 2: a mask file containing 0 where a mask is to be applied and 1 elsewhere<br>
Line 3: a file containing an NxN confusion matrix<br>
Line 4: a file containing the number of "counts"<br>
Line 5: a list of N files to contain the output means<br>
Line 6: a list of N files to contain the output stdvs<br>

Aside from the confusion matrix (which is an NxN ascii matrix) all the input files are plain, headerless
ascii text files containing a single matrix, each of which must have the same dimensions as the rest.


## Install

I have only ever tried installing on Linux. The only dependency is the GNU Scientific Library (GSL). As long as you have the GSL and standard dev tools installed (specifically gcc and make) it should just be a case of typing `make`.


