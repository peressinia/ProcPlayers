# ProcPlayer

R code for processing GSR time series data for use with SyncCalc (see <https://academic.mu.edu/peressini/synccalc/synccalc.htm>).

**Input:**  User chosen CSV file with time as first column (discarded), player GSR timeseries as next columns, and monster/control data as final column (discarded).

**Output:**  Six files whose names are the same as the input file's name appended as follows:

> 1.	-out.txt, raw output of analyses for debugging/verification
> 2. 	-lMat.txt, the linear sync matrix for use with SyncCalc
> 3.	-nl1Mat.txt, the Opt 1 nonlinear sync matrix for use with SyncCalc
> 4.	-nl2Mat.txt, the Opt 2 nonlinear sync matrix for use with SyncCalc
> 5.	-nl3Mat.txt, the Opt 3 nonlinear sync matrix for use with SyncCalc
> 6.	-sum.txt, which contains a summary of the results of the analysis.


The three nonlinear model options are:

> Option 1:	z<sub>2</sub> = &alpha;・e<sup>&beta;・z<sub>1</sup> + e<sup>&delta;・p<sub>1</sub></sup>.  
> Option 2:	z<sub>2</sub> = &alpha;・p<sub>1</sub>・z<sub>1</sub>・(1−z<sub>1</sub>).  
  > Option 3:	z<sub>2</sub> = &alpha;・p<sub>1</sub>・z<sub>1</sub></sup>・e<sup>&beta;・z<sub>1</sub></sup> + C.  

The non-linear matrices are generated with the nonlinear autocorrelation's (z<sub>2</sub> = A・e<sup>&beta;・z<sub>1</sub></sup>) R (square root of R<sup>2</sup>) on the diagonals (a[i,i]) and with the off-diagonals, a[i,j], i<>j, populated with:

> h<sub>i,j</sub> = &radic;|R'<sup>2</sup> - R<sup>2</sup>|,

where R'<sup>2</sup> is R<sup>2</sup> for the particular model (Option 1, 2, or 3) of player i's influence on player j, and R<sup>2</sup> is R<sup>2</sup> for the nonlinear autocorrelation (square of element on diagonal).


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

R programming language installation and GSR data in CVS file.

### Installing

Simply download and run the script in R.

## Running the tests

Test data included:  grp3_session3_game2.cvs


## Versioning

Version 2.0 employs different nonlinear model(s); see Sec. 3.1, p. 6-7 from Peressini & Guastello (2019) NSF proposal description.


## Authors

* **Anthony F. Peressini** - <https://github.com/peressinia>
* **Stephen J. Guastello**



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone we should
* Inspiration
* etc
