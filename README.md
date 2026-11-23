# ProcPlayer (version 2.1)

R code for processing GSR time series data for use with SyncCalc (see <https://academic.mu.edu/peressini/synccalc/synccalc.htm>).

**Input:** User-chosen CSV file with TIME as the first column, which is ignored by the R code. Participants' time series (GSR or other) appear in the next columns. If the user has time series that might be used in some analyses, but not in others, the optional time series should appear in the last columns on the right. When the script is run, a file-chooser dialog opens so you can select the CSV file.

There are two lines of code that offer some flexibility in the analysis. One defines the intended lag length (scroll down the script file to see it): `LAG <- 1` [Line 59]. The numeral 1 indicates a lag of one unit. It can be changed to 2, 3, etc., as needed.

The other is `COLUMNS_TO_IGNORE <- 0` [Line 61]. The numeral 0 tells the program not to exclude any of the columns in the CSV file (other than the TIME column) from the analysis. If 0 is changed to numeral 1, the program will ignore the last column in the spreadsheet. If 0 is changed to 2, it will ignore the last two columns, etc.

**Data preparation:** Before any model is fit, each participant's series is rescaled (by the function `tConvert` in the script) to (max − x)/SD. This reverses the direction of the series, sets its minimum to zero, and gives it unit variance. All models below are fit to the rescaled series.

**Output:** Six files whose names are the same as the input file's name (up to the first period in the name) are produced as follows. They appear in the same directory folder as the CSV file.

> 1. `-out.txt`, complete output of statistical analyses for debugging and verification
> 2. `-lMat.txt`, the linear sync matrix for use with SyncCalc
> 3. `-nl1Mat.txt`, the Opt 1 nonlinear sync matrix for use with SyncCalc
> 4. `-nl2Mat.txt`, the Opt 2 nonlinear sync matrix for use with SyncCalc
> 5. `-nl3Mat.txt`, the Opt 3 nonlinear sync matrix for use with SyncCalc
> 6. `-sum.txt`, which contains a summary of the results of the analysis.

In every matrix, entry a[i,j] (row i, column j) is player i's influence on player j.

The linear model for the autocorrelations is: X<sub>n</sub> = β<sub>0</sub> + β<sub>1</sub>X<sub>n-LAG</sub>. The diagonal entries of the linear matrix are R, the square root of the R<sup>2</sup> of this model. Off-diagonal entries are the regression weights β<sub>2</sub> in the transfer equation: X<sub>n</sub> = β<sub>0</sub> + β<sub>1</sub>X<sub>n-LAG</sub> + β<sub>2</sub>P<sub>n-LAG</sub>.

The three nonlinear model options are (z<sub>1</sub> is the player's own lagged value and p<sub>1</sub> is the other player's lagged value):

> Option 1: z<sub>2</sub> = A・e<sup>B・z<sub>1</sub></sup> + e<sup>D・p<sub>1</sub></sup>.  
> Option 2: z<sub>2</sub> = A・p<sub>1</sub>・z<sub>1</sub>・(1−z<sub>1</sub>).  
> Option 3: z<sub>2</sub> = A・p<sub>1</sub>・z<sub>1</sub>・e<sup>B・z<sub>1</sub></sup> + C.

Options 1 and 3 are fit by nonlinear least squares (Levenberg-Marquardt, via `nlsLM`). Option 2 is fit as an ordinary linear regression (with an intercept) of z<sub>2</sub> on p<sub>1</sub>・z<sub>1</sub> and p<sub>1</sub>・z<sub>1</sub><sup>2</sup>, which is the expanded form of the model above with the two coefficients left free (reported as B and Q in the `-sum.txt` file).

The nonlinear matrices for Options 1 and 3 are generated with the nonlinear autocorrelation's (z<sub>2</sub> = A・e<sup>B・z<sub>1</sub></sup>) R (square root of R<sup>2</sup>) on the diagonals (a[i,i]). The off-diagonals, a[i,j], i<>j, are populated with:

> h<sub>i,j</sub> = sgn(R'<sup>2</sup> − R<sup>2</sup>)・√|R'<sup>2</sup> − R<sup>2</sup>|,

where R'<sup>2</sup> is R<sup>2</sup> for the particular model (Option 1 or 3) of player i's influence on player j, and R<sup>2</sup> is the R<sup>2</sup> of the nonlinear autocorrelation. NOTE: The sign of the expression within the absolute value function is retained outside the square root function, so h is negative when the model with player i's influence fits worse than the autocorrelation alone. (If the nonlinear autocorrelation's R<sup>2</sup> comes out negative, it is set to zero.)

The Option 2 matrix is built the same way, except that the *linear* autocorrelation is used in place of the nonlinear one: its diagonal entries are the square root of the linear autocorrelation R<sup>2</sup> (the same diagonal as the linear matrix), and R<sup>2</sup> in the formula above is the linear autocorrelation R<sup>2</sup>. This follows Guastello & Peressini (2021a): because the polynomial logistic map model is fit by multiple linear regression, its transfer effect is assessed against the linear autocorrelation, whereas the exponential models, which are fit by nonlinear regression, are assessed against the nonlinear autocorrelation.

The options correspond to the models in that article as follows: Option 1 is Eq. 7 (double exponential), Option 2 is Eq. 11 (polynomial logistic map), and Option 3 is Eq. 12 (exponential logistic map).


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

R programming language installation (32 bit or 64 bit) and time series data in a CSV file. The script also requires the following R packages, which can be installed once with:

```r
install.packages(c("foreign", "stargazer", "minpack.lm"))
```

### Installing

Simply download and run the script in R.

## Running the tests

Test data set included: EmergencyResponseTeam.csv

Output files from test data set:

> EmergencyResponseTeam-lMat.txt  
> EmergencyResponseTeam-nl1Mat.txt  
> EmergencyResponseTeam-nl2Mat.txt  
> EmergencyResponseTeam-nl3Mat.txt  
> EmergencyResponseTeam-out.txt  
> EmergencyResponseTeam-sum.txt

## Versioning

Version 2.0 employs three different nonlinear models in addition to the linear model. See this article for theoretical background:

> Guastello, S. J., & Peressini, A. F. (2021a). A comparison of four dyadic synchronization models. *Nonlinear Dynamics, Psychology, and Life Sciences, 25*(1), 19-39. <http://www.societyforchaostheory.org/ndpls>

Version 2.1 is a slightly better internally documented version of 2.0.

## Authors

* **Anthony F. Peressini** - <https://github.com/peressinia>
* **Stephen J. Guastello**

## License

Copyright (c) 2016-2026 Anthony F. Peressini and Stephen J. Guastello

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Acknowledgments

* Portions of the work leading up to this program and `SyncCalc 1.1` were funded by a Way-Klingler fellowship grant from Marquette University (SJG, PI).
* Authors thank David Marra for his assistance with debugging this script.
