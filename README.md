
 NAME: -procPlayers-

 R code for processing GSR time series data for use with SyncCalc (see
   https://academic.mu.edu/peressini/synccalc/synccalc.htm)

      version 2.0
      Tony Peressini & Stephen Guastello
      November 15, 2019

   Function: To process player GSR sync data from Guastello lab ... (ver 1.0-1.1)
         - utilizing revised nonlinear model from 2019 NSF Proposal (ver 2.0-)

   Input:  User chosen CSV file with time as first column (discarded) and
           monster data as final column (discarded)

   Output:  Six files named the original filename appended with following
     (1) -out.txt,   raw output of analyses for debugging/verification
     (2) -lMat.txt,  the linear sync matrix for use with SyncCalc
     (3) -nl1Mat.txt, the Opt 1 nonlinear sync matrix for use with SyncCalc
     (4) -nl2Mat.txt, the Opt 2 nonlinear sync matrix for use with SyncCalc
     (5) -nl3Mat.txt, the Opt 3 nonlinear sync matrix for use with SyncCalc
     (6) -sum.txt,   which contains a summary of the results of the analysis

     Version 2.0 employs different nonlinear model(s); see Sec. 3.1, p. 6-7
     from Peressini & Guastello (2019) NSF proposal description.

     The three nonlinear model options are:
       Option 1:   z_2 = A*e^{B*z_1} + e^{D*p_1}
       Option 2:   z_2 = A*p_1*z_1*(1−z_1)
       Option 3:   z_2 = A*p_1e^{Bz_1}

       The non-linear matrices are generated with the nonlinear auto-
       correlation's (z_2 = A*e^{Bz_1}) R (square root of R^2) on the diagonals
       (a[i,i]) and with the off-diagonals, a[i,j], i<>j, populated with:
             h[i,j] = sqrt[ |R'^2 - R^2| ],
       where R'^2 is R^2 for model (Option 1, 2, or 3), and R^2 is R^2 for the
       nonlinear autocorrelation (square of element on diagonal).
