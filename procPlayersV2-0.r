################################################################################
#
# NAME: -procPlayers-
#
# R code for processing GSR time series data for use with SyncCalc (see
#   https://academic.mu.edu/peressini/synccalc/synccalc.htm)
#
#  		version 2.0
#  		Tony Peressini & Stephen Guastello
#  		November 15, 2019
#
# 	Function: To process player GSR sync data from Guastello lab ... (ver 1.0-1.1)
# 				- utilizing revised nonlinear model from 2019 NSF Proposal (ver 2.0-)
#
#   Input:  User chosen CSV file with time as first column (discarded) and
# 			    monster data as final column (discarded)
#
# 	Output:	Six files named the original filename appended with following
# 		(1) -out.txt,   raw output of analyses for debugging/verification
# 		(2) -lMat.txt,  the linear sync matrix for use with SyncCalc
# 		(3) -nl1Mat.txt, the Opt 1 nonlinear sync matrix for use with SyncCalc
#     		(4) -nl2Mat.txt, the Opt 2 nonlinear sync matrix for use with SyncCalc
#     		(5) -nl3Mat.txt, the Opt 3 nonlinear sync matrix for use with SyncCalc
# 		(6) -sum.txt,   which contains a summary of the results of the analysis
#
#     Version 2.0 employs different nonlinear model(s); see Sec. 3.1, p. 6-7
#     from Peressini & Guastello (2019) NSF proposal description.
#
#     The three nonlinear model options are:
#       Option 1:   z_2 = A*e^{B*z_1} + e^{D*p_1}
#       Option 2:   z_2 = A*p_1*z_1*(1−z_1)
#       Option 3:   z_2 = A*p_1*z_1*e^{Bz_1} + C
#
#       The non-linear matrices are generated with the nonlinear auto-
#       correlation's (z_2 = A*e^{Bz_1}) R (square root of R^2) on the diagonals
#       (a[i,i]) and with the off-diagonals, a[i,j], i<>j, populated with:
#             h[i,j] = sqrt[ |R'^2 - R^2| ],
#       where R'^2 is R^2 for model (Option 1, 2, or 3), and R^2 is R^2 for the
#       nonlinear autocorrelation (square of element on diagonal).   NOTE: The sign 
#	of the expression within the absolute value function is retained outside the
#	square root function.
#
################################################################################

# ######################
# Library Packages Used
# ######################
 library(foreign)								# to import `foreign' data files
 library(stargazer)								# for nice output
 library(minpack.lm)								# nlsLM routines (Levenberg-Marquardt)


# ############
# Constants
# ############
LAG <- 1                        # for the time series lag
PLEVEL <- 0.05                  # p level for Delta coef

# ######################################################################################
# Function to convert a column in dataframe x to `t-statistic' by centering and scaling
# ######################################################################################
 tConvert <- function (x){
	   ((x*-1000)-min(x*-1000))/sd(x*-1000)
 }
 #######################################################################################


 # ######################################################################################
 # Function to output results of nonLinear model
 # ######################################################################################
  nlmOut <- function (x){							                  # x is the result of nls call
	  lsep<-"---------------------------------------"			# r-squared was stored in x$control for convenience
	  outX<-capture.output(summary(x))
	  outXf<-c(lsep,outX[2],lsep,outX[3:11],sprintf("R-Squared for fit: %6.4f",x$control[1]),outX[12:length(outX)],lsep)
	  return(outXf)
  }
  #######################################################################################




#######################################################################################
# Get data file and initialize filenames and data structures
#######################################################################################
 dFile=file.choose()
 thefName<-basename(dFile)
 thePath<- dirname(dFile)
 outName<-unlist(strsplit(basename(dFile),"[.]"))[1]
 outFileName<-paste(outName, "-out.txt", sep = "")
 lMatrixFileName <- paste(outName, "-lMat.txt", sep = "")
 nlMatrixFileName1<- paste(outName, "-nl1Mat.txt", sep = "")
 nlMatrixFileName2<- paste(outName, "-nl2Mat.txt", sep = "")
 nlMatrixFileName3<- paste(outName, "-nl3Mat.txt", sep = "")
 summaryFileName<- paste(outName, "-sum.txt", sep = "")

 ds = read.csv(dFile)
 # head(ds)									# echo head of data to console for debugging

 dsz <- ds[ -c(1) ]								# remove first column - time
 N<-(dim(dsz)[2])-1								# set the number of players - ignoring the last monster column

 sink(paste(thePath, outFileName,sep="/"),split=FALSE)   			# to append use:  "sink(outFileName,append=TRUE)"
 cat("ProcPlayers (v2.0) output for PP run...\n")
 cat(sprintf("data file is: %s\n", dFile))
 cat(sprintf("data set has %d variables and %d rows.\n", dim(ds)[2],dim(ds)[1]))
 cat(sprintf("Lag = %d\n",LAG))
 cat("\n\n\n")

 invisible( dsz[1:dim(dsz)[2]] <- lapply(dsz[1:dim(dsz)[2]], tConvert) )  	# convert to T-stat each of remaining columns
 # head(dsz)									# echo head of data to out file


 #################
 # create matrices for sync data
 #################
 syncMatrixL<-matrix(rep(0.0,N*N),nrow=N,ncol=N)
 syncMatrixNL1<-matrix(rep(0.0,N*N),nrow=N,ncol=N)
 syncMatrixNL2<-matrix(rep(0.0,N*N),nrow=N,ncol=N)
 syncMatrixNL3<-matrix(rep(0.0,N*N),nrow=N,ncol=N)

 #################
 # create list to store summary statistics
 #################
 lAC<-c(rsq=0.0)					# r^2
 lPx<-c(rsq=0.0,b1=0.0,b2=0.0)   	# r^2 b1 b2
 nlAC<-c(rsq=0.0,B=0.0,A=0.0)				# r^2 B
 nlP1<-c(rsq=0.0,B=0.0,D=0.0,A=0.0)  # r^2 B D
 nlP2<-c(rsq=0.0,bi=0.0,qi=0.0)  # r^2 B D
 nlP3<-c(rsq=0.0,A=0.0,B=0.0,C=0.0)
 BlankPlayerData <-  list(LinearAutoCorr=lAC,LinearSync=rep(list(lPx),N),NonLinearAuto=nlAC,NonLinearSync1=rep(list(nlP1),N),
 	 						NonLinearSync2=rep(list(nlP2),N),NonLinearSync3=rep(list(nlP3),N))
 pList<- rep(list(BlankPlayerData), N)

 #################
 # do calculations
 #################
  for (i in 1:N ) {
	 cat(sprintf("=====================================================================|\n"))
	 cat(sprintf("          Analyzing Player %d as independent variable (I=%d)\n", i,i))
 	 cat(sprintf("=====================================================================|\n"))

	 #  Linear Autoregression
	 cat(sprintf("... Player %d LINEAR autocorrelation (I=%d)\n", i,i))
	 ded=ts.intersect(pI=ts(dsz[i]), pIl=lag( ts(dsz[i]) ,-LAG) )
	 theFit<-lm( pI ~ pIl, data=ded, na.action=NULL )
	 pList[[i]]$LinearAutoCorr["rsq"]<-summary(theFit)$r.squared

	 lPivotR2=summary(theFit)$r.squared
	 syncMatrixL[i,i]<-sqrt(lPivotR2)
	 syncMatrixNL2[i,i]<-sqrt(lPivotR2)

	 stargazer(theFit, type = "text")
	 cat("---------------------------------------------------------------------|\n\n")

	 #  nonLinear autoRegression Pi on lag(Pi)
	 cat(sprintf("... Player %d NONLINEAR autocorrelation (I=%d)\n", i,i))
	 pI<-dsz[[i]]
	 pIl<- c(rep(NA,times=LAG),dsz[[i]][1:(length(dsz[[i]])-LAG)])
	 theFit <- nlsLM(pI ~ alpha*exp(beta * pIl), start=c(alpha=0.5,beta=0.5), na.action = na.exclude )
	 rSq <- 1-(deviance(theFit)/sum((pI-mean(pI))^2))   # R^2 = 1 - [ (Residual Sum of Squares / Corrected Sum of Squares) ]
 	 pList[[i]]$NonLinearAuto["rsq"]<-rSq
   pList[[i]]$NonLinearAuto["B"]<-coef(theFit)["beta"]
   pList[[i]]$NonLinearAuto["A"]<-coef(theFit)["alpha"]
   theFit$control[1]<-rSq								# store r-squared in theFit$control for printing convenience in nlmOut

   if (rSq<0.0) rSq<-0.0								# Make R^2 zero if it is less than zero for Sync matrix
   nlPivotR2=rSq                        # Save autoRegression for h (aitch) calculations

   	 syncMatrixNL1[i,i]<-sqrt(rSq)
	 syncMatrixNL3[i,i]<-sqrt(rSq)

	 writeLines(nlmOut(theFit))
	 cat("---------------------------------------------------------------------|\n\n")

	 for (j in 1:N ) {

		 if (i==j)	{} 	# DO nothing - already handled the "diagonal" outside j-loop
		 else	{
			#  Linear Regression Pi on lag(Pi) + lag(Pj)
			cat(sprintf("... Player %d LINEAR regression on Player %d (I=%d, J=%d)\n", i,j,i,j))
		 	ded=ts.intersect(pI=ts(dsz[i]), pIl=lag( ts(dsz[i]) ,-LAG), pJl=lag( ts(dsz[j]) ,-LAG) )
			theFit<-lm( pI ~ pIl+pJl, data=ded, na.action=NULL )
		 	pList[[i]]$LinearSync[[j]]["rsq"]<-summary(theFit)["r.squared"]
		 	pList[[i]]$LinearSync[[j]]["b1"]<-coef(theFit)["pIl"]
		 	pList[[i]]$LinearSync[[j]]["b2"]<-coef(theFit)["pJl"]
			syncMatrixL[j,i]<-coef(theFit)["pJl"]
			stargazer(theFit, type = "text")
			cat("---------------------------------------------------------------------|\n\n")

      #
			#  Opt 1: nonLinear Regression:  z_2 = A*e^{B*z_1} + e^{D*p_1}
			#
			cat(sprintf("... Player %d NONLINEAR regression Opt 1 on Player %d (I=%d, J=%d)\n", i,j,i,j))
			pI<-dsz[[i]]
			pIl<- c(rep(NA,times=LAG),dsz[[i]][1:(length(dsz[[i]])-LAG)])
			pJl<-c(rep(NA,times=LAG),dsz[[j]][1:(length(dsz[[j]])-LAG)])
			theFit <- nlsLM(pI ~ alpha*exp(beta * pIl) + exp(delt * pJl), start=c(alpha=0.5,beta=0.5,delt=0.5), na.action = na.exclude )
		 	rSq <- 1-(deviance(theFit)/sum((pI-mean(pI))^2))   	# R^2 = 1 - [ (Residual Sum of Squares / Corrected Sum of Squares) ]
		 	theFit$control[1]<-rSq								# store r-squared in theFit$control for printing convenience in nlmOut

    	#if (summary(theFit)[["coefficients"]]["delt","Pr(>|t|)"]<PLEVEL) syncMatrixNL1[j,i]<-coef(theFit)["delt"] else syncMatrixNL1[j,i]<-0.0

      wun<-ifelse(rSq<nlPivotR2,-1,1)             # preserve sign
      aitch<-sqrt( abs(rSq-nlPivotR2) )*wun       # h = sqrt( this r^2 - auto corr r^2 of this row[i] NOTE: *wun preserves sign
      syncMatrixNL1[j,i]<-aitch
      pList[[i]]$NonLinearSync1[[j]]["rsq"]<-rSq
		 	pList[[i]]$NonLinearSync1[[j]]["B"]<-coef(theFit)["beta"]
		 	pList[[i]]$NonLinearSync1[[j]]["D"]<-coef(theFit)["delt"]		# it gets 0.0 if insignificant
      pList[[i]]$NonLinearSync1[[j]]["A"]<-coef(theFit)["alpha"]
			writeLines(nlmOut(theFit))
      cat(sprintf("h^2 = (R^2 model - R^2 NLAC) = (%6.4f - %6.4f) = %6.4f\n",rSq,nlPivotR2,aitch^2))
      cat(sprintf("h = %6.4f\n",aitch))
			cat("---------------------------------------------------------------------|\n\n")


      #
			#  Opt 2: nonLinear Regression:  X_2 =P_1*X_1*(1−X_1)
			#
		  cat(sprintf("... Player %d NONLINEAR regression Opt 2 on Player %d (I=%d, J=%d)\n", i,j,i,j))
      pIl<- c(rep(NA,times=LAG),dsz[[i]][1:(length(dsz[[i]])-LAG)])
      pJl<-c(rep(NA,times=LAG),dsz[[j]][1:(length(dsz[[j]])-LAG)])
      X1<-pIl*pJl
      X2<-pIl*pIl*pJl
      ded=ts.intersect( pI=ts(dsz[i]), X1, X2 )
      theFit<-lm( pI ~ X1 + X2, data=ded, na.action=na.exclude)
			rSq<-summary(theFit)[["r.squared"]]
		 	theFit$control[1]<-rSq								# store r-squared in theFit$control for printing convenience in nlmOut

      wun<-ifelse(rSq<lPivotR2,-1,1)             # preserve sign
      aitch<-sqrt( abs(rSq-lPivotR2) )*wun       # h = sqrt( this r^2 - auto corr r^2 of this row[i] NOTE: *wun preserves sign

      pList[[i]]$NonLinearSync2[[j]]["rsq"]<-rSq

			#	if (summary(theFit)[["coefficients"]]["delt","Pr(>|t|)"]<PLEVEL) syncMatrixNL2[j,i]<-coef(theFit)["delt"] else syncMatrixNL2[j,i]<-0.0

			syncMatrixNL2[j,i]<-aitch
			pList[[i]]$NonLinearSync2[[j]]["bi"]<-coef(theFit)["X1"]
		 	pList[[i]]$NonLinearSync2[[j]]["qi"]<-coef(theFit)["X2"]			# it gets 0.0 if insignificant
	    cat("---------------------------------------------------------------\n")
      cat("Formula:  pI ~ b1*X1 + b2*X2, where X1=pIl*pJl and X2=pIl^2*pJl\n")
  	  cat("---------------------------------------------------------------\n")
      stargazer(theFit, type = "text")
      cat(sprintf("h^2 = (R^2 model - R^2 LAC) = (%6.4f - %6.4f) = %6.4f\n",rSq,lPivotR2,aitch^2))
      cat(sprintf("h = %6.4f\n",aitch))
			cat("---------------------------------------------------------------------|\n\n")


      #
			#  Opt 3: nonLinear Regression:  z_2 = A*P_1*z_1*e^{B*z_1} + C
			#
			cat(sprintf("... Player %d NONLINEAR regression Opt 3 on Player %d (I=%d, J=%d)\n", i,j,i,j))
			pI<-dsz[[i]]
			pIl<- c(rep(NA,times=LAG),dsz[[i]][1:(length(dsz[[i]])-LAG)])
			pJl<-c(rep(NA,times=LAG),dsz[[j]][1:(length(dsz[[j]])-LAG)])
			theFit <- nlsLM(pI ~ alpha*pJl*pIl*exp(beta * pIl) + gamma, start=c(alpha=0.5,beta=0.5,gamma=0.5), na.action = na.exclude )
		 	rSq <- 1-(deviance(theFit)/sum((pI-mean(pI))^2))   	# R^2 = 1 - [ (Residual Sum of Squares / Corrected Sum of Squares) ]
		 	theFit$control[1]<-rSq								# store r-squared in theFit$control for printing convenience in nlmOut

      wun<-ifelse(rSq<nlPivotR2,-1,1)             # preserve sign
      aitch<-sqrt( abs(rSq-nlPivotR2) )*wun       # h = sqrt( this r^2 - auto corr r^2 of this row[i] NOTE: *wun preserves sign

      # if (summary(theFit)[["coefficients"]]["alpha","Pr(>|t|)"]<PLEVEL) syncMatrixNL3[j,i]<-coef(theFit)["alpha"] else syncMatrixNL3[j,i]<-0.0

      syncMatrixNL3[j,i]<-aitch
      pList[[i]]$NonLinearSync3[[j]]["rsq"]<-rSq
		 	pList[[i]]$NonLinearSync3[[j]]["B"]<-coef(theFit)["beta"]
		 	pList[[i]]$NonLinearSync3[[j]]["A"]<-coef(theFit)["alpha"]
      pList[[i]]$NonLinearSync3[[j]]["C"]<-coef(theFit)["gamma"]      			# ????? it gets 0.0 if insignificant
			writeLines(nlmOut(theFit))
      cat(sprintf("h^2 = (R^2 model - R^2 NLAC) = (%6.4f - %6.4f) = %6.4f\n",rSq,nlPivotR2,aitch^2))
      cat(sprintf("h = %6.4f\n",aitch))
			cat("---------------------------------------------------------------------|\n\n")

		 }
	 }
 }
  sink()
  #################
  # end of calculations
  #################




  #################
  # write linear sync matrix
  #################
   sink(paste(thePath, lMatrixFileName,sep="/"),split=FALSE)
   cat(outName,"Linear Sync Matrix","\n")
   cat(N)
   prmatrix(syncMatrixL, rowlab=rep("",N), collab=rep("",N))
   sink()

   #################
   # write nonlinear sync matrix - 1st based on  z_2 = Ae^{Bz_1} model with H
   #################
   sink(paste(thePath, nlMatrixFileName1,sep="/"),split=FALSE)
   cat(outName,"nonLinear Sync Matrix - Opt 1: z2 = A*e^{B*z_1} + e^{D*p_1}","\n")
   cat(N)
   prmatrix(syncMatrixNL1, rowlab=rep("",N), collab=rep("",N))
   sink()

   #################
   # write nonlinear sync matrix  - 2nd based on X_2 =P_1X_1(1−X_1) with H
   #################
   sink(paste(thePath, nlMatrixFileName2,sep="/"),split=FALSE)
   cat(outName,"nonLinear Sync Matrix - Opt 2: z2 = p1*z1*(1−z1)","\n")
   cat(N)
   prmatrix(syncMatrixNL2, rowlab=rep("",N), collab=rep("",N))
   sink()

   #################
   # write nonlinear sync matrix  - 3rd based on z_2 = AP_1z_1e^{Bz_1} + C model with H
   #################
   sink(paste(thePath, nlMatrixFileName3,sep="/"),split=FALSE)
   cat(outName,"nonLinear Sync Matrix - Opt 3: z2 = A*p1*z1*e^(B*z1) + C","\n")
   cat(N)
   prmatrix(syncMatrixNL3, rowlab=rep("",N), collab=rep("",N))
   sink()






 #################
 # write summary file
 #################
 sink(paste(thePath, summaryFileName,sep="/"),split=FALSE)
 cat(sprintf("========================================================\n"))
 cat("Data Results Summary:",outName," ( LAG=",LAG,")\n")
 cat(sprintf("========================================================\n\n"))
   for (i in 1:N ) {
   colHead<-paste("p",i,"w/p", 1:N,sep="")		# set up column headers
 	 colHead<-colHead[-c(i)]					# remove the ith column
 	 xSet<-1:N									# range to loop over for comparisons
 	 xSet<-xSet[-i]								# remove the "diagnoal" index

 	 cat(sprintf("--------------------------------------------------\n"))
 	 cat(sprintf("         Player %d Results Summary\n",i))
  	 cat(sprintf("---------------------------------------------------\n"))
 	 cat(sprintf("Linear Auto\n"))
 	 cat(sprintf("   r^2 = %9.6f\n\n",pList[[i]]$LinearAutoCorr["rsq"]))

 	 cat(sprintf("Linear Sync\n"))
 	 cat(sprintf("%16s",colHead),"\n")
 	 for (j in xSet ) { cat(sprintf("   R^2 = %9.6f",pList[[i]]$LinearSync[[j]]["rsq"])) }
 	 cat("\n")
 	 for (j in xSet ) {	cat(sprintf("    b1 = %9.6f",pList[[i]]$LinearSync[[j]]["b1"]))   }
 	 cat("\n")
 	 for (j in xSet ) {	cat(sprintf("    b2 = %9.6f",pList[[i]]$LinearSync[[j]]["b2"]))	 }
 	 cat("\n\n")

	 cat(sprintf("Nonlinear Auto - Model:  A*e^(B*z1)\n"))
 	 cat(sprintf("   r^2 = %9.6f\n",pList[[i]]$NonLinearAuto["rsq"]))
 	 cat(sprintf("     B = %9.6f\n",pList[[i]]$NonLinearAuto["B"]))
   cat(sprintf("     A = %9.6f\n\n",pList[[i]]$NonLinearAuto["A"]))

 	 cat(sprintf("Model:  A*e^(B*z1) + e^(D*p1) [Option 1]\n"))
 	 cat(sprintf("%16s",colHead),"\n")
 	 for (j in xSet ) {	cat(sprintf("   R^2 = %9.6f",pList[[i]]$NonLinearSync1[[j]]["rsq"])) }
 	 cat("\n")
 	 for (j in xSet ) {	cat(sprintf("     B = %9.6f",pList[[i]]$NonLinearSync1[[j]]["B"]))	}
 	 cat("\n")
 	 for (j in xSet ) {	cat(sprintf("     D = %9.6f",pList[[i]]$NonLinearSync1[[j]]["D"]))	}
   cat("\n")
   for (j in xSet ) {	cat(sprintf("     A = %9.6f",pList[[i]]$NonLinearSync1[[j]]["A"]))	}
 	 cat("\n\n")

 	 cat(sprintf("Model:  p1*z1*(1−z1) = B*(p1*z1) + Q*(p1+z1^2) [Option 2]\n"))
 	 cat(sprintf("%16s",colHead),"\n")
 	 for (j in xSet ) {	cat(sprintf("   R^2 = %9.6f",pList[[i]]$NonLinearSync2[[j]]["rsq"])) }
 	 cat("\n")
 	 for (j in xSet ) {	cat(sprintf("     B = %9.6f",pList[[i]]$NonLinearSync2[[j]]["bi"]))	}
 	 cat("\n")
 	 for (j in xSet ) {	cat(sprintf("     Q = %9.6f",pList[[i]]$NonLinearSync2[[j]]["qi"]))	}
 	 cat("\n\n")

 	 cat(sprintf("Model:  A*p1*z1*e^(B*z1) + C [Option 3]\n"))
 	 cat(sprintf("%16s",colHead),"\n")
 	 for (j in xSet ) {	cat(sprintf("   R^2 = %9.6f",pList[[i]]$NonLinearSync3[[j]]["rsq"])) }
 	 cat("\n")
 	 for (j in xSet ) {	cat(sprintf("     A = %9.6f",pList[[i]]$NonLinearSync3[[j]]["A"]))	}
 	 cat("\n")
 	 for (j in xSet ) {	cat(sprintf("     B = %9.6f",pList[[i]]$NonLinearSync3[[j]]["B"]))	}
   cat("\n")
   for (j in xSet ) {	cat(sprintf("     C = %9.6f",pList[[i]]$NonLinearSync3[[j]]["C"]))	}
 	 cat("\n\n")

	 }  # for (i in 1:N )
 sink()
