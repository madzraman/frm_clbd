*OPTIONS NOCENTER NoSource NoNotes NoSpool NoDate NoNumber;
*ods noresults; 
/* Generating multiple datasets of size 100 per group - power for time-varying group effect  */
/* Logistic mixed model for multiple observations withins subjects
/* and stratified conditional logisitic model 

/* adding in the BS and WS effects of the time-varying group 

/* this program gets a random GRPMEAN for each subject (which is correlated with the random effect)
/* and then samples GRP at each timepoint from Bernoulli with P=GRPMEAN

/* writing log and output info to external files */
*proc printto log  ='u:\DATA\POWER\BSWSsurvdat\RandomInt_MixedLogistic_CondLogistic_BSWS_N=100.log';
*proc printto print='u:\DATA\POWER\BSWSsurvdat\RandomInt_MixedLogistic_CondLogistic_BSWS_N=100.lst';
*run;

*DATA parms;



* ----> j loop: 3 DIFFERENT COV MATRICES: ; 



******************************************************************************************;
* NOTHING BELOW NEEDS TO BE SPECIFIED;
******************************************************************************************;



libname saveit '/home/u62617082/sasuser.v94'; 


proc iml;

covmatA = {0.8225 -0.2480, 
		  -0.2480 1.8664};
covmatB = {0.8225 0, 
		   0 1.3708};
covmatC = {0.8225 0.1820, 
		   0.1820 1.0068};
meanvec = {0 0};		  

all_covmats = ListCreate(3);
call ListSetItem(all_covmats, 1, covmatA);    /* add covmatA as 1st element to all_covmats */
call ListSetItem(all_covmats, 2, covmatB);    /* add covmatB as 2nd element to all_covmats */	  
call ListSetItem(all_covmats, 3, covmatC);    /* add covmatC as 3rd element to all_covmats */	  

*package load ListUtil;
*run struct(all_covmats);


************************************************************************************************;

* number of datasets ;
ndats = 100;

* number of subjects per group;
npergrp = 100;

* total number of subjects ;
nsubj = 2*npergrp;

* subject-specific beta for WS effect;
ssbetaWS = -.4;

* k loop: subject-specific beta for BS effect;
ssbetaBS = {-.4 0 .4};

* number of timepoints; 
ntime = 10;

***************************;

* c loop: ONE TABLE AT A TIME: covariance between GRPMEAN latent variable and random subject effect intercept ;
cov = {-.25 0 .25}; 

* ICC for random intercept;
iccINT = .2; * for all three cases, not actually used ;

* ICCs for random exposures;
iccSLOPE = {.362 .294 .234}; * not actually used ;

* correlation between random intercept and random exposure effects;
REcorr = {-.2 0 .2}; * used to make the 3 cov matrices respectively; * but not actually used ;

	  
/* can convert these matrices hard coded into loops later */

* calculate the subject variance FOR THE RANDOM INTERCEPT, given the ICC and error variance. One value.;
* or make this also a vector of three so we can basically vectorize the cov calculations. idk;

varerr = ((ATAN(1)*4)**2)/3;


* seed for random number generation;
seed   = 974657747;

c=1; k=1; j=1; ndat=0; n=0; id=0; subjINT=0; subjEXP=0; grpmeanZ=0; grpmean=0; nt=0; exp1=0; exp2=0; err=0; grp=0; grpdev=0; ystar=0; y=0;

************************************************************************************************;
xys = ndats || npergrp || nsubj || ssbetaWS || ssbetaBS[1]|| ssbetaBS[2] || ssbetaBS[3] || ntime || iccINT || iccSLOPE[1] || iccSLOPE[2] || iccSLOPE[3] || cov[1] || cov[2] || cov[3] || 
	  varerr || seed || c|| k || j || cov[c] || ssbetaBS[k] || REcorr[j] || ndat || n || id || subjINT || subjEXP || grpmeanZ || grpmean || nt || exp1 || exp2 || err || grp || grpdev || ystar || y;
cname = {"ndats" "npergrp" "nsubj" "ssbetaWS" "ssbetaBS[1]" "ssbetaBS[2]" "ssbetaBS[3]" "ntime" "iccINT" "iccSLOPE[1]" "iccSLOPE[2]" "iccSLOPE[3]" "cov[1]" "cov[2]" "cov[3]" "varerr" "seed" 
		 "c" "k" "j" "cov[c]" "ssbetaBS[k]" "REcorr[j]" "ndat" "n" "id" "subjINT" "subjEXP" "grpmeanZ" "grpmean" "nt" "exp1" "exp2" "err" "grp" "grpdev" "ystar" "y"};

create saveit.hetero_out_100d_n04_8 from xys [ colname=cname ];

do c = 1 to 3; * for each cov between ssbetaBS and random subject intercept;
	do k = 1 to 3; * for each ssBetaBS ;
		do j = 1 to 3; * ListLen(all_covmats); * for each covmatrix between v0 and v1 aka 3 of them ;
				
			covmatJ = ListGetItem(all_covmats, j); * get jth cov matrix of the 3;  	
		  	*print covmatJ;
		  	
		  	varsubINT = covmatJ[1,1];
		 	sdsubINT = SQRT(varsubINT);
		 	varsubEXP = covmatJ[2,2];
		 	sdsubEXP = SQRT(varsubEXP);
			covV0V1 = covmatJ[1,2]; * corr V0V1 is REcorr in tables;
		
			ndat=1;
			do while (ndat <= ndats);
			
			   	n=1;
			   	do while (n <= nsubj);
			    	id = n; 
			
				  	* correlate the random subject effect intercept and the grpmean, still want to do this for heterogeneous exposure models;
			      	*subjZINT = RANNOR(SEED);
			      	*subjZEXP = RANNOR(SEED);
			     	
			     	subj = RandNormal(1,meanvec,covmatJ);
				 	subjINT = subj[1,1];
			      	subjEXP = subj[1,2];
				  
				  	grpmeanZ = (cov[c]/sdsubINT)*RANNOR(SEED) + sqrt(1 - (cov[c]/sdsubINT)**2)*RANNOR(SEED); * N(0,1) distributed we know;
			      	grpmean = CDF('NORMAL',grpmeanZ,1);
			
			        nt = 1;
				    do while (nt <= ntime);
				        exp1 = ranexp(seed);
				        exp2 = ranexp(seed);
				        err  = log(exp1/exp2);
						grp  = RAND('BERNOULLI',grpmean);
						grpdev = grp-grpmean;
				        ystar = -.5 + grpmean*ssbetaBS[k] + grpdev*ssbetaWS + subjINT + subjEXP*grp + err; 
			            if ystar <= 0 then do;
			            	y = 0; 
			            	xys = ndats || npergrp || nsubj || ssbetaWS || ssbetaBS[1]|| ssbetaBS[2] || ssbetaBS[3] || ntime || iccINT || iccSLOPE[1] || iccSLOPE[2] || iccSLOPE[3] || cov[1] || cov[2] || cov[3] || 
	  							  varerr || seed || c|| k || j || cov[c] || ssbetaBS[k] || REcorr[j] || ndat || n || id || subjINT || subjEXP || grpmeanZ || grpmean || nt || exp1 || exp2 || err || grp || grpdev || ystar || y;
							cname = {"ndats" "npergrp" "nsubj" "ssbetaWS" "ssbetaBS[1]" "ssbetaBS[2]" "ssbetaBS[3]" "ntime" "iccINT" "iccSLOPE[1]" "iccSLOPE[2]" "iccSLOPE[3]" "cov[1]" "cov[2]" "cov[3]" "varerr" "seed" 
		 							 "c" "k" "j" "cov[c]" "ssbetaBS[k]" "REcorr[j]" "ndat" "n" "id" "subjINT" "subjEXP" "grpmeanZ" "grpmean" "nt" "exp1" "exp2" "err" "grp" "grpdev" "ystar" "y"};
			            	append from xys;
			            	nt = nt+1;
			            end;
			            if ystar > 0 then do;
			            	y = 1; 
			            	xys = ndats || npergrp || nsubj || ssbetaWS || ssbetaBS[1]|| ssbetaBS[2] || ssbetaBS[3] || ntime || iccINT || iccSLOPE[1] || iccSLOPE[2] || iccSLOPE[3] || cov[1] || cov[2] || cov[3] || 
	  							  varerr || seed || c|| k || j || cov[c] || ssbetaBS[k] || REcorr[j] || ndat || n || id || subjINT || subjEXP || grpmeanZ || grpmean || nt || exp1 || exp2 || err || grp || grpdev || ystar || y;
							cname = {"ndats" "npergrp" "nsubj" "ssbetaWS" "ssbetaBS[1]" "ssbetaBS[2]" "ssbetaBS[3]" "ntime" "iccINT" "iccSLOPE[1]" "iccSLOPE[2]" "iccSLOPE[3]" "cov[1]" "cov[2]" "cov[3]" "varerr" "seed" 
									 "c" "k" "j" "cov[c]" "ssbetaBS[k]" "REcorr[j]" "ndat" "n" "id" "subjINT" "subjEXP" "grpmeanZ" "grpmean" "nt" "exp1" "exp2" "err" "grp" "grpdev" "ystar" "y"};
			            	append from xys;
			            	nt = nt+1;
			            end;
					*print y;
					end; * end loop over total num timesteps;
					
				n=n+1;
				end; * end loop over total num subjects;
				
			ndat=ndat+1;
			end; * end loop over total num datasets;
		
		end; * end j each cov matrix (aka each REcorr between v0 and v1);
	end; * end k each ssBetaBS;
end; * end c each cov between ssBetaBS and v0;

close saveit.hetero_out_100d_n04_8;

proc print data=saveit.hetero_out_100d_n04_8; run;


