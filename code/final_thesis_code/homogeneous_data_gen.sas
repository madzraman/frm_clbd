 
/* Generating multiple datasets of size 100 per group - power for time-varying group effect  */
/* Logistic mixed model for multiple observations withins subjects
/* and stratified conditional logisitic model 

/* adding in the BS and WS effects of the time-varying group 

/* this program gets a random GRPMEAN for each subject (which is correlated with the random effect)
/* and then samples GRP at each timepoint from Bernoulli with P=GRPMEAN

/* writing log and output info to external files */


DATA parms;
* number of datasets ;
ndats = 250;
* number of subjects per group;
npergrp = 100;
* total number of subjects ;
nsubj = 2*npergrp;
* subject-specific beta for WS effect;
ssbeta = -.4;
* subject-specific beta for BS effect;
array ssbetaBS(3) ssbetaBS1-ssbetaBS3 (-.4 0 .4);

* number of timepoints; 
ntime = 10;
* ICC for timepoints;
icc =  .3;
* covariance between GRPMEAN latent variable and random subject effect intercept;
array cov(3) cov1-cov3 (-.25 0 .25); 

* calculate the subject variance, given the ICC and error variance;
varerr = ((ATAN(1)*4)**2)/3;
varsub = varerr*(icc/(1-icc));
sdsub  = sqrt(varsub);

* seed for random number generation;
seed   = 974657747;


******************************************************************************************;
* NOTHING BELOW NEEDS TO BE SPECIFIED;
******************************************************************************************;


libname saveit '/home/u62617082/sasuser.v94'; 


data saveit.homo_out_250d_n04; set parms;

array ssbetaBS(3) ssbetaBS1-ssbetaBS3; 
array cov(3) cov1-cov3;
do j = 1 to 3;
do k = 1 to 3;

ndat=1;
do while (ndat le ndats);

   n=1;
   do while (n le nsubj);
      id = n; 
	  * correlate the random subject effect and the grpmean;
      subjZ = RANNOR(SEED);
	   subj = SQRT(varsub)*subjZ;
	   grpmeanZ = (cov(j)/sqrt(varsub))*subjZ + sqrt(1 - (cov(j)/sqrt(varsub))**2)*RANNOR(SEED);
      grpmean = CDF('NORMAL',grpmeanz,1);

         nt = 1;
	      do while (nt le ntime);
            exp1 = ranexp(seed);
            exp2 = ranexp(seed);
            err  = log(exp1/exp2);
			grp  = RAND('BERNOULLI',grpmean);
			grpdev = grp-grpmean;
            ystar = -.5 + grpmean*ssbetaBS(k) + grpdev*ssbeta + subj + err;
            if ystar <= 0 then do;
               y = 0; output;
               nt = nt+1;
             end;
            if ystar > 0 then do;
               y = 1; output;
               nt = nt+1;
             end;
         end;

   n+1;
   end;
ndat+1;end;

end;end;

PROC PRINT data=saveit.homo_out_250d_n04; run;