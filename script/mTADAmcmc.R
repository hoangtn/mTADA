##Hoang Nguyen: Jan 12, 2017
DNandCC2traits <- "
data {
    int<lower=1> NN; //Number of genes
    int<lower=0> NCdn;
    int Ndn[NCdn]; //Number of trios
    int<lower=0> NC1dn;

    int dataDN[NN, NCdn]; //denovo data: Kdn classes
    real mutRate[NN, NCdn]; //mutation rates: Kdn classes

    int<lower=0> NCcc; //
    int<lower=0> NC1cc; //
    int Ncase[NCcc]; //Number of case samples
    int Ncontrol[NCcc]; //Number of control samples

    int dataCCcase[NN, NCcc]; //denovo data: Kdn classes
    int dataCCtotal[NN, NCcc]; //mutation rates: Kdn classes

    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    //These below parameters should be default
    real<lower=0> upperPi0;
    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks 
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaDN0[NCdn];
    real<lower=0> hyperBetaCC0[NCcc];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means

    real<lower=0> thetaH0[NCcc]; //Nca/(Nca + Ncn)
    real<lower=0> pi01;
    real<lower=0> pi02;
    real<lower=0> hyper2BetaDN[NCdn];
    real<lower=0> hyper2BetaCC[NCcc];
    }

parameters {
    real<lower=0,upper=fmin(pi01, pi02)> pi0; //Proportion of risk genes
    real<lower=lowerHyperGamma> hyperGammaMeanDN[NCdn]; //Hyper parameters for de novo relative risks
    real<lower=lowerGamma> gammaMeanDN[2, NCdn]; //parameters (in the sampling process) for de novo relative risks

    real<lower=lowerHyperGamma> hyperGammaMeanCC[NCcc]; //Hyper parameters for case-control relative risks
    real<lower=lowerGamma> gammaMeanCC[2, NCcc]; //parameters (in the sampling process) for case-control relative risks

}

transformed parameters {
    real hyperBetaDN[NCdn];
    real hyperBetaCC[NCcc];
    if (adjustHyperBeta != 0) {
      for (ii in 1:NCdn){
           hyperBetaDN[ii] = exp(betaPars[1]*hyperGammaMeanDN[ii]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);
      
       }
      for (ii in 1:NCcc){
           hyperBetaCC[ii] = exp(betaPars[1]*hyperGammaMeanCC[ii]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);
      
       }

   }
    else {
        hyperBetaDN = hyperBetaDN0;
        hyperBetaCC = hyperBetaCC0;
        }
    }

model {
//     pi0 ~ dirichlet(c(80, 3, 4, 5)); 
     real ps[4];
    pi0 ~ beta(1, 100);
  //De novo data: sample for hyper priors (NPdn populations and Kdn categories)
    for (ip in 1:NCdn){
         hyperGammaMeanDN[ip] ~ gamma(1, hyper2BetaDN[ip]);
         
     }
  //Case-control data: sample for hyper priors (NPcc populations and Kdn categories)
    for (ip in 1:NCcc){
         hyperGammaMeanCC[ip] ~ gamma(1, hyper2BetaCC[ip]);
         
     }
   //Start sampling for specific RRs

   for (ip in 1:NCdn){
     for (jj in 1:2){
          gammaMeanDN[jj, ip] ~ gamma(hyperGammaMeanDN[ip]*hyperBetaDN[ip], hyperBetaDN[ip]);
           }}
   for (ip in 1:NCcc){
     for (jj in 1:2){
          gammaMeanCC[jj, ip] ~ gamma(hyperGammaMeanCC[ip]*hyperBetaCC[ip], hyperBetaCC[ip]);
           }}

////Main program
//Loop through data points
////
     for (ii in 1:NN){

         ps[1] = log1m(pi01 + pi02 - pi0);
         ps[2] = log(pi0);
         ps[3] = log(pi01 - pi0);
         ps[4] = log(pi02 - pi0);

//prob0 <- c(0.88, 0.03, 0.04, 0.05) #Null, both, first, second
//// De novo data
       for (jj in 1:NCdn){
        ps[1] = ps[1] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]); //NULL
        ps[2] = ps[2] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]*gammaMeanDN[1, jj]); //Model 1: Both
        }
       for (jj in 1:NC1dn){
//1:NC1dn is for the 1st group; (NC1dn + 1):NCdn is for the second group
         ps[3] = ps[3] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]*gammaMeanDN[2, jj]); //Model 2: First trait
         ps[3] = ps[3] + poisson_lpmf(dataDN[ii, jj + NC1dn] | Ndn[jj + NC1dn]*2*mutRate[ii, jj + NC1dn]); //Model 2: Second trait
         ps[4] = ps[4] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]); //Model 4: First trait
         ps[4] = ps[4] + poisson_lpmf(dataDN[ii, jj + NC1dn] | Ndn[jj + NC1dn]*2*mutRate[ii, jj + NC1dn]*gammaMeanDN[2, jj + NC1dn]); //Model 4: Second trait

        }
//Case-control data
       for (jj in 1:NCcc){
        ps[1] = ps[1] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], thetaH0[jj]); //Add Null hypothesis
        ps[2] = ps[2] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj],
                         gammaMeanCC[1, jj]*Ncase[jj]/(gammaMeanCC[1, jj]*Ncase[jj] + Ncontrol[jj])); //Both
       }

       for (jj in 1:NC1cc){
        ps[3] = ps[3] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], gammaMeanCC[2, jj]*Ncase[jj]/(gammaMeanCC[2, jj]*Ncase[jj] + Ncontrol[jj])); //Model 3: first trait
        ps[3] = ps[3] + binomial_lpmf(dataCCcase[ii, jj + NC1cc] | dataCCtotal[ii, jj + NC1cc], thetaH0[jj + NC1cc]); //Model 3: second trait

        ps[4] = ps[4] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], thetaH0[jj]); //Model 3: first trait
        ps[4] = ps[4] + binomial_lpmf(dataCCcase[ii, jj + NC1cc] | dataCCtotal[ii, jj + NC1cc],
          gammaMeanCC[2, jj + NC1cc]*Ncase[jj + NC1cc]/(gammaMeanCC[2, jj + NC1cc]*Ncase[jj + NC1cc] + Ncontrol[jj + NC1cc])); //Model 4: second trait
       }

         target += log_sum_exp(ps);
         }
}
"
#######################
##Only DE NOVO
DN2traits <- "
data {
    int<lower=1> NN; //Number of genes
    int<lower=0> NCdn1;
    int Ndn1[NCdn1]; //Number of trios
    real<lower=0> hyperGammaMeanDN1[NCdn1];

    int<lower=0> NCdn2;
    int Ndn2[NCdn2]; //Number of trios
    real<lower=0> hyperGammaMeanDN2[NCdn2];

    int dataDN1[NN, NCdn1]; //denovo data: Kdn classes
    real mutRate1[NN, NCdn1]; //mutation rates: Kdn classes

    int dataDN2[NN, NCdn2]; //denovo data: Kdn classes
    real mutRate2[NN, NCdn2]; //mutation rates: Kdn classes


    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    //These below parameters should be default
    real<lower=0> lowerGamma; //Low limit for relative risks 
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaDN01[NCdn1];
    real<lower=0> hyperBetaDN02[NCdn2];

    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means

    real<lower=0> pi01;
    real<lower=0> pi02;
//    real UpperAlpha;
    }

parameters {
    real<lower=lowerGamma> gammaMeanDN1[NCdn1]; //parameters (in the sampling process) for de novo relative risks
    real<lower=lowerGamma> gammaMeanDN2[NCdn2]; //parameters (in the sampling process) for de novo relative risks

 //   real<upper=UpperAlpha> alpha0;
    real<lower=0,upper=fmin(pi01, pi02)> p12; //[NN]; //Proportion of risk genes
}

transformed parameters {
  //  real p12; //<lower=0,upper=fmin(pi01, pi02)> p12; //[NN]; //Proportion of risk genes
    real hyperBetaDN1[NCdn1];
    real hyperBetaDN2[NCdn2];

//       p12 = exp(alpha0)/(1 + exp(alpha0));


    if (adjustHyperBeta != 0) {
      for (i2i in 1:NCdn1){
            hyperBetaDN1[i2i] = exp(betaPars[1]*hyperGammaMeanDN1[i2i]^(betaPars[2]) + betaPars[3]); 
     
       }
      for (i2i in 1:NCdn2){
            hyperBetaDN2[i2i] = exp(betaPars[1]*hyperGammaMeanDN2[i2i]^(betaPars[2]) + betaPars[3]); 
     
       }

   }
    else {
        hyperBetaDN1 = hyperBetaDN01;
        hyperBetaDN2 = hyperBetaDN02;
        }
    }


model {

     real ps[4];
//   alpha0 ~ normal(0, 5);
//   p12 ~ 

   for (ip in 1:NCdn1){
          gammaMeanDN1[ip] ~ gamma(hyperGammaMeanDN1[ip]*hyperBetaDN1[ip], hyperBetaDN1[ip]);
           }

   for (ip in 1:NCdn2){
          gammaMeanDN2[ip] ~ gamma(hyperGammaMeanDN2[ip]*hyperBetaDN2[ip], hyperBetaDN2[ip]);
           }

////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log1m(pi01 + pi02 - p12);
         ps[2] = log(p12);
         ps[3] = log(pi01 - p12);
         ps[4] = log(pi02 - p12);
//prob0 <- c(0.88, 0.03, 0.04, 0.05) #Null, both, first, second
//// De novo data
       for (jj in 1:NCdn1){
        ps[1] = ps[1] + poisson_lpmf(dataDN1[ii, jj] | Ndn1[jj]*2*mutRate1[ii, jj]); //NULL: loop across two traits
        ps[2] = ps[2] + poisson_lpmf(dataDN1[ii, jj] | Ndn1[jj]*2*mutRate1[ii, jj]*gammaMeanDN1[jj]); //Model 1: Both traits/loop across two traits
        }
       for (jj in 1:NCdn2){
        ps[1] = ps[1] + poisson_lpmf(dataDN2[ii, jj] | Ndn2[jj]*2*mutRate2[ii, jj]); //NULL: loop across two traits
        ps[2] = ps[2] + poisson_lpmf(dataDN2[ii, jj] | Ndn2[jj]*2*mutRate2[ii, jj]*gammaMeanDN2[jj]); //Model 1: Both traits/loop across two traits
        }

  //Just loop a half of gamma for each trait: NC1dn = NCdn/2: Category number of de novo mutations 
       for (jj in 1:NCdn1){
         ps[3] = ps[3] + poisson_lpmf(dataDN1[ii, jj] | Ndn1[jj]*2*mutRate1[ii, jj]*gammaMeanDN1[jj]); //Model 2: First trait/loop for first trait 1:NC1dn
         ps[4] = ps[4] + poisson_lpmf(dataDN1[ii, jj] | Ndn1[jj]*2*mutRate1[ii, jj]); //Model 4: First trait/loop for first trait 1:NC1dn
}
       for (jj in 1:NCdn2){
         ps[3] = ps[3] + poisson_lpmf(dataDN2[ii, jj] | Ndn2[jj]*2*mutRate2[ii, jj]); //Model 2: Second trait/loop for second trait: jj + NC1dn
         ps[4] = ps[4] + poisson_lpmf(dataDN2[ii, jj] | Ndn2[jj]*2*mutRate2[ii, jj]*gammaMeanDN2[jj]); //Model 4: Second trait/loop for second trait: jj + NC1dn
        }

         target += log_sum_exp(ps);
         }
}
"
######################
#####Only CASE-CONTROL
CC2traits <- "
data {
    int<lower=1> NN; //Number of genes

    int<lower=0> NCcc; //
    int<lower=0> NC1cc; //
    int Ncase[NCcc]; //Number of case samples
    int Ncontrol[NCcc]; //Number of control samples

    int dataCCcase[NN, NCcc]; //denovo data: Kdn classes
    int dataCCtotal[NN, NCcc]; //mutation rates: Kdn classes

    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    //These below parameters should be default
    real<lower=0> upperPi0;
    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks 
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaCC0[NCcc];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means

    real<lower=0> thetaH0[NCcc]; //Nca/(Nca + Ncn)
    real<lower=0> pi01;
    real<lower=0> pi02;
    real<lower=0> hyper2BetaCC[NCcc];
    }

parameters {
    real<lower=0,upper=fmin(pi01, pi02)> pi0; //Proportion of risk genes
    real<lower=lowerHyperGamma> hyperGammaMeanCC[NCcc]; //Hyper parameters for case-control relative risks
    real<lower=lowerGamma> gammaMeanCC[2, NCcc]; //parameters (in the sampling process) for case-control relative risks

}

transformed parameters {
    real hyperBetaCC[NCcc];
    if (adjustHyperBeta != 0) {
      for (ii in 1:NCcc){
           hyperBetaCC[ii] = exp(betaPars[1]*hyperGammaMeanCC[ii]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);
      
       }

   }
    else {
        hyperBetaCC = hyperBetaCC0;
        }
    }

model {
     real ps[4];
//     pi0 ~ dirichlet(c(80, 3, 4, 5));
    pi0 ~ beta(1, 100);

  //Case-control data: sample for hyper priors (NPcc populations and Kdn categories)
    for (ip in 1:NCcc){
         hyperGammaMeanCC[ip] ~ gamma(1, hyper2BetaCC[ip]);
         
     }
   //Start sampling for specific RRs

   for (ip in 1:NCcc){
     for (jj in 1:2){
          gammaMeanCC[jj, ip] ~ gamma(hyperGammaMeanCC[ip]*hyperBetaCC[ip], hyperBetaCC[ip]);
           }}

////Main program
//Loop through data points
////
     for (ii in 1:NN){

         ps[1] = log1m(pi01 + pi02 - pi0);
         ps[2] = log(pi0);
         ps[3] = log(pi01 - pi0);
         ps[4] = log(pi02 - pi0);
//prob0 <- c(0.88, 0.03, 0.04, 0.05) #Null, both, first, second
//Case-control data
       for (jj in 1:NCcc){
        ps[1] = ps[1] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], thetaH0[jj]); //Add Null hypothesis
        ps[2] = ps[2] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj],
                         gammaMeanCC[1, jj]*Ncase[jj]/(gammaMeanCC[1, jj]*Ncase[jj] + Ncontrol[jj])); //Both
       }

       for (jj in 1:NC1cc){
        ps[3] = ps[3] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], gammaMeanCC[2, jj]*Ncase[jj]/(gammaMeanCC[2, jj]*Ncase[jj] + Ncontrol[jj])); //Model 3: first trait
        ps[3] = ps[3] + binomial_lpmf(dataCCcase[ii, jj + NC1cc] | dataCCtotal[ii, jj + NC1cc], thetaH0[jj + NC1cc]); //Model 3: second trait

        ps[4] = ps[4] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], thetaH0[jj]); //Model 3: first trait
        ps[4] = ps[4] + binomial_lpmf(dataCCcase[ii, jj + NC1cc] | dataCCtotal[ii, jj + NC1cc],
          gammaMeanCC[2, jj + NC1cc]*Ncase[jj + NC1cc]/(gammaMeanCC[2, jj + NC1cc]*Ncase[jj + NC1cc] + Ncontrol[jj + NC1cc])); //Model 4: second trait
       }

         target += log_sum_exp(ps);
         }
}
"


