#################################################################
# Some useful functions
#################################################################

# Fisher's method of combining p values
# Input: a vector of p-values
combpval <- function(p) {
  k <- length(p)
  T <- -2 * sum(log(p))
  return (1 - pchisq(T,2*k))
}

# Genomic control
# Input: a vector of p-values, the quantile (0.5 or 0.75)
genom.ctrl <- function(p, quant) {
  chisq.obs <- qchisq(1-p, 1) # convert p-values to chi-square statistics (dof 1)
  chisq.obs.quant <- quantile(chisq.obs, probs=quant, names=FALSE) # 0.75 if first quantile
  chisq.exp.quant <- qchisq(quant, 1)
  return (chisq.obs.quant / chisq.exp.quant)
}

# Bayesian FDR control (PMID:19822692, Section2.3)
# BF: a sorted vector of BFs (in decreasing order)
# pi0: the prior probability that the null model is true
# alpha: the FDR target
# Return: the q-value of each BF, and the number of findings with q below alpha. 
Bayesian.FDR <- function(BF, pi0, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(BF)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}

# draw QQ plot
# Input: a vector of p-values
plotQQ <- function(p.obs) {
  obs <- -log10(p.obs)
  theo <- -log10(ppoints(length(obs)))
  qqplot(theo, obs, xlab=expression("Theoretical " * -log[10](p)), ylab=expression("Observed "*-log[10](p)))
  abline(0,1,col='red')
}

# Similar, but show QQ plot in the original scale (not log. scale)
plotQQ.unif <- function(p.obs) {
  obs <- (p.obs)
  theo <- (ppoints(length(obs)))
  qqplot(theo, obs, xlab="Theoretical p-values", ylab="Observed p-values")
  abline(0,1,col='red')
}

#################################################################
# Bayes Factor Computation for a Single Gene
#################################################################

# model evidence of de novo data: P(x_d|H_0) 
# Input: the count data x, the sample size N, the mutation rate mu 
evidence.null.dn <- function(x, N, mu) {
  return (dpois(x, 2*N*mu))
}

# model evidence of de novo data: P(x_d|H_1) 
# Input: the count data x, the sample size N, the mutation rate mu and the parameters
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
evidence.alt.dn <- function(x, N, mu, gamma.mean, beta) {
  return (dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu)))
}


# model evidence of case-control data: P(x_1,x_0|H_0) 
# Input: the count data x, the sample size N and the parameters
# Prior distribution of q|H0: Gamma(rho0, nu0)
evidence.null.CC <- function(x, N, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x$cn, rho0, nu0/(nu0+N$cn)))
  marglik0.case.log <- log(dnbinom(x$ca, rho0+x$cn, (nu0+N$cn)/(nu0+N$cn+N$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  
  return (list(ctrl=exp(marglik0.ctrl.log), case=exp(marglik0.case.log), total=exp(marglik0.log)))
}

# model evidence of case-control data: P(x_1,x_0|H_1) 
# Input: the count data x, the sample size N and the parameters
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
# Prior distribution of q|H1: Gamma(rho1, nu1)
evidence.alt.CC <- function(x, N, gamma.mean, beta, rho1, nu1, q.lower=1e-8, q.upper=0.1, debug=FALSE) {
  integrand <- function(u) {
    q <- exp(u)
    return (dnbinom(x$ca, gamma.mean*beta, beta/(beta+N$ca*q)) * dgamma(q, rho1+x$cn, nu1+N$cn) * exp(u))
  }
  
  marglik1.ctrl <- dnbinom(x$cn, rho1, nu1/(nu1+N$cn))
  marglik1.case <- integrate(integrand, lower=log(q.lower), upper=log(q.upper))$value
  
  marglik1 <- marglik1.ctrl * marglik1.case
  
  #   return (exp(marglik1.ctrl.log+marglik1.case.log))
  return (list(ctrl=marglik1.ctrl, case=marglik1.case, total=marglik1))
}

# Bayes factor of the case-control data
# BF.cn and BF.ca: contribution from control and case data, respectively
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
bayes.factor.CC <- function(x, N, gamma.mean, beta, rho1, nu1, rho0, nu0) {
  marglik0.CC <- evidence.null.CC(x, N, rho0, nu0)
  marglik1.CC <- evidence.alt.CC(x, N, gamma.mean, beta, rho1, nu1)
  
  BF.cn <- marglik1.CC$ctrl / marglik0.CC$ctrl
  BF.ca <- marglik1.CC$case / marglik0.CC$case
  BF <- BF.cn * BF.ca
  
  return (list(BF=BF, BF.cn=BF.cn, BF.ca=BF.ca))
}

# Bayes factor of de novo counts of a gene 
# x: the de novo count
# N: the sample size (number of families)
# mu: the mutation rate (of this type of mutational events)
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
bayes.factor.denovo <- function(x, N, mu, gamma.mean, beta) {
  marg.lik0 <- dpois(x, 2*N*mu)
  marg.lik1 <- dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu))
  BF <- marg.lik1/marg.lik0
  
  return (BF)
}

# Bayes factor of the gene combining de novo and case-control
# x: a list of (dn, ca, cn), counts in de novo, cases and controls
# N: a list of (dn, ca, cn), sample sizes
# hyperpar: (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0)
# Prior distribution of RR in de novo: gamma.dn ~ Gamma(gamma.mean.dn*beta.dn, beta.dn)
# Prior distribution of RR in C/C data: gamma.CC ~ Gamma(gamma.mean.CC*beta.CC, beta.CC)
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
bayes.factor <- function(x, N, mu, hyperpar, debug=FALSE) {
  gamma.mean.dn <- hyperpar[1]
  beta.dn <- hyperpar[2]
  gamma.mean.CC <- hyperpar[3]
  beta.CC <- hyperpar[4]
  rho1 <- hyperpar[5]
  nu1 <- hyperpar[6]
  rho0 <- hyperpar[7]
  nu0 <- hyperpar[8]
  
  BF.dn <- bayes.factor.denovo(x$dn, N$dn, mu, gamma.mean.dn, beta.dn)
  x.CC <- list(ca=x$ca, cn=x$cn)
  N.CC <- list(ca=N$ca, cn=N$cn)
  BF.CC <- bayes.factor.CC(x.CC, N.CC, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0)$BF
  BF <- BF.dn * BF.CC
  
  if (debug) {
    cat("BF.dn = ", BF.dn, "\n")
    cat("BF.CC = ", BF.CC, "\n")
  }
  
  return (BF)
}

#################################################################
# TADA-Denovo: analysis of de novo data 
#################################################################

# Genome-wide application of TADA for denovo data
# Input: counts, N, mu, gamma.mean, beta 
# counts: m x K matrix, where m is the number of gene, and K is the number of mutational categories. counts[i,j] is the number of de novo mutation in the j-th category of the i-th gene. 
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# gamma.mean, beta: the parameters of RR. Vector (one value for each category)
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector)
TADA.denovo <- function(counts, N, mu, mu.frac, gamma.mean, beta) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  
  # Compute BFs
  for (i in 1:m) {
    for (j in 1:K)  BF[i,j] <- bayes.factor.denovo(counts[i,j], N, mu[i]*mu.frac[j], gamma.mean[j], beta[j])
  }
  
  # Total BF per gene
  BF.total <- exp(rowSums(log(BF)))
  
  return (list(BF=BF, BF.total=BF.total))
}

# Compute permutation BFs of one gene: de novo only
# mu.gene: the mutation rates of a gene (K-dim. vector)
# N: the sample size
# l: the number of permutations
# gamma.mean, beta: RR of de novo mutations (vectors)
# Output: BF - l BFs from permutation; sample - permutate data
permute.gene.denovo <- function(mu.gene, N, l, gamma.mean, beta, dn.max=5) {
  K <- length(mu.gene)
  BF.gene <- array(1, dim=c(l,K))
  BF.gene.total <- numeric(l)
  
  # permutation of l times
  count.gene <- array(0, dim=c(l,K))
  for (j in 1:K) {
    # pre-compute the de novo table for the j-th category
    table.dn <- numeric(dn.max+1)
    for (k in 0:dn.max) {
      table.dn[k+1] <- bayes.factor.denovo(k, N, mu.gene[j], gamma.mean[j], beta[j])
    }
    
    # permutation
    count.gene[,j] <- rpois(l, 2*N*mu.gene[j])
    for (i in 1:l) {
      x <- count.gene[i,j]
      cond.range.dn <- (x <= dn.max)
      if (cond.range.dn==TRUE) {
        BF.gene[i,j] <- table.dn[x+1]
      } else {
        BF.gene[i,j] <- bayes.factor.denovo(x, N, mu.gene[j], gamma.mean[j], beta[j]) 
      }
    }
  }
  
  BF.gene.total <- exp(rowSums(log(BF.gene)))
  
  return (list(BF=BF.gene.total, sample=count.gene))
}

# Genome-wide application of TADA for denovo data: the difference with TADA.denovo is the report of p-values. 
# Input: counts, N, mu, mu.frac, gamma.mena, beta 
# counts: m x K matrix, where m is the number of gene, and K is the number of mutational categories. counts[i,j] is the number of de novo mutation in the j-th category of the i-th gene. 
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# gamma.mean, beta: the parameters of RR. Vector (one value for each category)
# l: the number of permutations to obtain the null distribution
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector). pval: the p-values of all genes. BF.null: the null distribution of BFs. 
TADAp.denovo <- function(counts, N, mu, mu.frac, gamma.mean, beta, l=100, dn.max=5) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  pval <- numeric(m)
  
  # Compute BFs
  rs <- TADA.denovo(counts, N, mu, mu.frac, gamma.mean, beta)
  BF <- rs$BF
  BF.total <- rs$BF.total
  
  # Create the null distribution
  BF.null <- numeric(m*l)
  for (i in 1:m) {
    BF.null[((i-1)*l+1):(i*l)] <- permute.gene.denovo(mu[i]*mu.frac, N, l, gamma.mean, beta, dn.max=dn.max)$BF
  }
  
  # p-values of each gene
  BF.null.srt <- sort(BF.null, decreasing=TRUE)
  pval <- findInterval(-BF.total, -BF.null.srt)/length(BF.null.srt)
  pval[pval==0] <- 0.5/length(BF.null.srt)
  
  return (list(BF=BF, BF.total=BF.total, pval=pval, BF.null=BF.null))
}

#################################################################
# TADA: analysis of de novo and inherited data
#################################################################

# Genome-wide application of TADA 
# counts: m x 3K matrix, where m is the number of gene, and K is the number of mutational categories. Each category has three numbers: de novo, case and control. 
# N: sample size, three values for de novo, case and control
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# hyperpar: 8*K matrix, where each row is (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0), and each column corresponds to one mutation type
# denovo.only: whether using only de novo data (Boolean vector)
# pi.gene: for each gene, the estimated fractions of causal variants, one for each class of variants. These fractions will be used to set gene-specific RR (case-control)
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector)
TADA <- function(counts, N, mu, mu.frac, hyperpar, denovo.only=FALSE, pi.gene=1) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  
  if (length(denovo.only)==1) { denovo.only <- rep(denovo.only, m) }
  if (length(pi.gene)==1) { pi.gene <- array(1, dim=c(m,K))}
  
  gamma.mean.dn <- hyperpar[1,]
  beta.dn <- hyperpar[2,]
  gamma.mean.CC <- hyperpar[3,]
    
  # Compute BFs: BF[i,j] is the BF of the i-th gene in the j-th category
  for (i in 1:m) {
    if (denovo.only[i]==FALSE) {
      for (j in 1:K)  {
        # set hyperparameters
        hyperpar.gene <- hyperpar[,j]
        RR.product <- hyperpar.gene[3]*hyperpar.gene[4]
        hyperpar.gene[3] <- hyperpar.gene[3]*pi.gene[i,j] + (1-pi.gene[i,j])
        hyperpar.gene[4] <- RR.product/hyperpar.gene[3]
        
        # compute BF  
        start <- 3*(j-1)+1
        x <- list(dn=counts[i, start], ca=counts[i, start+1], cn=counts[i, start+2])
        BF[i,j] <- bayes.factor(x, N, mu[i]*mu.frac[j], hyperpar.gene)
      }
    } else {
      for (j in 1:K) {
        start <- 3*(j-1)+1
        x <- counts[i,start]
        BF[i,j] <- bayes.factor.denovo(x, N$dn, mu[i]*mu.frac[j], gamma.mean.dn[j], beta.dn[j])
      } 
    }
  }
  
  # Total BF per gene
  BF.total <- exp(rowSums(log(BF)))
  
  return (list(BF=BF, BF.total=BF.total))
}

# Compute permutation BFs of one gene
# mu.gene: the mutation rates of a gene (K-dim. vector), and the case-control counts (to be permuted): vectors (one value per category)
# N: sample size, three values for de novo, case and control
# l: number of permutations
# gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0: parameters. 
# table.CC: precomputed BFs. A 3-dim table, table.CC[i, j, k] stores the BF of (j-1, k-1) in the i-th category. 
# Output: BF - l BFs from permutation; sample.dn, sample.ca, sample.cn - permutate data
permute.gene <- function(mu.gene, count.ca, count.cn, N, l, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, table.CC, dn.max=5) {
  K <- length(mu.gene)
  BF.dn <- array(1, dim=c(l,K))
  BF.CC <- array(1, dim=c(l,K))
  BF.gene.total <- numeric(l)
  ca.max <- dim(table.CC)[1] - 1
  cn.max <- dim(table.CC)[2] - 1
  
  # permutation l times for each category
  sample.dn <- array(0, dim=c(l,K))
  sample.ca <- array(0, dim=c(l,K))
  sample.cn <- array(0, dim=c(l,K))
  for (j in 1:K) {
    # pre-compute the de novo table for the j-th category
    table.dn <- numeric(dn.max+1)
    for (k in 0:dn.max) {
      table.dn[k+1] <- bayes.factor.denovo(k, N$dn, mu.gene[j], gamma.mean.dn[j], beta.dn[j])
    }
    
    # generate permutation data
    sample.dn[,j] <- rpois(l, 2*N$dn*mu.gene[j])
    sample.ca[,j] <- rhyper(l, count.ca[j] + count.cn[j], N$ca+N$cn-count.ca[j]-count.cn[j], N$ca)
    sample.cn[,j] <- count.ca[j] + count.cn[j] - sample.ca[,j]
    
    # compute the BFs
    for (i in 1:l) {
      x <- list(dn=sample.dn[i,j], ca=sample.ca[i,j], cn=sample.cn[i,j])
      cond.range.dn <- (x$dn <= dn.max)
      if (cond.range.dn==TRUE) {
        BF.dn[i,j] <- table.dn[x$dn+1]
      } else {
        BF.dn[i,j] <- bayes.factor.denovo(x$dn, N$dn, mu.gene[j], gamma.mean.dn[j], beta.dn[j]) 
      }
      
      cond.range.CC <- (x$ca <= ca.max) & (x$cn <= cn.max)
      if ( cond.range.CC==TRUE ) {
        BF.CC[i,j] <- table.CC[j, x$ca + 1, x$cn + 1]
      } else {
        BF.CC[i,j] <- bayes.factor.CC(x, N, gamma.mean.CC[j], beta.CC[j], rho1[j], nu1[j], rho0[j], nu0[j])$BF
      }
    }
  }
  
  BF.gene <- BF.dn * BF.CC
  BF.gene.total <- exp(rowSums(log(BF.gene)))
  
  return (list(BF=BF.gene.total, sample.dn=sample.dn, sample.ca=sample.ca, sample.cn=sample.cn))
}

# Genome-wide application of TADA 
# counts: m x 3K matrix, where m is the number of gene, and K is the number of mutational categories. Each category has three numbers: de novo, case and control. 
# N: sample size, three values for de novo, case and control
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# hyperpar: 8*K matrix, where each row is (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0), and each column corresponds to one mutation type
# l: the number of permutations to obtain the null distribution
# dn.max, ca.max, cn.max: if counts are below these values, the BFs will be pre-computed. 
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector). pval: the p-values of all genes. BF.null: the null distribution of BFs.
TADAp <- function(counts, N, mu, mu.frac, hyperpar, l=100, dn.max=5, ca.max=10, cn.max=10) {
  gamma.mean.dn <- hyperpar[1,]
  beta.dn <- hyperpar[2,]
  gamma.mean.CC <- hyperpar[3,]
  beta.CC <- hyperpar[4,]
  rho1 <- hyperpar[5,]
  nu1 <- hyperpar[6,]
  rho0 <- hyperpar[7,]
  nu0 <- hyperpar[8,]
  
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  pval <- numeric(m)
  
  # Compute BFs of the observed data
  rs <- TADA(counts, N, mu, mu.frac, hyperpar)
  BF <- rs$BF
  BF.total <- rs$BF.total
  
  # Pre-compute the bayes-factors of the case-control data
  table.CC <- array(1, dim=c(K, (ca.max+1), (cn.max+1)))
  for (j in 1:K) {
    for (x1 in 0:ca.max) {
      for (x0 in 0:cn.max) {
        x <- list(ca=x1,cn=x0)
        table.CC[j, x1+1, x0+1] <- bayes.factor.CC(x, N, gamma.mean.CC[j], beta.CC[j], rho1[j], nu1[j], rho0[j], nu0[j])$BF
      }
    }
  }
  
  # Create the null distribution
  BF.null <- numeric(m*l)
  for (i in 1:m) {
#     print(i)
    BF.null[((i-1)*l+1):(i*l)] <- permute.gene(mu[i]*mu.frac, counts[i, 3*(1:K)-1], counts[i, 3*(1:K)], N, l, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, table.CC, dn.max=dn.max)$BF
  }
  
  # p-values of each gene
  BF.null.srt <- sort(BF.null, decreasing=TRUE)
  pval <- findInterval(-BF.total, -BF.null.srt)/length(BF.null.srt)
  pval[pval==0] <- 0.5/length(BF.null.srt)
  
  return (list(BF=BF, BF.total=BF.total, pval=pval, BF.null=BF.null))
}

#################################################################
# MOM estimation of hyperprior parameters from de novo data
#################################################################

# Prob. of having d or more de novo mutations under H1 
# Use simulation, but could also use analytic form 
multihit.prob <- function(N, mu, gamma.mean, beta, d=2, S=100) {
  p <- numeric(S)
  gamma <- rgamma(S, gamma.mean*beta, rate=beta)
  for (i in 1:S) {
    p[i] <- 1 - ppois(d-1, 2*N*mu*gamma[i])
  }
  return (mean(p))
}

# Estimate the number of multihit genes in a genome. 
# d: the parameter of the multiplicity test. 
# Returns: M0 - the number of multihit genes from non-risk genes; M1 - the number from risk genes. 
count.multihit <- function(N, mu, pi, gamma.mean, beta, d=c(2,3), S=2) {
  m <- length(mu)
  M0 <- numeric(length(d))
  M1 <- numeric(length(d))
  
  # M1: the number of causal genes having d or more de novo mutations
  p.alt <- array(0, dim=c(m, length(d)))  # p1[i,j] = P(X_i >= d_j|H1)
  for (i in 1:m) {
    for (j in 1:length(d)) {
      p.alt[i,j] <- multihit.prob(N, mu[i], gamma.mean, beta, d=d[j], S=S)
    }
  }
  for (j in 1:length(d)) { 
    M1[j] <- m * pi  * mean(p.alt[,j]) 
  }
  
  # M0: the number of non-causal genes having d or more de novo mutations
  p.null <- array(0, dim=c(m, length(d)))  # p0[i,j] = P(X_i >= d_j|H0)
  for (i in 1:m) {
    for (j in 1:length(d)) {
      p.null[i,j] <- 1 - ppois(d[j] - 1, 2*N*mu[i])
    }
  }
  for (j in 1:length(d)) { 
    M0[j] <- m * (1-pi) * mean(p.null[,j]) 
  }
  
  result <- data.frame(d=d, M0=M0, M1=M1)
  return (result)
}

# Estimating relative risk and the number of multiple hits from de novo data
# Input: sample size (N), mutation rates of all genes (mu), observed number of de novo events (C), beta (parameter of the prior distribution of gamma), k (number of disease genes)
# Output: the average relative risk (gamma.mean), the expected number of multi-hit genes (M)
denovo.MOM <- function(N, mu, C, beta, k) {
  m <- length(mu) # number of genes
  
  # enrichment of de novo events
  nu <- C / (2 * N * sum(mu))
  
  # MOM estimator of gamma.mean
  gamma.mean <- (nu-1)*m/k +1
  
  # expected M (choose d = 2)
  rs <- count.multihit(N, mu, k/m, gamma.mean, beta, d=2)    
  M <- sum(rs$M1) + sum(rs$M0)
  
  return (list(gamma.mean=gamma.mean, M=M))
}

#################################################################
# Empirical Bayes estimation of hyperprior parameters
#################################################################

# Evalaute the marginal log-likelihood at given parameters
# pi: the fraction of causal genes
# counts: the count data (of one mutational category), a date frame
# prior.weight: putting a prior so that q1 and q0 tend to be close (recommended value: 5000). Default 0. 
# Output: the negative log-likelihood, and posterior, BF of each gene
marginal <- function(hyperpar, pi, counts, N, mu, prior.weight=0, denovo.only=FALSE, debug=FALSE) {
  gamma.mean.dn <- hyperpar[1]
  beta.dn <- hyperpar[2]
  gamma.mean.CC <- hyperpar[3]
  beta.CC <- hyperpar[4]
  rho1 <- hyperpar[5]
  nu1 <- hyperpar[6]
  rho0 <- hyperpar[7]
  nu0 <- hyperpar[8]
  
  n <- nrow(counts)
  prob <- numeric(n)
  posterior <- numeric(n)
  bf <- numeric(n)
  for (i in 1:n)  {
    if (debug) cat("i = ", i, "\tdn = ", counts[i,]$dn, "\tca = ", counts[i,]$ca, "\tcn = ", counts[i,]$cn, "\n")
    if (denovo.only==TRUE) {
      prob.M1 <- evidence.alt.dn(counts[i,"dn"], N$dn, mu[i], gamma.mean.dn, beta.dn) 
      prob.M0 <- evidence.null.dn(counts[i,"dn"], N$dn, mu[i])
    } else {
      prob.M1 <-  evidence.alt(counts[i,],N,mu[i],gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1)
      prob.M0 <-  evidence.null(counts[i,],N,mu[i],rho0, nu0) 
    }
    if (debug) cat("prob.M1 = ", prob.M1,"\n")
    if (debug) cat("Prob.M0 = ", prob.M0, "\n")
    prob[i] <- pi * prob.M1 + (1 - pi) * prob.M0
    posterior[i] <- pi*prob.M1 / (prob[i])
    bf[i] <- prob.M1 / prob.M0
  }
  marg.nll <- -sum(log(prob[!is.na(prob)]))
  
  # add the prior
  if (prior.weight > 0) {
    u <- prior.weight
    q1.mean <- rho1/nu1
    q0.mean <- rho0/nu0
    marg.nll <- marg.nll - log(dgamma(q1.mean, q0.mean*u, u)) 
  }
  
  result <- list(prob=prob, marginal=marg.nll, posterior=posterior, bf=bf)
  return (result)
}

# Convert the full set of parameters (pi, hyperpar) to a subset based on the option
# est.option: a Boolean vector (size 8), one for each hyperpar. If FALSE, the corresponding parameter will be fixed. 
# est.pi: Boolean. 
fullpar2subpar <- function(hyperpar, pi, est.option, est.pi) {
  subpar <- NULL
  for (i in 1:length(est.option)) {
    if (est.option[i]) subpar <- c(subpar, hyperpar[i])
  }
  
  if (est.pi==TRUE)  subpar <- c(subpar, pi)
  
  return (subpar)
}

# Convert the subset of parameters (subpar) to the full set based on the option and intial values
# hyperpar.init, pi.init: used to set parameters not to be estimated
subpar2fullpar <- function(subpar, hyperpar.init, pi.init, est.option, est.pi) {
  hyperpar <- numeric(length(hyperpar.init))
  curr <- 1
  for (i in 1:length(est.option)) {
    if (est.option[i]) { hyperpar[i] <- subpar[curr]; curr <- curr+1 }
    else hyperpar[i] <- hyperpar.init[i]
  }
  if (est.pi==TRUE)  pi <- subpar[curr]
  else pi <- pi.init
  
  return (list(hyperpar=hyperpar, pi=pi))
}

# Empirical Bayes estimation of hyperparameters 
# counts: the count data - the number of mutations in cases, controls and de novo data respectively. 
# N: sample size (dn, ca, cn, respectively). 
# mu: vector of mutation rates
# lower, upper: the search range of the parameters 
# est.option: a Boolean vector specifying whether to estimate each parameter. If FALSE, use the corresponding parameter in hyperpar.init
# est.pi: whether to estimate pi. If FALSE, use the given value of pi; if TRUE, use the given value as initial
# prior.weight: putting a prior so that q1 and q0 tend to be close (recommended value: 5000). Default 0. 
# Output: the hyperparameters and the NLL (negative log likelihood) at the parameters. 
empBayes <- function(counts, N, mu, hyperpar.init, pi.init, lower, upper, lower.pi=1e-10, upper.pi=1, est.option=rep(FALSE, 8), est.pi=FALSE, prior.weight=0, debug=FALSE) {  
  # marginal likelihood function (parameters in log-scale)
  marginal.loglike <- function(subpar.log) {
    subpar <- exp(subpar.log)
    allpar <- subpar2fullpar(subpar, hyperpar.init, pi.init, est.option, est.pi)
    hyperpar <- allpar$hyperpar
    pi <- allpar$pi
    result <- marginal(hyperpar,  pi, counts, N, mu, prior.weight=prior.weight)$marginal  
    
    if (debug) {
      cat("Parameters = (", paste(subpar, collapse=", "), ")\tValue = ", result, "\n", sep="") 
    }
    
    return (result)
  }
  
  # no parameter to estimate
  if ( sum(est.option==TRUE) == 0 & est.pi == FALSE ) {
    value <- marginal(hyperpar.init, pi.init, counts, N, mu, prior.weight=prior.weight)$marginal
    return (list(hyperpar=hyperpar.init, value=value))
  }
  
  # initialize the sub. parameters
  subpar.init <- fullpar2subpar(hyperpar.init, pi.init, est.option, est.pi)
  sublower <- fullpar2subpar(lower, lower.pi, est.option, est.pi)
  subupper <- fullpar2subpar(upper, upper.pi, est.option, est.pi)
  
  # maximization of marginal likelihood
  like.optim <- optim(log(subpar.init), marginal.loglike, method="L-BFGS-B", lower=log(sublower), upper=log(subupper))
  subpar.est.log <- like.optim$par
  value.max <- like.optim$value
  allpar.est <- subpar2fullpar(exp(subpar.est.log), hyperpar.init, pi.init, est.option, est.pi)
  hyperpar.est <- allpar.est$hyperpar
  pi.est <- allpar.est$pi
  
  if (est.pi==TRUE) { return (list(pi=pi.est, hyperpar=hyperpar.est, value=value.max)) }
  else { return (list(hyperpar=hyperpar.est, value=value.max)) }
}

# Empirical Bayes estimation of hyperparameters (gamma.mean, beta and pi) for de novo data
# counts: count data - the number of de novo mutations. 
# N - sample size (trios). 
# mu: vector of mutation rate
# lower, upper: the search range of the parameters 
empBayes.denovo <- function(counts, N, mu, gamma.mean.init, beta.init, pi.init, lower, upper, lower.pi=1e-10, upper.pi=1, debug=FALSE) {  
  # marginal likelihood function (parameters in log-scale)
  marginal.loglike <- function(par.log) {
    par <- exp(par.log)
    gamma.mean <- par[1]
    beta <- par[2]
    pi <- par[3]
    hyperpar <- c(gamma.mean, beta, 0, 0, 0, 0, 0, 0)
    m <- length(counts)
    counts.full <- data.frame(dn=counts, ca=numeric(m), cn=numeric(m))
    N.full <- list(dn=N, ca=0, cn=0)
    result <- marginal(hyperpar, pi, counts.full, N.full, mu, denovo.only=TRUE)$marginal  
    
    if (debug) {
      cat("Parameters = (", paste(par, collapse=", "), ")\tValue = ", result, "\n", sep="") 
    }
    
    return (result)
  }
  
  # initialize the parameters
  par.init <- c(gamma.mean.init, beta.init, pi.init)
  par.lower <- c(lower, lower.pi)
  par.upper <- c(upper, upper.pi)
  
  # maximization of marginal likelihood
  like.optim <- optim(log(par.init), marginal.loglike, method="L-BFGS-B", lower=log(par.lower), upper=log(par.upper))
  par.est <- exp(like.optim$par)
  value.max <- like.optim$value
  gamma.mean.est <- par.est[1]
  beta.est <- par.est[2]
  pi.est <- par.est[3]
  
  return (list(gamma.mean=gamma.mean.est, beta=beta.est, pi=pi.est, value=value.max)) 
}

#################################################################
# Functions for simulation
#################################################################

# Generate simulation data of a set of genes (multiple mutational categories): de novo mutations only
# N: sample size (number of trios)
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean, beta: Relative risk of de novo mutation: gamma|M1 ~ Gamma(gamma.mean*beta, beta). Vectors (one per category) 
# Output: sample matrix (m by K), where m is the number of genes and K the number of variant categories. sample.info: more information of the samples, including the indicator (risk gene or not) and the RR. 
simulator.denovo <- function(N, mu, mu.frac, pi, gamma.mean, beta) {
  m <- length(mu) # number of genes
  K <- length(mu.frac) # number of mutational categories
  
  z <- rbinom(m, 1, pi)
  gamma <- array(1, dim=c(m,K))
  x <- array(0, dim=c(m,K))
  k <- sum(z==1)
  for (j in 1:K) {
    gamma[z==1, j] <- rgamma(k, gamma.mean[j]*beta[j], beta[j])
    x[,j] <- rpois(m, 2 * mu * mu.frac[j] * gamma[,j] * N)
  }
  
  sample.info <- cbind(mu, z, gamma, x)
  
  return (list(sample=x, sample.info=sample.info))
}

# Generate simulation data of a set of genes (multiple mutational categories)
# N: sample size (number of trios)
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean.dn, beta.dn: Relative risk of de novo mutation: gamma|M1 ~ Gamma(gamma.mean.dn*beta.dn, beta.dn). Vectors.
# gamma.mean.CC, beta.CC: Relative risk of inherited mutation (case/control): gamma.CC|M1 ~ Gamma(gamma.mean.CC*beta.CC, beta.CC). Vectors
# Frequency parameter of risk genes: q|M1 ~ Gamma(rho1, nu1)
# Frequency parameter of non-risk genes: q|M0 ~ Gamma(rho0, nu0)
# tradeoff option: if TRUE, implement q-gamma tradeoff (i.e. higher gamma means lower q). Suppose, gamma_i is the RR, then q_i is proportional to mu_i / gamma_i, where the constant is determined from the mean of q, mu and gamma. 
# Output: sample matrix (m by 3K), where m is the number of genes and K the number of variant categories.  sample.info: more information of the samples, including the indicator (risk gene or not) and the RR.
simulator <- function(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, tradeoff=FALSE) {
  m <- length(mu) # number of genes
  K <- length(mu.frac) # number of mutational categories
  
  # the tradeoff parameter (delta:=mu.mean/q.mean)
  delta <- mean(mu) * mu.frac / (rho0 / nu0)
  
  z <- rbinom(m, 1, pi)
  gamma.dn <- array(1, dim=c(m,K))
  gamma.CC <- array(1, dim=c(m,K))
  q <- array(0, dim=c(m,K))
  x <- array(0, dim=c(m,3*K))
  k <- sum(z==1)
  for (j in 1:K) {
    # sample de novo 
    gamma.dn[z==1, j] <- rgamma(k, gamma.mean.dn[j]*beta.dn[j], beta.dn[j])
    col <- 3*(j-1)+1
    x[,col] <- rpois(m, 2 * mu * mu.frac[j] * gamma.dn[,j] * N$dn)
    
    # sample case-control
    gamma.CC[z==1, j] <- rgamma(k, gamma.mean.CC[j]*beta.CC[j], beta.CC[j])
    q[z==0, j] <- rgamma(m-k, rho0[j], nu0[j])
    if (tradeoff==FALSE) {
      q[z==1, j] <- rgamma(k, rho1[j], nu1[j])
    } else {
      q[z==1, j] <- mu[z==1] * mu.frac[j] / (delta[j] * gamma.CC[z==1, j])
    }
    x[,col+1] <- rpois(m, q[,j] * gamma.CC[,j] *N$ca)
    x[,col+2] <- rpois(m, q[,j] * N$cn)
    
  }
  
  sample.info <- cbind(mu, z, gamma.dn, gamma.CC, q, x)
  
  return (list(sample=x, sample.info=sample.info))
}

# Evaluation of the power (FDR) of TADA.denovo
# N: sample size, the number of families of the de novo study
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean, beta: RR parameters for de novo mutations (vectors)
# gamma.mean.est, beta.est: the parameters used by TADA.denovo
# FDR: the FDR level to be controled
# Output: the number of discoveries at the specified FDR level
eval.TADA.denovo <- function(N, mu, mu.frac, pi, gamma.mean, beta, gamma.mean.est, best.est, FDR=0.1) {
  # sample the simulation data
  sample <- simulator.denovo(N, mu, mu.frac, pi, gamma.mean, beta)$sample
  
  # run TADA.denovo
  sampleBF <- TADA.denovo(sample, N, mu, mu.frac, gamma.mean.est, beta.est)$BF.total
  
  # Bayesian FDR control
  M1 <- Bayesian.FDR(sort(sampleBF, decreasing=TRUE), 1-pi, FDR)$ND
  
  return (M1)
}

# Evaluation of the power (FDR) of TADA
# N: sample size, the number of families of the de novo study
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean.dn, beta.dn: RR parameters for de novo mutations (vectors)
# gamma.mean.CC, beta.CC: RR parameters for inherited mutations (vectors)
# Frequency parameter of risk genes: q|M1 ~ Gamma(rho1, nu1)
# Frequency parameter of non-risk genes: q|M0 ~ Gamma(rho0, nu0)
# hyperpar.est: the parameters used by TADA
# FDR: the FDR level to be controled
# Output: the number of discoveries at the specified FDR level
eval.TADA <- function(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, hyperpar.est, FDR=0.1, tradeoff=FALSE) {
  # sample the simulation data
  sample <- simulator(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, tradeoff=tradeoff)$sample
  
  # run TADA
  sampleBF <- TADA(sample, N, mu, mu.frac, hyperpar.est)$BF.total
  
  # Bayesian FDR control
  M1 <- Bayesian.FDR(sort(sampleBF, decreasing=TRUE), 1-pi, FDR)$ND
  
  # sample information (counts and BFs)
  sim.results <- data.frame(counts=sample, BF=sampleBF)
  sim.results <- sim.results[order(-sim.results$BF),]
  
  return (list(M1=M1, sample=sim.results))
}
estimatePars <- function(pars, mcmcResult, nThin = NULL){
    mcmcDataFrame <- as.data.frame(mcmcResult)
        if (!is.null(nThin))
        mcmcDataFrame <- mcmcDataFrame[seq(1, dim(mcmcDataFrame)[1], by = nThin),]
    allPars <- colnames(mcmcDataFrame)

    if (is.null(pars)){
        message("\nYou are estimating all parameters, you can input a specific parameter(s) for pars\n")
        pars = allPars
    }

#    pars <- pars[grep("hyper|pi|alpha", pars)]
    message("====\nOnly pi, alpha and hyper parameters are estimated in this step\n",
            "The method does not calculate HPDs for hyper betas, just their medians\n===\n")

    if (length(pars[!is.element(pars, colnames(mcmcDataFrame))]) > 0)
        warning((pars[!is.element(pars, colnames(mcmcDataFrame))]), " is/are not in mcmc chains")
    pars <- pars[is.element(pars, colnames(mcmcDataFrame))]

    hpdList <- NULL
    for (iPar in pars) {
       # message("Estimating for ", iPar)
        xData <-  mcmcDataFrame[, iPar]
        if ((sd(xData) > 10^-10) & (length(grep("Beta", iPar)) == 0))
            hpdList <- rbind(hpdList, loc1stats(xData)[1:3])
        else
            hpdList <- rbind(hpdList, rep(median(xData), 3))
    }
    rownames(hpdList) <- pars

    colnames(hpdList) <- c("Mode", "lCI", "uCI")

return(hpdList)
}
mTADA <- function(geneName = NULL,
                  ntrio1 = NULL,
                  ntrio2 = NULL,
                  p1 = NULL,
                  p2 = NULL,
                  dataDN1 = NULL, #data.matrix(sDN[, 1:2]), #array(sDN, dim = c(NN, NCdn)), ##De novo data
                  dataDN2 = NULL, #data.matrix(sDN[, 3:4]), #array(sDN, dim = c(NN, NCdn)), ##De novo data
                      mutRate1 = NULL, #data.matrix(muAll[, 1:2]), ##Mutation rate data
                      mutRate2 = NULL, #data.matrix(muAll[, 3:4]), ##Mutation rate data
                  hyperGammaMeanDN1 = NULL, #c(gDN1),
                      hyperGammaMeanDN2 = NULL, #c(gDN2),
                  beta.dn = NULL,
                       adjustHyperBeta = 0, #adjustHyperBeta0,
                        betaPars = c(6.7771073, -1.7950864, -0.2168248),
                        upperPi0 = 0.5, lowerPi0 = 0,
                        lowerHyperGamma = 1,
                        lowerGamma = 1,
                        lowerBeta = 1,
                        hyperBetaDN01 = NULL, #array(c(1, 1)),
                        hyperBetaDN02 = NULL, #array(c(1, 1)),
                  nIteration = 5000,
                  outSample = 1000,
                  nChain = 1,
                  nCore = NULL,
                  printMessage = FALSE, #
                  useMCMC = FALSE,
                  iSeed = NULL
                  ){
    dataDN1 <- data.frame(dataDN1)
    dataDN2 <- data.frame(dataDN2)

    NCdn1 <- dim(dataDN1)[2]
    NCdn2 <- dim(dataDN2)[2]
    Ndn1 <- array(rep(ntrio1, NCdn1))
    Ndn2 <- array(rep(ntrio2, NCdn2))

    if ((is.null(nCore)) & (useMCMC)){
        message("No information for core numbers (nCore); therefore, nCore = nChain: ", nChain, " core(s) is/are used\n")   
        nCore = nChain
    }
    if (is.null(hyperBetaDN01))
        hyperBetaDN01 <- rep(1, NCdn1)
    if (is.null(hyperBetaDN02))
        hyperBetaDN02 <- rep(1, NCdn2)
                   
    mixDataKclasses <- list(NN = dim(dataDN1)[1], #Number of genes
                        NCdn1 = NCdn1, NCdn2 = NCdn2, 
                        Ndn1 = Ndn1, Ndn2 = Ndn2, 
                        dataDN1 = data.matrix(dataDN1),
                        dataDN2 = data.matrix(dataDN2),
                        hyperGammaMeanDN1 = array(hyperGammaMeanDN1),
                      hyperGammaMeanDN2 = array(hyperGammaMeanDN2), 
                      mutRate1 = data.matrix(mutRate1), 
                      mutRate2 = data.matrix(mutRate2), 
                        adjustHyperBeta = adjustHyperBeta,
                        betaPars = betaPars,
                        upperPi0 = upperPi0, lowerPi0 = lowerPi0,
                        lowerHyperGamma = lowerHyperGamma,
                        lowerGamma = lowerGamma,
                        lowerBeta = lowerBeta,
                        hyperBetaDN01 = array(hyperBetaDN01),
                        hyperBetaDN02 = array(hyperBetaDN02),
                        pi01 = p1,
                        pi02 = p2                      )

    library("rstan")
    library("locfit")
    if (is.null(iSeed))
      iSeed <- sample.int(.Machine$integer.max, 1)
    t1M <- min(p1, p2)

    initList = NULL
    for (ii in 1:nChain) {
        tempB <- runif(1, 0.01*t1M, 0.9*t1M)
        initList[[ii]] <- list(p12 = tempB) #alpha0 = log(tempB/(1 - tempB)))
        }
    initList

    message("===================\nBuilding the model\n=================\n")

    m1 = stan_model(model_code = DN2traits)

    if (useMCMC){
        message("\n=======Use MCMC===========\n")
        vb1 <- stan(model_code=DN2traits, data = mixDataKclasses,
                    init = initList, iter = nIteration, chains = nChain,
                    cores = nCore,
                    pars = c('p12', 'gammaMeanDN1'), 
                seed = iSeed)
    } else {
        message("\n================Not use MCMC===========\n")
        message("===================Running VB========================\n")
        
    vb1 <- vb(object = m1, data = mixDataKclasses, pars = c('p12', 'gammaMeanDN1'),
          init = initList[[1]], iter = nIteration, seed = iSeed)

    }

    par2 = estimatePars(pars = c("p12", "gammaMeanDN1[1]"),
             mcmcResult = vb1)

    if (as.numeric(par2[1, 1]) < 10^-4){
        message('pi3 is very small; therfore, mTADA considers no overlapping information and set pi3 =0\n')
        par2[1, 1] <- 0        
    }
    
    if (printMessage)
        print(par2)
 ########################Calculate PP
    g0 = rbind(rep(1, NCdn1 + NCdn2),
               c(hyperGammaMeanDN1, hyperGammaMeanDN2),
               c(hyperGammaMeanDN1, rep(1, NCdn2)),
               c(rep(1, NCdn1), hyperGammaMeanDN2))

    gamma.mean.dn = g0

    dnData = data.frame(dataDN1, dataDN2)
    muAll = data.frame(mutRate1, mutRate2)
    Ndn = c(Ndn1, Ndn2)

    if (is.null(beta.dn))
        beta.dn = matrix(rep(1, 4*(NCdn1 + NCdn2)), nrow = 4)
    par1 = par2[, 1]
    pBoth = par1['p12']
    p00 = c(1 -(p1 + p2 - pBoth), pBoth, p1 - pBoth, p2 - pBoth)
    prob0 = p00

#    if (!(identical(g0, beta.dn))){
 #       message("\nMean gammas and betas have different dimension\n")
#       print(g0)
 #       print(beta.dn)
  #  }
    
#    geneName = data[, 1]
    if (is.null(geneName))
        geneName = paste0("Gene", 1:dim(dnData)[1])
    x1 <- posProb.dn(dnData = dnData, muAll = muAll,
              gamma.mean.dn = gamma.mean.dn, Ndn = Ndn,
              prob0 = prob0, beta.dn = beta.dn)

    d11 <- data.frame(geneName = geneName, dnData, x1$PP) #, xGroup)

    pOut = p00[c(1, 3, 4, 2)]
    names(pOut) <- c("pNO", "pFIRST", "pSECOND", "pBOTH")
    return(list(probModel = pOut,
                pars = par2, data = d11, mcmcData = vb1))
    
}
##Hoang Nguyen: Jan 12, 2017
##If users have prior information for pi3 (p12 in this script) then they can set prior information for this parameter

##gO: gene-level genetic overl
##gO = p12/(pi01 + pi02 - p12)

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


#An example of gamma.mean.dn:
#1                   24.81551                         82.82128
#2                   24.81551                          1.00000
#3                    1.00000                         82.82128

#dataDN, muAll, Ndn should be across categories and 2 traits
#probs = NULL, BOTH, FIRST, SECOND

posProb.dn  <- function(dnData, muAll, gamma.mean.dn, Ndn, prob0, beta.dn){
    nGdn <- dim(dnData)[2]
    m <- dim(dnData)[1]
    bF <- array(0, dim = c(m, 4)) ##The first column is for NULL model
    ###############Calculate BF for all categories
    for (j in 1:nGdn){
#        message("j = ", j)
###NULL
        marg.lik0 <- dpois(dnData[, j], 2*Ndn[j]*muAll[, j])
#        marg.lik0 <- dnbinom(dnData[, j], beta.dn[i, j], beta.dn[i, j]/(beta.dn[i, j] + 2*Ndn[j]*muAll[, j]))
        ###Both traits, first trait and then second trait
        marg.lik1 <- array(1, dim = c(m, 4))
#        bF[, 1] <-  bF[, 1] + log(marg.lik0) #For NULL model
        for (i in 1:4){ ##Add NULL
            if (gamma.mean.dn[i, j] == 1)
                marg.lik1[, i] <- dpois(dnData[, j], 2*Ndn[j]*muAll[, j])
            else
                marg.lik1[, i] <- dnbinom(dnData[, j], gamma.mean.dn[i, j]*beta.dn[i, j], beta.dn[i, j]/(beta.dn[i, j] + 2*Ndn[j]*muAll[, j]))
            bF[, i] <- bF[, i ] + log(marg.lik1[, i])# - log(marg.lik0)
            ##bF = 1:4, but calculate for 2:4
        }
#        print(j)
        }
    ##########################
    bF <- exp(bF)
   
    PP <- t(apply(bF, 1, function(x)        x*prob0/sum(x*prob0)))
    colnames(PP) <- c("NO", "BOTH", "FIRST", "SECOND")
    return(list(BF = bF, PP = PP))

}
library(locfit)
loc2plot <- function(x,y,cprob=0.5, xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

# finds the mode for a bivariate density
loc2mode <- function(x,y,alpha=0.5,xlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	tt <- max(fitted(fit))
	wt <- fitted(fit) == tt
	c(x[wt][1],y[wt][1])
}

# this works for univariate data; gives a vector with
# mode (global), hpd1_low, hpd1_high, hpd2_low, hpd2_hi, etc
#The reason for multiple hpd limits is if you have multimodal data.
#prob = prob *outside* the limit; i.e for a normal distribution 0.05 is expected to give
#     the 0.025 and 0.975 quantiles.
# this won't work for weighted data, use loc1statsx instead.
# xlim is optional - use it to define the support of the density.
loc1stats <- function(x,xlim,prob=0.05,...)
{
	if(missing(xlim)){
		fit <- locfit(~x)
	}
	else {
		fit <- locfit(~x,xlim=xlim)
	}
	fx <- fitted(fit)
	x.modef <- max(fx)
	x.mode <- x[fx == x.modef]
	if(!missing(xlim)){
		if(predict(fit,xlim[1]) > x.modef){
			x.modef <- predict(fit,xlim[1])
			x.mode <- xlim[1]
		}
		if(predict(fit,xlim[2]) > x.modef){
			x.modef <- predict(fit,xlim[2])
			x.mode <- xlim[2]
		}
	}

	if(length(x.mode)>1)x.mode <- x.mode[1]
	lev <- sort(fx)[floor(prob*length(x))]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	indx <- order(x)
	ii <- 2
	flip <- TRUE
	for(j in 2:length(x)){
		if(flip && fx[indx[j]] > lev){
			l1[[ii]] <- x[indx[j-1]]
			if(j==2 && !missing(xlim)){
				if(predict(fit,xlim[1]) >= lev)l1[[ii]] <- xlim[1]
			}
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && fx[indx[j]] < lev){
			l1[[ii]] <- x[indx[j]]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip && !missing(xlim) && j == length(x)){
			l1[[ii]] <- xlim[2]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}








loc2plot.old <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}


#gives the HPD value for an observation px,py in a density constructed from x, y.
gethpdprob2 <- function(x,y,px,py,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
#	d1 <- (x-px)^2+(y-py)^2
#	best <- d1 == min(d1)
#	lev <- mean(fitted(fit)[best])
	lev <- predict.locfit(fit,list(px,py))
	slev <- sort(fitted(fit))
	indic <- slev <= lev
	sum(indic)/length(x)
}


plotParHeatmap1 <- function(pars, mcmcResult, nThin = NULL, color = "blue",
                           xLim = NULL, yLim = NULL,
                           changeXaxis = FALSE,
                           parXaxis = 19358,
                           mainLab = NULL, xLab = NULL, yLab = NULL,
                           maxk0 = 500, cprob = c(0.5, 0.05),
                           transFunction = NULL){
    mcmcDataFrame <- as.data.frame(mcmcResult)
    if (!is.null(nThin))
        mcmcDataFrame <- mcmcDataFrame[seq(1, dim(mcmcDataFrame)[1], by = nThin),]
    allPars <- colnames(mcmcDataFrame)


        
    x <- mcmcDataFrame[, pars[1]]
        if (pars[1] == 'alpha0')
            x <- exp(x)/(exp(x) + 1)
    if (changeXaxis){
        x <- x*parXaxis
        xLim = range(x)
    }
    y <- mcmcDataFrame[, pars[2]]
    sc1<-sd(x)
    sc2<-sd(y)
    fit <- locfit(~x+y,scale=c(sc1,sc2), maxk = maxk0)
    lev <- sort(fitted(fit))[floor(cprob*length(x))]
    max.i <- fitted(fit)==max(fitted(fit))
    mode.x<-x[max.i][[1]]
    mode.y<-y[max.i][[1]]
    if (is.null(yLim))
        yLim <- c(0, max(y))
    if (is.null(xLim))
        xLim <- c(0, max(x))
    if (is.null(xLab))
        xLab <- pars[1]
    if (is.null(yLab))
        yLab <- pars[2]

    plot(mode.x,mode.y,
     #xaxt="n",
         main = mainLab, #diseaseName,
         xlim = xLim, ylim = yLim,
         pch=-1, xlab = xLab, ylab = yLab)

my.color.palette=colorRampPalette(c("white", color), space = "Lab")
plot(fit, type="image", m=300, add=TRUE,
     col=my.color.palette(25)[-c(2,3,5,6,8,10,12)])

plot(fit,add=TRUE,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""), labcex=0.75,vfont=c("sans serif","bold"))
points(mode.x,mode.y,pch=3)

}
posProb.dnOneTrait <- function(dnData, muAll, gamma.mean.dn,
                                Ndn, prob0, beta.dn){
    nGdn <- dim(dnData)[2]
    m <- dim(dnData)[1]
    bF <- array(0, dim = c(m, 2))
    ###############Calculate BF for all categories
    for (j in 1:nGdn){
#        message("j = ", j)
###NULL
        marg.lik0 <- dpois(dnData[, j], 2*Ndn[j]*muAll[, j])
###Both traits, first trait and then second trait
        marg.lik1 <- dnbinom(dnData[, j], gamma.mean.dn[j]*beta.dn[j],
                                  beta.dn[j]/(beta.dn[j] + 2*Ndn[j]*muAll[, j]))
        bF[, j] <- log(marg.lik1) - log(marg.lik0)
        
        }        
    ##########################
    bF <- exp(bF)
    bF <- apply(bF, 1, prod)
     PP <- bF1*prob0/(1 - prob0 + prob0*bF1)
    
    return(list(BF = bF, PP = PP))

}
