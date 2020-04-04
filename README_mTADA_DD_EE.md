This notebook descibes steps used to jointly analyze two traits by
`mTADA`.

I. Introduction
---------------

**`mTADA` jointly analyze de novo mutations (DNMs) of two traits to 1)
estimate the gene-level genetic overlap of the two traits; 2) report
shared and specific risk genes; and 3) identify additional risk genes
for each analyzed trait.**

The method requires genetic parameters from single-trait analyses (the
third and fourth columns in Table 1 below). Users can obtain
single-trait parameters from `extTADA/TADA` methods.

**Table 1. `mTADA` model for one variant category at the
*i*<sup>*t**h*</sup> gene**.

<table>
<colgroup>
<col style="width: 23%" />
<col style="width: 21%" />
<col style="width: 25%" />
<col style="width: 29%" />
</colgroup>
<thead>
<tr class="header">
<th><strong>Hypothesis</strong></th>
<th>Proportion</th>
<th>First trait</th>
<th>Second trait</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><span class="math inline"><em>H</em><sub>0</sub></span></td>
<td><span class="math inline"><em>π</em><sub>0</sub></span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>1</sub></span> ~ <span class="math inline"><em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>1</sub><em>μ</em><sub><em>i</em></sub>)</span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>2</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>2</sub><em>μ</em><sub><em>i</em></sub>)</span></td>
</tr>
<tr class="even">
<td><span class="math inline"><em>H</em><sub>1</sub></span></td>
<td><span class="math inline"><em>π</em><sub>1</sub></span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>1</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>1</sub><em>γ</em><sub>1</sub><em>μ</em><sub><em>i</em></sub>)</span>; <span class="math inline"><em>γ</em><sub>1</sub></span> ~ <span class="math inline">$Gamma(\bar{\gamma_1}\beta_1,\beta_1)$</span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>2</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>2</sub><em>μ</em><sub><em>i</em></sub>)</span></td>
</tr>
<tr class="odd">
<td><span class="math inline"><em>H</em><sub>2</sub></span></td>
<td><span class="math inline"><em>π</em><sub>2</sub></span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>1</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>1</sub><em>μ</em><sub><em>i</em></sub>)</span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>2</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>2</sub><em>γ</em><sub>2</sub><em>μ</em><sub><em>i</em></sub>)</span>; <span class="math inline"><em>γ</em><sub>2</sub></span> ~ <span class="math inline">$Gamma(\bar{\gamma_2}\beta_2,\beta_2)$</span></td>
</tr>
<tr class="even">
<td><span class="math inline"><em>H</em><sub>3</sub></span></td>
<td><span class="math inline"><em>π</em><sub>3</sub></span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>1</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>1</sub><em>γ</em><sub>1</sub><em>μ</em><sub><em>i</em></sub>)</span>; <span class="math inline"><em>γ</em><sub>1</sub></span> ~ <span class="math inline"><em>G</em><em>a</em><em>m</em><em>m</em><em>a</em>(<em>γ̄</em><sub>1</sub><em>β</em><sub>1</sub>, <em>β</em><sub>1</sub>)</span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>2</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>2</sub><em>γ</em><sub>2</sub><em>μ</em><sub><em>i</em></sub>)</span>; <span class="math inline"><em>γ</em><sub>2</sub></span> ~ <span class="math inline"><em>G</em><em>a</em><em>m</em><em>m</em><em>a</em>(<em>γ̄</em><sub>2</sub><em>β</em><sub>2</sub>, <em>β</em><sub>2</sub>)</span></td>
</tr>
</tbody>
</table>

**Figure 1. `mTADA` framework.**

![](data/mTADAmethod.jpeg)

#### Data for reproducible analyses

Data used in the main manuscript are inside the folder [data](data):

1.  [FullDataSet\_DenovoMutations\_for\_mTADA.txt](data/FullDataSet_DenovoMutations_for_mTADA.txt):
    all gene-level de novo mutations.

2.  [SingleTrait\_Parameters.txt](data/SingleTrait_Parameters.txt): all
    single-trait parameters. We used `extTADA` to estimate these
    parameters.

    *Note*: Users can re-run all these single-trait analyses by
    following an example here:
    <a href="https://github.com/hoangtn/extTADA" class="uri">https://github.com/hoangtn/extTADA</a>.

II. Requirements
----------------

**`mTADA` is written in `R`**. Other `R` packages are required to run
`mTADA`:

-   `rstan`:
    <a href="https://mc-stan.org/rstan/" class="uri">https://mc-stan.org/rstan/</a>.

-   `locfit`:
    <a href="https://cran.r-project.org/web/packages/locfit/index.html" class="uri">https://cran.r-project.org/web/packages/locfit/index.html</a>.

III. An example: joint analysis of DD and EE DNMs.
--------------------------------------------------

Only one function `mTADA` (in the **Run `mTADA`** section) is used to
obtain results. However, some additional steps are described here.

### Load the source codes

``` r
dataDir <- "./data/"
sourceDir <- "./script/"
rFile <- dir(sourceDir, ".R$")
for (ir in 1:length(rFile)){
    source(paste0(sourceDir, rFile[ir]))
}
```

    ## locfit 1.5-9.1    2013-03-22

### Read the data and single-trait parameters

``` r
## De novo data
data <- read.table(paste0(dataDir, "FullDataSet_DenovoMutations_for_mTADA.txt"), header = TRUE, as.is = TRUE) 
## Single-trait parameters
sPar <- read.table(paste0(dataDir, "SingleTrait_Parameters.txt"), as.is = TRUE, header = TRUE)

trait1 = "DD"
trait2 = "EE"
##Take a quick look at the single-trait parameters of DD and EE
sPar[grep(trait1, sPar[, 1]), ] ##Trait 1
```

    ##                 Parameter EstimatedValue
    ## 8                DD_pi[1]     0.02936283
    ## 9  DD_hyperGammaMeanDN[1]    22.31762802
    ## 10 DD_hyperGammaMeanDN[2]    86.03966530
    ## 11      DD_hyperBetaDN[1]     0.82594514
    ## 12      DD_hyperBetaDN[2]     0.80689775

``` r
sPar[grep(trait2, sPar[, 1]), ] ##Trait 2
```

    ##                 Parameter EstimatedValue
    ## 18               EE_pi[1]     0.01548789
    ## 19 EE_hyperGammaMeanDN[1]    51.08181282
    ## 20 EE_hyperGammaMeanDN[2]    65.15189031
    ## 21      EE_hyperBetaDN[1]     0.80906448
    ## 22      EE_hyperBetaDN[2]     0.80774192

### Set parameters for two traits.

As described above, `mTADA` needs single-trait parameters:

-   the number of trios: *ntrio*;

-   the mean and disperson parameters of relative risks:
    $\\bar{\\gamma\_j}$ and *β*<sub>*j*</sub> (j=1, 2);

-   the proportion of risk genes: *π*<sub>1</sub><sup>*S*</sup> and
    *π*<sub>2</sub><sup>*S*</sup>.

**All these parameters are shown above.**

``` r
### Trait-1 INFORMATION
ntrio1 = 4293 #family numbers                                                                   
p1 = 0.02936283 #The proportion of risk genes, this is p1S                                
meanGamma1 = c(22.31762802,  86.03966530) #Mean Gamma of two categories                     
beta1 = c(0.82594514, 0.80689775) #Beta values inside the distribution RR ~ Gamma(meanRR*beta, beta)     
dataT1 <- data[, paste0(c("dn_damaging_", "dn_lof_"), trait1)] #De novo data              
muDataT1 <- data[, c("mut_damaging", "mut_lof")] #Mutation data of the first trait              
##########################################
### Trait-2 INFORMATION
ntrio2 = 356
p2 = 0.01548789 #This is p2S
meanGamma2 = c(51.08181282, 65.15189031)
beta2 = c(0.80906448, 0.80774192)
dataT2 <- data[, paste0(c("dn_damaging_", "dn_lof_"), trait2)]
muDataT2 <- muDataT1
```

### Run `mTADA`

In this example, we only use a small number of iterations and two MCMC
chains. However, users can change these parameters to obtain more
reliable results.

``` r
nIteration = 2000 #This should be higher to obtain better results.
nChain = 2 #The number of MCMC chains

##########MAIN ANALYSIS
mTADAresults <- mTADA(geneName = data[, 1],
    #######Trait-1 information
                  ntrio1 = ntrio1, # Trio number of Trait 1
                  p1 = p1, #Risk-gene proportion of Trait 1
                  dataDN1 = data.frame(dataT1), #De novo data of Trait 1
                  mutRate1 = data.frame(muDataT1), # Mutation rates of Trait 1
                  hyperGammaMeanDN1 = c(meanGamma1), # Mean relative risks of Trait 1
                  hyperBetaDN01 = beta1, #NULL, #array(c(1, 1)),                                        
    #######Trait-2 information
                  ntrio2 = ntrio2, # Trio number of Trait 2
                  p2 = p2, #Risk-gene proportion of Trait 2
                  dataDN2 = data.frame(dataT2), # De novo data of Trait 2
                  mutRate2 = data.frame(muDataT2), # Mutation rates of Trait 2
                  hyperGammaMeanDN2 = c(meanGamma2), # Mean relative risks of Trait 2
                  hyperBetaDN02 = beta2, #NULL, #array(c(1, 1)),                                    
    ####Other parameters               
                  nIteration = nIteration,
                  useMCMC = TRUE, #If FALSE, it will use the 'Variational Bayes' approach. 
                  nChain = nChain
                      )
```

    ## No information for core numbers (nCore); therefore, nCore = nChain: 2 core(s) is/are used

    ## Loading required package: ggplot2

    ## Loading required package: StanHeaders

    ## rstan (Version 2.18.2, GitRev: 2e1f913d3ca3)

    ## For execution on a local, multicore CPU with excess RAM we recommend calling
    ## options(mc.cores = parallel::detectCores()).
    ## To avoid recompilation of unchanged Stan programs, we recommend calling
    ## rstan_options(auto_write = TRUE)

    ## ===================
    ## Building the model
    ## =================

    ## 
    ## =======Use MCMC===========

    ## recompiling to avoid crashing R session

    ## ====
    ## Only pi, alpha and hyper parameters are estimated in this step
    ## The method does not calculate HPDs for hyper betas, just their medians
    ## ===

### Get results

`mTADA`’s output includes:

1.  `data`: main gene-level results (posterior probabilities for the
    four models as described in the main manuscript: PP0, PP1, PP2 and
    PP3).

2.  `probModel`: a vector of *π*<sub>*j*</sub>, (*j* = 0..3) in Table 1.

3.  `pars`: the estimated value and credible interval of *π*<sub>3</sub>
    (described as p12 in the our code).

4.  `mcmcData`: MCMC sampling results for *π*<sub>3</sub>.

The most important information is from `data`. **Users can use this
information to obtain top prioritized genes for downstream analyses
(e.g., top shared/specific genes, top genes for each trait)**. However,
we will also take a quick look at all these information.

#### Results for downstream analyses (gene-level posteior probabilities of four models)

``` r
fData <- mTADAresults$data ## Full analysis results of the two-trait analysis.
head(fData)
```

    ##   geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE        NO
    ## 1     A1BG              0         0              0         0 0.9785246
    ## 2 A1BG-AS1              0         0              0         0 0.9647711
    ## 3     A1CF              0         0              0         0 0.9894171
    ## 4      A2M              0         0              1         0 0.7720136
    ## 5  A2M-AS1              0         0              0         0 0.9637774
    ## 6    A2ML1              0         0              0         0 0.9920100
    ##           BOTH        FIRST      SECOND
    ## 1 0.0028679901 0.0102413410 0.008366030
    ## 2 0.0060546081 0.0204564780 0.008717788
    ## 3 0.0006461228 0.0026951985 0.007241535
    ## 4 0.0024959676 0.0002611957 0.225229188
    ## 5 0.0062901964 0.0212034136 0.008728950
    ## 6 0.0002087698 0.0009217631 0.006859482

##### Genes with PP3 &gt; 0.8 (Posterior probabilities of Model 3)

``` r
fData[fData$BOTH > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 2348   CACNA1A              5         0              2         0
    ## 3201      CHD2              0         6              0         1
    ## 6254    GABBR2              2         0              2         0
    ## 6265    GABRB3              2         0              2         0
    ## 6610     GNAO1              4         1              2         0
    ## 7165     HECW2              5         1              1         0
    ## 7426    HNRNPU              0         7              0         1
    ## 8283     KCNQ2              9         0              2         0
    ## 8284     KCNQ3              3         0              1         0
    ## 10146      MLL              1        26              1         0
    ## 12480     PHIP              1         2              0         1
    ## 14673    SCN2A              9         4              2         0
    ## 14681    SCN8A              6         0              2         0
    ## 16228   STXBP1              6         5              4         1
    ##                 NO      BOTH        FIRST       SECOND
    ## 2348  3.079347e-04 0.9934025 3.886573e-03 2.402945e-03
    ## 3201  4.562950e-10 0.9369294 6.307060e-02 2.069438e-10
    ## 6254  2.488554e-03 0.9526933 1.677815e-03 4.314029e-02
    ## 6265  9.116467e-04 0.9800866 1.564362e-03 1.743736e-02
    ## 6610  1.582325e-08 0.9984033 1.596416e-03 3.021221e-07
    ## 7165  1.913166e-06 0.8916632 1.083344e-01 4.807440e-07
    ## 7426  8.675178e-13 0.9362774 6.372262e-02 3.891492e-13
    ## 8283  3.207961e-13 0.9982139 1.786145e-03 5.473469e-12
    ## 8284  4.264476e-03 0.9130268 8.124566e-02 1.463106e-03
    ## 10146 1.465322e-48 0.8672381 1.327619e-01 2.922299e-49
    ## 12480 1.451496e-02 0.8909314 9.017537e-02 4.378234e-03
    ## 14673 3.079061e-18 0.9964297 3.570261e-03 2.623565e-17
    ## 14681 4.175364e-06 0.9959428 4.021457e-03 3.156982e-05
    ## 16228 7.565866e-24 1.0000000 9.837105e-09 2.348108e-17

##### Genes with PP1 &gt; 0.8 (Posterior probabilities of Model 1)

``` r
fData[fData$FIRST > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 347       ADNP              1        19              0         0
    ## 681    ANKRD11              0        32              0         0
    ## 1001    ARID1B              0        30              0         0
    ## 1002     ARID2              0         3              0         0
    ## 1153     ASXL1              0         4              0         0
    ## 1155     ASXL3              0        14              0         0
    ## 1317     AUTS2              0         4              0         0
    ## 1450    BCL11A              2         3              0         0
    ## 1630     BRPF1              0         4              0         0
    ## 2355   CACNA1E              2         2              0         0
    ## 2434    CAMTA1              1         2              0         0
    ## 3202      CHD3              3         1              0         0
    ## 3203      CHD4              5         1              0         0
    ## 3206      CHD7              2         2              0         0
    ## 3457      CLTC              2         3              0         0
    ## 3516     CNOT3              2         2              0         0
    ## 3599  COL4A3BP              4         0              0         0
    ## 3773    CREBBP              7         3              0         0
    ## 3876   CSNK2A1              4         0              0         0
    ## 3924      CTCF              5         0              0         0
    ## 3942    CTNNB1              0        11              0         0
    ## 4632    DNMT3A              4         1              0         0
    ## 4832    DYRK1A              4        14              0         0
    ## 4861      EBF3              2         3              0         0
    ## 4948    EFTUD2              3         2              0         0
    ## 4974     EHMT1              2         7              0         0
    ## 5157     EP300              3        12              0         0
    ## 6120     FOXP1              4         8              0         0
    ## 6121     FOXP2              1         2              0         0
    ## 7330    HIVEP2              2         2              0         0
    ## 7333       HK1              3         1              0         0
    ## 8168    KANSL1              0         8              0         0
    ## 8177     KAT6A              0         8              0         0
    ## 8178     KAT6B              0         8              0         0
    ## 8211     KCNB1              2         1              0         0
    ## 8228     KCNH1              4         0              0         0
    ## 8336     KDM5B              0         3              0         0
    ## 9618     LZTR1              2         1              0         0
    ## 9727    MAP4K4              3         2              0         0
    ## 9906    MED13L              5        13              0         0
    ## 9935     MEF2C              4         4              0         0
    ## 10670    MYT1L              2         2              0         0
    ## 10978     NFIX              1         4              0         0
    ## 11282     NSD1              1         7              0         0
    ## 12004    PACS1              8         0              0         0
    ## 12831     POGZ              0         6              0         0
    ## 12994    PPM1D              0         5              0         0
    ## 13062  PPP2R5D             12         0              0         0
    ## 13250  PRPF40A              1         2              0         0
    ## 13538    PUF60              0         3              0         0
    ## 13540     PUM2              1         2              0         0
    ## 13541     PURA              3         7              0         0
    ## 14894    SETD2              1         2              0         0
    ## 14897    SETD5              2        14              0         0
    ## 15074    SIN3A              1         3              0         0
    ## 15133  SLC12A2              2         1              0         0
    ## 15440   SLC6A1              6         2              0         0
    ## 15546  SMARCA2              9         0              0         0
    ## 15752      SON              0         3              0         0
    ## 15985    SRCAP              1         4              0         0
    ## 16337  SYNGAP1              0        13              0         0
    ## 16578    TCF12              1         2              0         0
    ## 16581    TCF20              0         5              0         0
    ## 16587     TCF4              4         9              0         0
    ## 17284    TNPO3              1         2              0         0
    ## 17548   TRIP12              2         3              0         0
    ## 18337    WDR26              1         2              0         0
    ## 18420    WHSC1              0         3              0         0
    ##                 NO       BOTH     FIRST       SECOND
    ## 347   4.081312e-37 0.19979514 0.8002049 3.111074e-39
    ## 681   2.296030e-60 0.13781713 0.8621829 1.120490e-62
    ## 1001  1.703280e-56 0.14546487 0.8545351 8.851994e-59
    ## 1002  2.192670e-03 0.17375369 0.8240395 1.411517e-05
    ## 1153  1.844496e-05 0.17170422 0.8282772 1.167373e-07
    ## 1155  2.521344e-25 0.19925379 0.8007462 1.915449e-27
    ## 1317  1.133388e-05 0.18373465 0.8162539 7.788809e-08
    ## 1450  7.454474e-07 0.18992799 0.8100713 5.335923e-09
    ## 1630  7.112469e-05 0.13193829 0.8679903 3.300679e-07
    ## 2355  4.910791e-02 0.08853795 0.8622002 1.539571e-04
    ## 2434  1.011637e-02 0.16928150 0.8205384 6.371803e-05
    ## 3202  5.972012e-02 0.11025244 0.8297852 2.422535e-04
    ## 3203  5.533953e-04 0.08705964 0.9123854 1.612133e-06
    ## 3206  4.700052e-02 0.10043834 0.8523921 1.690790e-04
    ## 3457  6.336230e-05 0.11705293 0.8828835 2.564701e-07
    ## 3516  1.533745e-04 0.17649399 0.8233516 1.003748e-06
    ## 3599  1.572764e-03 0.18437749 0.8140389 1.087559e-05
    ## 3773  1.301079e-10 0.09289069 0.9071093 4.067646e-13
    ## 3876  8.950877e-04 0.19630097 0.8027973 6.682039e-06
    ## 3924  3.491418e-05 0.19220150 0.8077633 2.536304e-07
    ## 3942  1.540409e-19 0.18126021 0.8187398 1.041165e-21
    ## 4632  3.850359e-05 0.17062349 0.8293378 2.418441e-07
    ## 4832  5.187261e-31 0.18689018 0.8131098 3.640007e-33
    ## 4861  8.610517e-07 0.18873681 0.8112623 6.115773e-09
    ## 4948  2.504174e-05 0.14848530 0.8514895 1.333201e-07
    ## 4974  2.139532e-13 0.14595440 0.8540456 1.116302e-15
    ## 5157  4.075793e-23 0.11063119 0.8893688 1.547871e-25
    ## 6120  1.000552e-18 0.17501722 0.8249828 6.480409e-21
    ## 6121  5.430066e-03 0.17543562 0.8190988 3.550696e-05
    ## 7330  1.638872e-03 0.14033801 0.8580149 8.183749e-06
    ## 7333  2.806476e-03 0.14326972 0.8539094 1.437576e-05
    ## 8168  1.947735e-13 0.18235870 0.8176413 1.326235e-15
    ## 8177  2.073457e-12 0.14973393 0.8502661 1.114775e-14
    ## 8178  1.000973e-12 0.15718066 0.8428193 5.699201e-15
    ## 8211  8.661089e-03 0.17635623 0.8149255 5.722321e-05
    ## 8228  2.402900e-03 0.16945402 0.8281281 1.501125e-05
    ## 8336  3.013121e-02 0.12229065 0.8474454 1.327471e-04
    ## 9618  2.093521e-02 0.16189829 0.8170399 1.266494e-04
    ## 9727  2.183341e-04 0.12614227 0.8736384 9.624482e-07
    ## 9906  1.525871e-28 0.12336520 0.8766348 6.555695e-31
    ## 9935  2.502455e-11 0.18020483 0.8197952 1.679400e-13
    ## 10670 5.476304e-04 0.15515371 0.8442956 3.072430e-06
    ## 10978 2.421798e-07 0.18195301 0.8180467 1.644546e-09
    ## 11282 2.438359e-11 0.13643772 0.8635623 1.176157e-13
    ## 12004 5.416622e-09 0.17368221 0.8263178 3.475871e-11
    ## 12831 1.677008e-09 0.19162401 0.8083760 1.213664e-11
    ## 12994 5.312755e-08 0.19344290 0.8065570 3.890132e-10
    ## 13062 5.255592e-15 0.17377213 0.8262279 3.374650e-17
    ## 13250 8.680420e-03 0.16339845 0.8278688 5.230630e-05
    ## 13538 2.907187e-04 0.19606753 0.8036396 2.165430e-06
    ## 13540 8.183432e-03 0.16058635 0.8311819 4.826973e-05
    ## 13541 2.498648e-16 0.19722763 0.8027724 1.874162e-18
    ## 14894 4.850155e-02 0.14938806 0.8018345 2.758758e-04
    ## 14897 1.794325e-27 0.16840092 0.8315991 1.109324e-29
    ## 15074 2.374885e-04 0.14759870 0.8521626 1.255826e-06
    ## 15133 5.043368e-02 0.14453631 0.8047535 2.765423e-04
    ## 15440 2.196144e-10 0.18559329 0.8144067 1.527948e-12
    ## 15546 1.417526e-08 0.11301504 0.8869849 5.514142e-11
    ## 15752 3.834163e-03 0.16609181 0.8300506 2.342292e-05
    ## 15985 3.770482e-05 0.12626269 0.8736994 1.663553e-07
    ## 16337 1.015931e-22 0.15150695 0.8484931 5.538286e-25
    ## 16578 3.103036e-03 0.17683118 0.8200454 2.042841e-05
    ## 16581 1.277927e-07 0.19048463 0.8095152 9.180530e-10
    ## 16587 1.135971e-20 0.16823416 0.8317658 7.014670e-23
    ## 17284 6.734956e-03 0.17313027 0.8200914 4.340826e-05
    ## 17548 4.282788e-05 0.13816533 0.8617916 2.096284e-07
    ## 18337 1.978841e-03 0.18601283 0.8119945 1.383973e-05
    ## 18420 3.983613e-03 0.14995575 0.8460391 2.155642e-05

##### Genes with PP2 &gt; 0.8 (Posterior probabilities of Model 2)

``` r
fData[fData$SECOND > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 14671    SCN1A              2         0              4         4
    ##                 NO     BOTH        FIRST   SECOND
    ## 14671 2.037003e-12 0.120551 8.524689e-15 0.879449

#### Use mTADA’s results for single-trait analyses.

We can obtain single-trait results by summing PP1 and PP3 (Trait 1) or
PP2 and PP3 (Trait 2).

##### Trait 1

``` r
fData[, 'pTrait1'] <- fData[, 'BOTH'] + fData[, 'FIRST']
fData1 <- fData[fData$pTrait1 > 0.8, ]
head(fData1[, c(1:5, 10)])
```

    ##      geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE   pTrait1
    ## 347      ADNP              1        19              0         0 1.0000000
    ## 447     AHDC1              0         8              0         0 1.0000000
    ## 681   ANKRD11              0        32              0         0 1.0000000
    ## 1000   ARID1A              1         2              0         0 0.9149099
    ## 1001   ARID1B              0        30              0         0 1.0000000
    ## 1002    ARID2              0         3              0         0 0.9977932

##### Trait 2

``` r
fData[, 'pTrait2'] <- fData[, 'BOTH'] + fData[, 'SECOND']
fData2 <- fData[fData$pTrait2 > 0.8, ]
head(fData2[, c(1:5, 11)])
```

    ##      geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE   pTrait2
    ## 2348  CACNA1A              5         0              2         0 0.9958055
    ## 3201     CHD2              0         6              0         1 0.9369294
    ## 6254   GABBR2              2         0              2         0 0.9958336
    ## 6265   GABRB3              2         0              2         0 0.9975240
    ## 6610    GNAO1              4         1              2         0 0.9984036
    ## 7165    HECW2              5         1              1         0 0.8916637

#### Other information

Some additional information can be obtained from mTADA’s results.

``` r
pCI <- mTADAresults$pars ## Genetic parameters
piValue <- mTADAresults$probModel ## Posterior probabilities of genes for four models
mcmcResult <- mTADAresults$mcmcData ##MCMC results
```

##### The proportions of risk genes

*piValue* is a vector of *π* values. In the result below, pNO, pFIRST,
pSECOND, and pBOTH are *π*<sub>0</sub>, *π*<sub>1</sub>, *π*<sub>2</sub>
and *π*<sub>3</sub> respectively in **Table 1**.

``` r
piValue
```

    ##        pNO     pFIRST    pSECOND      pBOTH 
    ## 0.96188882 0.02262329 0.00874835 0.00673954

##### Estimated information of *π*<sub>3</sub>.

Credible-interval information is from *p**C**I*.

``` r
pCI ## Mode: estimated values; CI: credible interval with low (l) and upper (u) values
```

    ##                        Mode          lCI         uCI
    ## p12              0.00673954  0.003683986  0.01016346
    ## gammaMeanDN1[1] 22.73366919 20.056821697 25.50755047

To check the convergent information of *π*<sub>3</sub>, we can visualize
MCMC results.

``` r
## p12 is pi3 in the model
plotParHeatmap1(mcmcResult = mcmcResult, pars = c('p12', 'gammaMeanDN1[1]'))
```

![](README_mTADA_DD_EE_files/figure-markdown_github/unnamed-chunk-14-1.png)

### Citation

**`mTADA`: a framework for identifying risk genes from de novo mutations
in multiple traits.** Hoang T. Nguyen, Amanda Dobbyn, Ruth C. Brown,
Brien P. Riley, Joseph Buxbaum, Dalila Pinto, Shaun M Purcell, Patrick F
Sullivan8, Xin He, Eli A. Stahl.
