-   [I. Introduction](#i.-introduction)
    -   [Data for reproducible
        analyses](#data-for-reproducible-analyses)
-   [II. Requirements](#ii.-requirements)
-   [III. An example: joint analysis of DD and EE
    DNMs.](#iii.-an-example-joint-analysis-of-dd-and-ee-dnms.)
    -   [Load the source codes](#load-the-source-codes)
    -   [Read the data and single-trait
        parameters](#read-the-data-and-single-trait-parameters)
    -   [Set parameters for two
        traits.](#set-parameters-for-two-traits.)
    -   [Run `mTADA`](#run-textttmtada)
    -   [Get results](#get-results)
-   [Citation](#citation)

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
<td><span class="math inline"><em>x</em><sub><em>i</em>2</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>2</sub><em>γ</em><sub>2</sub><em>μ</em><sub><em>i</em></sub>)</span>; <span class="math inline"><em>γ</em><sub>2</sub></span> ~ Gamma(<span class="math inline">$\bar{\gamma_2}$</span><span class="math inline"><em>β</em><sub>2</sub></span>,<span class="math inline"><em>β</em><sub>2</sub>)</span></td>
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
    all gene-level de novo mutations. These DNMs are used in the main
    manuscript.

2.  [SingleTrait\_Parameters.txt](data/SingleTrait_Parameters.txt): all
    single-trait parameters. We used `extTADA` to estimate these
    parameters from the DNMs above.

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

Software versions were used in our analyses: `R` version 3.5.2, `locfit`
version 1.5-9.1, and `rstan` version 2.18.2.

III. An example: joint analysis of DD and EE DNMs.
--------------------------------------------------

Only one function `mTADA` (in the **Run `mTADA`** section) is used to
obtain results. Therefore, users can go directly to the **Run `mTADA`**
section to run `mTADA`. However, some additional steps are described
here.

### Load the source codes

``` r
dataDir <- "./data/"
source("script/mTADA.R")
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
    ## 1     A1BG              0         0              0         0 0.9786222
    ## 2 A1BG-AS1              0         0              0         0 0.9648705
    ## 3     A1CF              0         0              0         0 0.9895017
    ## 4      A2M              0         0              1         0 0.7739798
    ## 5  A2M-AS1              0         0              0         0 0.9638768
    ## 6    A2ML1              0         0              0         0 0.9920889
    ##           BOTH        FIRST      SECOND
    ## 1 0.0029101744 0.0101964215 0.008271183
    ## 2 0.0061436834 0.0203668209 0.008618982
    ## 3 0.0006556171 0.0026833389 0.007159334
    ## 4 0.0025388770 0.0002606864 0.223220589
    ## 5 0.0063827379 0.0211104836 0.008630018
    ## 6 0.0002118363 0.0009177016 0.006781578

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
    ## 2348  3.035436e-04 0.9935409 3.813967e-03 2.341592e-03
    ## 3201  4.502601e-10 0.9380427 6.195728e-02 2.018715e-10
    ## 6254  2.455611e-03 0.9538138 1.648178e-03 4.208239e-02
    ## 6265  8.989622e-04 0.9805673 1.535676e-03 1.699810e-02
    ## 6610  1.559591e-08 0.9984333 1.566422e-03 2.943760e-07
    ## 7165  1.889474e-06 0.8934847 1.065129e-01 4.693610e-07
    ## 7426  8.560545e-13 0.9374014 6.259855e-02 3.796156e-13
    ## 8283  3.161881e-13 0.9982474 1.752593e-03 5.333154e-12
    ## 8284  4.209933e-03 0.9145154 7.984676e-02 1.427875e-03
    ## 10146 1.447843e-48 0.8694101 1.305899e-01 2.854420e-49
    ## 12480 1.433492e-02 0.8927332 8.865739e-02 4.274480e-03
    ## 14673 3.034935e-18 0.9964967 3.503311e-03 2.556395e-17
    ## 14681 4.115565e-06 0.9960190 3.946084e-03 3.076183e-05
    ## 16228 7.456939e-24 1.0000000 9.651991e-09 2.287837e-17

##### Genes with PP1 &gt; 0.8 (Posterior probabilities of Model 1)

``` r
fData[fData$FIRST > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 681    ANKRD11              0        32              0         0
    ## 1001    ARID1B              0        30              0         0
    ## 1002     ARID2              0         3              0         0
    ## 1153     ASXL1              0         4              0         0
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
    ## 681   2.300295e-60 0.14009002 0.8599100 1.109734e-62
    ## 1001  1.706194e-56 0.14784226 0.8521577 8.765737e-59
    ## 1002  2.195213e-03 0.17649621 0.8212946 1.396994e-05
    ## 1153  1.846725e-05 0.17442291 0.8255585 1.155418e-07
    ## 1317  1.134496e-05 0.18660092 0.8133877 7.707271e-08
    ## 1450  7.460883e-07 0.19286805 0.8071312 5.279439e-09
    ## 1630  7.126480e-05 0.13412927 0.8657991 3.269361e-07
    ## 2355  4.923472e-02 0.09006324 0.8605494 1.525895e-04
    ## 2434  1.012861e-02 0.17196206 0.8178463 6.306557e-05
    ## 3202  5.984664e-02 0.11209995 0.8278134 2.399905e-04
    ## 3203  5.549606e-04 0.08858122 0.9108622 1.598205e-06
    ## 3206  4.711161e-02 0.10214636 0.8505745 1.675405e-04
    ## 3457  6.350520e-05 0.11903062 0.8809056 2.541089e-07
    ## 3516  1.535457e-04 0.17927199 0.8205735 9.933766e-07
    ## 3599  1.574272e-03 0.18725018 0.8111648 1.076153e-05
    ## 3773  1.304617e-10 0.09450387 0.9054961 4.032063e-13
    ## 3924  3.494267e-05 0.19516824 0.8047966 2.509346e-07
    ## 3942  1.541989e-19 0.18409658 0.8159034 1.030314e-21
    ## 4632  3.855092e-05 0.17332863 0.8266326 2.393723e-07
    ## 4832  5.192023e-31 0.18979423 0.8102058 3.601684e-33
    ## 4861  8.618116e-07 0.19166280 0.8083363 6.051171e-09
    ## 4948  2.508314e-05 0.15090332 0.8490715 1.320133e-07
    ## 4974  2.143173e-13 0.14833840 0.8516616 1.105414e-15
    ## 5157  4.085488e-23 0.11251424 0.8874858 1.533809e-25
    ## 6120  1.001698e-18 0.17777712 0.8222229 6.413636e-21
    ## 6121  5.436110e-03 0.17819639 0.8163324 3.513998e-05
    ## 7330  1.641825e-03 0.14264455 0.8557055 8.104739e-06
    ## 7333  2.811361e-03 0.14561552 0.8515589 1.423610e-05
    ## 8168  1.949692e-13 0.18520838 0.8147916 1.312385e-15
    ## 8177  2.076835e-12 0.15216867 0.8478313 1.103822e-14
    ## 8178  1.002461e-12 0.15971374 0.8402863 5.642402e-15
    ## 8211  8.670451e-03 0.17912576 0.8121472 5.662997e-05
    ## 8228  2.405882e-03 0.17214265 0.8254366 1.485800e-05
    ## 8336  3.019208e-02 0.12432768 0.8453487 1.314941e-04
    ## 9618  2.096248e-02 0.16447724 0.8144349 1.253642e-04
    ## 9727  2.187883e-04 0.12825114 0.8715291 9.534212e-07
    ## 9906  1.529128e-28 0.12543443 0.8745656 6.494559e-31
    ## 9935  2.505072e-11 0.18302839 0.8169716 1.661931e-13
    ## 10670 5.484644e-04 0.15765985 0.8417886 3.041920e-06
    ## 10978 2.424250e-07 0.18479779 0.8152020 1.627385e-09
    ## 11282 2.442952e-11 0.13869152 0.8613085 1.164897e-13
    ## 12004 5.422963e-09 0.17642556 0.8235744 3.440144e-11
    ## 12831 1.678395e-09 0.19458402 0.8054160 1.200777e-11
    ## 12994 5.316966e-08 0.19642418 0.8035758 3.848694e-10
    ## 13062 5.261735e-15 0.17651660 0.8234834 3.339958e-17
    ## 13250 8.691954e-03 0.16600559 0.8252507 5.177685e-05
    ## 13538 2.909342e-04 0.19907902 0.8006279 2.142254e-06
    ## 13540 8.194765e-03 0.16315775 0.8285997 4.778381e-05
    ## 14897 1.796607e-27 0.17107812 0.8289219 1.098033e-29
    ## 15074 2.378849e-04 0.15000469 0.8497562 1.243537e-06
    ## 15133 5.050952e-02 0.14686814 0.8023485 2.737910e-04
    ## 15440 2.198215e-10 0.18848186 0.8115181 1.511899e-12
    ## 15546 1.420833e-08 0.11493342 0.8850666 5.463802e-11
    ## 15752 3.839143e-03 0.16873688 0.8274008 2.318514e-05
    ## 15985 3.778320e-05 0.12837338 0.8715887 1.647948e-07
    ## 16337 1.017552e-22 0.15396529 0.8460347 5.483686e-25
    ## 16578 3.106439e-03 0.17961098 0.8172624 2.021694e-05
    ## 16581 1.279012e-07 0.19343125 0.8065686 9.083251e-10
    ## 16587 1.137419e-20 0.17090926 0.8290907 6.943293e-23
    ## 17284 6.742710e-03 0.17586149 0.8173528 4.296126e-05
    ## 17548 4.290714e-05 0.14044301 0.8595139 2.076147e-07
    ## 18337 1.980674e-03 0.18890476 0.8091009 1.369413e-05
    ## 18420 3.990015e-03 0.15239075 0.8435979 2.134416e-05

##### Genes with PP2 &gt; 0.8 (Posterior probabilities of Model 2)

``` r
fData[fData$SECOND > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 14671    SCN1A              2         0              4         4
    ##                 NO      BOTH        FIRST    SECOND
    ## 14671 2.054044e-12 0.1233351 8.557447e-15 0.8766649

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
    ## 1000   ARID1A              1         2              0         0 0.9147526
    ## 1001   ARID1B              0        30              0         0 1.0000000
    ## 1002    ARID2              0         3              0         0 0.9977908

##### Trait 2

``` r
fData[, 'pTrait2'] <- fData[, 'BOTH'] + fData[, 'SECOND']
fData2 <- fData[fData$pTrait2 > 0.8, ]
head(fData2[, c(1:5, 11)])
```

    ##      geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE   pTrait2
    ## 2348  CACNA1A              5         0              2         0 0.9958825
    ## 3201     CHD2              0         6              0         1 0.9380427
    ## 6254   GABBR2              2         0              2         0 0.9958962
    ## 6265   GABRB3              2         0              2         0 0.9975654
    ## 6610    GNAO1              4         1              2         0 0.9984336
    ## 7165    HECW2              5         1              1         0 0.8934852

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

    ##         pNO      pFIRST     pSECOND       pBOTH 
    ## 0.961987972 0.022524138 0.008649198 0.006838692

##### Estimated information of *π*<sub>3</sub>.

Credible-interval information is from *p**C**I*.

``` r
pCI ## Mode: estimated values; CI: credible interval with low (l) and upper (u) values
```

    ##                         Mode          lCI         uCI
    ## p12              0.006838692  0.003715629  0.01008365
    ## gammaMeanDN1[1] 22.649020011 20.064977563 25.57691995

To check the convergent information of *π*<sub>3</sub>, we can visualize
MCMC results.

``` r
## p12 is pi3 in the model
plotParHeatmap1(mcmcResult = mcmcResult, pars = c('p12', 'gammaMeanDN1[1]'))
```

![](README_mTADA_DD_EE_files/figure-markdown_github/unnamed-chunk-14-1.png)

Citation
--------

**`mTADA`: a framework for identifying risk genes from de novo mutations
in multiple traits.** Hoang T. Nguyen, Amanda Dobbyn, Ruth C. Brown,
Brien P. Riley, Joseph Buxbaum, Dalila Pinto, Shaun M Purcell, Patrick F
Sullivan8, Xin He, Eli A. Stahl.
