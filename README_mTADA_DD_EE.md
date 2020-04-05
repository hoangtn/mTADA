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
    ## 1     A1BG              0         0              0         0 0.9783072
    ## 2 A1BG-AS1              0         0              0         0 0.9645497
    ## 3     A1CF              0         0              0         0 0.9892287
    ## 4      A2M              0         0              1         0 0.7676666
    ## 5  A2M-AS1              0         0              0         0 0.9635561
    ## 6    A2ML1              0         0              0         0 0.9918342
    ##           BOTH        FIRST      SECOND
    ## 1 0.0027739886 0.0103414379 0.008577383
    ## 2 0.0058561186 0.0206562639 0.008937962
    ## 3 0.0006249653 0.0027216273 0.007424715
    ## 4 0.0024010994 0.0002623218 0.229669985
    ## 5 0.0060839832 0.0214104926 0.008949406
    ## 6 0.0002019363 0.0009308141 0.007033091

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
    ## 2348  3.181953e-04 0.9930793 4.056230e-03 2.546308e-03
    ## 3201  4.703475e-10 0.9343370 6.566300e-02 2.187548e-10
    ## 6254  2.565267e-03 0.9500842 1.746828e-03 4.560374e-02
    ## 6265  9.412515e-04 0.9789649 1.631312e-03 1.846255e-02
    ## 6610  1.635466e-08 0.9983331 1.666528e-03 3.202287e-07
    ## 7165  1.968177e-06 0.8874338 1.125638e-01 5.071745e-07
    ## 7426  8.942092e-13 0.9336601 6.633992e-02 4.113474e-13
    ## 8283  3.315671e-13 0.9981354 1.864575e-03 5.801455e-12
    ## 8284  4.391295e-03 0.9095654 8.449825e-02 1.545022e-03
    ## 10146 1.505845e-48 0.8622025 1.377975e-01 3.079669e-49
    ## 12480 1.493303e-02 0.8867476 9.370025e-02 4.619161e-03
    ## 14673 3.182193e-18 0.9962733 3.726737e-03 2.780559e-17
    ## 14681 4.315122e-06 0.9957646 4.197617e-03 3.345821e-05
    ## 16228 7.820510e-24 1.0000000 1.026986e-08 2.489009e-17

##### Genes with PP1 &gt; 0.8 (Posterior probabilities of Model 1)

``` r
fData[fData$FIRST > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 347       ADNP              1        19              0         0
    ## 681    ANKRD11              0        32              0         0
    ## 1000    ARID1A              1         2              0         0
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
    ## 4903    EEF1A2              3         0              0         0
    ## 4948    EFTUD2              3         2              0         0
    ## 4974     EHMT1              2         7              0         0
    ## 5157     EP300              3        12              0         0
    ## 6120     FOXP1              4         8              0         0
    ## 6121     FOXP2              1         2              0         0
    ## 6351   GATAD2B              0         7              0         0
    ## 6606     GNAI1              5         0              0         0
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
    ## 9821      MBD5              0         3              0         0
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
    ## 13599   QRICH1              0         3              0         0
    ## 14605    SATB2              4         8              0         0
    ## 14894    SETD2              1         2              0         0
    ## 14897    SETD5              2        14              0         0
    ## 15010   SHANK3              0         3              0         0
    ## 15074    SIN3A              1         3              0         0
    ## 15133  SLC12A2              2         1              0         0
    ## 15440   SLC6A1              6         2              0         0
    ## 15546  SMARCA2              9         0              0         0
    ## 15752      SON              0         3              0         0
    ## 15985    SRCAP              1         4              0         0
    ## 16337  SYNGAP1              0        13              0         0
    ## 16537  TBL1XR1              2         4              0         0
    ## 16578    TCF12              1         2              0         0
    ## 16581    TCF20              0         5              0         0
    ## 16587     TCF4              4         9              0         0
    ## 17284    TNPO3              1         2              0         0
    ## 17548   TRIP12              2         3              0         0
    ## 18337    WDR26              1         2              0         0
    ## 18420    WHSC1              0         3              0         0
    ##                 NO       BOTH     FIRST       SECOND
    ## 347   4.075219e-37 0.19300103 0.8069990 3.185616e-39
    ## 681   2.286580e-60 0.13278089 0.8672191 1.144324e-62
    ## 1000  8.435859e-02 0.11162858 0.8036298 3.830051e-04
    ## 1001  1.696820e-56 0.14019460 0.8598054 9.043211e-59
    ## 1002  2.187024e-03 0.16766322 0.8301353 1.443771e-05
    ## 1153  1.839547e-05 0.16566762 0.8343139 1.193919e-07
    ## 1155  2.517522e-25 0.19247366 0.8075263 1.961299e-27
    ## 1317  1.130924e-05 0.17736565 0.8226230 7.969992e-08
    ## 1450  7.440225e-07 0.18339252 0.8166067 5.461483e-09
    ## 1630  7.081436e-05 0.12708530 0.8728435 3.370048e-07
    ## 2355  4.882756e-02 0.08516607 0.8658494 1.569803e-04
    ## 2434  1.008920e-02 0.16332956 0.8265161 6.516674e-05
    ## 3202  5.944000e-02 0.10616220 0.8341505 2.472635e-04
    ## 3203  5.499375e-04 0.08369862 0.9157498 1.642898e-06
    ## 3206  4.675476e-02 0.09665988 0.8564129 1.724824e-04
    ## 3457  6.304607e-05 0.11267638 0.8872603 2.616951e-07
    ## 3516  1.529943e-04 0.17032385 0.8295221 1.026783e-06
    ## 3599  1.569412e-03 0.17799381 0.8204256 1.112906e-05
    ## 3773  1.293262e-10 0.08932610 0.9106739 4.146271e-13
    ## 3876  8.936267e-04 0.18959920 0.8095003 6.841186e-06
    ## 3924  3.485082e-05 0.18560579 0.8143591 2.596237e-07
    ## 3942  1.536898e-19 0.17495858 0.8250414 1.065272e-21
    ## 4632  3.839853e-05 0.16461736 0.8353440 2.473322e-07
    ## 4832  5.176678e-31 0.18043596 0.8195640 3.725180e-33
    ## 4861  8.593624e-07 0.18223311 0.8177660 6.259367e-09
    ## 4903  1.656663e-02 0.17879743 0.8045156 1.203419e-04
    ## 4948  2.494996e-05 0.14312396 0.8568510 1.362175e-07
    ## 4974  2.131462e-13 0.14066931 0.8593307 1.140440e-15
    ## 5157  4.054346e-23 0.10646573 0.8935343 1.578975e-25
    ## 6120  9.980074e-19 0.16888786 0.8311121 6.628701e-21
    ## 6121  5.416644e-03 0.16930370 0.8252433 3.632206e-05
    ## 6351  4.184338e-12 0.19772811 0.8022719 3.370772e-14
    ## 6606  1.518119e-05 0.19663679 0.8033479 1.214570e-07
    ## 7330  1.632327e-03 0.13522628 0.8631330 8.358849e-06
    ## 7333  2.795649e-03 0.13806996 0.8591197 1.468534e-05
    ## 8168  1.943387e-13 0.17602709 0.8239729 1.357006e-15
    ## 8177  2.065967e-12 0.14433512 0.8556649 1.139062e-14
    ## 8178  9.976723e-13 0.15156121 0.8484388 5.825206e-15
    ## 8211  8.640292e-03 0.17020420 0.8210970 5.854099e-05
    ## 8228  2.396280e-03 0.16348476 0.8341036 1.535149e-05
    ## 8336  2.999639e-02 0.11777938 0.8520887 1.355218e-04
    ## 9618  2.087466e-02 0.15617370 0.8228221 1.295022e-04
    ## 9727  2.173284e-04 0.12147282 0.8783089 9.824358e-07
    ## 9821  3.825551e-04 0.19688392 0.8027305 3.066834e-06
    ## 9906  1.518661e-28 0.11878432 0.8812157 6.691039e-31
    ## 9935  2.496641e-11 0.17393210 0.8260679 1.718209e-13
    ## 10670 5.457804e-04 0.14959467 0.8498564 3.140106e-06
    ## 10978 2.416351e-07 0.17563247 0.8243673 1.682674e-09
    ## 11282 2.428181e-11 0.13144420 0.8685558 1.201104e-13
    ## 12004 5.402541e-09 0.16759011 0.8324099 3.555208e-11
    ## 12831 1.673923e-09 0.18504351 0.8149565 1.242312e-11
    ## 12994 5.303393e-08 0.18681438 0.8131856 3.982266e-10
    ## 13062 5.241949e-15 0.16767752 0.8323225 3.451691e-17
    ## 13250 8.654817e-03 0.15761178 0.8336799 5.348144e-05
    ## 13538 2.902396e-04 0.18937073 0.8103368 2.216970e-06
    ## 13540 8.158282e-03 0.15488003 0.8369123 4.934805e-05
    ## 13541 2.494646e-16 0.19050004 0.8095000 1.918858e-18
    ## 13599 3.046366e-04 0.19699891 0.8026940 2.443722e-06
    ## 14605 1.811653e-19 0.19372097 0.8062790 1.422728e-21
    ## 14894 4.834878e-02 0.14406860 0.8073006 2.820171e-04
    ## 14897 1.789259e-27 0.16245764 0.8375424 1.134390e-29
    ## 15010 5.824304e-04 0.19724989 0.8021630 4.681166e-06
    ## 15074 2.366097e-04 0.14226433 0.8574978 1.283073e-06
    ## 15133 5.026545e-02 0.13936361 0.8100883 2.826457e-04
    ## 15440 2.191543e-10 0.17917399 0.8208260 1.563614e-12
    ## 15546 1.410210e-08 0.10877080 0.8912292 5.625517e-11
    ## 15752 3.823108e-03 0.16022040 0.8359325 2.395074e-05
    ## 15985 3.753127e-05 0.12158918 0.8783731 1.698107e-07
    ## 16337 1.012337e-22 0.14605518 0.8539448 5.659371e-25
    ## 16537 5.958261e-09 0.19547396 0.8045260 4.731776e-11
    ## 16578 3.095478e-03 0.17065666 0.8262270 2.089813e-05
    ## 16581 1.275514e-07 0.18393435 0.8160655 9.396780e-10
    ## 16587 1.132756e-20 0.16229563 0.8377044 7.173123e-23
    ## 17284 6.717737e-03 0.16706473 0.8261731 4.440098e-05
    ## 17548 4.265225e-05 0.13311839 0.8668387 2.140906e-07
    ## 18337 1.974769e-03 0.17958572 0.8184253 1.416332e-05
    ## 18420 3.969415e-03 0.14455599 0.8514526 2.202713e-05

##### Genes with PP2 &gt; 0.8 (Posterior probabilities of Model 2)

``` r
fData[fData$SECOND > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 14671    SCN1A              2         0              4         4
    ##                 NO      BOTH       FIRST    SECOND
    ## 14671 2.000017e-12 0.1145081 8.45359e-15 0.8854919

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
    ## 1000   ARID1A              1         2              0         0 0.9152584
    ## 1001   ARID1B              0        30              0         0 1.0000000
    ## 1002    ARID2              0         3              0         0 0.9977985

##### Trait 2

``` r
fData[, 'pTrait2'] <- fData[, 'BOTH'] + fData[, 'SECOND']
fData2 <- fData[fData$pTrait2 > 0.8, ]
head(fData2[, c(1:5, 11)])
```

    ##      geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE   pTrait2
    ## 2348  CACNA1A              5         0              2         0 0.9956256
    ## 3201     CHD2              0         6              0         1 0.9343370
    ## 6254   GABBR2              2         0              2         0 0.9956879
    ## 6265   GABRB3              2         0              2         0 0.9974274
    ## 6610    GNAO1              4         1              2         0 0.9983335
    ## 7165    HECW2              5         1              1         0 0.8874343

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
    ## 0.961667875 0.022844235 0.008969295 0.006518595

##### Estimated information of *π*<sub>3</sub>.

Credible-interval information is from *p**C**I*.

``` r
pCI ## Mode: estimated values; CI: credible interval with low (l) and upper (u) values
```

    ##                         Mode          lCI          uCI
    ## p12              0.006518595  0.003932507  0.009906634
    ## gammaMeanDN1[1] 22.757501929 19.943026063 25.589192643

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
