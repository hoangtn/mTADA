---
output:
  html_document: default
  pdf_document: default
---

  `mTADA`
---------------

-   [I. Introduction](#i.-introduction)
    -   [Data for reproducible
        analyses](#data-for-reproducible-analyses)
-   [II. Requirements](#ii.-requirements)
-   [III. An example: joint analysis of DD and EE
    DNMs](#iii.-an-example-joint-analysis-of-dd-and-ee-dnms)
    -   [Load the source codes](#load-the-source-codes)
    -   [Read the data and single-trait
        parameters](#read-the-data-and-single-trait-parameters)
    -   [Set parameters for two traits](#set-parameters-for-two-traits)
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

**Table 1. `mTADA` model for one variant category at the**
*i*<sup>*th*</sup> **gene**.

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
<td><span class="math inline"><em>x</em><sub><em>i</em>1</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>1</sub><em>γ</em><sub>1</sub><em>μ</em><sub><em>i</em></sub>)</span>; <span class="math inline"><em>γ</em><sub>1</sub></span> ~ <span class="math inline"><em>G</em><em>a</em><em>m</em><em>m</em><em>a</em>(<em>γ̄</em><sub>1</sub><em>β</em><sub>1</sub>, <em>β</em><sub>1</sub>)</span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>2</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>2</sub><em>μ</em><sub><em>i</em></sub>)</span></td>
</tr>
<tr class="odd">
<td><span class="math inline"><em>H</em><sub>2</sub></span></td>
<td><span class="math inline"><em>π</em><sub>2</sub></span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>1</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>1</sub><em>μ</em><sub><em>i</em></sub>)</span></td>
<td><span class="math inline"><em>x</em><sub><em>i</em>2</sub> ∼ <em>P</em><em>o</em><em>i</em><em>s</em><em>s</em><em>o</em><em>n</em>(2<em>N</em><sub>2</sub><em>γ</em><sub>2</sub><em>μ</em><sub><em>i</em></sub>)</span>; <span class="math inline"><em>γ</em><sub>2</sub></span> ~ <span class="math inline"><em>G</em><em>a</em><em>m</em><em>m</em><em>a</em>(<em>γ̄</em><sub>2</sub><em>β</em><sub>2</sub>, <em>β</em><sub>2</sub>)</span></td>
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

III. An example: joint analysis of DD and EE DNMs
-------------------------------------------------

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

### Set parameters for two traits

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

#### Results for downstream analyses (gene-level posteior probabilities (PPs) of four models)

We demonstrate how to choose top proritized genes from `mTADA`’s results
using a PP threshold of 0.8. These genes can be shared genes, specific
genes; or genes for single traits.

``` r
fData <- mTADAresults$data ## Full analysis results of the two-trait analysis.
head(fData)
```

    ##   geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE        NO
    ## 1     A1BG              0         0              0         0 0.9783601
    ## 2 A1BG-AS1              0         0              0         0 0.9646035
    ## 3     A1CF              0         0              0         0 0.9892745
    ## 4      A2M              0         0              1         0 0.7687196
    ## 5  A2M-AS1              0         0              0         0 0.9636099
    ## 6    A2ML1              0         0              0         0 0.9918769
    ##           BOTH        FIRST      SECOND
    ## 1 0.0027968521 0.0103170918 0.008525977
    ## 2 0.0059043959 0.0206076714 0.008884411
    ## 3 0.0006301115 0.0027151990 0.007380160
    ## 4 0.0024240791 0.0002620490 0.228594301
    ## 5 0.0061341390 0.0213601262 0.008895786
    ## 6 0.0002035984 0.0009286126 0.006990863

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
    ## 2348  3.156368e-04 0.9931599 4.013926e-03 2.510560e-03
    ## 3201  4.668500e-10 0.9349822 6.501777e-02 2.158151e-10
    ## 6254  2.546173e-03 0.9507336 1.729651e-03 4.499059e-02
    ## 6265  9.338741e-04 0.9792444 1.614628e-03 1.820708e-02
    ## 6610  1.622213e-08 0.9983506 1.649042e-03 3.157130e-07
    ## 7165  1.954506e-06 0.8884849 1.115127e-01 5.006059e-07
    ## 7426  8.875660e-13 0.9343115 6.568851e-02 4.058226e-13
    ## 8283  3.288808e-13 0.9981550 1.845015e-03 5.719657e-12
    ## 8284  4.359755e-03 0.9104263 8.368933e-02 1.524650e-03
    ## 10146 1.495782e-48 0.8634529 1.365471e-01 3.040590e-49
    ## 12480 1.482913e-02 0.8877874 9.282421e-02 4.559283e-03
    ## 14673 3.156474e-18 0.9963123 3.687715e-03 2.741407e-17
    ## 14681 4.280269e-06 0.9958090 4.153687e-03 3.298729e-05
    ## 16228 7.756999e-24 1.0000000 1.016192e-08 2.453867e-17

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
    ## 347   4.076700e-37 0.19465193 0.8053481 3.167503e-39
    ## 681   2.288871e-60 0.13400221 0.8659978 1.138544e-62
    ## 1000  8.444503e-02 0.11265787 0.8025160 3.810791e-04
    ## 1001  1.698387e-56 0.14147299 0.8585270 8.996828e-59
    ## 1002  2.188395e-03 0.16914194 0.8286553 1.435940e-05
    ## 1153  1.840748e-05 0.16713313 0.8328483 1.187474e-07
    ## 1155  2.518451e-25 0.19412114 0.8058789 1.950158e-27
    ## 1317  1.131522e-05 0.17891246 0.8210761 7.925989e-08
    ## 1450  7.443686e-07 0.18498008 0.8150192 5.430983e-09
    ## 1630  7.088960e-05 0.12826197 0.8716668 3.353229e-07
    ## 2355  4.889546e-02 0.08598278 0.8649655 1.562480e-04
    ## 2434  1.009580e-02 0.16477452 0.8250649 6.481503e-05
    ## 3202  5.950790e-02 0.10715368 0.8330924 2.460491e-04
    ## 3203  5.507747e-04 0.08451237 0.9149352 1.635449e-06
    ## 3206  4.681431e-02 0.09757541 0.8554386 1.716577e-04
    ## 3457  6.312271e-05 0.11373702 0.8861996 2.604289e-07
    ## 3516  1.530866e-04 0.17182201 0.8280239 1.021189e-06
    ## 3599  1.570226e-03 0.17954423 0.8188745 1.106750e-05
    ## 3773  1.295155e-10 0.09018929 0.9098107 4.127231e-13
    ## 3876  8.939816e-04 0.19122750 0.8078717 6.802519e-06
    ## 3924  3.486621e-05 0.18720809 0.8127568 2.581678e-07
    ## 3942  1.537751e-19 0.17648890 0.8235111 1.059418e-21
    ## 4632  3.842403e-05 0.16607543 0.8338859 2.459999e-07
    ## 4832  5.179249e-31 0.18200362 0.8179964 3.704492e-33
    ## 4861  8.597727e-07 0.18381289 0.8161862 6.224487e-09
    ## 4903  1.657449e-02 0.18034688 0.8029590 1.196710e-04
    ## 4948  2.497223e-05 0.14442457 0.8555503 1.355146e-07
    ## 4974  2.133419e-13 0.14195132 0.8580487 1.134585e-15
    ## 5157  4.059543e-23 0.10747500 0.8925250 1.571439e-25
    ## 6120  9.986252e-19 0.17037605 0.8296239 6.592696e-21
    ## 6121  5.419903e-03 0.17079259 0.8237514 3.612415e-05
    ## 6351  4.185648e-12 0.19940945 0.8005905 3.351439e-14
    ## 6606  1.518612e-05 0.19831115 0.8016735 1.207618e-07
    ## 7330  1.633915e-03 0.13646603 0.8618917 8.316382e-06
    ## 7333  2.798275e-03 0.13933118 0.8578559 1.461025e-05
    ## 8168  1.944443e-13 0.17756475 0.8224352 1.349533e-15
    ## 8177  2.067784e-12 0.14564487 0.8543551 1.133170e-14
    ## 8178  9.984733e-13 0.15292481 0.8470752 5.794630e-15
    ## 8211  8.645342e-03 0.17169806 0.8195984 5.822101e-05
    ## 8228  2.397887e-03 0.16493385 0.8326530 1.526890e-05
    ## 8336  3.002908e-02 0.11887309 0.8509630 1.348491e-04
    ## 9618  2.088935e-02 0.15756325 0.8214186 1.288097e-04
    ## 9727  2.175722e-04 0.12260477 0.8771767 9.775905e-07
    ## 9821  3.826780e-04 0.19855968 0.8010546 3.049268e-06
    ## 9906  1.520409e-28 0.11989470 0.8801053 6.658233e-31
    ## 9935  2.498053e-11 0.17545536 0.8245446 1.708784e-13
    ## 10670 5.462293e-04 0.15094353 0.8485071 3.123685e-06
    ## 10978 2.417674e-07 0.17716742 0.8228323 1.673414e-09
    ## 11282 2.430649e-11 0.13265511 0.8673449 1.195054e-13
    ## 12004 5.405960e-09 0.16906919 0.8309308 3.535946e-11
    ## 12831 1.674672e-09 0.18664209 0.8133579 1.235353e-11
    ## 12994 5.305667e-08 0.18842472 0.8115752 3.959883e-10
    ## 13062 5.245261e-15 0.16915721 0.8308428 3.432986e-17
    ## 13250 8.661032e-03 0.15901633 0.8322694 5.319620e-05
    ## 13538 2.903560e-04 0.19099780 0.8087096 2.204448e-06
    ## 13540 8.164386e-03 0.15626495 0.8355216 4.908634e-05
    ## 13541 2.495618e-16 0.19213465 0.8078654 1.907999e-18
    ## 13599 3.047342e-04 0.19867545 0.8010174 2.429722e-06
    ## 14605 1.812297e-19 0.19537655 0.8046235 1.414628e-21
    ## 14894 4.838585e-02 0.14535955 0.8059741 2.805267e-04
    ## 14897 1.790489e-27 0.16390034 0.8360997 1.128305e-29
    ## 15010 5.826150e-04 0.19892791 0.8004848 4.654333e-06
    ## 15074 2.368229e-04 0.14355836 0.8562035 1.276464e-06
    ## 15133 5.030627e-02 0.14061877 0.8087938 2.811647e-04
    ## 15440 2.192661e-10 0.18073311 0.8192669 1.554952e-12
    ## 15546 1.411983e-08 0.10979924 0.8902007 5.598529e-11
    ## 15752 3.825791e-03 0.16164559 0.8345048 2.382262e-05
    ## 15985 3.757334e-05 0.12272212 0.8772401 1.689730e-07
    ## 16337 1.013209e-22 0.14737786 0.8526221 5.629995e-25
    ## 16537 5.960269e-09 0.19714085 0.8028591 4.704748e-11
    ## 16578 3.097313e-03 0.17215594 0.8247260 2.078407e-05
    ## 16581 1.276101e-07 0.18552553 0.8144743 9.344249e-10
    ## 16587 1.133537e-20 0.16373717 0.8362628 7.134659e-23
    ## 17284 6.721918e-03 0.16853741 0.8246965 4.415995e-05
    ## 17548 4.269484e-05 0.13434232 0.8656148 2.130085e-07
    ## 18337 1.975758e-03 0.18114678 0.8168634 1.408472e-05
    ## 18420 3.972860e-03 0.14586602 0.8501392 2.191293e-05

##### Genes with PP2 &gt; 0.8 (Posterior probabilities of Model 2)

``` r
fData[fData$SECOND > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 14671    SCN1A              2         0              4         4
    ##                NO      BOTH        FIRST    SECOND
    ## 14671 2.00889e-12 0.1159578 8.470647e-15 0.8840422

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
    ## 1000   ARID1A              1         2              0         0 0.9151739
    ## 1001   ARID1B              0        30              0         0 1.0000000
    ## 1002    ARID2              0         3              0         0 0.9977972

##### Trait 2

``` r
fData[, 'pTrait2'] <- fData[, 'BOTH'] + fData[, 'SECOND']
fData2 <- fData[fData$pTrait2 > 0.8, ]
head(fData2[, c(1:5, 11)])
```

    ##      geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE   pTrait2
    ## 2348  CACNA1A              5         0              2         0 0.9956704
    ## 3201     CHD2              0         6              0         1 0.9349822
    ## 6254   GABBR2              2         0              2         0 0.9957242
    ## 6265   GABRB3              2         0              2         0 0.9974515
    ## 6610    GNAO1              4         1              2         0 0.9983509
    ## 7165    HECW2              5         1              1         0 0.8884854

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
    ## 0.961721614 0.022790496 0.008915556 0.006572334

##### Estimated information of *π*<sub>3</sub>.

Credible-interval information is from *p**C**I*.

``` r
pCI ## Mode: estimated values; CI: credible interval with low (l) and upper (u) values
```

    ##                         Mode          lCI         uCI
    ## p12              0.006572334  0.003861267  0.01016222
    ## gammaMeanDN1[1] 22.628526092 20.062524125 25.39720292

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
Sullivan, Xin He, Eli A. Stahl.
