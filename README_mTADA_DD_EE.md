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
    ## 1     A1BG              0         0              0         0 0.9783450
    ## 2 A1BG-AS1              0         0              0         0 0.9645882
    ## 3     A1CF              0         0              0         0 0.9892615
    ## 4      A2M              0         0              1         0 0.7684200
    ## 5  A2M-AS1              0         0              0         0 0.9635947
    ## 6    A2ML1              0         0              0         0 0.9918648
    ##           BOTH        FIRST      SECOND
    ## 1 0.0027903541 0.0103240111 0.008540587
    ## 2 0.0058906752 0.0206214817 0.008899630
    ## 3 0.0006286489 0.0027170260 0.007392823
    ## 4 0.0024175419 0.0002621266 0.228900306
    ## 5 0.0061198844 0.0213744407 0.008911025
    ## 6 0.0002031260 0.0009292383 0.007002865

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
    ## 2348  3.163597e-04 0.9931371 4.025879e-03 2.520661e-03
    ## 3201  4.678387e-10 0.9347998 6.520017e-02 2.166461e-10
    ## 6254  2.551571e-03 0.9505500 1.734506e-03 4.516391e-02
    ## 6265  9.359589e-04 0.9791654 1.619343e-03 1.827927e-02
    ## 6610  1.625958e-08 0.9983457 1.653983e-03 3.169889e-07
    ## 7165  1.958372e-06 0.8881876 1.118099e-01 5.024634e-07
    ## 7426  8.894439e-13 0.9341274 6.587265e-02 4.073843e-13
    ## 8283  3.296398e-13 0.9981495 1.850541e-03 5.742768e-12
    ## 8284  4.368672e-03 0.9101829 8.391804e-02 1.530409e-03
    ## 10146 1.498629e-48 0.8630992 1.369008e-01 3.051643e-49
    ## 12480 1.485851e-02 0.8874933 9.307193e-02 4.576215e-03
    ## 14673 3.163741e-18 0.9963013 3.698741e-03 2.752469e-17
    ## 14681 4.290117e-06 0.9957965 4.166099e-03 3.312035e-05
    ## 16228 7.774944e-24 1.0000000 1.019242e-08 2.463796e-17

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
    ## 347   4.076279e-37 0.19418263 0.8058174 3.172652e-39
    ## 681   2.288220e-60 0.13365487 0.8663451 1.140188e-62
    ## 1000  8.442044e-02 0.11236514 0.8028328 3.816269e-04
    ## 1001  1.697941e-56 0.14110943 0.8588906 9.010018e-59
    ## 1002  2.188005e-03 0.16872150 0.8290761 1.438167e-05
    ## 1153  1.840407e-05 0.16671645 0.8332650 1.189307e-07
    ## 1155  2.518187e-25 0.19365281 0.8063472 1.953325e-27
    ## 1317  1.131352e-05 0.17847270 0.8215159 7.938500e-08
    ## 1450  7.442702e-07 0.18452875 0.8154705 5.439654e-09
    ## 1630  7.086820e-05 0.12792731 0.8720015 3.358012e-07
    ## 2355  4.887615e-02 0.08575045 0.8652169 1.564563e-04
    ## 2434  1.009392e-02 0.16436367 0.8254775 6.491503e-05
    ## 3202  5.948859e-02 0.10687167 0.8333933 2.463945e-04
    ## 3203  5.505365e-04 0.08428085 0.9151670 1.637569e-06
    ## 3206  4.679737e-02 0.09731498 0.8557158 1.718923e-04
    ## 3457  6.310091e-05 0.11343533 0.8865013 2.607891e-07
    ## 3516  1.530603e-04 0.17139606 0.8284499 1.022780e-06
    ## 3599  1.569994e-03 0.17910345 0.8193155 1.108500e-05
    ## 3773  1.294616e-10 0.08994372 0.9100563 4.132648e-13
    ## 3876  8.938807e-04 0.19076461 0.8083347 6.813511e-06
    ## 3924  3.486183e-05 0.18675258 0.8132123 2.585817e-07
    ## 3942  1.537509e-19 0.17605382 0.8239462 1.061082e-21
    ## 4632  3.841678e-05 0.16566085 0.8343005 2.463787e-07
    ## 4832  5.178518e-31 0.18155794 0.8184421 3.710373e-33
    ## 4861  8.596561e-07 0.18336377 0.8166354 6.234403e-09
    ## 4903  1.657226e-02 0.17990638 0.8034015 1.198617e-04
    ## 4948  2.496590e-05 0.14405471 0.8559202 1.357145e-07
    ## 4974  2.132863e-13 0.14158674 0.8584133 1.136250e-15
    ## 5157  4.058064e-23 0.10718791 0.8928121 1.573582e-25
    ## 6120  9.984496e-19 0.16995293 0.8300471 6.602933e-21
    ## 6121  5.418976e-03 0.17036927 0.8241756 3.618042e-05
    ## 6351  4.185276e-12 0.19893151 0.8010685 3.356934e-14
    ## 6606  1.518471e-05 0.19783519 0.8021495 1.209594e-07
    ## 7330  1.633463e-03 0.13611345 0.8622448 8.328459e-06
    ## 7333  2.797528e-03 0.13897251 0.8582153 1.463160e-05
    ## 8168  1.944143e-13 0.17712759 0.8228724 1.351658e-15
    ## 8177  2.067267e-12 0.14527241 0.8547276 1.134846e-14
    ## 8178  9.982456e-13 0.15253706 0.8474629 5.803324e-15
    ## 8211  8.643906e-03 0.17127333 0.8200244 5.831198e-05
    ## 8228  2.397430e-03 0.16452182 0.8330655 1.529238e-05
    ## 8336  3.001978e-02 0.11856202 0.8512832 1.350404e-04
    ## 9618  2.088517e-02 0.15716815 0.8218177 1.290066e-04
    ## 9727  2.175029e-04 0.12228281 0.8774987 9.789686e-07
    ## 9821  3.826430e-04 0.19808333 0.8015310 3.054262e-06
    ## 9906  1.519912e-28 0.11957888 0.8804211 6.667564e-31
    ## 9935  2.497651e-11 0.17502228 0.8249777 1.711464e-13
    ## 10670 5.461016e-04 0.15055996 0.8488908 3.128354e-06
    ## 10978 2.417298e-07 0.17673102 0.8232687 1.676047e-09
    ## 11282 2.429947e-11 0.13231072 0.8676893 1.196775e-13
    ## 12004 5.404988e-09 0.16864865 0.8313513 3.541423e-11
    ## 12831 1.674459e-09 0.18618763 0.8138124 1.237331e-11
    ## 12994 5.305021e-08 0.18796693 0.8120330 3.966246e-10
    ## 13062 5.244320e-15 0.16873650 0.8312635 3.438304e-17
    ## 13250 8.659265e-03 0.15861696 0.8326705 5.327731e-05
    ## 13538 2.903229e-04 0.19053526 0.8091722 2.208007e-06
    ## 13540 8.162650e-03 0.15587115 0.8359170 4.916076e-05
    ## 13541 2.495342e-16 0.19166997 0.8083300 1.911086e-18
    ## 13599 3.047064e-04 0.19819888 0.8014940 2.433702e-06
    ## 14605 1.812114e-19 0.19490592 0.8050941 1.416931e-21
    ## 14894 4.837531e-02 0.14499247 0.8063513 2.809505e-04
    ## 14897 1.790139e-27 0.16349013 0.8365099 1.130036e-29
    ## 15010 5.825625e-04 0.19845091 0.8009619 4.661961e-06
    ## 15074 2.367623e-04 0.14319037 0.8565716 1.278343e-06
    ## 15133 5.029466e-02 0.14026185 0.8091619 2.815859e-04
    ## 15440 2.192343e-10 0.18028985 0.8197102 1.557415e-12
    ## 15546 1.411478e-08 0.10950670 0.8904933 5.606206e-11
    ## 15752 3.825028e-03 0.16124036 0.8349108 2.385905e-05
    ## 15985 3.756137e-05 0.12239989 0.8775624 1.692113e-07
    ## 16337 1.012961e-22 0.14700173 0.8529983 5.638348e-25
    ## 16537 5.959698e-09 0.19666701 0.8033330 4.712431e-11
    ## 16578 3.096791e-03 0.17172967 0.8251527 2.081650e-05
    ## 16581 1.275934e-07 0.18507318 0.8149267 9.359183e-10
    ## 16587 1.133315e-20 0.16332729 0.8366727 7.145596e-23
    ## 17284 6.720729e-03 0.16811869 0.8251164 4.422848e-05
    ## 17548 4.268273e-05 0.13399424 0.8659629 2.133162e-07
    ## 18337 1.975477e-03 0.18070298 0.8173074 1.410707e-05
    ## 18420 3.971880e-03 0.14549348 0.8505127 2.194541e-05

##### Genes with PP2 &gt; 0.8 (Posterior probabilities of Model 2)

``` r
fData[fData$SECOND > 0.8, ]
```

    ##       geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE
    ## 14671    SCN1A              2         0              4         4
    ##                NO      BOTH        FIRST    SECOND
    ## 14671 2.00636e-12 0.1155445 8.465784e-15 0.8844555

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
    ## 1000   ARID1A              1         2              0         0 0.9151979
    ## 1001   ARID1B              0        30              0         0 1.0000000
    ## 1002    ARID2              0         3              0         0 0.9977976

##### Trait 2

``` r
fData[, 'pTrait2'] <- fData[, 'BOTH'] + fData[, 'SECOND']
fData2 <- fData[fData$pTrait2 > 0.8, ]
head(fData2[, c(1:5, 11)])
```

    ##      geneName dn_damaging_DD dn_lof_DD dn_damaging_EE dn_lof_EE   pTrait2
    ## 2348  CACNA1A              5         0              2         0 0.9956578
    ## 3201     CHD2              0         6              0         1 0.9347998
    ## 6254   GABBR2              2         0              2         0 0.9957139
    ## 6265   GABRB3              2         0              2         0 0.9974447
    ## 6610    GNAO1              4         1              2         0 0.9983460
    ## 7165    HECW2              5         1              1         0 0.8881881

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
    ## 0.961706341 0.022805769 0.008930829 0.006557061

##### Estimated information of *π*<sub>3</sub>.

Credible-interval information is from *p**C**I*.

``` r
pCI ## Mode: estimated values; CI: credible interval with low (l) and upper (u) values
```

    ##                         Mode          lCI         uCI
    ## p12              0.006557061  0.003891385  0.01000175
    ## gammaMeanDN1[1] 22.649165084 20.088115025 25.55547357

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
