# Analysis 200914_europe1

**DATE**: 2020-09-14

**LOCATION**: `/maceci/code/mt-analysis/200914_europe1`

**STATUS**: Finished

Analysis with same specifications than Sarah BDMM-prime analysis `2020-05-18_european_origins` (the one from the paper *Sarah Nadeau et al. doi:10.1101/2020.06.10.20127738*) but adding trajectory logger BDMM-Prime and an extended sequence dataset.


## Preprocessing

Preprocessing steps done by nextstrain pipeline:

- Filter sequences with bad quality: sequences less than 27,000 bases in length or with more than 3,000 N (unknown) bases are omitted from the analysis.
- Exclude sequences from exclude.txt file, this file is maintained by nextstrain with sequences that are from the same patient or cluster of sequences such as the Diamond Princess.
- Filter minimum and maximum date. Exclude incomplete dates.
- Mask 100 from beginning and 50 from end. Mask sites 13402, 24389 and 24390.
- Subsample sequences for each of the five demes, but including all the sequences in Sarah's original alignment.

## Data Alignment

File: [200914_europe1.fasta](data/200914_europe1.fasta)

Extended Sarah's alignment based on data available on GISAID as of 2020-09-11. 
Total number of sequences: 301
- France: 74
- Germany: 45
- Hubei: 31
- Italy: 71
- Other European: 80 


## XML

- BDMM-Prime, same specifications that in Sarah's analysis and 200910_ds_europe0.
- File: [200914_europe1.xml](analyses/200914_europe1.xml)

#### Tip Dates
Use tip dates auto-configure.


#### Site Model

| Gamma Site Model     |           | Modified     |
| -------------------- | :-------: | :----------: |
| **Parameter**        | **Value** | **estimate** |
| Substitution Rate    | 1.0       |              |
| Gamma Category Count | 4         |
| Shape                | 1.0       |  x           |
| Proportion Invariant | 0.0       |              |
| Subst Model          | HKY       | 
| Kappa                | 2.0       |  x           |
| Frequencies          | Empirical | 


#### Clock Model

| Strick Clock  |                      | Modified     |
| ------------- | :------------------: | :----------: |
| **Parameter** | **Value**            | **estimate** |
| Clock Rate    | 0.0008 (fixed) *[1]* |              |

*[1] Automatic set clock rate Mode disabled*

#### Priors

| Tree Prior                |                      |            | Modified     |
| ------------------------- | :------------------: | :---------:| :----------: |
| **Parameter**             | **Value**            | **scalar** | **estimate** |
| Tree                      | BDMMPrime            | 
| Parameterization          | Epi Parameterization |
| Set locations             | Autoconfigure        |
| **R0**                        
| Number of change times    | 0                    |
| Values F G H I O          | 1.0 1.0 1.0 1.0 1.0  |            | x            |
| **Become Uninfectious Rate**                        
| Number of change times    | 0                    |
| Values ALL                | 36.5  		       | x          |              |
| **Sampling proportion**                        
| Number of change times    | 2                    |
| Change times              | 0.12 0.205 *[2]*     |            |              |
| Time as ages              | true                 |
| Values  F G H I O         | E1: 1.0E-5 1.0E-5 0.0 1.0E-5 1.0E-5 <br> E2:  1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 <br> E3: 0.0 0.0 0.0 0.0 0.0 | | x |
| **Rho Sampling**                        
| Number of elements        | 0                    |
| **Removal Prob**                        
| Number of change times    | 0                    |
| Values ALL                | 1.0 			     | x          |              |
| **Migration Rte**                        
| Number of change times    | 0                    |
| Values  F G H I O         | -- 0.1 0.0 0.1 0.1 <br> 0.1 -- 0.0 0.1 0.1 <br> 0.1 0.1 -- 0.1 0.1 <br> 0.1 0.1 0.0 -- 0.1 <br> 0.1 0.1 0.0 0.1 -- <br> *[3]* | | x |
| **R0Among Demes**                        
| Number of change times    | 0                    |
| Values ALL                | 0.0 		           | x          |              |
| **Origin**                | 10.0                 |            | x            |
| **Final Sample Offset**   | ?                    |
| **Frequencies** F G H I O | 0.0 0.0 1.0 0.0 0.0  |
| Condition on Survival     | False?               |
| Use Analytical ST sol     | True?                |
| Rel Tolerance             | 1.0E-7               |
| Abs Tolerance             | 1.0E-100             |
| Parallelize               | True?                |
| Parallelization Factor    | 0.1                  |

*[1] 0.12 (just after Jan 23, when I stop sampling in Hubei) and 0.205 (~Dec 23, just before 1st sample that could have been selected)*

*[2] China assumed to be a source to all other demes and a sink to none*


| R0 Prior           |                                     | Modified     |
| ------------------ | :---------------------------------: | :----------: |
| **Parameter**      | **Value**                           | **estimate** |
| Log Normal         | [1.0, 1.0, 1.0, 1.0, 1.0][0.0, Inf] | x            |
| M                  | 0.8                                 |              |
| S                  | 0.5                                 |              |
| Mean in Real Space | False                               |
| Offset             | 0.0                                 | 
| Lower              | 0.0                                 |              
| Upper              | Inf                                 |    

| gammaShape Prior      |                     | Modified     |
| --------------------- | :-----------------: | :----------: |
| **Parameter**         | **Value**           | **estimate** |
| Exponential           | [1.0][-Inf, Inf]    | x            |
| Mean                  | 0.5                 |              |           
| Offset                | 0.0                 |
| Lower                 | -Inf                |              
| Upper                 | Inf                 |  

| kappa Prior        |                    |              |
| ------------------ | :----------------: | :----------: |
| **Parameter**      | **Value**          | **estimate** |
| Log Normal         | [2.0][0.0, Inf]    | x            |
| M                  | 1.0                |              |
| S                  | 1.25               |              |
| Mean in Real Space | False              |
| Offset             | 0.0                |
| Lower              | 0.0                |              
| Upper              | Inf                | 

| Migration Rate Prior |                    |              |
| -------------------- | :----------------: | :----------: |
| **Parameter**        | **Value**          | **estimate** |
| Log Normal           | [...][0.0, Inf]    | x            |
| M                    | 0.0                |              |
| S                    | 1.0                |              |
| Mean in Real Space   | False              |
| Offset               | 0.0                |
| Lower                | 0.0                |              
| Upper                | Inf                | 

| Origin Prior       |                    |              |
| ------------------ | :----------------: | :----------: |
| **Parameter**      | **Value**          | **estimate** |
| Log Normal         | [10.0][0.0, Inf]   | x            |
| M                  | -1.0               |              |
| S                  | 0.2                |              |
| Mean in Real Space | False              |
| Offset             | 0.0                |
| Lower              | 0.0                |              
| Upper              | Inf                | 

| samplingProportion Prior |                 | Modified     |
| ------------------------ | :-------------: | :----------: |
| **Parameter**            | **Value**       | **estimate** |
| Uniform                  | [...][0.0, 1.0] | x            |
| Lower                    | 0.0             |              
| Upper                    | 1.0             | 
| Offset                   | 0.0             |

The **sampling proportion prior** is manually modified afterwards:

- Upper bounds based on # seqs/confirmed cases according to WHO situation reports 
unless otherwise noted. This is done using function `feast.function.Slice` with a uniform distribution 
0-upperbound and removing the prior set in BEAUti. Lines 262-308 of xml.

| deme          | last_seq   | upper bound            |
| ------------- | ---------- | ---------------------- |
| France        | 08.03.2020 | 66/706 = 0.093         |
| Germany       | 03.03.2020 | 15/157 = 0.10          |
| Hubei         | 18.01.2020 | 10/66 = 0.15 *[4]*     |
| Italy	        | 04.03.2020 | 13/2502 = 0.005        |
| OtherEuropean	| 08.03.2020 | 41 / 714 = 0.057 *[5]* |

*[4] Source: https://www.statista.com/statistics/1103040/cumulative-coronavirus-covid19-cases-number-worldwide-by-day/*

*[5] Source: https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_Europe and linked country pages*


- Enforce that sampling proportion should be the same before and after the Jan 23 
 breakpoint for all demes except Hubei. We have to modify the sampling proportion 
operator scaler for each deme to be equal between Epochs 2 and 3. We do it with 
the operator `feast.operators.BlockScaleOperator`. Lines 349-368 of xml.

- Add trajectory logger (BDDM-Prime)

| Sampled Trajectory |           |              
| ------------------ | :-------: |
| **Parameter**      | **Value** |
| nParticles         | 1000      |
| useTauLeaping      | True      |              
| minLeapCount       | 100       |              
| epsilon            | 0.03      |  
| typeLabel          | type      |  
| logEvery           | 10000     |  


#### MCMC

Four chains run.

| Parameter              | Value    |
| ---------------------- | :------: | 
| Chain Length           | 10000000 |
| Store Every            | -1       |              
| Pre Burnin             | 0        | 
| Num Init Attempts      | 10       | 
| tracelog Every         | 1000     | 
| screenlog Every        | 1000     | 
| treelog Every          | 1000     | 
| typedTreeLog Every     | 1000     |     
| nodeTypedTreeLog Every | 1000     |        


## Analysis

- Analysis run using *BEAST v2.6.3* on Euler with BDMM-Prime.jar with snakemake workflow
```
bsub -W 120:00 snakemake --profile euler  -p -j 100
```
- Snakemake and config files in `workflow/Snakefile_europe12`, `workflow/config_europe1.yaml`, workflow scripts in `workflow/scripts`. Workflow dag `results/200914_europe1.dag.svg`
- Final chain length: 42121000, 10000000/run, other chains failed so chain with seed 4 was resume several times.
- Files raw results in 200914_europe1/results/eEruope: `200914_europe1.4.log, 200914_europe1.4.trees, 200914_europe1.4.typed.trees, 200914_europe1.4.typed.node.trees, 200914_europe1.4.TL.traj, 200914_europe1.4.xml.state`

## Results processing

- Summary log table 4% burnin [open](results/trace-tables/200914_europe1.summary.tsv)
- Summary tree, maximum clade credibility tree mean height 4% burnin [open](results/mcc-tree/200914_europe1.typed.node.tree)
- Trajectory figures 4% burnin 
 ![open](results/traj-figs) 


## Notes

In trajectory plots, the total cases estimates are higher than for the downsampled alignment (dsEurope0) showing that now we are able to infer the dynamics of a bigger outbreak as we expected. In comparison with the case counts, in all the European demes we can observe that our analysis estimates higher number of cases but quite close to the current knowledge. However, the estimates for Hubei are very low in comparison with the case counts that we have for this deme.

We think of several explanations for this underestimation of cases in Hubei:
- Our model assumes a sampling rate after 23 Jan in Hubei of 0 (we do not include sequences from Hubei after the lockdown). This should not affect the estimates, only the uncertainty about them (from 23 Jan we have less information about what happened in Hubei but the estimates will extrapolate the previous month behaviour).
- Our model also assumes a constant migration rate during the entire time of the analysis. This is a strong assumption that we now is not reflecting reality, since after the lockdown the migration from Hubei to other countries was drastically reduces and even before is very likely that there were a reduction of the travels from Hubei due to social behaviour.
- The Hubei sequences included in the analysis do not contain enough information to infer the dynamics of the complete epidemic in Hubei. Other works have tried to tackle this problem by considering as sequences from Hubei the ones from other countries where we have the model history in an attempt of increasing the diversity of Hubei sequences.

The first introduction by source plot show us that for all european demes, the first introduction was most likely from Hubei.

The births vs migration plots illustrates how in the beginning of the epidemic migrations have a principal role that decreases during the epidemics, having already in January for most demes more than half of the events as births. For Hubei we assume there are no migrations into Hubei. With respect to this plot I have my doubts, it is a really nice result that we could expect, but I wonder if we would obtain the same if we were not underestimating the number of cases in Hubei deme. Also out estimation for the migration rate is smaller for Hubei that for the other demes, and this could be partially due to our assumption of contant migration rates.

The source of the migrations is dominated at the beginning by Hubei, that become less important in time when the cases in the other demes started to increase. The destination of migrations from each demes remain more or less constant during the epidemic.




