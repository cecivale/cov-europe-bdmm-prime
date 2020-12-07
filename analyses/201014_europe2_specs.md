# Analysis 201014_europe2

**DATE**: 2020-10-14

**LOCATION**: `/maceci/code/mt-analysis/201014_europe2`

**STATUS**: Finished

Analysis with same specifications than Sarah BDMM-prime analysis `2020-05-18_european_origins` (the one from the paper *Sarah Nadeau et al. doi:10.1101/2020.06.10.20127738*) but adding a deme for Spain, trajectory logger BDMM-Prime and an extended sequence dataset.


## Preprocessing

Preprocessing steps done by nextstrain pipeline:

- Filter sequences with bad quality: sequences less than 27,000 bases in length or with more than 3,000 N (unknown) bases are omitted from the analysis.
- Exclude sequences from exclude.txt file, this file is maintained by nextstrain with sequences that are from the same patient or cluster of sequences such as the Diamond Princess.
- Filter minimum and maximum date. Exclude incomplete dates.
- Mask 100 from beginning and 50 from end. Mask sites 13402, 24389 and 24390.
- Subsample sequences uniformly at random for each of the five demes, **but including all the sequences in Sarah's original alignment.**

## Data Alignment

File: [201014_europe2.fasta](data/201014_europe2.fasta)

Extended Sarah's alignment based on data available on GISAID as of 2020-10-14. 
Total number of sequences: 322
- France: 75
- Germany: 51
- Hubei: 31
- Italy: 68
- Other European: 43
- Spain: 54


## XML

- BDMM-Prime, same specifications that in previous analysis.
- File: [201014_europe2.xml](analyses/201014_europe2.xml)

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
| **Origin**                | 1.0                  |            | x            |
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
unless otherwise noted. This is done using function `feast.function.Slice` with a uniform distribution 0-upperbound and removing the prior set in BEAUti. Lines 121-177 of xml. In this analysis the uuper bound was not modified exactly to the new number of sequences 

| deme          | until      | seqs/cases      | upper bound |
| ------------- | ---------- | --------------- | ----------- |
| France        | 08.03.2020 | 75/716 = 0.1    | 0.093       |
| Germany       | 08.03.2020 | 51/847 = 0.06   | 0.048       |
| Hubei         | 23.01.2020 | 31/623 = 0.05   | 0.06        |
| Italy	        | 08.03.2020 | 68/5883 = 0.011 | 0.011       |
| OtherEuropean	| 08.03.2020 | 43/1958 = 0.02  | 0.012       |
| Spain      	| 08.03.2020 | 54/1092 =0.049  | 0.046       |


*Case data info from ECDC*


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

Three succesful chains run with seeds 1,4,10.

| Parameter              | Value    |
| ---------------------- | :------: | 
| Chain Length           | 20000000 |
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
- Snakemake, config files and rule scripts in `workflow/`. Workflow dag `results/dags/201014_europe2.dag.svg`
- Final chain length: 54e6, 2e7/chain.
- Files raw results in 200914_europe1/results/europe1

## Results processing

- Summary log table 10% burnin [open](results/trace-tables/201014_europe2.summary.tsv)
- Summary tree, maximum clade credibility tree mean height 4% burnin [open](results/mcc-tree/201014_europe2.typed.node.tree)
- Trajectory figures 10% burnin 
 ![](results/traj-figs/201014_europe2_figtraj01.png) 



## Notes


TODO


