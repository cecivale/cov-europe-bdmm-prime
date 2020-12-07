# Analysis 201030_europe3

**DATE**: 2020-10-30

**LOCATION**: `/maceci/code/mt-analysis/201030_europe3`

**STATUS**: Finished

Analysis with same specifications than Sarah BDMM-prime analysis `2020-05-18_european_origins` (the one from the paper *Sarah Nadeau et al. doi:10.1101/2020.06.10.20127738*) but adding a deme for Spain, **using China instead of only Hubei sequences**,  trajectory logger BDMM-Prime and an extended sequence dataset **subsampled according to the number of covid cases or deaths in the deme (europeC).**

## Preprocessing

Preprocessing steps done by nextstrain pipeline:

- Filter sequences with bad quality: sequences less than 27,000 bases in length or with more than 3,000 N (unknown) bases are omitted from the analysis.
- Exclude sequences from exclude.txt file, this file is maintained by nextstrain with sequences that are from the same patient or cluster of sequences such as the Diamond Princess.
- Filter minimum and maximum date. Exclude incomplete dates.
- Mask 100 from beginning and 50 from end. Mask sites 13402, 24389 and 24390.
- Subsampling scheme ssprop_China (340 seqs and new include.txt commit e1abd47f):

## Data Alignment

Initially, we tried three different alignments:
- Europe 3: same time intervals than original analysis
- Europe 3b: Hubei samples till end of analysis (2020-03-08)
- Europe 3c: same time intervals but China samples instead of Hubei, build europeC used also in europeGLM analysis

And select 3c since we obtain a better representation of the epidemic in China (more diversity in the sequences, higher R0 and infered number of cases).

File: [201014_europe2.fasta](data/201030_europe3.fasta)

New alignment based on data available on GISAID as of 2020-11-13. 
Total number of sequences: 340
 ssprop_China:
    France:
      country: "France"
      region: "Europe"
      min_date: "2020-01-23"
      max_date: "2020-03-08"
      seq_per_deme: 60
      prob: "cases"
    Germany:
      country: "Germany"
      region: "Europe"
      min_date: "2020-01-28"
      max_date: "2020-03-08"
      seq_per_deme: 50
      prob: "cases"
    Hubei-China:
      country: "China"
      region: "Asia"
      min_date: "2019-12-24"
      max_date: "2020-01-23"
      seq_per_deme: 60
      prob: "deaths"
    Italy:
      country: "Italy"
      region: "Europe"
      min_date: "2020-01-29"
      max_date: "2020-03-08"
      seq_per_deme: 60
      prob: "cases"
    OtherEuropean:
      region: "Europe"
      exclude_country: ["France","Germany","Italy","Spain"]
      min_date: "2020-01-29"
      max_date: "2020-03-08"
      seq_per_deme: 50
      prob: "cases"
    Spain:
      country: "Spain"
      region: "Europe"
      min_date: "2020-02-24"
      max_date: "2020-03-08"
      seq_per_deme: 60
      prob: "cases"


## XML

- BDMM-Prime, same specifications that in previous analysis.
- File: [201030_europe3.xml](analyses/201030_europe3.xml)
- We tried also [201113_europe3v2.xml](analyses/201113_europe3v2.xml) with uninformative priors for samplin proportion (uniform 0-1) but results in imposible sampling proportions.
- I used Hubei-China as deme name so the alphabetical order of the parameters does not change (todo: change to China)
- Stochastic mapping done after mcmc with [201030_europe3.trajectoryMapper.xml](analyses/201030_europe3.trajectoryMapper.xml)


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
| Values F G H I O  S       | 1.0 1.0 1.0 1.0 1.0 1.0 |            | x            |
| **Become Uninfectious Rate**                        
| Number of change times    | 0                    |
| Values ALL                | 36.5  		       | x          |              |
| **Sampling proportion**                        
| Number of change times    | 2                    |
| Change times              | 0.12 0.205 *[2]*     |            |              |
| Time as ages              | true                 |
| Values  F G H I O S       | E1: 1.0E-5 1.1E-5 0.0 1.3E-5 1.4E-5 1.5E-5 <br> E2:  1.0E-5 1.1E-5 1.2E-5 1.3E-5 1.4E-5 1.5E-5 <br> E3: 0.0 0.0 0.0 0.0 0.0 0.0 | | x |
| **Rho Sampling**                        
| Number of elements        | 0                    |
| **Removal Prob**                        
| Number of change times    | 0                    |
| Values ALL                | 1.0 			     | x          |              |
| **Migration Rte**                        
| Number of change times    | 0                    |
| Values  F G H I O S       | -- 0.1 0.0 0.1 0.1 0.1 <br> 0.1 -- 0.0 0.1 0.1 0.1 <br> 0.1 0.1 -- 0.1 0.1 0.1 <br> 0.1 0.1 0.0 -- 0.1 0.1 <br> 0.1 0.1 0.0 0.1 --  0.1 <br> 0.1 0.1 0.0 0.1 0.1 -- <br> *[3]* | | x |
| **R0Among Demes**                        
| Number of change times    | 0                    |
| Values ALL                | 0.0 		           | x          |              |
| **Origin**                | 1.0                  |            | x            |
| **Final Sample Offset**   | ?                    |
| **Frequencies** F G H I O S| 0.0 0.0 1.0 0.0 0.0 0.0 |
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
| Log Normal         | [1.0, 1.0, 1.0, 1.0, 1.0, 1.0][0.0, Inf] | x            |
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
| Log Normal         | [1.0][0.0, Inf]   | x            |
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

| deme          | until      | seqs/cases      | 
| ------------- | ---------- | --------------- | 
| France        | 08.03.2020 | 60/716 = 0.1    | 
| Germany       | 08.03.2020 | 50/847 = 0.06   | 
| Hubei-China   | 23.01.2020 | 60/623 = 0.097  |
| Italy	        | 08.03.2020 | 60/5883 = 0.011 |
| OtherEuropean	| 08.03.2020 | 50/1958 = 0.026 | 
| Spain      	| 08.03.2020 | 60/1092 = 0.055 | 


*Case data info from ECDC*

- Enforce that sampling proportion should be the same before and after the Jan 23 
 breakpoint for all demes except Hubei. Done with SmartScaleOperator and initial values.

- Stochastic mapping done after mcmc with [201030_europe3.trajectoryMapper.xml](analyses/201030_europe3.trajectoryMapper.xml)


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

Three succesful chains run with seeds 1,2,6.

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
- Snakemake, config files and rule scripts in `workflow/`. Workflow dag `results/dags/201030_europe3.dag.svg`
- Final chain length: 27030000 states
- Files raw results in 201030_europe3/results/europe3c

## Results processing

- Summary log table 10% burnin [open](results/trace-tables/201030_europe3.summary.tsv)
- Summary tree, maximum clade credibility tree mean height 4% burnin [open](results/mcc-tree/201030_europe3.typed.node.tree)
- Trajectory figures 10% burnin 
<img src="results/traj-figs/201030_europe3_figtraj01.png"  width="500">


## Notes

Higher number of inferred cases for all demes but bigger uncertainty. Estimation for China closer to ECDC reportetd cases.


