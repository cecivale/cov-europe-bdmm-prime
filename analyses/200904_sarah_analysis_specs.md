# Analysis 200904_sarah

**DATE**: 2020-09-04

**LOCATION**: `/maceci/code/mt-analysis/200904_sarah`

**STATUS**: In progress

Rerun Sarah BDMM-prime analysis `2020-05-18_european_origins` (the one from the paper *Sarah Nadeau et al. doi:10.1101/2020.06.10.20127738*) with same `.xml` file.

All the files are found in the `2019-nCov Data` gitlab repository. Analysis notes in `2020-05-18_european_origins/analysis_notes.txt`.


## Data Alignment

File: `2019-ncov-data/analyses/2020-04-19_europe_bdmm_prime/europe_demes.fasta` 
Sarah's alignment based on data available on GISAID as of 2020-04-01. 


## Preprocessing

Preprocessing steps already done by Sarah in `analyses/2020-05-04_bdmm_european_origins/april_19_redo_downsample_based_on_mar_8_deaths/assign_demes_and_downsample_alignment.R`

Assign demes, filter to date of Lombardy lockdown and downsample within demes.
 - Only European samples collected on or beforde Mar 8 (Lombardy lockdown).
 - Only Hubei samples before Jan 23 (Wuhan lockdown)
 - Downsampled OtherEuropean by taking # seqs = # deaths (or 1 seq where no dealths)
 - Downsample valid sequences so that BDMM runs in reasonable time


## XML

- Created using BEAUti and modify afterwards
- Template: BEAUti BDMM-prime?
- File: `europe_demes.xml`, Remake `200904_sarah_europe_demes.xml`, not used.

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

| Strick Clock  |                | Modified     |
| ------------- | :-------:      | :----------: |
| **Parameter** | **Value**      | **estimate** |
| Clock Rate    | 0.0008 (fixed) |              |


#### Priors

| Tree Prior              |                      |            | Modified     |
| ----------------------- | :------------------: | :---------:| :----------: |
| **Parameter**           | **Value**            | **scalar** | **estimate** |
| Tree                    | BDMMPrime            | 
| Parameterization        | Epi Parameterization |
| Set locations           | Autoconfigure        |
| **R0**                        
| Number of change times  | 0                    |
| Values F G H I O        | 1.0 1.0 1.0 1.0 1.0  |            | x            |
| **Become Uninfectious Rate**                        
| Number of change times  | 0                    |
| Values ALL              | 1.0 				    | x          |              |
| **Sampling proportion**                        
| Number of change times  | 2                    |
| Change times            | 0.12 0.205 *[1]*     |            |              |
| Time as ages            | true                 |
| Values  F G H I O       | E1: 1.0E-5 1.0E-5 0.0 1.0E-5 1.0E-5 <br> E2:  1.0E-5 1.0E-5 1.0E-5 1.0E-5 1.0E-5 <br> E3: 0.0 0.0 0.0 0.0 0.0 | | x |
| **Rho Sampling**                        
| Number of elements      | 0                    |
| **Removal Prob**                        
| Number of change times  | 0                    |
| Values ALL              | 1.0 				    | x          |              |
| **Migration Rte**                        
| Number of change times  | 0                    |
| Values  F G H I O       | -- 0.1 0.0 0.1 0.1 <br> 0.1 -- 0.0 0.1 0.1 <br> 0.1 0.1 -- 0.1 0.1 <br> 0.1 0.1 0.0 -- 0.1 <br> 0.1 0.1 0.0 0.1 -- <br> *[2]* | | x |
| **R0Among Demes**                        
| Number of change times  | 0                    |
| Values ALL              | 1.0 				    | x          |              |
| **Origin**              | 10.0                 |            | x            |
| **Final Sample Offset** | ?                    |
| **Frequencies**         | 0.0 0.0 1.0 0.0 0.0  |
| Condition on Survival   | False?               |
| Use Analytical ST sol   | True?                |
| Rel Tolerance           | 1.0E-7               |
| Abs Tolerance           | 1.0E-100             |
| Parallelize             | True?                |
| Parallelization Factor  | 0.1                  |

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


#### MCMC

| Parameter              | Value    |
| ---------------------- | :------: | 
| Chain Length           | 10000000 |
| Store Every            | -1       |              
| Pre Burnin             | 0        | 
| Num Init Attempts      | 10       | 
| tracelog Every         | 1000     | 
| screenlog Every        | 1000     | 
| treelog Every          | 1000    | 
| typedTreeLog Every     | 1000    |     
| nodeTypedTreeLog Every | 1000    |        


## Analysis

- Analysis run using *BEAST v2.6.3* in local BDMM-Prime IntelliJ project.
- Time: - seconds
- Files raw results: `europe_demes.log, europe_demes.trees, europe_demes.europe_demes.typed.trees, europe_demes.europe_demes.typed.node.trees, europe_demes.xml.state`

TODO