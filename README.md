# A comprehensive study of the phylodynamics of SARS-CoV-2 in Europe
#### Cecilia Valenzuela Agüí
##### Supervised by Dr. Timothy Vaughan, Prof. Dr. Tanja Stadler, Sarah Nadeau. Computational Biology Group 
##### CBB Master’s Thesis Project Outline September 1, 2020 - March 16, 2021


This project will focus on the spatial dynamics of the early spread of SARS-CoV-2 in Europe. We will apply a novel approach based on the Multi-type Birth Death phylodynamic model to infer structured population dynamics jointly with between-subpopulation transmission rates from viral genome sequences. The inferred epidemic trajectories for the combined outbreak responsible for the observed sequence data will allow us to better understand the entry into and early spread of SARS-CoV-2 in Europe.

For the analysis, we will use the software package BEAST2, a tool for Bayesian evolutionary analysis of molecular sequences using MCMC. The sequence data will be gathered from GISAID publicly available sequences.

###### i. Main Tasks

1. Create a current and representative sequence alignment for the early phase of COVID-19 pandemic from available SARS-CoV sequences. This includes sequence selection, curation and alignment.
2. Use travel and geographical data as the basis for informative priors on migration rates, this could significantly improve the precision and accuracy of the results.
3. Reconstruction of a partial transmission tree of the early pandemic, including inferred geographic location of ancestral lineages.
4. Inference of multi-type epidemic trajectories, allowing us to directly infer the number SARS-CoV-2 of introductions into a country.

###### ii. Tools

* Bayesian Analysis - `BEAST2 v2.6.3` and related software
* IDE - `IntelliJ IDEA 2020.2.1 (Community Edition)`
* Pre and post analysis - `R`?
* Reference management - `endNote X9.3.3`
* Version Control - `Git and Gitlab SARS-CoV-2 EU phylodynamics repo`
* Mackbook Pro 2.2 GHz 6-Core Intel Core i7
* ETH Euler cluster

###### iii. Workflow 

* **Analysis**

    - XML definition with BEAUti inside IntelliJ 
    - Run analysis using BEAST2 on Euler with `.jar` file from IntelliJ project
    - Analyze results with Tracer, R...?


* **Documentation**

    - _Gitlab repository_: Analysis specs .md file for each analysis in analysis folder. README files for each folder with file information.

    - _Wiki Gitlab Repository_: Research Plan, Sources, Analysis Specs Template, useful guidelines and theoretical topics description.


* **Reference management** with endNote X9, annotating and research notes in the pdf and in the reference itself. endNote compressed library backup in gitlab repo.

* **SARS-CoV-2 EU phylodynamics Gitlab repository**

| Folder          | Description |
| --------------- | ----------- |
| */alignments*   | Sequence data and alignments used in the analysis |
| */analysis*     | XML files from the analysis and Analysis specs .md file with all the information about the analysis steps|
| */results*      | Processed results from the analysis |
| */workflow*     | Snakefile and scripts for workflow |
| */reports*      | PDFs documents and reports |


* **Local File system**

| Folder             | Description |
| ------------------ | ----------- |
| */Master's Thesis* | Local copy of gitlab repo and wiki, endNote library and other organisational files |
| */mt-analysis*     | All files from analysis, raw and processed results. One subfolder for each analysis. Important for reproducibility files are copied to gitlab repository. |
| */mt-beast-dev*    | BEAST2, BDMM-Prime and other packages for IntelliJ project.|


* **Backups**

Weekly, every Friday. Backup of Master's Thesis, mt-analysis, mt-beast-dev local folders to Lacie hard drive and important results to my bs-home folder.
