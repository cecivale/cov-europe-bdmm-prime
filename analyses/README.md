# Analyses

| Nr | Analysis name      | Description |  Files | MCMC | ESS | Results | Notes |
| -- | ------------------ | ----------- | ------ | ---- | --- | ------- | ----- | 
| 1  | **[200904_sarah](analyses/200904_sarah_analysis_specs.md)**  | Sarah's analysis from *[1]* repeated, same XML and sequence data. | * Seq. data:  [sarah_europe<br>_demes.fasta](data/sarah_europe_demes.fasta) <br>* XML: [200904_sarah<br>_europe_demes.xml](analyses/200904_sarah_europe_demes.xml)  | 1e7 | ~ 300 |  | Similar results to original analysis. |
| 2  | **[200907_sarah_traj](analyses/200907_sarah_traj_specs.md)** | Sarah's analysis but adding trajectory logger BDMM-Prime |  * Seq. data:  [sarah_europe<br>_demes.fasta](data/sarah_europe_demes.fasta) <br>* XML: [200907_sarah<br>_traj.xml](analyses/200904_sarah_europe_demes.xml)  | 12000 | <<< 300 |  | Too slow, downsample dataset |
| 3  | **[200910_ds_europe0](analyses/200910_ds_europe0_specs.md)** | Downsample dataset, same analysis specs than for analysis 2 | * Seq. data:  [200910_ds_<br>europe0.fasta](data/200910_ds_europe0.fasta) <br>* XML: [200910_ds<br>_europe0.xml](analyses/200910_ds_europe0.xml)  | | |  | |

*[1] Nadeau, S. A., et al. (2020). "The origin and early spread of SARS-CoV-2 in Europe." medRxiv: 2020.2006.2010.20127738.*