![image](/pics/covflo_logo.png)

## Introduction

Nextflow pipeline for generation of phylogenetic trees to be visualized with Auspice. 
COVFLO is written by JMC and adapted from snakefile, R, Python scripts written by Kimia Kamelian which generates 
phylogenies using tools like FastTree, Augur bioinformatic toolkit (MAFFT, IQ-TREE, TimeTree..), Goalign, cov2clusters, and TreeCluster that can be visualized in Auspice from [Nextstrain](https://docs.nextstrain.org/projects/auspice/en/stable/index.html). This pipeline consolidates the environments and scripts previously used for routine phylogentic analysis of SARS-CoV-2 sequences at the BCCDC into a portable, version-controlled command line tool.

## Table of Contents

- [Introduction](#introduction)
- [Quick-Start Guide](#quick-start%guide)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Input](#input)
- [Output](#output)
- [Workflow](#workflow)
- [References](#references)

## Quick-Start Guide

Run covflo pipeline:
```
nextflow run j3551ca/covflo -profile conda --conda_cache /path/to/caches --dir /home/user/sarscov2/input_data -r main
```
For details on available arguments, enter:
```
nextflow run j3551ca/covflo -r main --help
```

## Dependencies

[Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is required to build an environment with required workflow dependencies.

This bioinformatic pipeline requires Nextflow:
```
conda install -c bioconda nextflow
```
or download and add the nextflow executable to a location in your user $PATH variable:
```
curl -fsSL get.nextflow.io | bash
mv nextflow ~/bin/
```
Nextflow requires Java v8.0+, so check that it is installed:
```
java -version
```
The OS-independent conda environment activated upon running covflo is specified in the
```environment.yml``` file of the project directory and is built when 
```-profile conda``` is included in the command line. Nextflow will save
the environment to the project directory by default. Alternatively, the 
necessary conda environment can be saved to a different shared location 
accesible to compute nodes by adding ```--conda_cache /path/to/new/location/```.

## Installation

To copy the program into a directory of your choice, from desired directory run:
```
git clone https://github.com/j3551ca/covflo.git
cd covflo
nextflow run main.nf -profile conda --dir /home/user/sarscov2/input_data/
```
or run directly from Github using:
```
nextflow run j3551ca/covflo -profile conda --dir /home/user/sarscov2/input_data
```

## Input

The pipeline requires the following files which should be present in the config
and data folders of the directory containing sequences to be analyzed. These
are named the same within different directories - the only thing that needs to be changed
each run is the input directory, which can be specified with the --dir flag on the
command line.

- Multi-fasta file containing consensus sequences of interest [./data/sequences.fasta]
- Reference genome used to align reads to during guided assembly [./config/Ref.gb]
- File containing metadata for sequences under analysis [./data/metadata.csv]
- Excluded strains/ samples [./config/dropped_strains.txt]
- Strains/ samples to ensure are included [./config/included_strains.txt]
- Genomic cluster text files from previous build to be used as input for current [./config/SARS-CoV-2_{0.8,0.9}\_GenomicClusters.txt]
- Pairwise transmission probabilities (>0.8 or 0.9) between samples text files from previous build [./config/SARS-CoV-2_{0.8,0.9}\_TransProbs.txt]
- Colors used in final auspice visualization [./config/colors.csv]
- Sample latitudes and longitudes [./config/lat_longs.csv]
- Specifications for visualization in auspice (ex. title) [./config/auspice_config.json]

## Output

The output directories are 'results', 'auspice', and 'reports'.

results:
- filtered.fasta
- removedpercent.fasta
- replaced.fasta
- informative.fasta (no log)
- names.dedup
- deduped.fasta
- compressed.fasta
- weights
- fasttree.nwk
- resolvedtree.nwk
- blscaled.raxml{\*.startTree, \*.rba, \*.log, \*.bestTreeCollapsed, \*.bestTree, \*.bestModel}
- brlen_round.nwk
- collapse_length.nwk
- repopulate.nwk
- order.nwk
- tree.nwk
- branch_lengths.json
- nt_muts.json
- aa_muts.json
- SARS-CoV-2_{0.8,0.9}\_TransProbs.txt (pairwise transmission probabilities used as input for next tree build)
- SARS-CoV-2_{0.8,0.9}\_GenomicClusters.txt (genomic clusters used as input for next tree build)
- SARS-CoV-2_{0.8,0.9}\_ClustersSummary.csv
- tree_collapse_snp.nwk
- tc_cluster.tsv


*NOTE: the above files are listed in order of appearance in the 'main.nf' script, where process used to generate them as well as short description of process can be found in 'tag' directive.

auspice:
- ncov_na.json (final tree)

reports:
- covflo_usage.html
- covflo_timeline.html
- covflo_dag.html


## Workflow

![image](/pics/covflo_workflow.png)

## References

1. Hadfield, J. et al. NextStrain: Real-time tracking of pathogen evolution. Bioinformatics 34, 4121–3 (2018).

2. Huddleston J, Hadfield J, Sibley TR, Lee J, Fay K, Ilcisin M, Harkins E, Bedford T, Neher RA, Hodcroft EB, (2021). Augur: a bioinformatics toolkit for phylogenetic analyses of human pathogens. Journal of Open Source Software, 6(57), 2906, https://doi.org/10.21105/joss.02906

3. Katoh, K., Misawa, K., Kuma, K., & Miyata, T. (2002). MAFFT: a novel method for rapid
multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research,
30(14), 3059–3066. https://doi.org/10.1093/nar/gkf436

4. Nguyen, L.-T., Schmidt, H. A., Haeseler, A. von, & Minh, B. Q. (2014). IQ-TREE: A Fast and
Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Molecular Biology and Evolution, 32(1), 268–274. https://doi.org/10.1093/molbev/msu300

5. Sagulenko, P., Puller, V., & Neher, R. A. (2018). TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution, 4(1). https://doi.org/10.1093/ve/vex042

6. Lemoine, F., Gascuel, O. (2021). Gotree/Goalign: toolkit and Go API to facilitate the development of phylogenetic workflows,
NAR Genomics and Bioinformatics, 3(3), lqab075, https://doi.org/10.1093/nargab/lqab075

7. Steenwyk J.L., Buida III T.J., Li Y., Shen X-X., Rokas A. (2020) Clipkit: A multiple sequence alignment trimming software for accurate phylogenomic inference. PLOS Biology Available at: https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3001007. 

8. Price, M.N., Dehal, P.S., Arkin, A.P. (2009). FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix, Molecular Biology and Evolution, 26(7), 1641–1650, https://doi.org/10.1093/molbev/msp077

9. Kozlov, A.M., Darriba, D., Flouri, T., Morel, B., Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference, Bioinformatics, 35(21), 4453–4455, https://doi.org/10.1093/bioinformatics/btz305

10. Sobkowiak, B., Kamelian, K., Zlosnik, J. E. A., Tyson, J., Silva, A. G. D., Hoang, L. M. N., Prystajecky, N., & Colijn, C. (2022). Cov2clusters: genomic clustering of SARS-CoV-2 sequences. BMC genomics, 23(1), 710. https://doi.org/10.1186/s12864-022-08936-4

11. Balaban, M., Moshiri, N., Mai, U., Jia, X., Mirarab, S. (2019). "TreeCluster: Clustering biological sequences using phylogenetic trees." PLoS ONE. 14(8):e0221068. doi:10.1371/journal.pone.0221068
