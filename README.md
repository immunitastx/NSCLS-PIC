---
editor_options: 
  markdown: 
    wrap: 72
---

## Tim README

This repository is required to download and preprocess the PIC-seq data
from the paper referenced below. The code currently contained will only
help in analyses related to the human PIC-seq data available at
GSE160903.

### Files of interest

1.  scripts/get_geo_files.R

The purpose of this file is to download the the GEO files and prepare
them for preprocessing, one issue however is you need to copy and paste
the shell commands into the terminal, I was unable to use `system()` to
get them working properly.

1.  scripts/human_data_selector.R

This script processes the data as in the paper, creates several outputs
files, and sends them to the cloud.

### File output
#### Non-metacell information
-   _output/immunitas_output/singlet_myeloid_pic.mtx_: A sparse matrix of
    the umi counts for singlets sorted on myeloid markers
-   _output/immunitas_output/singlet_cd3_pic.mtx_: A sparse matrix of the
    umi counts for singlets sorted on CD3+ cells
-   _output/immunitas_output/doublet_cd3_myeloid_pic.mtx_: A sparse matrix
    of the umi counts for doublets sorted on CD3+ and myeloid markers
#### Metacell information
-   _output/immunitas_output/metacell_singlet_average_umi_exprs.tsv_:
    A tsv file of the average umi expression across singlets contained in a
    metacell. From the sin_cl object in the @e_gc slot.
-   _output/immunitas_output/metacell_singlet_proportion_umi_exprs.tsv_:
    A tsv file of the proportional umi expression across singlets contained in a
    metacell. From the sin_cl object in the @e_cov slot.


## Original README below

This repository contains the PIC-seq all the needed code and metadata
files to analyze and generate figures for Cohen M. & Giladi A. et al.
Nature Cancer 2022

In order to run the scripts, download processed data from GSE160903 to
the folder output/umi.tab

Unzip and change file names by running these shell commands (in
output/umi.tab): gzip -d \* ls -1 \| awk -F'\_' '{print \$0,\$2}' \|
xargs -n 2 mv

Download published data into their respective folders: GSE123139
output/published_data/melanoma/ EGAD00001006608
output/published_data/breast/ GSE135382 output/published_data/pic-seq/

To start analysis, run from the root directory: Rscript scripts/run.r

Please send questions to Amir Giladi:
[aygoldberg\@gmail.com](mailto:aygoldberg@gmail.com)

Required R packages:

Package version glmnet 2.0-16 foreach 1.4.4 Matrix 1.2-18 compositions
1.40-2 bayesm 3.1-0.1 energy 1.7-5 robustbase 0.93-3 tensorA 0.36.1
gplots 3.0.1.1 plotrix 3.7-4 plyr 1.8.4 RANN 2.6.1 reshape2 1.4.3
KernSmooth 2.23-15 dendextend 1.9.0 Hmisc 4.2-0 Formula 1.2-3 survival
3.2-3 lattice 0.20-38 tglkmeans 0.2.0 ggrepel 0.8.1 ggplot2 3.3.2.9000
ape 5.2 scales 1.0.0 metacell 0.3.41 tgstat 2.3.5 misha 4.0.10
