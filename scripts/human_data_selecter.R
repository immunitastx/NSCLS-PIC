set.seed(1111)
library(metacell)
library(scales)
library(ape)
library(ggrepel)
library(tglkmeans)
library(Hmisc)
library(dendextend)
require("KernSmooth")
require("reshape2")
require("RANN")
require("plyr")
require("plotrix")
require("gplots")
require("parallel")
library("compositions")
require(glmnet)

source("scripts/metacell_functions.r")
source("scripts/pic_parser.r")

dir.create("saved_work")
dir.create("figures")
scdb_init("saved_work", force_reinit=T)

##################################
# Import the entire count matrix
message ("Preprocessing")

metadata = read.delim("annotations/TableS1.txt", stringsAsFactor=F)
write.table(metadata[ metadata$organism == "Human",], quote=F, sep = "\t", row.names=F, file = "annotations/metadata_human.txt")
#write.table(metadata[ metadata$organism == "Mouse",], quote=F, sep = "\t", row.names=F, file = "annotations/metadata_mouse.txt")

## Loads the data from the MARS-seq bumti-batch dataset
mcell_import_multi_mars("all_human", "annotations/metadata_human.txt", "output/umi.tab/", force = T)
#mcell_import_multi_mars("all_mouse", "annotations/metadata_mouse.txt", "output/umi.tab/", force = T)

human_mat = scdb_mat("all_human")
human_stats = human_mat@cell_metadata
singlets = rownames(human_stats)[ human_stats$PIC == "Singlets"]
doublets = rownames(human_stats)[ human_stats$PIC == "PIC"]
### Save the singlet and pic data for matt B
myeloid_singlets_ind <- human_stats$sorting.scheme == "Myelo+"
myeloid_matrix <- human_mat@mat[,myeloid_singlets_ind]
tcell_singlets_ind <- human_stats$sorting.scheme == "CD3+"
cd3_matrix <- human_mat@mat[,tcell_singlets_ind]
PIC_doublets_ind <- human_stats$sorting.scheme == "CD3+Myelo+"
pic_matrix <- human_mat@mat[,PIC_doublets_ind]
### 

mcell_mat_ignore_cells("human_singlets", "all_human", singlets, reverse=T)
mcell_mat_ignore_cells("human_PIC", "all_human", doublets, reverse=T)

mat = scdb_mat("all_human")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^MT-", nms, v=T)))

import_metacell_structure("human_singlets", "import/human/", "human_singlets", bad_genes)

####### End of preprocessing section #######


##### Begin figure 1 generation code ######

