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

pic_bucket <- pic_bucket

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

Matrix::writeMM(myeloid_matrix,file = "output/immunitas_output/singlet_myeloid_pic.mtx")
Matrix::writeMM(cd3_matrix,file = "output/immunitas_output/singlet_cd3_pic.mtx")
Matrix::writeMM(pic_matrix,file = "output/immunitas_output/doublet_cd3_myeloid_pic.mtx")

ImmunitasR::SendToBucket(local_file_path = "output/immunitas_output/singlet_myeloid_pic.mtx",
                         bucket_location = pic_bucket)
ImmunitasR::SendToBucket(local_file_path = "output/immunitas_output/singlet_cd3_pic.mtx",
                         bucket_location = pic_bucket)
ImmunitasR::SendToBucket(local_file_path = "output/immunitas_output/doublet_cd3_myeloid_pic.mtx",
                         bucket_location = pic_bucket)

ImmunitasR::SendToBucket(local_file_path = "annotations/metadata_human.txt",
                         bucket_location = pic_bucket)

### 

mcell_mat_ignore_cells("human_singlets", "all_human", singlets, reverse=T)
mcell_mat_ignore_cells("human_PIC", "all_human", doublets, reverse=T)

mat = scdb_mat("all_human")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^MT-", nms, v=T)))

import_metacell_structure("human_singlets", "import/human/", "human_singlets", bad_genes)

####### End of preprocessing section #######


##### Begin figure 1 generation code ######
set.seed(1111)

id = "human_singlets"
## Below they load in all the metacell data
sin_2d = scdb_mc2d(id); sin_cl = scdb_mc(id); sin_mat = scdb_mat(id)

cells = names(sin_cl@mc)
cell_stats = sin_mat@cell_metadata[cells,]
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)

outdir = paste0("figures/figure1")
supdir = paste0("figures/figureS1")
dir.create(outdir)
dir.create(supdir)

sin_stats = sin_mat@cell_metadata[names(sin_cl@mc),]
sin_umis = read_large_umis(id, cells = cells)
sin_n = sweep(sin_umis,2,colSums(sin_umis),"/") * 1000

### Here is how you can see their annotation selection for single cells
color_scheme = sin_cl@color_key
color2name = as.vector(color_scheme$group); names(color2name) = color_scheme$color
name2color = as.vector(color_scheme$color); names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)

sin_comb = with(sin_stats, paste0(sorting.scheme, ":", PIC, "@", organ, "@", site, "@", patient))
names(sin_comb) = rownames(sin_stats)

anno_colors = as.matrix(read.delim("annotations/human_lung_annotations.txt", stringsAsFactor=F, row.names=1))
anno_names = anno_colors[,2]
lin_ord = names(anno_names)
umicount = colSums(sin_umis)

bad_clusts = colnames(lfp)[ lfp["TRBC2",] > -1 & lfp["CST3",] > 2]
## Separate population
t_pops = names(which(anno_colors[,3] == "T"))
dc_pops = names(which(anno_colors[,3] == "Myeloid"))

t_clusts = setdiff(which(color2name[ sin_cl@colors] %in% t_pops), bad_clusts)
t_cells = setdiff(names(sin_cl@mc)[ sin_cl@mc %in% t_clusts], names(which(umicount > 1e4)))

dc_clusts =  setdiff(which(color2name[ sin_cl@colors] %in% dc_pops), bad_clusts)
dc_cells = setdiff(names(sin_cl@mc)[ sin_cl@mc %in% dc_clusts], names(which(umicount > 1e4)))

cells = union(t_cells, dc_cells)
#############
### In this section I will try to create two annotated metacell files one for 
### e_gc another for cov_gc for DC and tcell
dim(sin_cl@e_gc)
sin_cl@color_key
init_df <- as.data.frame(sin_cl@colors)
gene_to_colors <- init_df %>%
  rename(color = `sin_cl@colors`) %>% 
  full_join(sin_cl@color_key %>% select(-gene))
ind <- match(sin_cl@colors, sin_cl@color_key$color) 
meta_cell_anno <- cbind(sin_cl@colors, sin_cl@color_key$group[ind]) |>
  as.data.frame() |>
  select(cell_type = V2,color = V1 )

meta_cell_anno

meta_cell_anno$metacell <- paste0("metacell_", seq(nrow(meta_cell_anno )))

## Get the singlets that make up each metacell
singlet_metacell_anno <- sin_cl@mc %>% as.data.frame() %>%
  rownames_to_column("singlet") %>%
  rename(., metacell = `.`) %>%
  mutate(metacell = paste0("metacell_", metacell))

singlet_metacell_anno <- left_join(singlet_metacell_anno,
          meta_cell_anno)
write_tsv(singlet_metacell_anno,
          file = "output/immunitas_output/metacell_annos_singlets.tsv")
write_tsv(meta_cell_anno,
          file = "output/immunitas_output/metacell_annos_metacells.tsv")
## Extract metacell matrices
metacell_colnames <- paste0("metacell_", colnames(sin_cl@e_gc))
umi_average_mat <- sin_cl@e_gc
umi_percent_exprs_mat <- sin_cl@cov_gc
colnames(umi_average_mat) <- metacell_colnames
colnames(umi_percent_exprs_mat) <- metacell_colnames
write_tsv(as.data.frame(umi_average_mat) %>%rownames_to_column(var = "genes") ,
          file = "output/immunitas_output/metacell_singlet_average_umi_exprs.tsv")
write_tsv(as.data.frame(umi_percent_exprs_mat) %>% rownames_to_column(var = "genes") ,
          file = "output/immunitas_output/metacell_singlet_proportion_umi_exprs.tsv")


ImmunitasR::SendToBucket(local_file_path = "output/immunitas_output/metacell_annos_singlets.tsv",
                         bucket_location = pic_bucket)
ImmunitasR::SendToBucket(local_file_path = "output/immunitas_output/metacell_annos_metacells.tsv",
                         bucket_location = pic_bucket)
ImmunitasR::SendToBucket(local_file_path = "output/immunitas_output/metacell_singlet_average_umi_exprs.tsv",
                         bucket_location = pic_bucket)
ImmunitasR::SendToBucket(local_file_path = "output/immunitas_output/metacell_singlet_proportion_umi_exprs.tsv",
                         bucket_location = pic_bucket)


### Check that this is correct by checking expression of BIRC3 and LAMP3
# cbind(meta_cell_anno, t(sin_cl@cov_gc[rownames(sin_cl@cov_gc) %in% c("BIRC3", "LAMP3"),])) %>%
#   as.data.frame() %>% filter(cell_type == "MigDC")

##############

lateral_genes = as.matrix(read.delim("annotations/lateral_genes.txt", stringsAsFactor=F, row.names=1))[,1]
human_lateral = names(which(lateral_genes == "human_lateral"))
bad_genes = grep("^RPL|^RPS|^MT[A-Z]|^HLA-|^AC[0-9]|^SNOR|^HIST[0-9]|^HSP[0-9]|^MIR[0-9]", rownames(sin_umis), v=T)
bad_markers = c("Metazoa_SRP", "MTND1P23", "MALAT1",bad_genes, human_lateral)

lr_features =  choose_lr_features(id, t_cells, dc_cells, bad_genes, must_haves = names(scdb_gset(id)@gene_set))
mle_features = choose_mle_features(id, id, t_cells, dc_cells, union(bad_markers, bad_genes), existing_list= names(scdb_gset(id)@gene_set))

mle_features = union(mle_features, "CXCL13")

############

mc = names(table(sin_cl@mc[cells]))
mc = mc[ order(factor(color2name[ sin_cl@colors[ as.numeric(mc)]], levels = lin_ord))]
patient = as.vector(sin_stats$patient); names(patient) = rownames(sin_stats)
sample_dist = table(sin_cl@mc[cells], patient[cells])
sample_dist = t(apply(sample_dist[mc,],1,sort,T))
dist_n = sample_dist / rowSums(sample_dist)
