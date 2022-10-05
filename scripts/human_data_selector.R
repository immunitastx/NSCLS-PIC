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

library(tidyverse)

source("scripts/metacell_functions.r")
source("scripts/pic_parser.r")

dir.create("saved_work")
dir.create("figures")
scdb_init("saved_work", force_reinit=T)

pic_bucket <- "gs://immunitas_public_datasets/NSCLS_PIC_Cohen"

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
# NOW WE GET THE PIC-SEQ DATA
#############
outdir = paste0("figures/figure2")
dir.create(outdir)
supdir = paste0("figures/figureS2")
dir.create(supdir)
##############

lateral_genes = as.matrix(read.delim("annotations/lateral_genes.txt", stringsAsFactor=F, row.names=1))[,1]
human_lateral = names(which(lateral_genes == "human_lateral"))
bad_genes = grep("^RPL|^RPS|^MT[A-Z]|^HLA-|^AC[0-9]|^SNOR|^HIST[0-9]|^HSP[0-9]|^MIR[0-9]", rownames(sin_umis), v=T)
bad_markers = c("Metazoa_SRP", "MTND1P23", "MALAT1",bad_genes, human_lateral)

lr_features =  choose_lr_features(id, t_cells, dc_cells, bad_genes, must_haves = names(scdb_gset(id)@gene_set))
mle_features = choose_mle_features(id, id, t_cells, dc_cells, union(bad_markers, bad_genes), existing_list= names(scdb_gset(id)@gene_set))

mle_features = union(mle_features, "CXCL13")

############

id_d = "human_PIC" 
db_mat = scdb_mat(id_d)
db_umis = read_large_umis(id_d)
db_cells = db_mat@cells
umis = cbind(sin_umis, db_umis[ rownames(sin_umis),])
numis=800
ds = .downsamp(umis, numis)
comb = with(cell_stats, paste0(organ, "@", site, "@", patient)); names(comb) = rownames(cell_stats)
X = table(comb[ cells], cells %in% t_cells)
good_combs = names(which(rowSums(X < 20) == 0))
good_cells = cells[ comb[cells] %in% good_combs]
## Save time by reloading mle_res file
mle_res_file <- "output/immunitas_output/mle_res.rds"
if (!file.exists(mle_res_file)) {

mle_res = run_pic_seq(id, id, ds, intersect(good_cells, t_cells), intersect(good_cells, dc_cells), lr_features, mle_features, paste0(outdir, "/alpha_estimation.png"),
                      numis, comb = comb, reg = 1e-4)
mle_res$well = rownames(mle_res)

cell_stats = scdb_mat("all_human")@cell_metadata[ union(cells, db_cells),]
mle_res$type = as.vector(cell_stats[rownames(mle_res), "PIC"])
mle_res[ rownames(mle_res) %in% t_cells, "type"] = "T"
mle_res[ rownames(mle_res) %in% dc_cells, "type"] = "DC"
mle_res$sin_alpha = with(mle_res, ifelse(type == "PIC", 1 * (alpha >= 0.5), 1 * (well %in% t_cells)))
mle_res$sin_ll = with(mle_res, pic_ll_to_pair(id, id, ds[,well], a_mc, b_mc, sin_alpha, reg = 1e-4, markers = mle_features))
mle_res$diff = with(mle_res, ll - sin_ll)

saveRDS(mle_res, mle_res_file)
} else{
  mle_res <- readRDS(mle_res_file)
}

#######
bad_cells = names(sin_cl@mc)[ sin_cl@mc %in% bad_clusts]

gate = as.vector(cell_stats$sorting.scheme); names(gate) = rownames(cell_stats)
gate = c("T", "PIC", "APC")[ as.numeric(factor(gate))]; names(gate) = rownames(cell_stats)

nk_pic = rownames(mle_res)[ color2name[ sin_cl@colors[ mle_res$b_mc]] %in% c("NK", "NK_CX3CR1")]
good_pics = setdiff(rownames(mle_res)[ mle_res$type == "PIC" & mle_res$diff > 0 & mle_res$alpha > 0 & mle_res$alpha < 1], c(NA, nk_pic))

alpha = mle_res[ good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc[ good_pics]]]; names(parser_t) = good_pics
parser_dc = color2name[ sin_cl@colors[ dc_mc[ good_pics]]]; names(parser_dc) = good_pics

sin_cells = setdiff(intersect(cells, names(gate)[ gate %in% c("APC", "T")]), 
                    c(names(sin_names)[ sin_names %in% c("NK", "NK_CX3CR1", "Baso")], bad_cells))

parsed_pics_df <- cbind(parser_t, parser_dc) %>%
  as.data.frame() %>%
  rownames_to_column(var = "pic_ids") %>%
  rename(pic_tcell = parser_t, pic_myeloid = parser_dc)
write_tsv(parsed_pics_df,file = "output/immunitas_output/parsed_pics_cell_types.tsv")
ImmunitasR::SendToBucket(local_file_path = "output/immunitas_output/parsed_pics_cell_types.tsv",
                         bucket_location = pic_bucket)

##############

t_col = "limegreen"; dc_col =  "firebrick3"; db_col = "orange2"
cols = c(t_col, dc_col, db_col); names(cols) = c("T", "DC", "PIC")
pdf(paste0(supdir, "/FigS2c.pdf"), useDingbats=F)
plot(1,1,xlim = quantile(mle_res$diff, c(0,1), na.rm=T), ylim = c(0,1), type="n", axes=F, xlab="", ylab="")
abline(v=0, lwd=3); grid(col="black")
axis(1); axis(2)
with(mle_res, sapply(names(cols), function(x) lines(ecdf(diff[type == x]), col = cols[x], lwd=4)))
dev.off()

all_pics = db_mat@cells
patient_dist = table(cell_stats[ all_pics, "patient"], all_pics %in% good_pics)
patient_dist = patient_dist[ rowSums(patient_dist) > 0,]
rownames(patient_dist) = seq_len(nrow(patient_dist))
dist_n = patient_dist / rowSums(patient_dist)
ylim = c(0, max(dist_n[,2]) + 0.05)
pdf(paste0(supdir, "/FigS2d.pdf"), useDingbats=F)
X = barplot(dist_n[,2], ylim = ylim)
text(X, dist_n[,2] + 0.02, patient_dist[,2])
dev.off()

##############

cells = union(t_cells, dc_cells)
t_nms = rev(read.table(paste0(supdir, "/FigS2f-h_bottom.txt"), stringsAsFactor=F)[[1]])
dc_nms = rev(read.table(paste0(supdir, "/FigS2f-h_top.txt"), stringsAsFactor=F)[[1]])

a_lfp = log2(scdb_mc(paste0(id, "_a"))@mc_fp)
t_nms = t_nms[ order(max.col(a_lfp[t_nms, as.numeric(factor(intersect(clust_ord, t_clusts)))]))]

b_lfp = log2(scdb_mc(paste0(id, "_b"))@mc_fp)
dc_nms = dc_nms[ order(max.col(b_lfp[dc_nms, as.numeric(factor(intersect(clust_ord, dc_clusts)))]))]

nms = union(t_nms, dc_nms)
bad_cells = names(sin_cl@mc)[ sin_cl@mc %in% bad_clusts]

all_pics = good_pics


#### Fig 2b
cell_ord = sample(good_pics[ order(factor(parser_t, levels = lin_ord))])
disp_genes = rev(read.table( paste0(outdir, "/Fig2b.txt"), stringsAsFactor=F)[[1]])

IM = log(1 + 7 * sweep(ds[disp_genes,cell_ord],2,alpha[cell_ord],"*"))
vct = factor(parser_t, levels = lin_ord[1:11]); names(vct) = good_pics
IM2 = log(1 + 7 * ds[disp_genes, sample(intersect(colnames(ds), t_cells))])
vct2 = factor(sin_names[t_cells], levels = lin_ord[1:11]); names(vct) = good_pics
grad = colorRampPalette(c("white", "orange", "tomato", "mediumorchid4", "midnightblue"))(1000)
pdf(paste0(outdir, "/Fig2b.pdf"), useDingbats=F, height=5, width=10)
par(fig = c(0,0.5,0.1,1), mar = c(0.5,5,0.5,0.5))
image.2(IM2, col = grad, vct = vct2[colnames(IM2)], annotate="rows", zlim = c(0, max(c(IM, IM2)))); box(lwd=1)
par(fig = c(0,0.5,0,0.1), new=T)
image(matrix(seq_len(ncol(IM2))), axes=F, col = name2color[ sin_names[ colnames(IM2)[ order(vct2[ colnames(IM2)])]]]); box(lwd=1)
par(fig = c(0.5,1,0.1,1), new=T) #, mar = rep(0.5,4))
image.2(IM, col = grad, vct = vct[colnames(IM)], annotate="none", zlim = c(0, max(c(IM, IM2)))); box(lwd=1)
par(fig = c(0.5,1,0,0.1), new=T)
image(matrix(seq_len(ncol(IM))), axes=F, col = name2color[ parser_t[ colnames(IM)[ order(vct[ colnames(IM)])]]]); box(lwd=1)
dev.off()

#### Fig 2e