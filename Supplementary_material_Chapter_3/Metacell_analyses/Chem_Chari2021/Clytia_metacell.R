######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Chem single cell dataset obtained from here: https://data.caltech.edu/records/pg2v4-0mm09
# Corresponding research paper: https://www.science.org/doi/10.1126/sciadv.abh1683


# Set library path for R on ALICE
.libPaths('/data/evassvis/software/R/4.0')

#Check current directory
getwd()

# set path
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab')
print(getwd())

# Load Seurat and related packages
library(Seurat)
#library(patchwork)
#library(dplyr)
library(SingleCellExperiment)
library(scater)
library(loomR)
library(SeuratDisk)
library(SeuratData)
# load metacell
library(metacell)

# make directories for metacell pipeline
if(!dir.exists("chem")) dir.create("chem")
# initialise the scdb (necessary for metacell pipeline)
scdb_init("chem/", force_reinit=T)

## chem dataset by Chari et al 2021 was in .h5ad format. 
## Converting from AnnData to Seurat via h5Seurat:

# convert the AnnData file to an h5Seurat
#Convert("bus_fs_combo_raw.h5ad", dest = "h5seurat", overwrite = TRUE) # this can be done only the first time.
# load the h5Seurat file into a Seurat object
chem <- LoadH5Seurat("bus_fs_combo_raw.h5seurat") # this has to be loaded everytime.
chem # summary of object

# View metadata data frame, stored in object@meta.data
chem[[]]

# write a table with umi counts (can be done only the first time!!)
#write.table(chem[["RNA"]]@counts, file = 'chem_TableFile_CountsPerCell.tsv', quote = FALSE, sep='\t', col.names = TRUE)
# and then import the table into metacell (has to be imported everytime)
mcell_import_scmat_tsv(mat_nm="chem", dset_nm="chem",
                       fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/chem_TableFile_CountsPerCell.tsv')

mat = scdb_mat("chem")
# the current size of the matrix is:
print(dim(mat@mat))
#[1] 46716 13673 (genes,cells)

# summary info:
print(mat)
#An object of class tgScMat, stat type umi.
#13673 cells by 46716 genes. median cell content 1802.

#### metacell pipeline

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

### Exploring and filtering the UMI matrix

# remove any potential ercc standards
erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="chem_noercc",mat_id="chem",ig_genes=erccs)
mat <- scdb_mat("chem_noercc")
print(dim(mat@mat))
#[1] 46716 13673

# Check umi distribution and choose cut-off. 
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
# Check plot to choose range:
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 100)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 200)
mcell_plot_umis_per_cell("chem",min_umis_cutoff = 225) # use as minimum
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 250)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 1000)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 15000)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 17000)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 20000)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 22000)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 24000)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 26000)
#mcell_plot_umis_per_cell("chem",min_umis_cutoff = 28000)
mcell_plot_umis_per_cell("chem",min_umis_cutoff = 30000) # use as maximum

cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>30000))
small_cells <- names(which(cell_sizes<225))
mcell_mat_ignore_cells("chem_filt","chem_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("chem_filt")
print(dim(mat@mat))
#[1] 46716 13410

### Selecting gene markers

# calculate statistics.
mcell_add_gene_stat(gstat_id="chem", mat_id="chem_filt", force=T)

# we skip the blacklist filtering, keep rest as tutorial
# Also we adjust the minimum tot UMIs to 225.
mcell_gset_filter_multi(gset_id = "chem_feats",gstat_id="chem",T_tot=225,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#Selected 1370 markers

mcell_plot_gstats(gstat_id="chem", gset_id="chem_feats")

### Building the balanced cell graph

mcell_add_cgraph_from_mat_bknn(mat_id="chem_filt",
                               gset_id = "chem_feats",
                               graph_id="chem_graph",
                               K=150,
                               dsamp=F)

### Resampling and generating the co-clustering graph

mcell_coclust_from_graph_resamp(
  coc_id="chem_coc1000",
  graph_id="chem_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)

print('Done resampling.')


# Try with alpha = 3 (instead of default 2)
mcell_mc_from_coclust_balanced(
  coc_id="chem_coc1000",
  mat_id= "chem_filt",
  mc_id= "chem_mc",
  K=20, min_mc_size=20, alpha=3)
# messages

print('Done coclustering.')

### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="chem_mc", mat_id = "chem_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="chem_mc_f",
                    mc_id="chem_mc",
                    mat_id="chem_filt",
                    T_lfc=3, plot_mats=F)
#messages


### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("chem_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("chem_mc_f",mc_f)
mc_f <- scdb_mc("chem_mc_f")

mcell_gset_from_mc_markers(gset_id="chem_markers", mc_id="chem_mc_f")

mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="chem_markers", mat_id="chem_filt",plot_cells = F)


### Make barplots for each gene of interest (how much expression in each metacell):

# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Chem_genes-and-gene-types_best_correct_codes.txt', sep='') # Check filename is correct!!
if(!dir.exists(barplots_path)) dir.create(barplots_path)
if(!dir.exists(metacells_path)) dir.create(metacells_path)
lfp = log2(mc_f@mc_fp)
genedata <- read.delim(genedata_path, sep=",", header = FALSE)
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  geneID = paste(genetype, '_', genename,'.png', sep = '')
  figname = paste(barplots_path, geneID, sep = '')
  if (genename %in% rownames(lfp)){  #create a plot only if the gene is in the metacells
    png(figname,h=400,w=1500);barplot(lfp[genename,],col=mc_f@colors,las=2,main=geneID,cex.main=3,cex.axis=1,ylab="log2FC",xlab="metacells");dev.off()
  }
}

print(paste('Done creating bar plots. Please see', barplots_path, '.'))


### Projecting metacells and cells in 2D

# Use the same pre-made config file as Amphimedon tutorial.
download.file("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/config.yaml","config.yaml")

tgconfig::override_params("config.yaml","metacell")
mcell_mc2d_force_knn(mc2d_id="chem_2dproj",mc_id="chem_mc_f", graph_id="chem_graph")

mcell_mc2d_plot(mc2d_id="chem_2dproj")


## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="chem_2dproj",  # change name here to appropriate dataset.
                         gene = genename, show_mc_ids = F, 
                         show_legend = T, color_cells=T, mat_ds = NULL, 
                         zero_sc_v = 0, one_sc_v =1, 
                         base_dir = basedir)  
  }
}

### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="chem_mc_f", graph_id="chem_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="chem_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="chem_mc_f",
                        graph_id="chem_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)

print('Done confusion matrix.')


#### Personalising heatmaps for genes of interest

## Make custom gset for Chem opsins
#sets = c('opsin','opsin')
#names(sets) = c('XLOC-037999','XLOC-041051') # recovered from reconciliation(s) results.
#gs = gset_new_gset(sets, 'Chem opsin set')
#scdb_add_gset("Chem_opsin_set", gs)

## plot heatmap for opsins using metacells
#mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_opsin_set", mat_id="chem_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/#Chem_heatmap_opsins_metacells.png')
## plot heatmap for opsins using cells
#mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_opsin_set", mat_id="chem_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_opsins_cells.png')


# Make custom gset for all good Chem phototr markers.
# Selected based on reconciliations.
# Both Rhabdomeric and ciliary.
sets = c('CNG_a_and_b','NCKX','NCKX','NCKX','NCKX','NCKX','NCKX','PDE6_a_and_b','RGS9','GB5','GC_2DEF','GCAP','G_aC','G_aC','opsins','opsins','G_b','calmodulin','calmodulin','arrestin','GRK','actin','actin','CamkII','DAGL','G_aR','INAD','IP3R','NINAC','PKC','PLC','rdgC')
names(sets) = c('XLOC-045142','TRINITY-DN14902-c0-g1-i6','XLOC-008369','XLOC-035197','XLOC-029206','TRINITY-DN9887-c0-g1-i3','TRINITY-DN9887-c0-g1-i4','XLOC-013466','XLOC-030008','XLOC-014565','XLOC-013404','XLOC-043843','XLOC-035588','XLOC-004433','XLOC-037999','XLOC-041051','XLOC-001915','XLOC-001548','XLOC-030364','TRINITY-DN4403-c0-g1-i6','XLOC-024731','XLOC-040487','XLOC-021750','XLOC-007220','TRINITY-DN8636-c0-g1-i3','XLOC-034293','TRINITY-DN1492-c0-g1-i4','XLOC-029997','XLOC-042993','TRINITY-DN2344-c0-g2-i7','XLOC-036877','XLOC-033745')
gs = gset_new_gset(sets, 'Chem phototr markers set')
scdb_add_gset("Chem_phototr_marks_set", gs)

# plot heatmap for good phototr markers using metacells
mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_phototr_marks_set", mat_id="chem_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_phototr_marks_best_metacells.png')
# plot heatmap for good phototr markers using cells
mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_phototr_marks_set", mat_id="chem_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_phototr_marks_best_cells.png')



# Make custom gset for Chem phototransduction genes (all)
#sets = c('CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','G_aC','G_aC','GB5','GC_2DEF','GCAP','NCKX','NCKX','NCKX','NCKX','NCKX','NCKX','NCKX','NCKX','PDE6_a_and_b','PDE6_a_and_b','PDE6_a_and_b','PDE6_a_and_b','RGS9','arrestin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','G_b','G_b','opsin','opsin','RK','CamkII','CamkII','DAGL','G_aR','INAD','INAD','INAD','INAD','IP3R','NINAC','PKC','PKC','PKC','PKC','PLC','PLC','PLC','PLC','PLC','PLC','rdgC')
#names(sets) = c('XLOC_010579','XLOC_010579','XLOC_045142','XLOC_035588','XLOC_004433','XLOC_014565','XLOC_013404','XLOC_043843','TRINITY_DN14902_c0_g1_i6','XLOC_008369','XLOC_035197','XLOC_029206','XLOC_035197','XLOC_035286','TRINITY_DN9887_c0_g1_i3','TRINITY_DN9887_c0_g1_i4','XLOC_013466','XLOC_040048','XLOC_044146','XLOC_011509','XLOC_030008','TRINITY_DN4403_c0_g1_i6','XLOC_001548','XLOC_001566','XLOC_001856','XLOC_001932','XLOC_002284','XLOC_004852','XLOC_009018','XLOC_009932','XLOC_012373','XLOC_012381','XLOC_012988','XLOC_033007','XLOC_014720','XLOC_017415','TRINITY_DN19575_c0_g1_i2','XLOC_030364','XLOC_031677','TRINITY_DN14543_c0_g1_i8','XLOC_036368','XLOC_036846','XLOC_037566','XLOC_039176','XLOC_039445','XLOC_040704','XLOC_041059','XLOC_043394','XLOC_043888','XLOC_044235','XLOC_045454','TRINITY_DN8170_c0_g1_i6','XLOC_001915','TRINITY_DN20474_c0_g1_i2','XLOC_037999','XLOC_041051','XLOC_024731','XLOC_007220','XLOC_013565','TRINITY_DN8636_c0_g1_i3','XLOC_034293','XLOC_001006','XLOC_005484','TRINITY_DN1492_c0_g1_i4','TRINITY_DN1492_c0_g1_i4','XLOC_029997','XLOC_042993','XLOC_003040','XLOC_020605','XLOC_030361','TRINITY_DN2344_c0_g2_i7','XLOC_003332','XLOC_014839','XLOC_043829','XLOC_033790','TRINITY_DN14106_c0_g1_i3','XLOC_036877','XLOC_033745') 
#gs = gset_new_gset(sets, 'Chem phototr set')
#scdb_add_gset("Chem_phototr_set", gs)

# plot heatmap for phototr genes (all) using metacells
#mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_phototr_set", mat_id="chem_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_phototr_metacells.png')
# plot heatmap for phototr genes (all) using cells
#mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_phototr_set", mat_id="chem_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_phototr_cells.png')


# Make custom gset for Chem Ciliary genes (including common)
#sets = c('CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','G_aC','G_aC','GB5','GC_2DEF','GCAP','NCKX','NCKX','NCKX','NCKX','NCKX','NCKX','NCKX','NCKX','PDE6_a_and_b','PDE6_a_and_b','PDE6_a_and_b','PDE6_a_and_b','RGS9','arrestin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','G_b','G_b','opsin','opsin','RK')
#names(sets) = c('XLOC_010579','XLOC_010579','XLOC_045142','XLOC_035588','XLOC_004433','XLOC_014565','XLOC_013404','XLOC_043843','TRINITY_DN14902_c0_g1_i6','XLOC_008369','XLOC_035197','XLOC_029206','XLOC_035197','XLOC_035286','TRINITY_DN9887_c0_g1_i3','TRINITY_DN9887_c0_g1_i4','XLOC_013466','XLOC_040048','XLOC_044146','XLOC_011509','XLOC_030008','TRINITY_DN4403_c0_g1_i6','XLOC_001548','XLOC_001566','XLOC_001856','XLOC_001932','XLOC_002284','XLOC_004852','XLOC_009018','XLOC_009932','XLOC_012373','XLOC_012381','XLOC_012988','XLOC_033007','XLOC_014720','XLOC_017415','TRINITY_DN19575_c0_g1_i2','XLOC_030364','XLOC_031677','TRINITY_DN14543_c0_g1_i8','XLOC_036368','XLOC_036846','XLOC_037566','XLOC_039176','XLOC_039445','XLOC_040704','XLOC_041059','XLOC_043394','XLOC_043888','XLOC_044235','XLOC_045454','TRINITY_DN8170_c0_g1_i6','XLOC_001915','TRINITY_DN20474_c0_g1_i2','XLOC_037999','XLOC_041051','XLOC_024731') 
#gs = gset_new_gset(sets, 'Chem ciliary set')
#scdb_add_gset("Chem_ciliary_set", gs)

# plot heatmap for ciliary genes (all) using metacells
#mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_ciliary_set", mat_id="chem_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_ciliary_metacells.png')
# plot heatmap for ciliary genes (all) using cells
#mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_ciliary_set", mat_id="chem_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_ciliary_cells.png')


# Make custom gset for Chem Rhabdomeric genes (including common)
#sets = c('arrestin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','G_b','G_b','opsin','opsin','RK','CamkII','CamkII','DAGL','G_aR','INAD','INAD','INAD','INAD','IP3R','NINAC','PKC','PKC','PKC','PKC','PLC','PLC','PLC','PLC','PLC','PLC','rdgC')
#names(sets) = c('TRINITY_DN4403_c0_g1_i6','XLOC_001548','XLOC_001566','XLOC_001856','XLOC_001932','XLOC_002284','XLOC_004852','XLOC_009018','XLOC_009932','XLOC_012373','XLOC_012381','XLOC_012988','XLOC_033007','XLOC_014720','XLOC_017415','TRINITY_DN19575_c0_g1_i2','XLOC_030364','XLOC_031677','TRINITY_DN14543_c0_g1_i8','XLOC_036368','XLOC_036846','XLOC_037566','XLOC_039176','XLOC_039445','XLOC_040704','XLOC_041059','XLOC_043394','XLOC_043888','XLOC_044235','XLOC_045454','TRINITY_DN8170_c0_g1_i6','XLOC_001915','TRINITY_DN20474_c0_g1_i2','XLOC_037999','XLOC_041051','XLOC_024731','XLOC_007220','XLOC_013565','TRINITY_DN8636_c0_g1_i3','XLOC_034293','XLOC_001006','XLOC_005484','TRINITY_DN1492_c0_g1_i4','TRINITY_DN1492_c0_g1_i4','XLOC_029997','XLOC_042993','XLOC_003040','XLOC_020605','XLOC_030361','TRINITY_DN2344_c0_g2_i7','XLOC_003332','XLOC_014839','XLOC_043829','XLOC_033790','TRINITY_DN14106_c0_g1_i3','XLOC_036877','XLOC_033745') 
#gs = gset_new_gset(sets, 'Chem rhabdomeric set')
#scdb_add_gset("Chem_rhabdomeric_set", gs)

# plot heatmap for rhabdomeric genes (all) using metacells
#mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_rhabdomeric_set", mat_id="chem_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_rhabdomeric_metacells.png')
# plot heatmap for rhabdomeric genes (all) using cells
#mcell_mc_plot_marks(mc_id="chem_mc_f", gset_id="Chem_rhabdomeric_set", mat_id="chem_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/figs/Chem_heatmap_rhabdomeric_cells.png')



### Finding genes expressed in a metacell of interest

lfp = log2(mc_f@mc_fp) 

# writing the metacells genes, a single file per metacell in the metacells_path folder
for (metacell in colnames(lfp)){
  lines = ''
  for (gene in rownames(lfp)){
    newline = paste(gene,lfp[gene,metacell],sep = ',')
    lines <- paste(lines,newline,'\n',sep='')
  }
  write(lines, file=paste(metacells_path,'metacell-',metacell,'.csv',sep=''))
}
print(paste('Done writing the metacell genes in files. Please see', metacells_path, '.'))


# We can extract the top X genes (by lfp) into file for each metacell:
# Extract top 100 for each metacell:
x <- 100  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Clytia_Chari_et_al_2021_PatcherLab/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
}


# writing the metacells genes with lfp > 0, a single file per metacell in the metacells_path folder
for (metacell in colnames(lfp)){
  lines = ''
  for (gene in rownames(lfp)){
    if (lfp[gene,metacell] > 0){
      newline = paste(gene,lfp[gene,metacell],sep = ',')
      lines <- paste(lines,newline,'\n',sep='')
    }
  }
  write(lines, file=paste('metacell-',metacell,'-genes-with-positive-lfp.csv',sep=''))
}
print(paste('Done writing the metacell genes in files. Please see', metacells_path, '.'))

