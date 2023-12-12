######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Dmel single cell dataset obtained from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142787
# Corresponding research paper: https://www.nature.com/articles/s41586-020-2879-3


# Set library path for R on ALICE
.libPaths('/data/evassvis/software/R/4.0')

#Check current directory
getwd()


# load metacell
library(metacell)


# set path
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel')
print(getwd())

# make directories
if(!dir.exists("dmel")) dir.create("dmel/")
# initialise the scdb
scdb_init("dmel/", force_reinit=T)

### Conversion of Ozel dataset into dataset that can be used with metacell

## Load Seurat and related packages
#library(Seurat)
#library(patchwork)
#library(dplyr)

## choose the file (seurat object)
#filename <- 'GSE142787_Adult.rds'

## Give custom name to the file. Here calling it "Ozel".
#Ozel <- readRDS(filename)

## View metadata data frame, stored in object@meta.data
#Ozel[[]]

# Write a table with umi counts (to be done only once!)
#write.table(Ozel[["RNA"]]@counts, file = 'Ozel_TableFile_CountsPerCell.tsv', quote = FALSE, sep='\t', col.names = TRUE)
# and then import the table into metacell (has to be loaded everytime)
mcell_import_scmat_tsv(mat_nm="dmel", dset_nm="dmel",
                       fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/Ozel_TableFile_CountsPerCell.tsv')

#### metacell pipeline

mat = scdb_mat("dmel")
print(dim(mat@mat))
#[1]  12028 109743

# the current size of the matrix is:
#print(dim(mat@mat))
#[1]  12028 109743

print(mat)
#An object of class tgScMat, stat type umi.
#109743 cells by 12028 genes. median cell content 1805.


if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

### Exploring and filtering the UMI matrix

erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="dmel_noercc",mat_id="dmel",ig_genes=erccs)
mat <- scdb_mat("dmel_noercc")
print(dim(mat@mat))
#[1]  12028 109743

# Check umi distribution and choose cut-off.
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
# Check plot to choose range:
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 100)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 200)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 250)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 800)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 900)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 950)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 975)# keep as minimum
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 1000)
mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 4000)# keep as maximum.
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 5000)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 5500)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 6000)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 7000)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 7500) 
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 7750) 
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 8000) 
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 8500) 
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 9000)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 10000) 
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 11000)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 12000)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 13000)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 14000)
#mcell_plot_umis_per_cell("dmel",min_umis_cutoff = 15000) 


cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>4000))
small_cells <- names(which(cell_sizes<975))
mcell_mat_ignore_cells("dmel_filt","dmel_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("dmel_filt")
print(dim(mat@mat))
#[1]


### Selecting gene markers

# calculate statistics.
print('Doing mcell_add_gene_stat step')
mcell_add_gene_stat(gstat_id="dmel", mat_id="dmel_filt", force=T)

# we skip the blacklist filtering, keep rest as tutorial
# Also we adjust the minimum tot UMIs to 975.
print('Doing mcell_gset_filter_multi step')
mcell_gset_filter_multi(gset_id = "dmel_feats",gstat_id="dmel",T_tot=975,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#If filtering is done at min 250, max 7500: Selected 442 markers.
#If filtering is done at min 950, max 5500: Selected 385 markers.
#Selected ? markers

print('Doing mcell_plot_gstats step')
mcell_plot_gstats(gstat_id="dmel", gset_id="dmel_feats")

### Building the balanced cell graph

print('Doing mcell_add_cgraph_from_mat_bknn step')
mcell_add_cgraph_from_mat_bknn(mat_id="dmel_filt",
                               gset_id = "dmel_feats",
                               graph_id="dmel_graph",
                               K=150,
                               dsamp=F)

### Resampling and generating the co-clustering graph

print('Doing Resampling step')
#In tutorial: n_resamp = 1000
#In tutorial: p_resamp = 0.75
#If graphs are very large, we can decrease resampling. Here keep 1000.
mcell_coclust_from_graph_resamp(
  coc_id="dmel_coc1000",
  graph_id="dmel_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)

print('Doing co-clustering step')
#Here we have "dmel_coc1000".
# We are also changing alpha. Default is 2. We will put alpha = 50. Less harsh filtering of edges.
mcell_mc_from_coclust_balanced(
  coc_id="dmel_coc1000",
  mat_id= "dmel_filt",
  mc_id= "dmel_mc",
  K=20, min_mc_size=20, alpha=50)
# messages

### Removing outlier cells

print('Doing mcell_plot_outlier_heatmap step')
mcell_plot_outlier_heatmap(mc_id="dmel_mc", mat_id = "dmel_filt", T_lfc=3)

print('Doing mcell_mc_split_filt step')
mcell_mc_split_filt(new_mc_id="dmel_mc_f",
                    mc_id="dmel_mc",
                    mat_id="dmel_filt",
                    T_lfc=3, plot_mats=F)
#messages

### Creating heatmaps of metacells and genes

print('Creating heatmaps of metacells and genes')
mc_f<- scdb_mc("dmel_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("dmel_mc_f",mc_f)
mc_f <- scdb_mc("dmel_mc_f")

mcell_gset_from_mc_markers(gset_id="dmel_markers", mc_id="dmel_mc_f")

mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="dmel_markers", mat_id="dmel_filt",plot_cells = F)


### Make barplots for each gene of interest (how much expression in each metacell):

print('Make barplots for each gene of interest')
# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Dmel_genes-and-gene-types.txt', sep='') # Check filename is correct!!
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

print('Projecting metacells and cells in 2D')
# Use the same config file as Amphimedon tutorial.
download.file("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/config.yaml","config.yaml")

tgconfig::override_params("config.yaml","metacell")
mcell_mc2d_force_knn(mc2d_id="dmel_2dproj",mc_id="dmel_mc_f", graph_id="dmel_graph")

mcell_mc2d_plot(mc2d_id="dmel_2dproj")

## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="dmel_2dproj",  # change name here to appropriate dataset.
                         gene = genename, show_mc_ids = F,
                         show_legend = T, color_cells=T, mat_ds = NULL,
                         zero_sc_v = 0, one_sc_v =1,
                         base_dir = basedir)
  }
}

### Visualizing the metacell confusion matrix

print('Visualization metacell confusion matrix')
mc_hc <- mcell_mc_hclust_confu(mc_id="dmel_mc_f", graph_id="dmel_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="dmel_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="dmel_mc_f",
                        graph_id="dmel_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)

#### Personalising heatmaps for genes of interest

print('Making heatmaps for genes of interest')
## Note: In case of Drosophila, the components of the rhabdomeric PRC are well known and sure!
## We will anyway look at potential "ciliary" components.

# Make custom gset for Dmel opsins
sets = c('opsin_ninaE','opsin_Rh2','opsin_Rh3','opsin_Rh4','opsin_Rh5','opsin_Rh6','opsin_Rh7')
names(sets) = c('ninaE','Rh2','Rh3','Rh4','Rh5','Rh6','Rh7')
gs = gset_new_gset(sets, 'Dmel opsin set')
scdb_add_gset("Dmel_opsin_set", gs)

# plot heatmap for opsins using metacells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_opsin_set", mat_id="dmel_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_opsins_metacells.png')
# plot heatmap for opsins using cells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_opsin_set", mat_id="dmel_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_opsins_cells.png')


# Make custom gset for ALL components of Dmel PRCs (carefully selected, checked with flybase)
sets = c('Actin_Act5C','CamKII_CaMKII','DAGL_inaE','INAD_inaD','IP3R_Itp-r83A','NINAC_ninaC','PKC_inaC','PLC_norpA','rdgC_rdgC','TRP_trp','TRP_gamma_Trpgamma','TRPL_trpl','Arrestin_Arr2','Arrestin_Arr1','calmodulin_Cam','G_alpha_CG30054','G_alpha_q','G_beta_Gbeta76C','G_gamma_Ggamma30A','GRK_Gprk1','GRK_Gprk2','opsin_ninaE','opsin_Rh2','opsin_Rh3','opsin_Rh4','opsin_Rh5','opsin_Rh6','opsin_Rh7')
names(sets) = c('Act5C','CaMKII','inaE','inaD','Itp-r83A','ninaC','inaC','norpA','rdgC','trp','Trpgamma','trpl','Arr2','Arr1','Cam','CG30054','Galphaq','Gbeta76C','Ggamma30A','Gprk1','Gprk2','ninaE','Rh2','Rh3','Rh4','Rh5','Rh6','Rh7')
gs = gset_new_gset(sets, 'Dmel all PRC components set')
scdb_add_gset("Dmel_PRC_set", gs)

# plot heatmap for all dmel PRC components using metacells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_PRC_set", mat_id="dmel_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_all_PRC_components_metacells.png')
# plot heatmap for all dmel PRC components using cells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_PRC_set", mat_id="dmel_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_all_PRC_components_cells.png')


# Make custom gset for components of a potential "ciliary" PRC in dmel
# We will do 3 versions of this:

# a) Any potential ciliary gene: Using both the genes that are "common" but are surely used in true dmel
# rhabdomeric PRC and any additional "common" genes, plus all the candidate ciliary components.
# Basically, everything that is not "rhabdomeric only".

# b) Using improved list of genes based on reconciliation results.

# c) An even stricter list: best candidates based on homology with Hsap ciliary genes as found by reconciliations.

# a) CILIARY ALL
sets = c('arrestin','arrestin','Arrestin_Arr1','Arrestin_Arr2','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin_Cam','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','G_aC','G_aC','G_alpha_CG30054','G_alpha_q','G_alpha_q','G_b','G_beta_Gbeta76C','G_g','G_g','G_gamma_Ggamma30A','GB5','GC_2DEF','NCKX','NCKX','NCKX','NCKX','NCKX','opsin_ninaE','opsin_Rh2','opsin_Rh3','opsin_Rh4','opsin_Rh5','opsin_Rh6','opsin_Rh7','PDE6_a_and_b','PDE6_a_and_b','RGS9','RGS9','GRK','GRK_Gprk1')
names(sets) = c('krz','CG32683','Arr1','Arr2','CG11638','CG13526','CG13898','Eip63F-1','CG31960','CG17770','azot','CG30378','CG5024','Acam','CG17272','Cam','CngB','CngA','Cngl','CG42260','Galphai','Galphao','CG30054','CG17760','Galphaq','Gbeta13F','Gbeta76C','Ggamma1','CG43324','Ggamma30A','Gbeta5','CG34357','CG17167','zyd','Nckx30C','CG12061','CG1090','ninaE','Rh2','Rh3','Rh4','Rh5','Rh6','Rh7','Pde6','Pde11','RSG7','CG42450','Gprk2','Gprk1')
gs = gset_new_gset(sets, 'Dmel potential ciliary set')
scdb_add_gset("Dmel_ciliary_set", gs)

# plot heatmap for all dmel ciliary using metacells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_ciliary_set", mat_id="dmel_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_ciliary_ALL_metacells.png')
# plot heatmap for all dmel ciliary using cells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_ciliary_set", mat_id="dmel_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_ciliary_ALL_cells.png')


# b) CILIARY BETTER
sets = c('arrestin','arrestin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin','calmodulin_Cam','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','G_aC','G_aC','G_b','G_g','G_g','G_gamma_Ggamma30A','GB5','GC_2DEF','NCKX','NCKX','NCKX','NCKX','NCKX','opsin_ninaE','opsin_Rh2','opsin_Rh3','opsin_Rh4','opsin_Rh5','opsin_Rh6','opsin_Rh7','PDE6_a_and_b','RGS9','GRK')
names(sets) = c('krz','CG32683','CG13898','CG31960','azot','CG30378','Acam','CG17272','Cam','CngB','CngA','Cngl','CG42260','Galphai','Galphao','Gbeta13F','Ggamma1','CG43324','Ggamma30A','Gbeta5','CG34357','CG17167','zyd','Nckx30C','CG12061','CG1090','ninaE','Rh2','Rh3','Rh4','Rh5','Rh6','Rh7','Pde6','CG42450','Gprk2')
gs = gset_new_gset(sets, 'Dmel better ciliary set')
scdb_add_gset("Dmel_ciliary_set_better", gs)

# plot heatmap for all dmel ciliary using metacells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_ciliary_set_better", mat_id="dmel_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_ciliary_BETTER_metacells.png')
# plot heatmap for all dmel ciliary using cells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_ciliary_set_better", mat_id="dmel_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_ciliary_BETTER_cells.png')


# c) CILIARY BEST
sets = c('arrestin','arrestin','calmodulin_Cam','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','G_aC','G_b','G_gamma_Ggamma30A','GB5','GC_2DEF','NCKX','NCKX','NCKX','NCKX','NCKX','opsin_ninaE','opsin_Rh2','opsin_Rh3','opsin_Rh4','opsin_Rh5','opsin_Rh6','opsin_Rh7','PDE6_a_and_b','RGS9','GRK')
names(sets) = c('krz','CG32683','Cam','CngB','CngA','Cngl','CG42260','Galphai','Gbeta13F','Ggamma30A','Gbeta5','CG34357','CG17167','zyd','Nckx30C','CG12061','CG1090','ninaE','Rh2','Rh3','Rh4','Rh5','Rh6','Rh7','Pde6','CG42450','Gprk2')
gs = gset_new_gset(sets, 'Dmel best ciliary set')
scdb_add_gset("Dmel_ciliary_set_best", gs)

# plot heatmap for all dmel ciliary using metacells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_ciliary_set_best", mat_id="dmel_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_ciliary_BEST_metacells.png')
# plot heatmap for all dmel ciliary using cells
mcell_mc_plot_marks(mc_id="dmel_mc_f", gset_id="Dmel_ciliary_set_best", mat_id="dmel_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/figs/Dmel_heatmap_ciliary_BEST_cells.png')



### Finding genes expressed in metacells:

lfp = log2(mc_f@mc_fp)


print('Writing metacell genes: a file for each metacell')
## writing the metacells genes, a single file per metacell in the metacells_path folder
for (metacell in colnames(lfp)){
  lines = ''
  for (gene in rownames(lfp)){
    newline = paste(gene,lfp[gene,metacell],sep = ',')
    lines <- paste(lines,newline,'\n',sep='')
  }
  write(lines, file=paste(metacells_path,'metacell-',metacell,'.csv',sep=''))
}
print(paste('Done writing the metacell genes in files. Please see', metacells_path, '.'))


## Extract top 200 for each metacell:
print('Extracting top genes per metacell (lfp)')
x <- 200  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/MydataDmel/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
}



## writing the metacells genes with lfp > 0, a single file per metacell in the metacells_path folder
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
