######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Cint single cell dataset obtained from here: 
# Corresponding research paper: https://www.sciencedirect.com/science/article/pii/S0012160618300137


# set library path
.libPaths('/data/evassvis/software/R/4.0')

# set Working Directory here: https://data.mendeley.com/datasets/r7vzs2gt7x/1
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Cint_Sharma_et_al_2019')

library(metacell)

# make directories
if(!dir.exists("cint")) dir.create("cint/")
scdb_init("cint/", force_reinit=T)

# import Cint single cell data set
# Sharma et al 2019 data is 10X, so load 10X matrix:
mcell_import_scmat_10x(
  "cint",
  base_dir = "/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Cint_Sharma_et_al_2019/Sharma2019_dataset",
)

mat = scdb_mat("cint")
print(dim(mat@mat))
# [1] 15232  2737

# Note, extract from paper:
# " After quality control, the final dataset comprised of 2607 cells 
# with 13,435 genes detected across them, out of 15,288 total genes 
# annotated in the C. robusta genome. "

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")


### Exploring and filtering the UMI matrix

erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="cint_noercc",mat_id="cint",ig_genes=erccs)
mat <- scdb_mat("cint_noercc")

print(dim(mat@mat))
#[1] 15232  2737


# Check umi distribution and choose cut-off. 
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 180)
mcell_plot_umis_per_cell("cint",min_umis_cutoff = 300) # choose as min.
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 750)
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 800)
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 1000)
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 7500)
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 8000)
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 12000)
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 14000)
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 16000) 
#mcell_plot_umis_per_cell("cint",min_umis_cutoff = 17000)
mcell_plot_umis_per_cell("cint",min_umis_cutoff = 18000) # choose as max.

# For Cint use lower cut-off 300 and higher cut-off 18000
cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>18000))
small_cells <- names(which(cell_sizes<300))
mcell_mat_ignore_cells("cint_filt","cint_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("cint_filt")
print(dim(mat@mat))
#[1] 15232  2725


### Selecting gene markers

# calculate statistics.
mcell_add_gene_stat(gstat_id="cint", mat_id="cint_filt", force=T)

# we skip the blacklist filtering, keep rest as tutorial
# Also we ajust the minimum tot UMIs to 300.
mcell_gset_filter_multi(gset_id = "cint_feats",gstat_id="cint",T_tot=300,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#Selected 1011 markers

mcell_plot_gstats(gstat_id="cint", gset_id="cint_feats")


### Building the balanced cell graph

mcell_add_cgraph_from_mat_bknn(mat_id="cint_filt",
                               gset_id = "cint_feats",
                               graph_id="cint_graph",
                               K=150,
                               dsamp=F)
#will build balanced knn graph on 2725 cells and 1011 genes, this can be a bit heavy for >20,000 cells

### Resampling and generating the co-clustering graph

mcell_coclust_from_graph_resamp(
  coc_id="cint_coc1000",
  graph_id="cint_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)
print("Resampling done.")

# Default alpha is 2.
# To be more stringent with filtering, will try alpha=1.5.
mcell_mc_from_coclust_balanced(
  coc_id="cint_coc1000",
  mat_id= "cint_filt",
  mc_id= "cint_mc",
  K=20, min_mc_size=20, alpha=1.5)
# messages:
#filtered 445257 left with 143252 based on co-cluster imbalance
#building metacell object, #mc 33
#add batch counts
#compute footprints
#compute absolute ps
#compute coverage ps
#reordering metacells by hclust and most variable two markers
#reorder on KH2013:KH.C6.197 vs KH2013:KH.C5.298

print("Co-clustering graph done.")


### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="cint_mc", mat_id = "cint_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="cint_mc_f",
                    mc_id="cint_mc",
                    mat_id="cint_filt",
                    T_lfc=3, plot_mats=F)
#Messages:
#starting split outliers 
#splitting metacell 11
#splitting metacell 27
#splitting metacell 33
#add batch counts
#compute footprints
#compute absolute ps
#compute coverage ps

print("Finished removing outliers.")


### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("cint_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("cint_mc_f",mc_f)
mc_f <- scdb_mc("cint_mc_f")

mcell_gset_from_mc_markers(gset_id="cint_markers", mc_id="cint_mc_f")

mcell_mc_plot_marks(mc_id="cint_mc_f", gset_id="cint_markers", mat_id="cint_filt",plot_cells = F)


### Make barplots for each gene of interest (how much expression in each metacell):

# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Cint_Sharma_et_al_2019/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Cint_genes-and-gene-types.txt', sep='') # check filename is correct!
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

# As a test use the same config file as Amphimedon tutorial.
download.file("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/config.yaml","config.yaml")

tgconfig::override_params("config.yaml","metacell")
mcell_mc2d_force_knn(mc2d_id="cint_2dproj",mc_id="cint_mc_f", graph_id="cint_graph")

mcell_mc2d_plot(mc2d_id="cint_2dproj")


## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene

for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="cint_2dproj", 
                         gene = genename, show_mc_ids = F, 
                         show_legend = T, color_cells=T, mat_ds = NULL, 
                         zero_sc_v = 0, one_sc_v =1, 
                         base_dir = basedir)  
  }
}




### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="cint_mc_f", graph_id="cint_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="cint_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="cint_mc_f",
                        graph_id="cint_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)


#### Personalising heatmaps for genes of interest

# Make custom gset for opsins
sets = c('opsin','opsin','opsin','opsin')
names(sets) = c('KH2013:KH.L38.6','KH2013:KH.L171.13','KH2013:KH.C4.578','KH2013:KH.C7.385') # four opsins chosen from reconciliation(s) results.
gs = gset_new_gset(sets, 'Cint opsin set')
scdb_add_gset("Cint_opsin_set", gs)


# plot heatmap for opsins using metacells
mcell_mc_plot_marks(mc_id="cint_mc_f", gset_id="Cint_opsin_set", mat_id="cint_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Cint_Sharma_et_al_2019/figs/Cint_heatmap_opsins_metacells.png')
# plot heatmap for opsins using cells
mcell_mc_plot_marks(mc_id="cint_mc_f", gset_id="Cint_opsin_set", mat_id="cint_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Cint_Sharma_et_al_2019/figs/Cint_heatmap_opsins_cells.png')


# Making gset for all phototr genes:
sets = c('CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','NCKX','PDE6_a_and_b','RGS9','GB5','GC_2DEF','recoverin','G_aC','G_aC','G_aC','opsin','opsin','opsin','opsin','G_b','calmodulin','arrestin','RK','RK','RK','actin','actin','actin','actin','actin','CamkII','DAGL','G_aR','G_aR','INAD','INAD','IP3R','PKC','PLC','PLC','PLC','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL')
names(sets) = c('KH2013:KH.C2.249','KH2013:KH.C7.605','KH2013:KH.L42.6','KH2013:KH.L96.97','KH2013:KH.L50.5','KH2013:KH.L18.65','KH2013:KH.C2.958','KH2013:KH.L96.22','KH2013:KH.L18.91','KH2013:KH.C7.343','KH2013:KH.C1.612','KH2013:KH.C1.612','KH2013:KH.C9.789','KH2013:KH.L38.6','KH2013:KH.L171.13','KH2013:KH.C4.578','KH2013:KH.C7.385','KH2013:KH.C7.88','KH2013:KH.C3.968','KH2013:KH.C1.1125','KH2013:KH.C8.374','KH2013:KH.L112.32','KH2013:KH.L56.1','KH2013:KH.S1440.1','KH2013:KH.C8.649','KH2013:KH.S1440.1','KH2013:KH.C1.120','KH2013:KH.S1440.1','KH2013:KH.C2.758','KH2013:KH.C2.705','KH2013:KH.L154.32','KH2013:KH.L44.22','KH2013:KH.C5.57','KH2013:KH.C5.57','KH2013:KH.C8.70','KH2013:KH.C3.45','KH2013:KH.C4.381','KH2013:KH.L138.1','KH2013:KH.S2332.1','KH2013:KH.C10.392','KH2013:KH.C10.619','KH2013:KH.C2.1064','KH2013:KH.C2.757','KH2013:KH.C4.585','KH2013:KH.C8.233','KH2013:KH.C8.571','KH2013:KH.C9.379') 
gs = gset_new_gset(sets, 'Cint phototr set')
scdb_add_gset("Cint_phototr_set", gs)

# plot heatmap for phototr using metacells
mcell_mc_plot_marks(mc_id="cint_mc_f", gset_id="Cint_phototr_set", mat_id="cint_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Cint_Sharma_et_al_2019/figs/Cint_heatmap_phototr_metacells.png')
# plot heatmap for phototr using cells
mcell_mc_plot_marks(mc_id="cint_mc_f", gset_id="Cint_phototr_set", mat_id="cint_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Cint_Sharma_et_al_2019/figs/Cint_heatmap_phototr_cells.png')



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


# Extract top 100 genes
x <- 100  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Cint_Sharma_et_al_2019/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
}

# Extract all genes with lfp >0 for each metacell
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

