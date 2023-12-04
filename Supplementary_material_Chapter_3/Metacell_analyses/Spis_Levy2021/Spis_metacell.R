######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Spis single cell dataset obtained from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166901
# Corresponding research paper: https://www.sciencedirect.com/science/article/pii/S0092867421004402


# set library path
.libPaths('/data/evassvis/software/R/4.0')

# set Working Directory here:
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Stylopora_Levy_2021')

library(metacell)

# make directories
if(!dir.exists("spis")) dir.create("spis/")
scdb_init("spis/", force_reinit=T)

# import Spis single cell data set
# Levy et al 2021 data of adult:
# Note: use correct formated files.
mcell_import_scmat_10x(
  "spis",
  base_dir = "/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Stylopora_Levy_2021/Files_from_NCBI/adult/Renamed/change",
  matrix_fn = "matrix.mtx",
  genes_fn = "genes.tsv",
  cells_fn = "barcodes.tsv",
)

mat = scdb_mat("spis")
print(dim(mat@mat))
# [1] 37380 26850


if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")


### Exploring and filtering the UMI matrix

erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="spis_noercc",mat_id="spis",ig_genes=erccs)
mat <- scdb_mat("spis_noercc")

print(dim(mat@mat))
#[1] 37380 26850

# Check umi distribution and choose cut-off.
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
# Note: doesn't make sense to keep a cell with less than 80 counts!
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 80) 
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 85) # choose as min
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 90)
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 100)
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 150)
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 180)
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 1000)
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 7500)
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 8000)
mcell_plot_umis_per_cell("spis",min_umis_cutoff = 9000) # choose as max.
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 10000)
#mcell_plot_umis_per_cell("spis",min_umis_cutoff = 12000) 


# For Spis use lower cut-off 85 and higher cut-off 9000
cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>9000))
small_cells <- names(which(cell_sizes<85))
mcell_mat_ignore_cells("spis_filt","spis_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("spis_filt")
print(dim(mat@mat))
#[1] 37380 25567


### Selecting gene markers

# calculate statistics.
mcell_add_gene_stat(gstat_id="spis", mat_id="spis_filt", force=T)

# we skip the blacklist filtering, keep rest as tutorial
# Also we ajust the minimum tot UMIs to 85.
mcell_gset_filter_multi(gset_id = "spis_feats",gstat_id="spis",T_tot=85,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#Selected 225 markers

mcell_plot_gstats(gstat_id="spis", gset_id="spis_feats")


### Building the balanced cell graph

mcell_add_cgraph_from_mat_bknn(mat_id="spis_filt",
                               gset_id = "spis_feats",
                               graph_id="spis_graph",
                               K=150,
                               dsamp=F)


### Resampling and generating the co-clustering graph

mcell_coclust_from_graph_resamp(
  coc_id="spis_coc1000",
  graph_id="spis_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)
print("Resampling done.")

# Use alpha of 4 instead of default of 2.
mcell_mc_from_coclust_balanced(
  coc_id="spis_coc1000",
  mat_id= "spis_filt",
  mc_id= "spis_mc",
  K=20, min_mc_size=20, alpha=4)
# messages:


print("Co-clustering graph done.")


### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="spis_mc", mat_id = "spis_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="spis_mc_f",
                    mc_id="spis_mc",
                    mat_id="spis_filt",
                    T_lfc=3, plot_mats=F)
#Messages:


print("Finished removing outliers.")


### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("spis_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("spis_mc_f",mc_f)
mc_f <- scdb_mc("spis_mc_f")

mcell_gset_from_mc_markers(gset_id="spis_markers", mc_id="spis_mc_f")

mcell_mc_plot_marks(mc_id="spis_mc_f", gset_id="spis_markers", mat_id="spis_filt",plot_cells = F)



### Make barplots for each gene of interest (how much expression in each metacell):


# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Stylopora_Levy_2021/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Spis_genes-and-gene-types.txt', sep='') # check filename is correct!
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
mcell_mc2d_force_knn(mc2d_id="spis_2dproj",mc_id="spis_mc_f", graph_id="spis_graph")

mcell_mc2d_plot(mc2d_id="spis_2dproj")


## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene

for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="spis_2dproj", 
                         gene = genename, show_mc_ids = F, 
                         show_legend = T, color_cells=T, mat_ds = NULL, 
                         zero_sc_v = 0, one_sc_v =1, 
                         base_dir = basedir)  
  }
}



### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="spis_mc_f", graph_id="spis_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="spis_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="spis_mc_f",
                        graph_id="spis_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)


#### Personalising heatmaps for genes of interest


# Make custom gset for opsins
sets = c('opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin')
names(sets) = c('Spis_XP_022791651_1','Spis_XP_022782127_1','Spis_XP_022809903_1','Spis_XP_022793729_1','Spis_XP_022793717_1','Spis_XP_022803935_1','Spis_XP_022803931_1','Spis_XP_022803925_1','Spis_XP_022798035_1','Spis_XP_022792203_1') # opsins chosen from reconciliation(s) results.
gs = gset_new_gset(sets, 'Spis opsin set')
scdb_add_gset("Spis_opsin_set", gs)


# plot heatmap for opsins using metacells
mcell_mc_plot_marks(mc_id="spis_mc_f", gset_id="Spis_opsin_set", mat_id="spis_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Stylopora_Levy_2021/figs/Spis_heatmap_opsins_metacells.png')
# plot heatmap for opsins using cells
mcell_mc_plot_marks(mc_id="spis_mc_f", gset_id="Spis_opsin_set", mat_id="spis_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Stylopora_Levy_2021/figs/Spis_heatmap_opsins_cells.png')


# Make custom gset for all phototr genes
sets = c('CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','NCKX','NCKX','NCKX','NCKX','NCKX','PDE6_a_and_b','RGS9','GB5','G_aC','G_aC','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','G_b','G_g','calmodulin','calmodulin','arrestin','RK','RK','actin','actin','actin','actin','actin','CamkII','CamkII','DAGL','DAGL','G_aR','INAD','INAD','NINAC','PKC','PLC','PLC','rdgC','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRP')
names(sets) = c('Spis_XP_022800306_1','Spis_XP_022794597_1','Spis_XP_022779110_1','Spis_XP_022780572_1','Spis_XP_022779111_1','Spis_XP_022806875_1','Spis_XP_022806881_1','Spis_XP_022803387_1','Spis_XP_022800270_1','Spis_XP_022784892_1','Spis_XP_022796614_1','Spis_XP_022787547_1','Spis_XP_022781927_1','Spis_XP_022806086_1','Spis_XP_022788667_1','Spis_XP_022791651_1','Spis_XP_022782127_1','Spis_XP_022809903_1','Spis_XP_022793729_1','Spis_XP_022793717_1','Spis_XP_022803935_1','Spis_XP_022803931_1','Spis_XP_022803925_1','Spis_XP_022798035_1','Spis_XP_022792203_1','Spis_XP_022793483_1','Spis_XP_022799130_1','Spis_XP_022785429_1','Spis_XP_022785402_1','Spis_XP_022792530_1','Spis_XP_022792971_1','Spis_XP_022808146_1','Spis_XP_022780033_1','Spis_XP_022777933_1','Spis_XP_022795448_1','Spis_XP_022795461_1','Spis_XP_022779968_1','Spis_XP_022783573_1','Spis_XP_022779812_1','Spis_XP_022792240_1','Spis_XP_022780293_1','Spis_XP_022796132_1','Spis_XP_022779762_1','Spis_XP_022779808_1','Spis_XP_022785172_1','Spis_XP_022794886_1','Spis_XP_022785315_1','Spis_XP_022785277_1','Spis_XP_022781094_1','Spis_XP_022804911_1','Spis_XP_022798091_1','Spis_XP_022791467_1','Spis_XP_022789137_1','Spis_XP_022789667_1','Spis_XP_022789140_1','Spis_XP_022789647_1') 
gs = gset_new_gset(sets, 'Spis phototr set')
scdb_add_gset("Spis_phototr_set", gs)


# plot heatmap for phototr genes using metacells
mcell_mc_plot_marks(mc_id="spis_mc_f", gset_id="Spis_phototr_set", mat_id="spis_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Stylopora_Levy_2021/figs/Spis_heatmap_phototr_metacells.png')
# plot heatmap for phototr genes using cells
mcell_mc_plot_marks(mc_id="spis_mc_f", gset_id="Spis_phototr_set", mat_id="spis_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Stylopora_Levy_2021/figs/Spis_heatmap_phototr_cells.png')



# Extract top 100 genes
x <- 100  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Stylopora_Levy_2021/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
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
  write(lines, file=paste(metacells_path,'metacell-',metacell,'-genes-with-positive-lfp.csv',sep=''))
}
print(paste('Done writing the metacell genes in files. Please see', metacells_path, '.'))



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

