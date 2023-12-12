######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Hsap single cell dataset obtained from here: https://zenodo.org/record/5515631#.YWVQBWLMKUl
# Corresponding research paper: https://www.embopress.org/doi/full/10.15252/embj.2018100811


# set library path
.libPaths('/data/evassvis/software/R/4.0')

# set Working Directory here:
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Homo_Lukowski2019/10Xdata_Lukowski2019/')
print(getwd())

# load metacell
library("metacell")

# make directories for metacell analysis
if(!dir.exists("hsap")) dir.create("hsap/")
# initialise the scdb
scdb_init("hsap/", force_reinit=T)

### Input

# Using raw counts matrix available here: https://zenodo.org/record/5515631#.YWVQBWLMKUl
# That is in .csv format. We need .tsv.
# See script preprocessing_code.py.

# Dataset is converted to .tsv format. Import as tsv:
mcell_import_scmat_tsv(mat_nm="hsap", dset_nm="hsap",
                       fn='cleaned_matrix.tsv')

mat = scdb_mat("hsap")
print(dim(mat@mat))
#[1] 21989 20009
print(mat)
#An object of class tgScMat, stat type umi.
#20009 cells by 21989 genes. median cell content 2038.795.

# Make directory for figures
if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

### Exploring and filtering the UMI matrix

# Filter "ERCC" genes in this dataset.
erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
#Check ERCC genes found:
erccs
# Going to filter out also "MT-" (mitochondrial) genes.
mt <- rownames(mat@mat)[grepl("MT-",rownames(mat@mat))]
#Check MT genes found:
mt
# Merge two lists:
erccsmt <- append(erccs, mt)
#Check
erccsmt
#Now filter out all genes:
mcell_mat_ignore_genes(new_mat_id="hsap_noerccsmt",mat_id="hsap",ig_genes=erccsmt)
mat <- scdb_mat("hsap_noerccsmt")
print(dim(mat@mat))
#21967 20009

# Check umi distribution and choose cut-off. 
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
# Check plot to choose range:
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 100)
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 300)
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 400) # choose as low cut-off
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 500)
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 800)
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 900) 
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 1000)
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 3500) 
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 5000)
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 7000)
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 7500)
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 8000)
mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 8500) # choose as high cut-off
#mcell_plot_umis_per_cell("hsap",min_umis_cutoff = 10000)


cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>8500))
small_cells <- names(which(cell_sizes<400))
mcell_mat_ignore_cells("hsap_filt","hsap_noerccsmt",ig_cells=c(small_cells,large_cells)) 
mat <- scdb_mat("hsap_filt")
print(dim(mat@mat))
#[1] 21967 19626

### Selecting gene markers

# calculate statistics.
print('Doing mcell_add_gene_stat step')
mcell_add_gene_stat(gstat_id="hsap", mat_id="hsap_filt", force=T)
#Calculating gene statistics... will downsamp
#done downsamp
#will gen mat_n
#done gen mat_n
#done computing basic gstat, will compute trends
#..done


# we skip the blacklist filtering, keep rest as tutorial
# Also we adjust the minimum tot UMIs to 400.
print('Doing mcell_gset_filter_multi step')
mcell_gset_filter_multi(gset_id = "hsap_feats",gstat_id="hsap",T_tot=400,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#Selected 217 markers

print('Doing mcell_plot_gstats step')
mcell_plot_gstats(gstat_id="hsap", gset_id="hsap_feats")


### Building the balanced cell graph

print('Doing mcell_add_cgraph_from_mat_bknn step')
mcell_add_cgraph_from_mat_bknn(mat_id="hsap_filt",
                               gset_id = "hsap_feats",
                               graph_id="hsap_graph",
                               K=150,
                               dsamp=F)
#messages


### Resampling and generating the co-clustering graph

print('Doing Resampling step')
#In tutorial: n_resamp = 1000
#In tutorial: p_resamp = 0.75
#Here keeping default parameters.
mcell_coclust_from_graph_resamp(
  coc_id="hsap_coc1000",
  graph_id="hsap_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)


print('Doing co-clustering step')
#In tutorial: alpha = 2
#Here increase alpha to 5, less stringent filtering. So, less metacells.
mcell_mc_from_coclust_balanced(
  coc_id="hsap_coc1000",
  mat_id= "hsap_filt",
  mc_id= "hsap_mc",
  K=20, min_mc_size=20, alpha=5)
# messages

print("Done co-clustering.")


### Removing outlier cells

print('Doing mcell_plot_outlier_heatmap step')
mcell_plot_outlier_heatmap(mc_id="hsap_mc", mat_id = "hsap_filt", T_lfc=3)

print('Doing mcell_mc_split_filt step')
mcell_mc_split_filt(new_mc_id="hsap_mc_f",
                    mc_id="hsap_mc",
                    mat_id="hsap_filt",
                    T_lfc=3, plot_mats=F)
#messages


### Creating heatmaps of metacells and genes

print('Creating heatmaps of metacells and genes')
mc_f<- scdb_mc("hsap_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("hsap_mc_f",mc_f)
mc_f <- scdb_mc("hsap_mc_f")

mcell_gset_from_mc_markers(gset_id="hsap_markers", mc_id="hsap_mc_f")

mcell_mc_plot_marks(mc_id="hsap_mc_f", gset_id="hsap_markers", mat_id="hsap_filt",plot_cells = F)


### Make barplots for each gene of interest (how much expression in each metacell):

print('Make barplots for each gene of interest')
# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Homo_Lukowski2019/10Xdata_Lukowski2019/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Hsap_genes-and-gene-types.txt', sep='') # Check filename is correct!!
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
# As a test use the same config file as Amphimedon tutorial.
download.file("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/config.yaml","config.yaml")

tgconfig::override_params("config.yaml","metacell")
mcell_mc2d_force_knn(mc2d_id="hsap_2dproj",mc_id="hsap_mc_f", graph_id="hsap_graph")

mcell_mc2d_plot(mc2d_id="hsap_2dproj")


## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="hsap_2dproj",  # change name here to appropriate dataset.
                         gene = genename, show_mc_ids = F,
                         show_legend = T, color_cells=T, mat_ds = NULL,
                         zero_sc_v = 0, one_sc_v =1,
                         base_dir = basedir)
  }
}


### Visualizing the metacell confusion matrix

print('Visualization metacell confusion matrix')
mc_hc <- mcell_mc_hclust_confu(mc_id="hsap_mc_f", graph_id="hsap_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="hsap_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="hsap_mc_f",
                        graph_id="hsap_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)



#### Personalising heatmaps for genes of interest

print('Making heatmaps for genes of interest')
## Note: In case of Human, the components of the cililary PRC are well known and sure!
## We will anyway look at potential "rhabdomeric" components.

# Make custom gset for Hsap opsins
# Do only for opsins expressed in PRC:
sets = c('Opsin_OPN1LW','Opsin_OPN1MW3','Opsin_OPN1SW','Opsin_OPN3','Opsin_RHO','Opsin_OPN4')
names(sets) = c('OPN1LW','OPN1MW3','OPN1SW','OPN3','RHO','OPN4')
gs = gset_new_gset(sets, 'Hsap opsin set')
scdb_add_gset("Hsap_opsin_set", gs)

# plot heatmap for opsins using metacells
mcell_mc_plot_marks(mc_id="hsap_mc_f", gset_id="Hsap_opsin_set", mat_id="hsap_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Homo_Lukowski2019/10Xdata_Lukowski2019/figs/Hsap_heatmap_opsins_metacells.png')
# plot heatmap for opsins using cells
mcell_mc_plot_marks(mc_id="hsap_mc_f", gset_id="Hsap_opsin_set", mat_id="hsap_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Homo_Lukowski2019/10Xdata_Lukowski2019/figs/Hsap_heatmap_opsins_cells.png')


# Make custom gset for ALL components of Hsap PRCs (carefully selected, checked with genecards) 
# and candidate rhabdomeric based on reconciliations.
sets = c('CNG_CNGA1','CNG_CNGA2','CNG_CNGA3','CNG_CNGA4','CNG_CNGB1','CNG_CNGB3','GCAP_GUCA1A','GCAP_GUCA1B','GCAP_GUCA1C','GC_GUCY2D','GC_GUCY2F','PDE6_PDE6A','PDE6_PDE6B','PDE6_PDE6C','PDE6_PDE6G','PDE6_PDE6H','Recoverin_RCVRN','RGS9_RGS9','R9AP_RGS9BP','GB5_GNB5','NCKX_SLC24A1','NCKX_SLC24A2','NCKX_SLC24A4','G_alpha_GNAT1','G_alpha_GNAT2','G_beta_GNB1','G_beta_GNB2','G_beta_GNB3','G_beta_GNB4','G_gamma_GNGT1','G_gamma_GNGT2','G_gamma','Arrestin_ARR3','Arrestin_SAG','Calmodulin_CALM1','RK_GRK1','RK_GRK7','Opsin_OPN1LW','Opsin_OPN1MW3','Opsin_OPN1SW','Opsin_OPN3','Opsin_RHO','Opsin_OPN4','Arrestin','Arrestin','G_gamma','RK','RK','RK','RK','RK','INAD','INAD','TRP','TRP','TRP','NINAC','NINAC','PLC','PLC','IP3R','IP3R','IP3R','Ga_R','Ga_R','Ga_R','Actin','Actin','Actin','Actin','DAGL','PKC','PKC','PKC','rdgC','rdgC','CamKII','CamKII','CamKII','CamKII')
names(sets) = c('CNGA1','CNGA2','CNGA3','CNGA4','CNGB1','CNGB3','GUCA1A','GUCA1B','GUCA1C','GUCY2D','GUCY2F','PDE6A','PDE6B','PDE6C','PDE6G','PDE6H','RCVRN','RGS9','RGS9BP','GNB5','SLC24A1','SLC24A2','SLC24A4','GNAT1','GNAT2','GNB1','GNB2','GNB3','GNB4','GNGT1','GNGT2','GNG11','ARR3','SAG','CALM1','GRK1','GRK7','OPN1LW','OPN1MW3','OPN1SW','OPN3','RHO','OPN4','ARRB1','ARRB2','GNG13','GRK2','GRK3','GRK4','GRK5','GRK6','MPDZ','PATJ','TRPC1','TRPC4','TRPC5','MYO3A','MYO3B','PLCB4','PLCB1','ITPR2','ITPR3','ITPR1','GNA14','GNA11','GNAQ','ACTB','ACTBL2','ACTG1','POTEI','DAGLA','PRKCG','PRKCB','PRKCA','PPEF1','PPEF2','CAMK2B','CAMK2G','CAMK2D','CAMK2A')
gs = gset_new_gset(sets, 'Hsap all PRC components set')
scdb_add_gset("Hsap_PRC_set", gs)

# plot heatmap for all hsap PRC components using metacells
mcell_mc_plot_marks(mc_id="hsap_mc_f", gset_id="Hsap_PRC_set", mat_id="hsap_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Homo_Lukowski2019/10Xdata_Lukowski2019/figs/Hsap_heatmap_all_PRC_components_metacells.png')
# plot heatmap for all dmel PRC components using cells
mcell_mc_plot_marks(mc_id="hsap_mc_f", gset_id="Hsap_PRC_set", mat_id="hsap_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Homo_Lukowski2019/10Xdata_Lukowski2019/figs/Hsap_heatmap_all_PRC_components_cells.png')



### Finding genes expressed in metacells:

lfp = log2(mc_f@mc_fp)

print('Writing metacell genes: a file for each metacell')
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


# Extract top 200 for each metacell:
print('Extracting top genes per metacell (lfp)')
x <- 200  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Homo_Lukowski2019/10Xdata_Lukowski2019/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
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
