######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Mmus single cell dataset obtained from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63472
# Corresponding research paper: https://www.sciencedirect.com/science/article/pii/S0092867415005498?via%3Dihub


# set library path
.libPaths('/home/vboxuser/R/x86_64-pc-linux-gnu-library/4.0')

# set Working Directory here:
setwd('/home/vboxuser/Documents/Metacell_R.4.0.0/Mmus_Macosko2015/')
print(getwd())

# load metacell
library("metacell")

# make directories for metacell analysis
if(!dir.exists("mmusp")) dir.create("mmus/")
# initialise the scdb
scdb_init("mmus/", force_reinit=T)

# Dataset is already in .tsv format. Import as tsv:
mcell_import_scmat_tsv(mat_nm="mmus", dset_nm="mmus",
                       fn='GSE63472_P14Retina_merged_digital_expression.txt')
#[1] TRUE

mat = scdb_mat("mmus")
print(dim(mat@mat))
#[1] 24658 49300

print(mat)
#An object of class tgScMat, stat type umi.
#49300 cells by 24658 genes. median cell content 744.

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
mcell_mat_ignore_genes(new_mat_id="mmus_noerccsmt",mat_id="mmus",ig_genes=erccsmt)
mat <- scdb_mat("mmus_noerccsmt")
print(dim(mat@mat))
#[1] 24618 49300
# makes the "mat.mmus_noerccsmt.Rda" file.

# Check umi distribution and choose cut-off. 
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
# Check plot to choose range:
#mcell_plot_umis_per_cell("mmus",min_umis_cutoff = 100)
#mcell_plot_umis_per_cell("mmus",min_umis_cutoff = 200) #min
#mcell_plot_umis_per_cell("mmus",min_umis_cutoff = 225)
#mcell_plot_umis_per_cell("mmus",min_umis_cutoff = 250)
#mcell_plot_umis_per_cell("mmus",min_umis_cutoff = 1000)
#mcell_plot_umis_per_cell("mmus",min_umis_cutoff = 10000)
#mcell_plot_umis_per_cell("mmus",min_umis_cutoff = 12000)
mcell_plot_umis_per_cell("mmus",min_umis_cutoff = 14000) #max

cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>14000))
small_cells <- names(which(cell_sizes<200))
mcell_mat_ignore_cells("mmus_filt","mmus_noerccsmt",ig_cells=c(small_cells,large_cells)) 
mat <- scdb_mat("mmus_filt")
print(dim(mat@mat))
#[1] 24618 49152
# makes the "mat.mmus_filt.Rda" file.


### Selecting gene markers

# calculate statistics.
print('Doing mcell_add_gene_stat step')
#mcell_add_gene_stat(gstat_id="mmus", mat_id="mmus_filt", force=T)
# Load pre-computed results
load("mmus/gstat.mmus.Rda")


# We skip the blacklist filtering, keep rest as tutorial
# Also we adjust the minimum tot UMIs to 200.
print('Doing mcell_gset_filter_multi step')
#mcell_gset_filter_multi(gset_id = "mmus_feats",gstat_id="mmus",T_tot=200,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#Selected 536 markers
# Load pre-computed results
load("mmus/gset.mmus_feats.Rda")


print('Doing mcell_plot_gstats step')
mcell_plot_gstats(gstat_id="mmus", gset_id="mmus_feats")


### Building the balanced cell graph

print('Doing mcell_add_cgraph_from_mat_bknn step')
#mcell_add_cgraph_from_mat_bknn(mat_id="mmus_filt",
#                               gset_id = "mmus_feats",
#                               graph_id="mmus_graph",
#                               K=150,
#                               dsamp=F)
#
# Load pre-computed results
load("mmus/cgraph.mmus_graph.Rda")

### Resampling and generating the co-clustering graph

print('Doing Resampling step')
#In tutorial: n_resamp = 1000
#In tutorial: p_resamp = 0.75
#Here keeping default parameters.
#mcell_coclust_from_graph_resamp(
# coc_id="mmus_coc1000",
#  graph_id="mmus_graph",
#  min_mc_size=20,
#  p_resamp=0.75, n_resamp=1000)
# Load pre-computed results
load("mmus/coclust.mmus_coc1000.Rda")

print('Doing co-clustering step')
#In tutorial: alpha = 2. Here we use aplpha = 4.
#
#mcell_mc_from_coclust_balanced(
#  coc_id="mmus_coc1000",
#  mat_id= "mmus_filt",
#  mc_id= "mmus_mc",
#  K=20, min_mc_size=20, alpha=4)
#
# Load pre-computed results
load("mmus/mc.mmus_mc.Rda")

print("Done co-clustering.")


### Removing outlier cells

print('Doing mcell_plot_outlier_heatmap step')
mcell_plot_outlier_heatmap(mc_id="mmus_mc", mat_id = "mmus_filt", T_lfc=3)

print('Doing mcell_mc_split_filt step')
#mcell_mc_split_filt(new_mc_id="mmus_mc_f",
#                    mc_id="mmus_mc",
#                    mat_id="mmus_filt",
#                    T_lfc=3, plot_mats=F)
#
# Load pre-computed results
load("mmus/mc.mmus_mc_f.Rda")

### Creating heatmaps of metacells and genes

print('Creating heatmaps of metacells and genes')
mc_f<- scdb_mc("mmus_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("mmus_mc_f",mc_f)
mc_f <- scdb_mc("mmus_mc_f")

mcell_gset_from_mc_markers(gset_id="mmus_markers", mc_id="mmus_mc_f")

mcell_mc_plot_marks(mc_id="mmus_mc_f", gset_id="mmus_markers", mat_id="mmus_filt",plot_cells = F)




### Make barplots for each gene of interest (how much expression in each metacell):

# change the working directory for the dataset at hand:
mypath = '/home/vboxuser/Documents/Metacell_R.4.0.0/Mmus_Macosko2015/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Mmus_genes-and-gene-types.txt', sep='') # Check filename is correct!!
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
mcell_mc2d_force_knn(mc2d_id="mmus_2dproj",mc_id="mmus_mc_f", graph_id="mmus_graph")

mcell_mc2d_plot(mc2d_id="mmus_2dproj")


## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="mmus_2dproj",  # change name here to appropriate dataset.
                         gene = genename, show_mc_ids = F, 
                         show_legend = T, color_cells=T, mat_ds = NULL, 
                         zero_sc_v = 0, one_sc_v =1, 
                         base_dir = basedir)  
  }
}



### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="mmus_mc_f", graph_id="mmus_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="mmus_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="mmus_mc_f",
                        graph_id="mmus_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)

print('Done confusion matrix.')




#### Personalising heatmaps for genes of interest


## Make custom gset for Mmus opsins
sets = c('opsin','opsin','opsin','opsin','opsin')
names(sets) = c('OPN3','RHO','OPN1SW','OPN1MW','OPN4') # recovered from reconciliation(s) results.
gs = gset_new_gset(sets, 'Mmus opsin set')
scdb_add_gset("Mmus_opsin_set", gs)

# plot heatmap for opsins using metacells
mcell_mc_plot_marks(mc_id="mmus_mc_f", gset_id="Mmus_opsin_set", mat_id="mmus_filt",plot_cells = F, fig_fn='/home/vboxuser/Documents/Metacell_R.4.0.0/Mmus_Macosko2015/figs/Mmus_heatmap_opsins_metacells.png')
# plot heatmap for opsins using cells
mcell_mc_plot_marks(mc_id="mmus_mc_f", gset_id="Mmus_opsin_set", mat_id="mmus_filt",plot_cells = T, fig_fn='/home/vboxuser/Documents/Metacell_R.4.0.0/Mmus_Macosko2015/figs/Mmus_heatmap_opsins_cells.png')



# Make custom gset for all good Mmus phototr markers.
# Selected based on reconciliations. Kept best.
# Both Rhabdomeric and ciliary.
sets = c('CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','NCKX','NCKX','NCKX','NCKX','GC2','GC2','GC2','GC2','GCAP','GCAP','PDE6_a_and_b','PDE6_a_and_b','PDE6_a_and_b','PDE6_g','PDE6_g','recoverin','RGS9','RGS9','GB5','R9AP','G_aC','G_aC','G_aC','G_aC','G_aC','arrestin','arrestin','arrestin','arrestin','calmodulin','calmodulin','calmodulin','G_b','G_b','G_g','G_g','G_g','G_g','G_g','G_g','G_g','G_g','G_g','G_g','RK','RK','RK','RK','RK','RK','opsin','opsin','opsin','opsin','opsin','actin','actin','CamKII','CamKII','CamKII','DAGL','DAGL','G_aR','G_aR','G_aR','G_aR','INAD','INAD','IP3R','IP3R','IP3R','NINAC','NINAC','PKC','PKC','PKC','PLC','PLC','PLC','PLC','rdgC','rdgC','TRPC','TRPC','TRPC','TRPC','TRPC','TRPC','TRPC')
names(sets) = c('CNGA2','CNGA1','CNGA3','CNGA4','CNGB3','CNGB1','SLC24A2','SLC24A1','SLC24A4','SLC24A3','GUCY2F','GUCY2E','GUCY2D','GUCY2G','GUCA1A','GUCA1B','PDE6C','PDE6B','PDE6A','PDE6H','PDE6G','RCVRN','RGS11','RGS9','GNB5','RGS9BP','GNAT2','GNAT1','GNAT3','GNAI3','GNAI2','ARR3','SAG','ARRB2','ARRB1','CALM1','CALM1','CALM1','GNB4','GNB3','GM5741','GNG12','GNG7','GNG8','GNG4','GNG10','GNG5','GNG13','GNGT2','GNGT1','GRK2','GRK3','GRK6','GRK4','GRK5','GRK1','OPN3','RHO','OPN1SW','OPN1MW','OPN4','ACTB','ACTBL2','CAMK2A','CAMK2G','CAMK2B','DAGLB','DAGLA','GNA14','GNA11','GNAQ','GNA15','MPDZ','PATJ','ITPR1','ITPR3','ITPR2','MYO3A','MYO3B','PRKCA','PRKCB','PRKCG','PLCB4','PLCB2','PLCB3','PLCB1','PPEF2','PPEF1','TRPC5','TRPC4','TRPC1','XNTRPC','TRPC7','TRPC3','TRPC6')
gs = gset_new_gset(sets, 'Mmus phototr markers set')
scdb_add_gset("Mmus_phototr_marks_set", gs)

# plot heatmap for good phototr markers using metacells
mcell_mc_plot_marks(mc_id="mmus_mc_f", gset_id="Mmus_phototr_marks_set", mat_id="mmus_filt",plot_cells = F, fig_fn='/home/vboxuser/Documents/Metacell_R.4.0.0/Mmus_Macosko2015/figs/Mmus_heatmap_phototr_marks_best_metacells.png')
# plot heatmap for good phototr markers using cells
mcell_mc_plot_marks(mc_id="mmus_mc_f", gset_id="Mmus_phototr_marks_set", mat_id="mmus_filt",plot_cells = T, fig_fn='/home/vboxuser/Documents/Metacell_R.4.0.0/Mmus_Macosko2015/figs/Mmus_heatmap_phototr_marks_best_cells.png')





# Extract top 100 for each metacell:
x <- 100  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste(metacells_path,'metacell-',metacell,'_top',x,'.csv',sep=''))
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






