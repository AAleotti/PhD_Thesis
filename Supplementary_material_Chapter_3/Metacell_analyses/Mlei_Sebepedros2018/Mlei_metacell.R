######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Mlei single cell dataset obtained from here: https://www.wisdom.weizmann.ac.il/~arnau/
# Corresponding research paper: https://www.nature.com/articles/s41559-018-0575-6


# set library path
.libPaths('/data/evassvis/software/R/4.0')

# set Working Directory here:
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018')

library(metacell)

# make directories
if(!dir.exists("mlei")) dir.create("mlei/")
scdb_init("mlei/", force_reinit=T)

# import UMIs for Mnemiopsis
mcell_import_multi_mars("mlei", "MARS_Batches.txt", base_dir="/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/raw_umi_tables/")

mat = scdb_mat("mlei")
print(dim(mat@mat))
#[1] 21622  6144

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

### Exploring and filtering the UMI matrix

erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="mlei_noercc",mat_id="mlei",ig_genes=erccs)
mat <- scdb_mat("mlei_noercc")
print(dim(mat@mat))
#[1] 21530  6144

# Check umi distribution and choose cut-off. 
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
# Do a few tests with Mlei dataset:
mcell_plot_umis_per_cell("mlei",min_umis_cutoff = 100)
mcell_plot_umis_per_cell("mlei",min_umis_cutoff = 12000)

# For Mnemiopsis use lower cut-off 100 and higher cut-off 12000
cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>12000))
small_cells <- names(which(cell_sizes<100))
mcell_mat_ignore_cells("mlei_filt","mlei_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("mlei_filt")
print(dim(mat@mat))
#[1] 21530  5492

### Selecting gene markers

mcell_add_gene_stat(gstat_id="mlei", mat_id="mlei_filt", force=T)

# we skip the blacklist filtering, keep rest as tutorial
# Also we ajust the minimum tot UMIs to 100.
mcell_gset_filter_multi(gset_id = "mlei_feats",gstat_id="mlei",T_tot=100,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#> Selected 959 markers

mcell_plot_gstats(gstat_id="mlei", gset_id="mlei_feats")

### Building the balanced cell graph

mcell_add_cgraph_from_mat_bknn(mat_id="mlei_filt",
                               gset_id = "mlei_feats",
                               graph_id="mlei_graph",
                               K=150,
                               dsamp=F)

### Resampling and generating the co-clustering graph

mcell_coclust_from_graph_resamp(
  coc_id="mlei_coc1000",
  graph_id="mlei_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)

mcell_mc_from_coclust_balanced(
  coc_id="mlei_coc1000",
  mat_id= "mlei_filt",
  mc_id= "mlei_mc",
  K=20, min_mc_size=20, alpha=2)
# messages:
#filtered 1314362 left with 380500 based on co-cluster imbalance
#building metacell object, #mc 51
#add batch counts
#compute footprints
#compute absolute ps
#compute coverage ps
#reordering metacells by hclust and most variable two markers
#reorder on ML07216a vs ML21632a

### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="mlei_mc", mat_id = "mlei_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="mlei_mc_f",
                    mc_id="mlei_mc",
                    mat_id="mlei_filt",
                    T_lfc=3, plot_mats=F)

### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("mlei_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("mlei_mc_f",mc_f)
mc_f <- scdb_mc("mlei_mc_f")

mcell_gset_from_mc_markers(gset_id="mlei_markers", mc_id="mlei_mc_f")

mcell_mc_plot_marks(mc_id="mlei_mc_f", gset_id="mlei_markers", mat_id="mlei_filt",plot_cells = F)



## Make barplots for each gene of interest:

# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Mlei_genes-and-gene-types.txt', sep='') # check filename is correct!
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

# As a test use the same config file as Amphimedon tutorial. (but then verify if should be specific to Mlei)
# Or could it make sense to skip the configuration file?
download.file("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/config.yaml","config.yaml")

tgconfig::override_params("config.yaml","metacell")
mcell_mc2d_force_knn(mc2d_id="mlei_2dproj",mc_id="mlei_mc_f", graph_id="mlei_graph")

mcell_mc2d_plot(mc2d_id="mlei_2dproj")

## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="mlei_2dproj", 
                         gene = genename, show_mc_ids = F, 
                         show_legend = T, color_cells=T, mat_ds = NULL, 
                         zero_sc_v = 0, one_sc_v =1, 
                         base_dir = basedir)  
  }
}


### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="mlei_mc_f", graph_id="mlei_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="mlei_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="mlei_mc_f",
                        graph_id="mlei_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)




#### Personalising heatmap


# Make custom gset for opsins
sets = c('opsin','opsin')
names(sets) = c('ML12047a', 'ML13055a') # two opsins recovered from reconciliation(s) results.
gs = gset_new_gset(sets, 'Mlei opsin set')
scdb_add_gset("Mlei_opsin_set", gs)

# check how gset looks
print(gs)

# plot heatmap for opsins using metacells
mcell_mc_plot_marks(mc_id="mlei_mc_f", gset_id="Mlei_opsin_set", mat_id="mlei_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/figs/Mlei_heatmap_opsins_metacells.png')
# plot heatmap for opsins using cells
mcell_mc_plot_marks(mc_id="mlei_mc_f", gset_id="Mlei_opsin_set", mat_id="mlei_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/figs/Mlei_heatmap_opsins_cells.png')


# Making gset for all phototr genes:
sets = c('CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','NCKX','G_aC','G_aC','opsin','opsin','G_b','G_b','G_b','calmodulin','calmodulin','calmodulin','arrestin','GRK','GRK','GRK','actin','actin','actin','actin','actin','CamkII','DAGL','G_aR','G_aR','G_aR','INAD','IP3R','PKC','PLC','TRP_TRPL','TRP_TRPL','TRP_TRPL')
names(sets) = c('ML054419a','ML08605a','ML30617a','ML261714a','ML156513a','ML156514a','ML12047a','ML13055a','ML10611a','ML093041a','ML02234a','ML149632a','ML149631a','ML104636a','ML047926a','ML02651a','ML009124a','ML04904a','ML174735a','ML20265a','ML322211a','ML35651a','ML35935a','ML35309a','ML12863a','ML003265a','ML009153a','ML154511a','ML14122a','ML07011a','ML13931a','ML04921a','ML03701a','ML181719a','ML234550a')
gs = gset_new_gset(sets, 'Mlei all phototr')
scdb_add_gset("Mlei_all_phototr", gs)

# plot heatmap for all phototr using metacells
mcell_mc_plot_marks(mc_id="mlei_mc_f", gset_id="Mlei_all_phototr", mat_id="mlei_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/figs/Mlei_heatmap_phototr_metacells.png')
# plot heatmap for all phototr using cells
mcell_mc_plot_marks(mc_id="mlei_mc_f", gset_id="Mlei_all_phototr", mat_id="mlei_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/figs/Mlei_heatmap_phototr_cells.png')


# Make custom gset for subset of photoreceptor genes.
sets = c('CNG_a_and_b','NCKX','opsin','opsin','PKC','PLC')
names(sets) = c('ML08605a','ML261714a','ML12047a','ML13055a','ML13931a','ML04921a') 
gs = gset_new_gset(sets, 'Mlei subset set')
scdb_add_gset("Mlei_subset_set", gs)

# plot heatmap for subset using metacells
mcell_mc_plot_marks(mc_id="mlei_mc_f", gset_id="Mlei_subset_set", mat_id="mlei_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/figs/Mlei_heatmap_subset_metacells.png')
# plot heatmap for subset using cells
mcell_mc_plot_marks(mc_id="mlei_mc_f", gset_id="Mlei_subset_set", mat_id="mlei_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/figs/Mlei_heatmap_subset_cells.png')



### Finding genes expressed in a metacells:


# Writing the metacells genes, a single file per metacell in the metacells_path folder
for (metacell in colnames(lfp)){
  lines = ''
  for (gene in rownames(lfp)){
    newline = paste(gene,lfp[gene,metacell],sep = ',')
    lines <- paste(lines,newline,'\n',sep='')
  }
  write(lines, file=paste(metacells_path,'metacell-',metacell,'.csv',sep=''))
}
print(paste('Done writing the metacell genes in files. Please see', metacells_path, '.'))


# If you want to output results for only a subset of metacells:
# e.g. if I am interested in metacell 39:
#for (metacell in list(39)){
#  lines = ''
#  for (gene in rownames(lfp)){
#    newline = paste(gene,lfp[gene,metacell],sep = ',')
#    lines <- paste(lines,newline,'\n',sep='')
#  }
#  write(lines, file=paste(metacells_path,'metacell-',metacell,'.csv',sep=''))
#}




# Extract top 100 for each metacell:
x <- 100  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Mnemiopsis_sebepedros2018/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
}



# Writing the metacells genes with lfp > 0, a single file per metacell in the metacells_path folder
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

