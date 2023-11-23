######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Spur single cell dataset obtained from here: https://datadryad.org/stash/dataset/doi:10.5061/dryad.n5tb2rbvz
# Corresponding research paper: https://elifesciences.org/articles/70416


# set library path
.libPaths('/data/evassvis/software/R/4.0')

# set Working Directory here:
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Spur_Paganos_et_al_2021/Spur_metacell_analysis_multi-batch')

library(metacell)

# make directories
if(!dir.exists("spur")) dir.create("spur/")
scdb_init("spur/", force_reinit=T)

# import Spur single cell data sets
# Paganos et al 2021 data is 10X, so load 10X multi-batch data:
mcell_import_multi_scmat_10x("spur", "spur_dataset_table_fn.txt", base_dir="/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Spur_Paganos_et_al_2021/Spur_metacell_analysis_multi-batch")

mat = scdb_mat("spur")
print(dim(mat@mat))
#[1] 21090 29130 [29130 cells corresponds to number reported by authors in paper]

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")


# check if code format for genes is correct.
WHLcode <- rownames(mat@mat)[grepl("WHL22.92304",rownames(mat@mat))]
#Check if gene is found:
WHLcode
#[1] "WHL22.92304"  #correct!


# In some cases also the SPU format is found:
SPUcode <- rownames(mat@mat)[grepl("SPU_006683",rownames(mat@mat))]
#Check MT genes found:
SPUcode


### Exploring and filtering the UMI matrix

erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="spur_noercc",mat_id="spur",ig_genes=erccs)
mat <- scdb_mat("spur_noercc")

# Check umi distribution and choose cut-off. 
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
mcell_plot_umis_per_cell("spur",min_umis_cutoff = 180) 
mcell_plot_umis_per_cell("spur",min_umis_cutoff = 7500)

# Extract from paper: 
# "Genes that are transcribed in less than three cells and cells that have 
# less than a minimum of 200 transcribed genes were excluded from the analysis."

# For Spur use lower cut-off 180 and higher cut-off 7500
cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>7500))
small_cells <- names(which(cell_sizes<180))
mcell_mat_ignore_cells("spur_filt","spur_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("spur_filt")
print(dim(mat@mat))
# [1] 21090 29097


### Selecting gene markers

# calculate statistics.
mcell_add_gene_stat(gstat_id="spur", mat_id="spur_filt", force=T)

# we skip the blacklist filtering, keep rest as tutorial
# Also we adjust the minimum tot UMIs to 180.
mcell_gset_filter_multi(gset_id = "spur_feats",gstat_id="spur",T_tot=180,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#Selected 305 markers

mcell_plot_gstats(gstat_id="spur", gset_id="spur_feats")


### Building the balanced cell graph

mcell_add_cgraph_from_mat_bknn(mat_id="spur_filt",
                               gset_id = "spur_feats",
                               graph_id="spur_graph",
                               K=150,
                               dsamp=F)

### Resampling and generating the co-clustering graph

mcell_coclust_from_graph_resamp(
  coc_id="spur_coc1000",
  graph_id="spur_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)
print("Resampling done.")

# To decrease number of metacells, increase alpha from default of 2 to 5.
mcell_mc_from_coclust_balanced(
  coc_id="spur_coc1000",
  mat_id= "spur_filt",
  mc_id= "spur_mc",
  K=20, min_mc_size=20, alpha=5)
# messages:
print("Co-clustering graph done.")


### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="spur_mc", mat_id = "spur_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="spur_mc_f",
                    mc_id="spur_mc",
                    mat_id="spur_filt",
                    T_lfc=3, plot_mats=F)
#messages:
print("Finished removing outliers.")


### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("spur_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("spur_mc_f",mc_f)
mc_f <- scdb_mc("spur_mc_f")

mcell_gset_from_mc_markers(gset_id="spur_markers", mc_id="spur_mc_f")

mcell_mc_plot_marks(mc_id="spur_mc_f", gset_id="spur_markers", mat_id="spur_filt",plot_cells = F)


### Make barplots for each gene of interest (how much expression in each metacell):

# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Spur_Paganos_et_al_2021/Spur_metacell_analysis_multi-batch/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Spur_genes-and-gene-types.txt', sep='') # check filename is correct!
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

# Use the same config file as Amphimedon tutorial.
download.file("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/config.yaml","config.yaml")

tgconfig::override_params("config.yaml","metacell")
mcell_mc2d_force_knn(mc2d_id="spur_2dproj",mc_id="spur_mc_f", graph_id="spur_graph")

mcell_mc2d_plot(mc2d_id="spur_2dproj")



## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene


for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="spur_2dproj", 
                         gene = genename, show_mc_ids = F, 
                         show_legend = T, color_cells=T, mat_ds = NULL, 
                         zero_sc_v = 0, one_sc_v =1, 
                         base_dir = basedir)  
  }
}



### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="spur_mc_f", graph_id="spur_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="spur_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="spur_mc_f",
                        graph_id="spur_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)


#### Personalising heatmaps for genes of interest


# Make custom gset for opsins
sets = c('opsin','opsin','opsin')
names(sets) = c('WHL22.272775','WHL22.338995','WHL22.290080') # opsins chosen from reconciliation(s) results.
gs = gset_new_gset(sets, 'Spur opsin set')
scdb_add_gset("Spur_opsin_set", gs)


# plot heatmap for opsins using metacells
mcell_mc_plot_marks(mc_id="spur_mc_f", gset_id="Spur_opsin_set", mat_id="spur_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Spur_Paganos_et_al_2021/Spur_metacell_analysis_multi-batch/figs/Spur_heatmap_opsins_metacells.png')
# plot heatmap for opsins using cells
mcell_mc_plot_marks(mc_id="spur_mc_f", gset_id="Spur_opsin_set", mat_id="spur_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Spur_Paganos_et_al_2021/Spur_metacell_analysis_multi-batch/figs/Spur_heatmap_opsins_cells.png')


# Make custom gset for all phototr genes
sets = c('CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','CNG_a_and_b','NCKX','NCKX','PDE6_a_and_b','R9AP','GB5','RGS9','GC_2DEF','G_aC','G_aC','opsin','opsin','opsin','G_b','G_g','G_g','calmodulin','calmodulin','calmodulin','calmodulin','arrestin','RK','RK','actin','actin','CamkII','DAGL','DAGL','DAGL','INAD','INAD','IP3R','IP3R','NINAC','NINAC','NINAC','NINAC','PKC','PKC','PLC','PLC','rdgC','TRP_TRPL','TRP_TRPL','TRP_TRPL','TRP_TRPL')
names(sets) = c('WHL22.92304','WHL22.751707','SPU_006683','WHL22.652873','WHL22.133623','WHL22.133648','WHL22.754407','WHL22.668748','WHL22.554631','WHL22.692880','WHL22.51293','WHL22.488878','WHL22.75475','WHL22.731493','WHL22.762062','WHL22.272775','WHL22.338995','WHL22.290080','WHL22.101602','WHL22.557017','WHL22.253891','WHL22.521535','WHL22.521547','WHL22.4013','WHL22.4059','WHL22.709217','WHL22.64904','WHL22.421157','WHL22.408109','WHL22.608502','WHL22.677837','WHL22.458877','WHL22.324483','WHL22.697968','WHL22.761541','WHL22.761563','WHL22.635333','WHL22.635570','WHL22.309651','WHL22.72626','WHL22.550412','WHL22.324194','WHL22.343472','WHL22.343484','WHL22.499001','WHL22.169178','WHL22.469735','WHL22.3956','WHL22.34179','WHL22.316252','WHL22.77417') 
gs = gset_new_gset(sets, 'Spur phototr set')
scdb_add_gset("Spur_phototr_set", gs)


# plot heatmap for all phototr using metacells
mcell_mc_plot_marks(mc_id="spur_mc_f", gset_id="Spur_phototr_set", mat_id="spur_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Spur_Paganos_et_al_2021/Spur_metacell_analysis_multi-batch/figs/Spur_heatmap_phototr_metacells.png')
# plot heatmap for all phototr using cells
mcell_mc_plot_marks(mc_id="spur_mc_f", gset_id="Spur_phototr_set", mat_id="spur_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Spur_Paganos_et_al_2021/Spur_metacell_analysis_multi-batch/figs/Spur_heatmap_phototr_cells.png')



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
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Spur_Paganos_et_al_2021/Spur_metacell_analysis_multi-batch/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
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

