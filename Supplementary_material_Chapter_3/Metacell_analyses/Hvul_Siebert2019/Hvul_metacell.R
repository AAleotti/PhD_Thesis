######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Hvul single cell dataset obtained from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121617
# Corresponding research paper: https://www.science.org/doi/10.1126/science.aav9314


# Set library path for R on ALICE
.libPaths('/data/evassvis/software/R/4.0')


#Check current directory
getwd()

# set path
setwd('/scratch/evassvis/aa1176/Hydra/')
print(getwd())

# load metacell
library(metacell)


# make directories for metacell pipeline
if(!dir.exists("hvul")) dir.create("hvul")
# initialise the scdb (necessary for metacell pipeline)
scdb_init("hvul/", force_reinit=T)


# Import the table into metacell (has to be imported everytime)
mcell_import_scmat_tsv(mat_nm="hvul", dset_nm="hvul",
                       fn='/scratch/evassvis/aa1176/Hydra/Hydra_DS_transcriptome_UMICounts.txt')



mat = scdb_mat("hvul")
# the current size of the matrix is:
print(dim(mat@mat))
#[1] 37114 27992

# summary info:
print(mat)
#An object of class tgScMat, stat type umi.
#27992 cells by 37114 genes. median cell content 5622.


#### metacell pipeline

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")



### Exploring and filtering the UMI matrix

# remove any potential ercc standards
erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="hvul_noercc",mat_id="hvul",ig_genes=erccs)
mat <- scdb_mat("hvul_noercc")
print(dim(mat@mat))
#[1] 37107 27992


# Check umi distribution and choose cut-off.
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
# Check plot to choose range:
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 100)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 250)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 300)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 350) # choose as min.
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 400)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 500)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 1000)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 10000)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 15000)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 20000)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 30000)
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 35000) # choose as max.
#mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 40000)
mcell_plot_umis_per_cell("hvul",min_umis_cutoff = 50000)


cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>35000))
small_cells <- names(which(cell_sizes<350))
mcell_mat_ignore_cells("hvul_filt","hvul_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("hvul_filt")
print(dim(mat@mat))
#[1] 37107 27026


### Selecting gene markers

# calculate statistics.
mcell_add_gene_stat(gstat_id="hvul", mat_id="hvul_filt", force=T)
print("Done calculate gene statistics.")

# we skip the blacklist filtering, keep rest as tutorial
# Also we adjust the minimum tot UMIs to 350.
mcell_gset_filter_multi(gset_id = "hvul_feats",gstat_id="hvul",T_tot=350,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#

mcell_plot_gstats(gstat_id="hvul", gset_id="hvul_feats")


### Building the balanced cell graph

print("Doing step build balanced cell grapph.")
mcell_add_cgraph_from_mat_bknn(mat_id="hvul_filt",
                               gset_id = "hvul_feats",
                               graph_id="hvul_graph",
                               K=150,
                               dsamp=F)



### Resampling and generating the co-clustering graph

mcell_coclust_from_graph_resamp(
  coc_id="hvul_coc1000",
  graph_id="hvul_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)

print('Done resampling.')


# alpha default is 2. Try with alpha =3.
mcell_mc_from_coclust_balanced(
  coc_id="hvul_coc1000",
  mat_id= "hvul_filt",
  mc_id= "hvul_mc",
  K=20, min_mc_size=20, alpha=3)
# messages

print('Done coclustering.')


### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="hvul_mc", mat_id = "hvul_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="hvul_mc_f",
                    mc_id="hvul_mc",
                    mat_id="hvul_filt",
                    T_lfc=3, plot_mats=F)
#messages



### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("hvul_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("hvul_mc_f",mc_f)
mc_f <- scdb_mc("hvul_mc_f")

mcell_gset_from_mc_markers(gset_id="hvul_markers", mc_id="hvul_mc_f")

mcell_mc_plot_marks(mc_id="hvul_mc_f", gset_id="hvul_markers", mat_id="hvul_filt",plot_cells = F)



### Make barplots for each gene of interest (how much expression in each metacell):

# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/Hydra/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Hvul_genes-and-gene-types.txt', sep='') # Check filename is correct!!
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
mcell_mc2d_force_knn(mc2d_id="hvul_2dproj",mc_id="hvul_mc_f", graph_id="hvul_graph")

mcell_mc2d_plot(mc2d_id="hvul_2dproj")



## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="hvul_2dproj",  # change name here to appropriate dataset.
                         gene = genename, show_mc_ids = F,
                         show_legend = T, color_cells=T, mat_ds = NULL,
                         zero_sc_v = 0, one_sc_v =1,
                         base_dir = basedir)
  }
}



### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="hvul_mc_f", graph_id="hvul_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="hvul_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="hvul_mc_f",
                        graph_id="hvul_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)

print('Done confusion matrix.')



#### Personalising heatmaps for genes of interest


## Make custom gset for Hvul opsins
sets = c('opsin','opsin','opsin','opsin','opsin','opsin','opsin')
names(sets) = c('t24564aep|OPSR_CAPHI','t10575aep|OPSR_BOVIN','t24564aep|OPSR_CAPHI','t26793aep|OPSO_RUTRU','t24044aep|OPSP_ICTPU','t27688aep|OPSP_COLLI','t29959aep|OPN4B_XENLA') # recovered from reconciliation(s) results.
gs = gset_new_gset(sets, 'Hvul opsin set')
scdb_add_gset("Hvul_opsin_set", gs)

# plot heatmap for opsins using metacells
mcell_mc_plot_marks(mc_id="hvul_mc_f", gset_id="Hvul_opsin_set", mat_id="hvul_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/Hydra/figs/Hvul_heatmap_opsins_metacells.png')
# plot heatmap for opsins using cells
mcell_mc_plot_marks(mc_id="hvul_mc_f", gset_id="Hvul_opsin_set", mat_id="hvul_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/Hydra/figs/Hvul_heatmap_opsins_cells.png')



# Make custom gset for all good Hvul phototr markers.
# Selected based on reconciliations. Kept best.
# Both Rhabdomeric and ciliary.
sets = c('CNG_a_and_b','G_aC','G_aC','GB5','GC2','GC2','NCKX','NCKX','PDE6_a_and_b','RGS9','arrestin','calmodulin','calmodulin','G_b','G_g','RK','RK','opsin','opsin','opsin','opsin','opsin','opsin','opsin','actin','CamKII','DAGL','G_aR','INAD','IP3R','NINAC','PKC','PLC','PLC','TRP_TRPL')
names(sets) = c('t27655aep|CNGA3_HUMAN','t12031aep|GNAI_PATPE','t15699aep|GNAO_BOVIN','t32166aep|GBB5_PONAB','t26339aep|ANPRA_HUMAN','t31398aep|ANPRB_HUMAN','t27685aep|NCKX5_DANRE','t1695aep|NCKX5_HUMAN','t38455aep|PDE11_DROME','t20170aep|RGS7_HUMAN','t14420aep|ARRB1_MACFA','t28436aep|CALM_METSE','t25988aep|CALM_HALOK','t19432aep|GBB_PINFU','t26652aep|GBG7_HUMAN','t17074aep|ARBK2_BOVIN','t13671aep|GRK5_BOVIN','t24564aep|OPSR_CAPHI','t10575aep|OPSR_BOVIN','t24564aep|OPSR_CAPHI','t26793aep|OPSO_RUTRU','t24044aep|OPSP_ICTPU','t27688aep|OPSP_COLLI','t29959aep|OPN4B_XENLA','t11116aep|ACT_HYDVU','t11598aep|KCC2A_DROME','t30622aep|DGLA_RAT','t612aep|GNAQ_MIZYE','t22223aep|MPDZ_HUMAN','t21920aep|ITPR1_BOVIN','t27676aep|MYO3A_HUMAN','t31809aep|KPCA_RAT','t29781aep|PIP1_DROME','t25874aep|PLCB4_RAT','t31083aep|TRPM3_HUMAN')
gs = gset_new_gset(sets, 'Hvul phototr markers set')
scdb_add_gset("Hvul_phototr_marks_set", gs)

# plot heatmap for good phototr markers using metacells
mcell_mc_plot_marks(mc_id="hvul_mc_f", gset_id="Hvul_phototr_marks_set", mat_id="hvul_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/Hydra/figs/Hvul_heatmap_phototr_marks_best_metacells.png')
# plot heatmap for good phototr markers using cells
mcell_mc_plot_marks(mc_id="hvul_mc_f", gset_id="Hvul_phototr_marks_set", mat_id="hvul_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/Hydra/figs/Hvul_heatmap_phototr_marks_best_cells.png')




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
