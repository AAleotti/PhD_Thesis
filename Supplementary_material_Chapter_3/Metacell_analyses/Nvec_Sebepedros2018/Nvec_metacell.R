######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Nvec single cell dataset obtained from here: https://www.wisdom.weizmann.ac.il/~arnau/
# Corresponding research paper: https://www.cell.com/cell/fulltext/S0092-8674(18)30596-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418305968%3Fshowall%3Dtrue


# set library path
.libPaths('/data/evassvis/software/R/4.0')

# set Working Directory here:
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella')

library(metacell)

# make directories
if(!dir.exists("nvec")) dir.create("nvec/")
scdb_init("nvec/", force_reinit=T)

# import UMIs for Adults only
mcell_import_multi_mars("nvec", "MARS_Batches_Adult_Only.txt", base_dir="/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella/UMI_tables")

mat = scdb_mat("nvec")
print(dim(mat@mat))
#[1] 32390 14208

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

### Exploring and filtering the UMI matrix

erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="nvec_noercc",mat_id="nvec",ig_genes=erccs)
mat <- scdb_mat("nvec_noercc")

# in tutorial cut-off: 200
# mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 200)

# nematostella paper cut-off: 100
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 100)
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 50)
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 80)
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 70)
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 60)
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 11000)
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 9000)
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 8000)
#mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 7000)
mcell_plot_umis_per_cell("nvec",min_umis_cutoff = 6000)

# We choose cut-off
cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>6000))
small_cells <- names(which(cell_sizes<60))
mcell_mat_ignore_cells("nvec_filt","nvec_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("nvec_filt")
print(dim(mat@mat))
#[1] 32298 12964 #how big the matrix is now.


### Selecting gene markers

mcell_add_gene_stat(gstat_id="nvec", mat_id="nvec_filt", force=T)

# we skip the blacklist filtering, keep rest as tutorial
mcell_gset_filter_multi(gset_id = "nvec_feats",gstat_id="nvec",T_tot=60,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#> Selected 338 markers

mcell_plot_gstats(gstat_id="nvec", gset_id="nvec_feats")

### Building the balanced cell graph

mcell_add_cgraph_from_mat_bknn(mat_id="nvec_filt",
                               gset_id = "nvec_feats",
                               graph_id="nvec_graph",
                               K=150,
                               dsamp=F)


### Resampling and generating the co-clustering graph

print('Resampling and co-clustering')
mcell_coclust_from_graph_resamp(
  coc_id="nvec_coc1000",
  graph_id="nvec_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)
print('Done resampling.')

# Try less stringent: alpha = 3 (instead of default 2).
mcell_mc_from_coclust_balanced(
  coc_id="nvec_coc1000",
  mat_id= "nvec_filt",
  mc_id= "nvec_mc",
  K=20, min_mc_size=20, alpha=3)
# messages:

print('Done co-clustering graph.')

### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="nvec_mc", mat_id = "nvec_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="nvec_mc_f",
                    mc_id="nvec_mc",
                    mat_id="nvec_filt",
                    T_lfc=3, plot_mats=F)



### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("nvec_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("nvec_mc_f",mc_f)
mc_f <- scdb_mc("nvec_mc_f")

mcell_gset_from_mc_markers(gset_id="nvec_markers", mc_id="nvec_mc_f")

mcell_mc_plot_marks(mc_id="nvec_mc_f", gset_id="nvec_markers", mat_id="nvec_filt",plot_cells = F)



## Make barplots for each gene of interest:

# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Nvec_genes-and-gene-types_best.txt', sep='') # check filename is correct!
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
mcell_mc2d_force_knn(mc2d_id="nvec_2dproj",mc_id="nvec_mc_f", graph_id="nvec_graph")

mcell_mc2d_plot(mc2d_id="nvec_2dproj")

## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="nvec_2dproj", 
                         gene = genename, show_mc_ids = F,
                         show_legend = T, color_cells=T, mat_ds = NULL,
                         zero_sc_v = 0, one_sc_v =1,
                         base_dir = basedir)
  }
}



### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="nvec_mc_f", graph_id="nvec_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="nvec_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="nvec_mc_f",
                        graph_id="nvec_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)


#### Personalising heatmap


# making gene set for Nematostella opsins that were recovered from reconciliation(s) results
sets = c('opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin')
names(sets) = c('v1g13116','v1g33918','v1g208768','v1g214774','v1g214775','v1g214772','v1g214773','v1g199627','v1g102435','v1g123451','v1g94740','v1g202741','v1g131013','v1g156694','v1g85309','v1g130042','v1g219988','v1g96290','v1g123690','v1g197433','v1g136537','v1g95791','v1g25025','v1g24997')
gs = gset_new_gset(sets, 'Nvec all opsins')
scdb_add_gset("Nvec_all_opsins", gs)

# opsin heatmap plot with metacells
mcell_mc_plot_marks(mc_id="nvec_mc_f", gset_id="Nvec_all_opsins", mat_id="nvec_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella/figs/nvec_heatmap_all_opsins_metacells.png')

# opsin heatmap plot with cells
mcell_mc_plot_marks(mc_id="nvec_mc_f", gset_id="Nvec_all_opsins", mat_id="nvec_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella/figs/nvec_heatmap_all_opsins_cells.png')


# Making gene set with all Nematostella phototr genes
sets = c('CNG','CNG','CNG','CNG','CNG','NCKX','NCKX','NCKX','NCKX','NCKX','NCKX','PDE6_A_B','RGS9','GB5','GC2','GC2','GC2','G_alpha_i_o','G_alpha_i_o','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','opsin','G_beta','G_gamma','G_gamma','calmodulin','calmodulin','calmodulin','arrestin','GRK','GRK','Actin','Actin','Actin','Actin','Actin','CamKII','CamKII','CamKII','DAGL','DAGL','DAGL','G_alpha_q','INAD','IP3R','IP3R','NINAC','PKC','PLC','PLC','rdgC','TRPC','TRPC','TRPC','TRPC','TRPC')
names(sets) = c('v1g122927','v1g125284','v1g20862','v1g21809','v1g81021','v1g134893','v1g198113','v1g198469','v1g206699','v1g83881','v1g85050','v1g164927','v1g230662','v1g244849','v1g120467','v1g196879','v1g93709','v1g177378','v1g237532','v1g13116','v1g33918','v1g208768','v1g214774','v1g214775','v1g214772','v1g214773','v1g199627','v1g102435','v1g123451','v1g94740','v1g202741','v1g131013','v1g156694','v1g85309','v1g130042','v1g219988','v1g96290','v1g123690','v1g197433','v1g136537','v1g95791','v1g25025','v1g24997','v1g105798','v1g238752','v1g212694','v1g88970','v1g188289','v1g239788','v1g97737','v1g179295','v1g246204','v1g159077','v1g93634','v1g93180','v1g157748','v1g186733','v1g153185','v1g157808','v1g220432','v1g202361','v1g60191','v1g93878','v1g242737','v1g216282','v1g143444','v1g244707','v1g164057','v1g93932','v1g101739','v1g174','v1g239552','v1g120174','v1g138017','v1g31537','v1g5948','v1g92648')
gs = gset_new_gset(sets, 'Nvec all phototr genes')
scdb_add_gset("Nvec_all_phototr",gs)

# all phototr heatmap plot with metacells
mcell_mc_plot_marks(mc_id="nvec_mc_f", gset_id="Nvec_all_phototr", mat_id="nvec_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella/figs/nvec_heatmap_all_phototr_metacells.png')

# all phototr heatmap plot with cells
mcell_mc_plot_marks(mc_id="nvec_mc_f", gset_id="Nvec_all_phototr", mat_id="nvec_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella/figs/nvec_heatmap_all_phototr_cells.png')





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
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
}

# Extract top 200 for each metacell:
x <- 200  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Nematostella_sebepedros2018/Arnau_datasets_Nematostella/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
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
