######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Aque single cell dataset obtained from here: https://www.wisdom.weizmann.ac.il/~arnau/
# Corresponding research paper: https://www.nature.com/articles/s41559-018-0575-6


# set library path
.libPaths('/data/evassvis/software/R/4.0')

# set Working Directory here:
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Amphimedon_sebepedros2018/Amphimedon')

library(metacell)

# make directories
if(!dir.exists("aque")) dir.create("aque/")
scdb_init("aque/", force_reinit=T)

# import UMIs for Adults
mcell_import_multi_mars("aque", "MARS_Batches.txt", base_dir="/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Amphimedon_sebepedros2018/Amphimedon/raw_umi_tables/Adult")

mat = scdb_mat("aque")
print(dim(mat@mat))
#[1] 49905  4992

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

### Exploring and filtering the UMI matrix

erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="aque_noercc",mat_id="aque",ig_genes=erccs)
mat <- scdb_mat("aque_noercc")

print(dim(mat@mat))
#[1] 49813  4992

# check the umis distribution
#mcell_plot_umis_per_cell("aque",min_umis_cutoff = 150)
#mcell_plot_umis_per_cell("aque",min_umis_cutoff = 200)
#mcell_plot_umis_per_cell("aque",min_umis_cutoff = 12000)
#mcell_plot_umis_per_cell("aque",min_umis_cutoff = 13000)
#mcell_plot_umis_per_cell("aque",min_umis_cutoff = 14000)
#mcell_plot_umis_per_cell("aque",min_umis_cutoff = 15000)
#mcell_plot_umis_per_cell("aque",min_umis_cutoff = 16000)

# We choose cut-off
cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>16000))
small_cells <- names(which(cell_sizes<150))
mcell_mat_ignore_cells("aque_filt","aque_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("aque_filt")
print(dim(mat@mat))
#[1] 49813  4883

### Selecting gene markers

mcell_add_gene_stat(gstat_id="aque", mat_id="aque_filt", force=T)

# This is same dataset as tutorial, so we can do exact same filtering including the blacklist filtering
bl<- scan("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/bl_genes",what="")
mcell_gset_filter_multi(gset_id = "aque_feats",gstat_id="aque",T_tot=200,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T,blacklist=bl)
#> Selected 1006 markers

mcell_plot_gstats(gstat_id="aque", gset_id="aque_feats")

### Building the balanced cell graph

mcell_add_cgraph_from_mat_bknn(mat_id="aque_filt",
                               gset_id = "aque_feats",
                               graph_id="aque_graph",
                               K=150,
                               dsamp=F)
#will build balanced knn graph on 4883 cells and 1037 genes, this can be a bit heavy for >20,000 cells
#sim graph is missing 10 nodes, out of 4883


### Resampling and generating the co-clustering graph

mcell_coclust_from_graph_resamp(
  coc_id="aque_coc1000",
  graph_id="aque_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)

mcell_mc_from_coclust_balanced(
  coc_id="aque_coc1000",
  mat_id= "aque_filt",
  mc_id= "aque_mc",
  K=20, min_mc_size=20, alpha=2)
# messages:
#filtered 1264488 left with 374744 based on co-cluster imbalance
#building metacell object, #mc 47
#add batch counts
#compute footprints
#compute absolute ps
#compute coverage ps
#reordering metacells by hclust and most variable two markers
#reorder on Contig13211_51996 vs Aqu2.1.43479_001

### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="aque_mc", mat_id = "aque_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="aque_mc_f",
                    mc_id="aque_mc",
                    mat_id="aque_filt",
                    T_lfc=3, plot_mats=F)
# messages:
#starting split outliers 
#splitting metacell 26
#splitting metacell 44
#add batch counts
#compute footprints
#compute absolute ps
#compute coverage ps

### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("aque_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("aque_mc_f",mc_f)
mc_f <- scdb_mc("aque_mc_f")

mcell_gset_from_mc_markers(gset_id="aque_markers", mc_id="aque_mc_f")

mcell_mc_plot_marks(mc_id="aque_mc_f", gset_id="aque_markers", mat_id="aque_filt",plot_cells = F)



## Make barplots for each gene of interest:

# change the working directory to yours here
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Amphimedon_sebepedros2018/Amphimedon/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Aque_genes-and-gene-types_best_markers.txt', sep='') # change here with correct name.
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

download.file("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/config.yaml","config.yaml")

tgconfig::override_params("config.yaml","metacell")
mcell_mc2d_force_knn(mc2d_id="aque_2dproj",mc_id="aque_mc_f", graph_id="aque_graph")

mcell_mc2d_plot(mc2d_id="aque_2dproj")


# Highlight specific genes in 2D projection graph

for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="aque_2dproj", 
                         gene = genename, show_mc_ids = F, 
                         show_legend = T, color_cells=T, mat_ds = NULL, 
                         zero_sc_v = 0, one_sc_v =1, 
                         base_dir = basedir)  
  }
}
print(paste('Done creating 2D plots. Please see', paste(mypath, 'figs/', sep='') , '.'))



### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="aque_mc_f", graph_id="aque_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="aque_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="aque_mc_f",
                        graph_id="aque_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)


#### Personalising heatmaps

# Aque doesn't have opsins, so will do heatmap for all other phototr genes directly.

# making gene set for phototr genes (best markers)
sets = c('CNG','CNG','PDE6 α_β','RGS9','GB5','G β','G γ','calmodulin','calmodulin','calmodulin','calmodulin','arrestin','GRK','GRK','Actin','CamKII','CamKII','DAGL','DAGL','G α q','G α q','INAD','INAD','IP3R','NINAC','PKC','PLC','PLC')
names(sets) = c('Aqu2.1.34432_001','Aqu2.1.40012_001','Aqu2.1.34466_001','Aqu2.1.36017_001','Aqu2.1.44326_001','Aqu2.1.25927_001','Aqu2.1.39671_001','Aqu2.1.37080_001','Aqu2.1.37081_001','Aqu2.1.37083_001','Aqu2.1.37082_001','Aqu2.1.41113_001','Aqu2.1.38104_001','Aqu2.1.39733_001','Aqu2.1.43436_001','Aqu2.1.27833_001','Aqu2.1.41582_001','Aqu2.1.38628_001','Aqu2.1.39940_001','Aqu2.1.43540_001','Aqu2.1.43541_001','Aqu2.1.02748_001','Aqu2.1.24128_001','Aqu2.1.42503_001','Aqu2.1.26962_001','Aqu2.1.43329_001','Aqu2.1.27482_001','Aqu2.1.32747_001')
gs = gset_new_gset(sets, 'Aque all phototr')
scdb_add_gset("Aque_all_phototr", gs)

# all phototr heatmap plot with metacells
mcell_mc_plot_marks(mc_id="aque_mc_f", gset_id="Aque_all_phototr", mat_id="aque_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Amphimedon_sebepedros2018/Amphimedon/figs/aque_heatmap_all_phototr_metacells.png')

# all phototr heatmap plot with cells
mcell_mc_plot_marks(mc_id="aque_mc_f", gset_id="Aque_all_phototr", mat_id="aque_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Amphimedon_sebepedros2018/Amphimedon/figs/aque_heatmap_all_phototr_cells.png')




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




# Extract top 100 for each metacell:
x <- 100  #set the number of top expressed genes here
for (metacell in colnames(lfp)){
  lines = ''
  top_x <- tail(sort(lfp[,metacell]), x)
  for (gene in names(top_x)){
    newline = paste(gene, top_x[gene][[1]], sep = ',')
    lines <- paste(lines, newline, '\n', sep = '')
  }
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Amphimedon_sebepedros2018/Amphimedon/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
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

