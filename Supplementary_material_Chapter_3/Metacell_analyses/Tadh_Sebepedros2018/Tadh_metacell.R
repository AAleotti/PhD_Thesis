######## Metacell Analysis #######


# Based on the tutorial: https://tanaylab.github.io/metacell/articles/d-amphimedon.html
# Modified by: Alessandra Aleotti, Maryam Ghaffari Saadat

# Tadh single cell dataset obtained from here: https://www.wisdom.weizmann.ac.il/~arnau/
# Corresponding research paper: https://www.nature.com/articles/s41559-018-0575-6


# set library path
.libPaths('/data/evassvis/software/R/4.0')

# set Working Directory here:
setwd('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018')

library(metacell)

# make directories
if(!dir.exists("tadh")) dir.create("tadh/")
scdb_init("tadh/", force_reinit=T)

# import UMIs for Trichoplax
mcell_import_multi_mars("tadh", "MARS_Batches.txt", base_dir="/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/raw_umi_tables")

mat = scdb_mat("tadh")
print(dim(mat@mat))
#[1] 16723  4608

if(!dir.exists("figs")) dir.create("figs/")
scfigs_init("figs/")

### Exploring and filtering the UMI matrix

erccs <- rownames(mat@mat)[grepl("ERCC",rownames(mat@mat))]
mcell_mat_ignore_genes(new_mat_id="tadh_noercc",mat_id="tadh",ig_genes=erccs)
mat <- scdb_mat("tadh_noercc")

# Check umi distribution and choose cut-off. 
# In tutorial of Amphimedon, the lower cut-off is 200 and higher one is 12000.
mcell_plot_umis_per_cell("tadh",min_umis_cutoff = 80) 
mcell_plot_umis_per_cell("tadh",min_umis_cutoff = 7000)

# For Trichoplax use lower cut-off 80 and higher cut-off 7000
cell_sizes <- colSums(as.matrix(mat@mat))
large_cells <- names(which(cell_sizes>7000))
small_cells <- names(which(cell_sizes<80))
mcell_mat_ignore_cells("tadh_filt","tadh_noercc",ig_cells=c(small_cells,large_cells))
mat <- scdb_mat("tadh_filt")
print(dim(mat@mat))
#[1] 16631  3658

### Selecting gene markers

# calculate statistics.
mcell_add_gene_stat(gstat_id="tadh", mat_id="tadh_filt", force=T)

# Skip the blacklist filtering, keep rest as tutorial
# Also we adjust the minimum tot UMIs to 80.
mcell_gset_filter_multi(gset_id = "tadh_feats",gstat_id="tadh",T_tot=80,T_top3=2,T_szcor=-0.1,T_niche=0.08,force_new=T)
#Selected 258 markers

mcell_plot_gstats(gstat_id="tadh", gset_id="tadh_feats")

### Building the balanced cell graph

mcell_add_cgraph_from_mat_bknn(mat_id="tadh_filt",
                               gset_id = "tadh_feats",
                               graph_id="tadh_graph",
                               K=150,
                               dsamp=F)

### Resampling and generating the co-clustering graph

mcell_coclust_from_graph_resamp(
  coc_id="tadh_coc1000",
  graph_id="tadh_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=1000)

mcell_mc_from_coclust_balanced(
  coc_id="tadh_coc1000",
  mat_id= "tadh_filt",
  mc_id= "tadh_mc",
  K=20, min_mc_size=20, alpha=2)
# messages:
#filtered 749139 left with 248761 based on co-cluster imbalance
#building metacell object, #mc 35
#add batch counts
#compute footprints
#compute absolute ps
#compute coverage ps
#reordering metacells by hclust and most variable two markers
#reorder on Tadh_P64165 vs bin_44398

### Removing outlier cells

mcell_plot_outlier_heatmap(mc_id="tadh_mc", mat_id = "tadh_filt", T_lfc=3)

mcell_mc_split_filt(new_mc_id="tadh_mc_f",
                    mc_id="tadh_mc",
                    mat_id="tadh_filt",
                    T_lfc=3, plot_mats=F)
#messages:
#starting split outliers 
#splitting metacell 29
#splitting metacell 29
#add batch counts
#compute footprints
#compute absolute ps
#compute coverage ps

### Creating heatmaps of metacells and genes

mc_f<- scdb_mc("tadh_mc_f")
mc_f@colors <- colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(ncol(mc_f@mc_fp))
scdb_add_mc("tadh_mc_f",mc_f)
mc_f <- scdb_mc("tadh_mc_f")

mcell_gset_from_mc_markers(gset_id="tadh_markers", mc_id="tadh_mc_f")

mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="tadh_markers", mat_id="tadh_filt",plot_cells = F)

### Make barplots for each gene of interest (how much expression in each metacell):

# change the working directory for the dataset at hand:
mypath = '/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/'
barplots_path = paste(mypath, 'figs/barplots/', sep='')
metacells_path = paste(mypath, 'metacells/', sep='')
genedata_path = paste(mypath, 'Tadh_genes-and-gene-types_best.txt', sep='') # Check filename is correct!!
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
mcell_mc2d_force_knn(mc2d_id="tadh_2dproj",mc_id="tadh_mc_f", graph_id="tadh_graph")

mcell_mc2d_plot(mc2d_id="tadh_2dproj")

## 2d graphs showing expression for each gene of interest:
# will create a folder for each gene where inside there are the 2d graphs of each sequence for that gene
for (i in seq(1:nrow(genedata))){
  genename = genedata[i,1]
  genetype = genedata[i,2]
  basedir = paste(mypath,'figs/',genetype,'/', sep='')
  if(!dir.exists(basedir)) dir.create(basedir)
  if(genename %in% rownames(lfp)){ #create a plot only if the gene is in the metacells
    mcell_mc2d_plot_gene(mc2d_id="tadh_2dproj",  # change name here to appropriate dataset.
                         gene = genename, show_mc_ids = F, 
                         show_legend = T, color_cells=T, mat_ds = NULL, 
                         zero_sc_v = 0, one_sc_v =1, 
                         base_dir = basedir)  
  }
}


### Visualizing the metacell confusion matrix

mc_hc <- mcell_mc_hclust_confu(mc_id="tadh_mc_f", graph_id="tadh_graph")

mc_sup <- mcell_mc_hierarchy(mc_id="tadh_mc_f",mc_hc=mc_hc, T_gap=0.04)

mcell_mc_plot_hierarchy(mc_id="tadh_mc_f",
                        graph_id="tadh_graph",
                        mc_order=mc_hc$order,
                        sup_mc = mc_sup,
                        width=2800,
                        height=2000,
                        min_nmc=2)


#### Personalising heatmaps for genes of interest


# Make custom gset for Tadh placopsins
sets = c('placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin')
names(sets) = c('Tadh_P54494','Tadh_P61884','Tadh_P58589','Tadh_P58559','Tadh_P5735','Tadh_P58557','Tadh_P58576','Tadh_P58578','Tadh_P4017','Tadh_P15129','Tadh_P28334','Tadh_P53608','Tadh_P58590','Tadh_P951','Tadh_P28157','Tadh_P4346') # placopsins recovered from reconciliation(s) results.
gs = gset_new_gset(sets, 'Tadh placopsin set')
scdb_add_gset("Tadh_placopsin_set", gs)

# plot heatmap for placopsins using metacells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_placopsin_set", mat_id="tadh_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_placopsins_metacells.png')
# plot heatmap for placopsins using cells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_placopsin_set", mat_id="tadh_filt",plot_cells = T, fig_fn='//scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_placopsins_cells.png')


# Making gset for all phototr genes:
sets = c('CNG_a_and_b','G_aC','G_aC','G_aC','G_aC','GB5','NCKX','NCKX','NCKX','PDE6_a_and_b','PDE6_a_and_b','PDE6_a_and_b','RGS9','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','arrestin','calmodulin','calmodulin','G_b','G_g','RK','RK','CamkII','CamkII','DAGL','DAGL','DAGL','G_aR','INAD','IP3R','IP3R','IP3R','NINAC','PKC','PKC','PKC','PLC','PLC','PLC','rdgC','TRP_TRPL')
names(sets) = c('Tadh_P13081','Tadh_P24489','Tadh_P30457','Tadh_P50796','Tadh_P50797','Tadh_P49804','Tadh_P14805','Tadh_P33154','Tadh_P63252','Tadh_P21752','Tadh_P28362','Tadh_P50691','Tadh_P52866','Tadh_P54494','Tadh_P61884','Tadh_P58589','Tadh_P58559','Tadh_P5735','Tadh_P58557','Tadh_P58576','Tadh_P58578','Tadh_P4017','Tadh_P15129','Tadh_P28334','Tadh_P53608','Tadh_P58590','Tadh_P951','Tadh_P28157','Tadh_P4346','Tadh_P64255','Tadh_P24722','Tadh_P37105','Tadh_P64010','Tadh_P37130','Tadh_P28308','Tadh_P59012','Tadh_P31065','Tadh_P50756','Tadh_P53544','Tadh_P53545','Tadh_P53546','Tadh_P64258','Tadh_P27973','Tadh_P55241','Tadh_P56365','Tadh_P60448','Tadh_P22946','Tadh_P20496','Tadh_P30067','Tadh_P57608','Tadh_P29799','Tadh_P61893','Tadh_P61894','Tadh_P56332','Tadh_P60588')
gs = gset_new_gset(sets, 'Tadh all phototr')
scdb_add_gset("Tadh_all_phototr", gs)

# plot heatmap for all phototr using metacells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_all_phototr", mat_id="tadh_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_phototr_metacells.png')
# plot heatmap for all phototr using cells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_all_phototr", mat_id="tadh_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_phototr_cells.png')


# Making heatmap for all ciliary (including common) genes:
sets = c('CNG_a_and_b','G_aC','G_aC','G_aC','G_aC','GB5','NCKX','NCKX','NCKX','PDE6_a_and_b','PDE6_a_and_b','PDE6_a_and_b','RGS9','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','arrestin','calmodulin','calmodulin','G_b','G_g','RK','RK')
names(sets) = c('Tadh_P13081','Tadh_P24489','Tadh_P30457','Tadh_P50796','Tadh_P50797','Tadh_P49804','Tadh_P14805','Tadh_P33154','Tadh_P63252','Tadh_P21752','Tadh_P28362','Tadh_P50691','Tadh_P52866','Tadh_P54494','Tadh_P61884','Tadh_P58589','Tadh_P58559','Tadh_P5735','Tadh_P58557','Tadh_P58576','Tadh_P58578','Tadh_P4017','Tadh_P15129','Tadh_P28334','Tadh_P53608','Tadh_P58590','Tadh_P951','Tadh_P28157','Tadh_P4346','Tadh_P64255','Tadh_P24722','Tadh_P37105','Tadh_P64010','Tadh_P37130','Tadh_P28308','Tadh_P59012')
gs = gset_new_gset(sets, 'Tadh ciliary')
scdb_add_gset("Tadh_ciliary", gs)

# plot heatmap for all ciliary genes using metacells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_ciliary", mat_id="tadh_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_all_ciliary_metacells.png')
# plot heatmap for all ciliary genes using cells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_ciliary", mat_id="tadh_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_all_ciliary_cells.png')


# Making heatmap for only ciliary (NO common) genes:
sets = c('CNG_a_and_b','G_aC','G_aC','G_aC','G_aC','GB5','NCKX','NCKX','NCKX','PDE6_a_and_b','PDE6_a_and_b','PDE6_a_and_b','RGS9')
names(sets) = c('Tadh_P13081','Tadh_P24489','Tadh_P30457','Tadh_P50796','Tadh_P50797','Tadh_P49804','Tadh_P14805','Tadh_P33154','Tadh_P63252','Tadh_P21752','Tadh_P28362','Tadh_P50691','Tadh_P52866')
gs = gset_new_gset(sets, 'Tadh only ciliary')
scdb_add_gset("Tadh_only_ciliary", gs)

# plot heatmap for only ciliary genes using metacells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_only_ciliary", mat_id="tadh_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_ciliary_no_common_metacells.png')
# plot heatmap for only ciliary genes using cells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_only_ciliary", mat_id="tadh_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_ciliary_no_common_cells.png')


# Making heatmap for only rhabdomeric (NO common) genes:
sets = c('CamkII','CamkII','DAGL','DAGL','DAGL','G_aR','INAD','IP3R','IP3R','IP3R','NINAC','PKC','PKC','PKC','PLC','PLC','PLC','rdgC','TRP_TRPL')
names(sets) = c('Tadh_P31065','Tadh_P50756','Tadh_P53544','Tadh_P53545','Tadh_P53546','Tadh_P64258','Tadh_P27973','Tadh_P55241','Tadh_P56365','Tadh_P60448','Tadh_P22946','Tadh_P20496','Tadh_P30067','Tadh_P57608','Tadh_P29799','Tadh_P61893','Tadh_P61894','Tadh_P56332','Tadh_P60588')
gs = gset_new_gset(sets, 'Tadh only rhabdomeric')
scdb_add_gset("Tadh_only_rhabdomeric", gs)

# plot heatmap for only rhabdomeric genes using metacells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_only_rhabdomeric", mat_id="tadh_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_rhabdomeric_no_common_metacells.png')
# plot heatmap for only ciliary genes using cells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_only_rhabdomeric", mat_id="tadh_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_rhabdomeric_no_common_cells.png')


# Making heatmap for all rhabdomeric (including common) genes:
sets = c('placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','arrestin','calmodulin','calmodulin','G_b','G_g','RK','RK','CamkII','CamkII','DAGL','DAGL','DAGL','G_aR','INAD','IP3R','IP3R','IP3R','NINAC','PKC','PKC','PKC','PLC','PLC','PLC','rdgC','TRP_TRPL')
names(sets) = c('Tadh_P54494','Tadh_P61884','Tadh_P58589','Tadh_P58559','Tadh_P5735','Tadh_P58557','Tadh_P58576','Tadh_P58578','Tadh_P4017','Tadh_P15129','Tadh_P28334','Tadh_P53608','Tadh_P58590','Tadh_P951','Tadh_P28157','Tadh_P4346','Tadh_P64255','Tadh_P24722','Tadh_P37105','Tadh_P64010','Tadh_P37130','Tadh_P28308','Tadh_P59012','Tadh_P31065','Tadh_P50756','Tadh_P53544','Tadh_P53545','Tadh_P53546','Tadh_P64258','Tadh_P27973','Tadh_P55241','Tadh_P56365','Tadh_P60448','Tadh_P22946','Tadh_P20496','Tadh_P30067','Tadh_P57608','Tadh_P29799','Tadh_P61893','Tadh_P61894','Tadh_P56332','Tadh_P60588')
gs = gset_new_gset(sets, 'Tadh rhabdomeric')
scdb_add_gset("Tadh_rhabdomeric", gs)

# plot heatmap for all rhabdomeric genes using metacells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_rhabdomeric", mat_id="tadh_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_all_rhabdomeric_metacells.png')
# plot heatmap for rhabdomeric genes using cells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_rhabdomeric", mat_id="tadh_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_all_rhabdomeric_cells.png')


# Making heatmap for common genes:
sets = c('placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','placopsin','arrestin','calmodulin','calmodulin','G_b','G_g','RK','RK')
names(sets) = c('Tadh_P54494','Tadh_P61884','Tadh_P58589','Tadh_P58559','Tadh_P5735','Tadh_P58557','Tadh_P58576','Tadh_P58578','Tadh_P4017','Tadh_P15129','Tadh_P28334','Tadh_P53608','Tadh_P58590','Tadh_P951','Tadh_P28157','Tadh_P4346','Tadh_P64255','Tadh_P24722','Tadh_P37105','Tadh_P64010','Tadh_P37130','Tadh_P28308','Tadh_P59012')
gs = gset_new_gset(sets, 'Tadh common')
scdb_add_gset("Tadh_common", gs)

# plot heatmap for common genes using metacells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_common", mat_id="tadh_filt",plot_cells = F, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_common_metacells.png')
# plot heatmap for common genes using cells
mcell_mc_plot_marks(mc_id="tadh_mc_f", gset_id="Tadh_common", mat_id="tadh_filt",plot_cells = T, fig_fn='/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/figs/Tadh_heatmap_common_cells.png')



### Finding genes expressed in a metacell of interest:

lfp = log2(mc_f@mc_fp) 

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
  write(lines, file=paste('/scratch/evassvis/aa1176/scRNAseq_analysis_working_dir/Trichoplax_sebepedros2018/metacells/metacell-',metacell,'_top',x,'.csv',sep=''))
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