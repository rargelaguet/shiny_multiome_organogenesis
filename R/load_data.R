if (Sys.info()[['nodename']]=="BI2404M") {
  data_folder <- "/Users/argelagr/data/gastrulation_multiome_10x/test/shiny"
} else if (Sys.info()[['nodename']]=="CAM-JTYCW2FJ91") {
  data_folder <- "/Users/rargelaguet/data/gastrulation_multiome_10x/shiny"
} 

###########################
## Load global variables ##
###########################

# source("R/global.R")

stages <- fread(paste0(data_folder,"/stages.txt"), header=F)[[1]]
samples <- fread(paste0(data_folder,"/samples.txt"), header=F)[[1]]
celltypes <- fread(paste0(data_folder,"/celltypes.txt"), header=F)[[1]]
TFs_all <- fread(paste0(data_folder,"/TFs.txt"), header=F)[[1]]
rna_genes <- fread(paste0(data_folder,"/rna_expression/genes.txt"), header=F)[[1]]
atac_genes <- fread(paste0(data_folder,"/chromatin_accessibility/atac_genes.txt"), header=F)[[1]]
atac_peaks <- fread(paste0(data_folder,"/chromatin_accessibility/atac_peaks.txt"), header=F)[[1]]

TFs <- list("CISBP"=fread(paste0(data_folder,"/tfs_cisbp.txt"))[[1]], "JASPAR"=fread(paste0(data_folder,"/tfs_jaspar.txt"))[[1]])
  
########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(file.path(data_folder,"metadata/cell_metadata.txt.gz")) %>% 
  .[,celltype := factor(celltype, levels = names(celltype_colours), ordered = TRUE)] %>%
  .[,stage := factor(stage, levels = names(stage_colours), ordered = TRUE)] %>%
  # .[,sample = factor(sample, levels = names(samples), ordered = TRUE)
  setkey(cell)

metacell_metadata.dt <- fread(file.path(data_folder,"metadata/metacell_metadata.txt.gz")) %>% 
  .[,celltype := factor(celltype, levels = names(celltype_colours), ordered = TRUE)] %>%
  .[,stage := factor(stage, levels = names(stage_colours), ordered = TRUE)] %>%
  # .[,sample = factor(sample, levels = names(samples), ordered = TRUE)
  setkey(metacell)

#########################
## Load RNA expression ##
#########################

# cells
cells <- fread(paste0(data_folder,"/rna_expression/rna_expr_cells.txt"), header=F)[[1]]
rna_expr_cells.array = HDF5Array(file = paste0(data_folder,"/rna_expression/rna_expr_cells.hdf5"), name = "rna_expr_logcounts")
colnames(rna_expr_cells.array) <- cells
rownames(rna_expr_cells.array) <- rna_genes

# metacells
metacells <- fread(paste0(data_folder,"/rna_expression/rna_expr_metacells.txt"), header=F)[[1]]
rna_expr_metacells.array = HDF5Array(file = paste0(data_folder,"/rna_expression/rna_expr_metacells.hdf5"), name = "rna_expr_logcounts")
colnames(rna_expr_metacells.array) <- metacells
rownames(rna_expr_metacells.array) <- rna_genes

# pseudobulk
rna_expr_pseudobulk.mtx <- fread(paste0(data_folder,"/rna_expression/rna_expr_pseudobulk.txt.gz"), sep = ",") %>% matrix.please

# pseudobulk with replicates
rna_expr_pseudobulk_replicates.mtx <- fread(paste0(data_folder,"/rna_expression/rna_expr_pseudobulk_replicates.txt.gz"), sep = ",") %>% matrix.please

##################################
## Load chromatin accessibility ##
##################################

metacell_metadata_atac.dt <- fread(paste0(data_folder,"/chromatin_accessibility/metacells_metaadata.txt"))

atac_cells <- fread(paste0(data_folder,"/chromatin_accessibility/atac_cells.txt"), header=F)[[1]]
atac_metacells <- fread(paste0(data_folder,"/chromatin_accessibility/atac_metacells.txt"), header=F)[[1]]

# Gene scores (cells)
atac_genes_cells.array = HDF5Array(file = paste0(data_folder,"/chromatin_accessibility/atac_genes_cells.hdf5"), name = "atac_genes")
colnames(atac_genes_cells.array) <- atac_cells
rownames(atac_genes_cells.array) <- atac_genes

# Gene scores (metacells)
atac_genes_metacells.array = HDF5Array(file = paste0(data_folder,"/chromatin_accessibility/atac_genes_metacells.hdf5"), name = "atac_genes_logcounts")
colnames(atac_genes_metacells.array) <- atac_metacells
rownames(atac_genes_metacells.array) <- atac_genes

# Gene scores (pseudobulk)
atac_genes_pseudobulk.mtx = fread(paste0(data_folder,"/chromatin_accessibility/atac_genes_pseudobulk.txt.gz"), sep = ",") %>% matrix.please

# Gene scores (pseudobulk with replicates)
atac_genes_pseudobulk_replicates.mtx = fread(paste0(data_folder,"/chromatin_accessibility/atac_genes_pseudobulk_replicates.txt.gz"), sep = ",") %>% matrix.please

# ATAC peaks (cells)
atac_peaks_cells.array = HDF5Array(file = paste0(data_folder,"/chromatin_accessibility/atac_peaks_cells.hdf5"), name = "atac_peaks")
colnames(atac_peaks_cells.array) <- atac_cells
rownames(atac_peaks_cells.array) <- atac_peaks

# ATAC peaks (metacells)
atac_peaks_metacells.array = HDF5Array(file = paste0(data_folder,"/chromatin_accessibility/atac_peaks_metacells.hdf5"), name = "atac_peaks_logcounts")
colnames(atac_peaks_metacells.array) <- atac_metacells
rownames(atac_peaks_metacells.array) <- atac_peaks

# ATAC peaks (pseudobulk)
atac_peaks_pseudobulk.mtx = fread(paste0(data_folder,"/chromatin_accessibility/atac_peaks_pseudobulk.txt.gz"), sep = ",") %>% matrix.please

# ATAC peaks (pseudobulk with replicates)
atac_peaks_pseudobulk_replicates.mtx = fread(paste0(data_folder,"/chromatin_accessibility/atac_peaks_pseudobulk_replicates.txt.gz"), sep = ",") %>% matrix.please

#####################
## Load PAGA graph ##
#####################

connectivity.mtx <- fread(file.path(data_folder,"paga/paga_connectivity.csv")) %>%
  matrix.please %>% .[celltypes,celltypes]

coordinates.mtx <- fread(file.path(data_folder,"paga/paga_coordinates.csv")) %>% 
  matrix.please %>% .[celltypes,]

# Parse data
connectivity.mtx[connectivity.mtx<0.20] <- 0
connectivity.mtx[connectivity.mtx>=0.20] <- 1

# Create igraph object
igraph.paga <- graph_from_adjacency_matrix(connectivity.mtx, mode = "undirected")

# Create tbl_graph object
igraph.paga.tbl <- as_tbl_graph(igraph.paga) %>%
  activate(nodes) %>%
  mutate(celltype=rownames(connectivity.mtx)) %>%
  mutate(x=coordinates.mtx[,1]) %>% mutate(y=coordinates.mtx[,2])

# Create network object
net.paga = network(connectivity.mtx)
net.paga %v% "x" = coordinates.mtx[,1]
net.paga %v% "y" = coordinates.mtx[,2]

###############
## TF motifs ##
###############

motif2gene.list <- list()
motif2gene.list[["JASPAR"]] <- fread(file.path(data_folder,"tf_motifs/JASPAR_motif2gene.txt.gz"))[,c("motif","gene")]
motif2gene.list[["CISBP"]] <- fread(file.path(data_folder,"tf_motifs/CISBP_motif2gene.txt.gz"))[,c("motif","gene")]

tf_motifs_pwm.list <- readRDS(file.path(data_folder,"tf_motifs/tf_motifs.rds"))

# motifs_all <- unique(unlist(sapply(motif2gene.list, function(x) x$motif)))
# write.table(motifs_all, paste0(data_folder,"/tf_motifs/motif_names.txt.gz"), row.names = F, col.names = F)
motifs_all <- fread(paste0(data_folder,"/tf_motifs/motif_names.txt.gz"), header=F)[[1]]

######################
## TF marker scores ##
######################

tf_markers_rna.dt <- fread(paste0(data_folder,"/tf_markers/CISBP/tf_markers_rna.txt.gz"))
tf_markers_chromvar.dt <- fread(paste0(data_folder,"/tf_markers/CISBP/tf_markers_chromvar.txt.gz"))

#########################
## Load reference UMAP ##
#########################

umap_reference.dt <- fread(paste0(data_folder,"/mapping/umap_coordinates.txt.gz"))

######################################
## Load differential RNA expression ##
######################################

# diff_rna_cells.dt <- fread(paste0(data_folder,"/differential/diff_rna_cells.txt.gz"))
# diff_rna_metacells.dt <- fread(paste0(data_folder,"/differential/diff_rna_metacells.txt.gz"))
# diff_rna_pseudobulk.dt <- fread(paste0(data_folder,"/differential/diff_rna_pseudobulk.txt.gz"))

###############################################
## Load differential chromatin accessibility ##
###############################################

# diff_rna_cells.dt <- fread(paste0(data_folder,"/differential/diff_rna_cells.txt.gz"))
# diff_rna_metacells.dt <- fread(paste0(data_folder,"/differential/diff_rna_metacells.txt.gz"))
# diff_rna_pseudobulk.dt <- fread(paste0(data_folder,"/differential/diff_rna_pseudobulk.txt.gz"))

##############################
## Load NMP differentiation ##
##############################

nmp_metadata.dt <- fread(paste0(data_folder,"/nmp_differentiation/metadata.txt.gz"))

nmp_metacells.rna <- fread(paste0(data_folder,"/nmp_differentiation/rna_metacells.txt"), header=F)[[1]]
nmp_genes <- fread(paste0(data_folder,"/nmp_differentiation/rna_genes.txt"), header=F)[[1]]
nmp_rna_expr.array = HDF5Array(paste0(data_folder,"/nmp_differentiation/rna_expr_metacells.hdf5"), name = "rna_expr_logcounts")
colnames(nmp_rna_expr.array) <- nmp_metacells.rna
rownames(nmp_rna_expr.array) <- nmp_genes

nmp_metacells.atac <- fread(paste0(data_folder,"/nmp_differentiation/atac_metacells.txt"), header=F)[[1]]
nmp_atac_peaks <- fread(paste0(data_folder,"/nmp_differentiation/atac_peaks.txt"), header=F)[[1]]
nmp_atac_peaks.array = HDF5Array(paste0(data_folder,"/nmp_differentiation/atac_peaks_metacells.hdf5"), name = "atac_peaks_logcounts")
colnames(nmp_atac_peaks.array) <- nmp_metacells.atac
rownames(nmp_atac_peaks.array) <- nmp_atac_peaks

# Load GRN
nmp_grn.igraph <- readRDS(paste0(data_folder,"/nmp_differentiation/network_igraph.rds"))

#####################################
## Load in silico ChIP-seq results ##
#####################################

insilico_chip_stats.dt <- fread(paste0(data_folder,"/insilico_chip/insilico_chip_stats.txt.gz"))

