############################################
## Define computer (for testing purposes) ##
############################################

if (Sys.info()[['nodename']]=="BI2404M") {
  data_folder <- "/Users/argelagr/data/gastrulation_multiome_10x/shiny"
} else if (Sys.info()[['nodename']]=="BAI-JTYCW2FJ91") {
  data_folder <- "/Users/rargelaguet/data/gastrulation_multiome_10x/shiny"
}

###############
## libraries ##
###############

library(R.utils)
library(HDF5Array)
library(data.table)
library(purrr)
library(GGally)
library(DT)

# general viz
library(cowplot)
library(ggrepel)
library(ggplot2)
require(patchwork)
require(ggpubr) # to do: remove this dependency?

# tf motifs
library(TFBSTools)
library(ggseqlogo)

# shiny
library(shiny)
library(shinyFiles)
library(shinythemes)
library(ggiraph)

# graph viz
require(visNetwork)
# library(sna)
require(igraph)
library(network)
library(ggraph)
library(tidygraph)

######################
## Global variables ##
######################

nmp_celltypes <- c("Spinal_cord","NMP","Somitic_mesoderm")
chr_mm10 <- paste0("chr",c(1:19,"X","Y"))

#####################
## Colour palettes ##
#####################

celltype_colours <- c(
  "Epiblast" = "#635547",
  "Primitive_Streak" = "#DABE99",
  "Caudal_epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Anterior_Primitive_Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def._endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent_mesoderm" = "#C594BF",
  "Mixed_mesoderm" = "#DFCDE4",
  "Intermediate_mesoderm" = "#139992",
  "Caudal_Mesoderm" = "#3F84AA",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE_mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial_progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  # "Blood_progenitors" = "#c9a997",
  "Blood_progenitors_1" = "#f9decf",
  "Blood_progenitors_2" = "#c9a997",
  # "Erythroid" = "#EF4E22",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  # "Neurectoderm" = "#65A83E",
  "Rostral_neurectoderm" = "#65A83E",
  # "Caudal_neurectoderm" = "#354E23",
  "Neural_crest" = "#C3C388",
  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "Visceral_endoderm" = "#F6BFCB",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#989898",
  "Parietal_endoderm" = "#1A1A1A"
)

celltype_palette = scale_color_manual(values = celltype_colours, name = "", drop=TRUE)
celltype_palette_fill = scale_fill_manual(values = celltype_colours, name = "", drop=TRUE)
sample_palette <- scale_color_brewer(palette="Dark2")
rna_palette <- scale_color_gradient(low = "gray80", high = "red")
atac_palette <- scale_color_gradient(low = "gray80", high = "blue")


stage_colours <- c(
  "E7.5" = "#edf8b1",
  "E7.75" = "#c7e9b4",
  "E8.0" = "#7fcdbb",
  "E8.25" = "#41b6c4",
  "E8.5" = "#1d91c0",
  "E8.75" = "#225ea8"
)

stage_palette = scale_color_manual(values = stage_colours, name = "stage")
stage_palette_fill = scale_fill_manual(values = stage_colours, name = "stage")


###############
## Functions ##
###############

minmax.normalisation <- function(x) {
  return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]

give.n <- function(x) { return(c(y = max(x), label = length(x))) }
