
#####################
## Define settings ##
#####################

# if (Sys.info()[['nodename']]=="BI2404M") {
#   data_folder <- "/Users/argelagr/data/gastrulation_multiome_10x/test/shiny"
# } else if (grepl("CAM-",Sys.info()[['nodename']])) {
#   data_folder <- "/Users/rargelaguet/data/gastrulation_multiome_10x/shiny"
# }

# for testing
# shiny::loadSupport()

###########################
## Load global variables ##
###########################

chr_mm10 <- paste0("chr",c(1:19,"X","Y")); names(chr_mm10) <- chr_mm10
stages <- fread(paste0(data_folder,"/stages.txt"), header=F)[[1]]
samples <- fread(paste0(data_folder,"/samples.txt"), header=F)[[1]]
celltypes <- fread(paste0(data_folder,"/celltypes.txt"), header=F)[[1]]
TFs_all <- fread(paste0(data_folder,"/TFs.txt"), header=F)[[1]]
rna_genes <- fread(paste0(data_folder,"/rna_expression/genes.txt"), header=F)[[1]]
# atac_peaks <- fread(paste0(data_folder,"/chromatin_accessibility/atac_peaks.txt"), header=F)[[1]]
# atac_genes <- fread(paste0(data_folder,"/chromatin_accessibility/atac_genes.txt"), header=F)[[1]]
TFs <- list("CISBP"=fread(paste0(data_folder,"/tfs_cisbp.txt"))[[1]], "JASPAR"=fread(paste0(data_folder,"/tfs_jaspar.txt"))[[1]])
nmp_celltypes <- c("Spinal_cord","NMP","Somitic_mesoderm")

##################
## Shiny app UI ##
##################

ui <- shinyUI(
  fluidPage(
    navbarPage(
        title = "Decoding gene regulation in the mouse embryo using single-cell multi-omics",
        theme = shinytheme("spacelab"),
        
        ##########
        ## UMAP ##
        ##########
        
        tabPanel(
          title = "UMAPs", id = "umap",
          sidebarPanel(width=3,
            selectInput(inputId = "umap_modality", label = "Modality", choices = c("ATAC" = "atac", "RNA" = "rna", "RNA+ATAC (MOFA)"="rna_atac"), selected = "rna_atac"),
            selectInput(inputId = "umap_stage", label = "Stages", choices = c("All stages"="all", "E7.5"="E7.5", "E7.75"="E7.75", "E8.5"="E8.5", "E8.75"="E8.75"), selected = "all"),
            # selectInput(inputId = "umap_sample", label = "Samples", choices = samples_choice, selected = "all"),
            selectInput(inputId = "umap_colourby", label = "Plot colour", choices = c("Cell type"="celltype", "Stage"="stage", "Sample"="sample", "Gene expression"="gene_expression"), selected = "celltype"),
            checkboxInput("umap_subset_cells", "Subset cells (faster plotting)", value=TRUE),
            conditionalPanel(
              condition = "input.umap_colourby == 'gene_expression'",
              selectizeInput("umap_gene_rna", "Select gene", choices = NULL, selected = "T")
            ),
            conditionalPanel(
              condition = "input.umap_colourby == 'chromatin_accessibility'",
              selectizeInput("umap_gene_atac", "Select gene", choices = NULL, selected = "T")
            )
          ),
          mainPanel(
            girafeOutput("umap", width = "1100px", height = "800px")
            # plotOutput("stage_contribution", width = big_plot_width)
          )
        ),
        
        #############
        ## Mapping ##
        #############
        
        tabPanel(
          title = "Mapping", id = "mapping",
          sidebarPanel(width=3, 
                       selectizeInput("mapping_sample", "Select sample", choices = samples, selected = "E8.5_rep1"),
                       checkboxInput("mapping_subset_cells", "Subset atlas cells", value=TRUE)
                       # selectInput("mapping_sample", "Select sample", choices = samples, selected = samples, multiple = TRUE),
          ),
          mainPanel(
            girafeOutput("mapping_plot", width = "1000px", height = "600px")
          )
        ),
        
        #####################
        ## Gene expression ##
        #####################
        
        tabPanel(
          title = "Gene expression", id = "gene_expr",
          sidebarPanel(width=3,
                       selectizeInput("gene_expr_gene", "Select gene", choices = NULL),
                       selectInput("gene_expr_celltypes", "Select celltype", choices = celltypes, selected = celltypes, multiple = TRUE)
                       # checkboxGroupInput("gene_expr_stages", label = "Stages", choices = list("E7.5"="E7.5", "E7.75"="E7.75", "E8.0"="E8.0", "E8.5"="E8.5", "E8.75"="E8.75"), selected = c("E7.5","E7.75","E8.0","E8.5","E8.75"))
          ),
          mainPanel(
            girafeOutput("gene_expr", width = "700px", height = "800px")
            # plotOutput("stage_contribution", width = big_plot_width)
          )
        ),
        
        #############################
        ## Chromatin accessibility ##
        #############################
        
        tabPanel(
          title = "Chromatin accessibility", id = "chr_acc",
          sidebarPanel(width=3,
                       selectInput("chr_acc_feature_type", "Select feature type", choices = c("Genes (TSS)"="tss", "ATAC peaks"="peaks"), selected = "tss"),
                       conditionalPanel(
                         condition = "input.chr_acc_feature_type == 'peaks'",
                         selectInput("chr_acc_feature_chr", "Select chromosome", choices = chr_mm10, selected = "chr1"),
                       ),
                       selectizeInput("chr_acc_feature", "Select feature", choices = NULL),
                       selectInput("chr_acc_celltypes", "Select celltype", choices = celltypes, selected = celltypes, multiple = TRUE),
          ),
          mainPanel(
            girafeOutput("chr_acc", width = "700px", height = "800px")
          )
        ),
        
        #####################################
        ## TF activities (RNA vs chromVAR) ##
        #####################################
        
        tabPanel(
          title = "TF activities", id = "rna_vs_chromvar",
          sidebarPanel(width=3,
            selectInput("rna_vs_chromvar_motif_annotation", "Motif annotation", choices = c("JASPAR","CISBP"), selected = "CISBP"),
            # selectInput("chromvar_algorithm", "chromVAR algorithm", choices = c("chromVAR (standard)"="chromvar","chromVAR-Multiome (+ binding)"="chromvar_multiome_positive","chromVAR-Multiome (all binding)"="chromvar_multiome_all"), selected = "chromVAR (standard)"),
            selectInput("chromvar_algorithm", "chromVAR algorithm", choices = c("chromVAR (standard)"="chromvar","chromVAR-Multiome (+)"="chromvar_multiome_positive"), selected = "chromVAR (standard)"),
            selectizeInput(inputId = "rna_vs_chromvar_tf", label = "Select transcription Factor", choices=NULL),
            selectizeInput(inputId = "rna_vs_chromvar_motif", label = "Select motif", choices=NULL),
            plotOutput("rna_vs_chromvar_plot_motif", height="100px", width="300px"),
            checkboxInput("rna_vs_chromvar_scatterplot_add_text", "Add cell type labels to the scatterplot", value = T),
            sliderInput("rna_vs_chromvar_cor_range", label = "RNA vs chromVAR corr. interval", min=-1, max=1, step=0.01, value = c(-1,1))
          ),
          mainPanel(
            h6("chromVAR standard scores are calculated using all ATAC peaks that contain the TF motif instance"),
            # br(),
            h6("chromVAR-Multiome (+) scores are calculated using all putative binding sites that TF RNA expr. and chr. accessibility correlate positively"),
            # conditionalPanel(
            #   condition = "chromvar_algorithm == 'chromvar'",
            #   HTML("standard chromVAR scores are calculated using all ATAC peaks that contain the TF motif instance")
            # ),
            girafeOutput("rna_vs_chromvar_plot", width = "900px", height = "600px")
            # verbatimTextOutput('rna_vs_chromvar_print')
          )
        ),

        ###########################
        ## Differential analysis ##
        ###########################
        
        tabPanel(
          title = "Differential analysis", id = "diff",
          sidebarPanel(width=3,
            selectInput(inputId = "diff_celltypeA", label = "Select celltype A", choices = celltypes, selected = "Gut"),
            selectInput(inputId = "diff_celltypeB", label = "Select celltype B", choices = celltypes, selected = "Neural_crest"),
            selectInput(inputId = "diff_resolution", label = "Select resolution", choices = c("Cells","Metacells","Pseudobulk"), selected = "Metacells"),
            selectInput(inputId = "diff_modality", label = "Select data modality", choices = c("RNA","ATAC"), selected = "RNA"),
            conditionalPanel(
              condition = "input.diff_modality == 'ATAC'",
              selectInput("diff_atac_chr", "Select chromosome", choices = chr_mm10, selected = "chr1")
              ),
            selectizeInput("diff_feature", "Select feature", choices = NULL),
            sliderInput("min_log_pval", label = "Minimum log pvalue", min=0, max=100, step=1, value = 25),
            sliderInput("diff_range", label = "Range of differential values", min=-8, max=8, step=0.5, value = c(-8,8)),
            # sliderInput("highlight_top_n_genes", label = "Highlight top N genes", min=0, max=100, step=1, value = 0)
          ),
          mainPanel(
            girafeOutput("diff_plot", width = "900px", height = "600px"),
            HTML("For visualisation efficiency, there is a maximum of 10,000 features to be displayed")
          )
        ),
        
        
        #########################
        ## Celltype TF markers ##
        #########################
        
        tabPanel(
          title = "Celltype TF markers", id = "celltype_tf_markers",
          sidebarPanel(width=3,
            # selectInput(inputId = "RNA_vs_chromVAR_pseudobulk_per_celltype_colour_by", label = "Colour by", choices = c("Activator/Repressor" = "activator_repressor", "Marker strength" = "marker_strength")),
            selectInput("celltype_tf_markers_celltype", "Select celltype", choices = celltypes, selected = "Gut"),
            selectInput("celltype_tf_markers_tf", "Select TF", choices = TFs_all, selected = "FOXA2"),
            # selectInput("celltype_tf_markers_motif_annotation", "Motif annotation", choices = c("JASPAR","CISBP"), selected = "CISBP"),
            selectInput("celltype_tf_markers_modality", "Select modality", choices = c("RNA"="rna","chromVAR-Multiome"="chromvar"), selected = "RNA")
            # sliderInput("celltype_tf_markers_range", label = "TF marker score range", min=0, max=1, step=0.01, value = c(0,1))
          ),
          mainPanel(
            girafeOutput("celltype_tf_markers_plot", width = "1000px", height = "600px")
          )
        ),
        
        #########################
        ## NMP differentiation ##
        #########################
        
        tabPanel(
          title = "NMP differentiation", id = "nmp_diff",
          sidebarPanel(width=3,
            selectInput(inputId = "nmp_diff_colourby", label = "Trajectory colour", choices = c("Cell type"="celltype", "Gene expression"="gene_expression", "Chromatin accessibility"="chromatin_accessibility"), selected = "celltype"),
            conditionalPanel(
              condition = "input.nmp_diff_colourby == 'chromatin_accessibility'",
              selectInput("nmp_diff_atac_chr", "Select chromosome", choices = chr_mm10, selected = "chr1")
            ),
            conditionalPanel(
              condition = "input.nmp_diff_colourby == 'chromatin_accessibility' || input.nmp_diff_colourby == 'gene_expression'",
              selectizeInput("nmp_diff_feature", "Select feature", choices = NULL),
            ),
            selectInput(inputId = "nmp_grn_colourby", label = "GRN colour", choices = c("Cell type"="celltype", "Cell type-specific expression"="celltype_expression", "Eigenvalue centrality"="eigenvalue_centrality", "Degree centrality"="degree_centrality"), selected = "celltype"),
            conditionalPanel(
              condition = "input.nmp_grn_colourby == 'celltype_expression'",
              selectizeInput("nmp_grn_colourby_celltype", "Select cell type", choices = nmp_celltypes, selected = "NMP")
            ),
          ),
          mainPanel(
            fluidRow(
              splitLayout(
                cellWidths = c("57%", "43%"), 
                girafeOutput("nmp_diff_plot", width = "430px", height = "350px"), 
                visNetworkOutput("nmp_grn_plot", width = "420px", height = "350px")
              )
            )
          )
        ),
        
        ########################
        ## in silico ChIP-seq ##
        ########################
        
        tabPanel(
          title = "In silico ChIP-seq", id = "insilico_chipseq",
          sidebarPanel(width=3,
            selectInput("insilico_chip_tf", "Transcription Factor", choices = TFs_all, selected = c("FOXA2","T","CDX2","RUNX1","TFAP2A"), multiple = TRUE),
            selectInput("insilico_chip_motif_annotation", "Motif annotation", choices = c("JASPAR","CISBP"), selected = "CISBP"),
            # selectInput("celltype_chip", "Celltype", choices = celltypes, selected = "All celltypes"),
            sliderInput("insilico_chip_min_score", label = "Minimum score threshold", min = 0, max = 1, value = 0.15),
            shinyFilesButton('insilico_chip_files', label='File select', title='Please select a file', multiple=T) ,
            downloadButton("insilico_chip_download", "Download"),
            verbatimTextOutput('insilico_chip_download_print')
          ),
          mainPanel(
            girafeOutput("insilico_chip_plot"),
          )
        ),
        
        ##################
        ## Brachyury KO ##
        ##################
        
        #################
        ## IGV session ##
        #################
        
        tabPanel(
          title = "Genome browser", id = "genome_browser",
          mainPanel(
            # HTML(sprintf("%s: %s you can download a precomputed session to interactively explore the ATAC data on the %s (see snapshot below). Watch %s for instructions", shiny::tags$b("Download IGV session \n"), a("Here", href="igv_session_template.xml"), a("IGV browser", href="https://software.broadinstitute.org/software/igv"), a("this video", href="XXX"))),
            # HTML(sprintf("%s you can find the instructions to download a precomputed session to interactively explore the ATAC data on the %s (see snapshot below). You can also watch %s.", a("Here", href="https://github.com/rargelaguet/mouse_organogenesis_10x_multiome_publication#igv-genome-browser-session"), a("IGV browser", href="https://software.broadinstitute.org/software/igv"), a("this videotutorial", href="XXX"))),
            HTML(sprintf("%s you can find the instructions to download a precomputed session to interactively explore the ATAC data on the %s (see snapshot below)", a("Here", href="https://github.com/rargelaguet/mouse_organogenesis_10x_multiome_publication#igv-genome-browser-session"), a("IGV browser", href="https://software.broadinstitute.org/software/igv"))),
            br(), br(),
            # shiny::img(src = "igv_session_overview.png", height = 450, width = 900),
            HTML('<center><img src="igv_session_overview.png", height="450px", width="900px"></center>'),
            # br(), br(),
            # HTML(sprintf("%s: %s you can download a precomputed session to interactively explore the ATAC data on the %s (see snapshot below). Watch %s for instructions", shiny::tags$b("Download UCSC session \n"), a("Here", href="XXX"), a("UCSC browser", href="https://software.broadinstitute.org/software/igv"), a("this video", href="XXX")))
          )
        ),
        shiny::tags$style(HTML(".navbar-header { width:100% } .navbar-brand { width: 100%; text-align: center; font-size: 27px }")) # center text
        
      )
    )
  )
