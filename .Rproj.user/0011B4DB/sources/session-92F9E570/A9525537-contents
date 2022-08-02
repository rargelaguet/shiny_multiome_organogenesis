
#####################
## Define settings ##
#####################

# if (Sys.info()[['nodename']]=="BI2404M") {
#   setwd("/Users/argelagr/gastrulation_multiome_10x/shiny")
#   data_folder <- "/Users/argelagr/data/gastrulation_multiome_10x/test/shiny"
# } else if (Sys.info()[['nodename']]=="rargelaguet.local") {
#   setwd("/Users/rargelaguet/gastrulation_multiome_10x/shiny")
#   data_folder <- "/Users/rargelaguet/data/gastrulation_multiome_10x/shiny"
# } 

# for testing
# shiny::loadSupport()
# source("load_data.R")

###############
## Shiny app ##
###############

server <- function(input, output, session) {
  
  ###########
  ## UMAPs ##
  ###########
  
  updateSelectizeInput(session = session, inputId = 'umap_gene_rna', choices = rna_genes, server = TRUE, selected = "T")
  
  plot_UMAP <- reactive({
      
    ## START TEST ##
    # input <- list()
    # input$umap_modality <- "rna_atac"
    # input$umap_gene_rna <- "T"
    # input$umap_stage <- "all"
    # input$umap_colourby <- "gene_expression"
    # input$umap_subset_cells <- TRUE
    ## END TEST ##
    
    ## Fetch data ##
    
    if(input$umap_modality == "rna") {
      umap.dt <- fread(file.path(data_folder,sprintf("dimensionality_reduction/rna/%s/umap_cells.txt.gz",input$umap_stage))) %>% 
        setnames(c("X","Y","cell")) %>% setkey(cell)
    } else if(input$umap_modality == "atac") {
      umap.dt <- fread(file.path(data_folder,sprintf("dimensionality_reduction/atac/%s/umap_cells.txt.gz",input$umap_stage))) %>%
        setnames(c("cell","X","Y")) %>% setkey(cell)
    } else if(input$umap_modality == "rna_atac") {
      umap.dt <- fread(file.path(data_folder,sprintf("dimensionality_reduction/rna_atac/%s/umap_cells.txt.gz",input$umap_stage))) %>% 
        setnames(c("cell","X","Y")) %>% setkey(cell)
    }
    
    to.plot <- umap.dt %>% 
      merge(cell_metadata.dt[sample!="E8.5_CRISPR_T_KO",c("cell","celltype","stage","sample")], by=c("cell"))
    
    if (input$umap_colourby == "gene_expression") {
      tmp <- data.table(
        cell = colnames(rna_expr_cells.array),
        color = as.numeric(rna_expr_cells.array[input$umap_gene_rna,])
      )
      to.plot <- to.plot %>% merge(tmp,by="cell") %>% setorder(color)
    } else {
      to.plot$color <- to.plot[[input$umap_colourby]]
    }
    
    ## Plot PAGA ##
    
    celltypes <- sapply(net.paga$val,"[[","vertex.names")
    alphas <- rep(1.0,length(celltypes)); names(alphas) <- celltypes
    sizes <- rep(8,length(celltypes)); names(sizes) <- celltypes
    
    # p1 <- ggnet2(
    p.paga <- ggnet2_interactive(
      net = net.paga,
      mode = c("x", "y"),
      color = celltype_colours[celltypes],
      node.alpha = alphas,
      node.size = sizes,    
      edge.size = 0.15,
      edge.color = "grey",
      label = TRUE,
      label.size = 2.5
    )
    
    ## Plot UMAP ##
    
    # Subset number of cells
    if (input$umap_subset_cells) {
      to.plot <- to.plot[sample(1:.N, size = min(2e4,.N))]
      dot_size <- 1
    } else {
      dot_size <- 0.5
      
    }
    
    p.umap <- ggplot(to.plot, aes(x = X, y = Y, color = color)) +
      geom_point_interactive(aes(tooltip = celltype, data_id = celltype), size=dot_size, alpha=0.9) +
      coord_fixed(ratio = 0.9) +
      theme_classic() +
      ggplot_theme_NoAxes()
    
    if (input$umap_colourby==c("celltype")) {
      p.umap <- p.umap + theme(legend.position = "none")
    }
      
    # Modify legends    
    if (input$umap_colourby%in%c("celltype","stage","sample")) {
      p.umap <- p.umap + 
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))
    } else {
      p.umap <- p.umap +
        theme(
          legend.title = element_blank()
        )
    }
    
    # Define palette
    palette <- switch(input$umap_colourby, 
                      "celltype" = celltype_palette, 
                      "stage" = stage_palette, 
                      "sample" = sample_palette, 
                      "gene_expression" = rna_palette, 
                      "chromatin_accessibility" = atac_palette
    )
    p.umap <- p.umap + palette
    
    ## Barplot/Violin plot with statistics ##
    
    if (input$umap_colourby=="celltype") {
      # tab <- table(to.plot$celltype, droplevels(to.plot$stage))
      tab <- table(cell_metadata.dt$celltype, cell_metadata.dt$stage)
      fractions <- sweep(tab, 1, rowSums(tab), "/")
      means <- apply(fractions, 1, function(x) sum(x * 1:length(x)))
      to.plot.barplot <- fractions %>% as.data.table(keep.rownames = T) %>% setnames(c("celltype","stage","value"))
      
      p.barplot <- ggplot(to.plot.barplot, aes(x = factor(celltype, levels = names(means)[order(means)]), y = value, fill = stage)) +
        # geom_bar(stat = "identity") +
        geom_bar_interactive(aes(tooltip=celltype, data_id=celltype), stat = "identity") +
        labs(y = "Fraction of cells") +
        stage_palette_fill +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=rel(0.80)),
          axis.text.y = element_text(color="black", size=rel(0.8)),
          axis.title.x = element_blank()
        )
    } else if (input$umap_colourby%in%c("stage","sample")) {
      
      to.plot.barplot <- to.plot[,.N,by=c(input$umap_colourby,"celltype")]
      
      p.barplot <- ggplot(to.plot.barplot, aes_string(x = input$umap_colourby, y = "N", fill = "celltype")) +
        # geom_bar(stat = "identity") +
        geom_bar_interactive(aes(tooltip=celltype, data_id=celltype), stat = "identity") +
        labs(y = "Number of cells") +
        celltype_palette_fill +
        theme_classic() +
        theme(
          legend.position = "none",
          # axis.text.x = element_text(size=rel(1.25), angle = 90, vjust = 0.5, hjust=1, color="black"),
          axis.text.x = element_text(size=rel(1.25), color="black"),
          axis.text.y = element_text(color="black", size=rel(0.8)),
          axis.title.x = element_blank()
        )
      
    } else if (input$umap_colourby%in%c("gene_expression")) {
      
      clust.sizes <- table(to.plot$celltype)
      
      p.barplot <- ggplot(to.plot, aes(x = celltype, y = color, fill = celltype)) +
        geom_violin_interactive(aes(tooltip=celltype, data_id=celltype), scale = "width", alpha=0.8) +
        geom_boxplot_interactive(aes(tooltip=celltype, data_id=celltype), width=0.5, outlier.shape=NA, alpha=0.8) +
        labs(y = "log2 normalised counts") +
        annotate(
          geom = "text",
          x = names(clust.sizes),
          y = rep_len(c(max(to.plot$color)*1.1, max(to.plot$color) * 1.2), length.out = length(clust.sizes)),
          label = as.vector(clust.sizes),
          size = 3
        ) +
        celltype_palette_fill +
        theme_classic() +
        theme(
          axis.title = element_text(size = 11, color="black"),
          axis.text.y = element_text(size = rel(1.0), color="black"),
          axis.text.x = element_text(size = rel(1.0), color="black", angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none",
          axis.title.x = element_blank()
        )
    }
    
    layout <- "
    ABB
    ABB
    ABB
    ABB
    CCC
    "
    girafe(
      # code = print((p.paga+p.umap)/p2 + plot_layout(widths = c(1,10), heights=c(4,1))),
      code = print(p.paga+p.umap+p.barplot + plot_layout(design = layout)),
      width_svg = 14, height_svg = 9,
      options = list( 
        opts_sizing(rescale = FALSE),
        # opts_selection(type = "single", css = "cursor:pointer;fill:magenta;color:magenta"),
        opts_selection(type = "single", css = ""),
        # opts_hover_inv(css = "opacity:0.45;"),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
    
  })
  
  output$umap = renderGirafe({
    shiny::validate(need(input$umap_gene_rna%in%rna_genes, "" ))
    plot_UMAP()
  })
  
  ########################################################
  ## gene expression (cells vs metacells vs pseudobulk) ##
  ########################################################

  updateSelectizeInput(session = session, inputId = 'gene_expr_gene', choices = rna_genes, server = TRUE, selected = "T")
  
  plot_gene_expr <- reactive({
    
    ## START TEST ##
    # input$gene_expr_celltypes <- celltypes
    # input$gene_expr_gene <- "T"
    ## END TEST ##
    
    # Plot single cells 
    cells_expr.dt <- data.table(
      expr = as.numeric(rna_expr_cells.array[input$gene_expr_gene,]),
      celltype = cell_metadata.dt[colnames(rna_expr_cells.array),celltype],
      class = "cells"
    ) %>% .[celltype%in%input$gene_expr_celltypes] %>% .[,celltype:=factor(celltype,levels=input$gene_expr_celltypes)]
    
    p.cells <- ggplot(cells_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
      geom_violin_interactive(aes(tooltip=celltype, data_id=celltype), scale="width", alpha=0.75) +
      geom_boxplot_interactive(aes(tooltip=celltype, data_id=celltype), width=0.5, outlier.shape=NA, alpha=0.75) +
      stat_summary(fun.data = function(x) { return(c(y = max(cells_expr.dt$expr)+0.5, label = length(x)))}, geom = "text", size=2) +
      scale_fill_manual(values=celltype_colours[input$gene_expr_celltypes]) +
      labs(x="",y="", title=sprintf("%s expression (cells)",input$gene_expr_gene)) +
      guides(x = guide_axis(angle = 90)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    # Plot metacells
    metacells_expr.dt <- data.table(
      expr = as.numeric(rna_expr_metacells.array[input$gene_expr_gene,]),
      celltype = metacell_metadata.dt[colnames(rna_expr_metacells.array),celltype],
      class = "metacells"
    ) %>% .[celltype%in%input$gene_expr_celltypes] %>% .[,celltype:=factor(celltype,levels=input$gene_expr_celltypes)]
    
    p.metacells <- ggplot(metacells_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
      geom_violin(scale="width", alpha=0.75) +
      geom_boxplot_interactive(aes(tooltip=celltype, data_id=celltype), width=0.5, outlier.shape=NA, alpha=0.75) +
      stat_summary(fun.data = function(x) { return(c(y = max(metacells_expr.dt$expr)+1, label = length(x)))}, geom = "text", size=2.5) +
      scale_fill_manual(values=celltype_colours[input$gene_expr_celltypes]) +
      labs(x="",y="", title=sprintf("%s expression (metacells)",input$gene_expr_gene)) +
      guides(x = guide_axis(angle = 90)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    pseudobulk_expr.dt <- data.table(
      expr = rna_expr_pseudobulk_replicates.mtx[input$gene_expr_gene,],
      sample = colnames(rna_expr_pseudobulk_replicates.mtx),
      class = "pseudobulk"
    ) %>% 
      .[,celltype:=strsplit(sample,"-") %>% map_chr(1)] %>%
      .[celltype%in%input$gene_expr_celltypes] %>% .[,celltype:=factor(celltype,levels=input$gene_expr_celltypes)]
    
    tmp <- pseudobulk_expr.dt[,.(expr=mean(expr), sd=sd(expr)), by=c("celltype")]
    
    p.pseudobulk <- ggplot(pseudobulk_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
      geom_bar_interactive(aes(tooltip=celltype, data_id=celltype), stat="identity", color="black", alpha=0.9, data=tmp) +
      geom_jitter(size=1.5, alpha=0.9, width=0.08, shape=21) +
      geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), width=0.25, alpha=1, size=0.6, data=tmp) +
      stat_summary(fun.data = function(x) { return(c(y = max(pseudobulk_expr.dt$expr)+1, label = length(x)))}, geom = "text", size=2.5) +
      scale_fill_manual(values=celltype_colours[input$gene_expr_celltypes]) +
      labs(x="",y="", title=sprintf("%s expression (pseudobulk with bootstrapped replicates)",input$gene_expr_gene)) +
      guides(x = guide_axis(angle = 90)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(colour="black",size=rel(0.9)),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    girafe(
      code = print(p.cells+p.metacells+p.pseudobulk + plot_layout(nrow=3)),
      width_svg = 8, height_svg = 9,
      options = list(
        # opts_zoom(min = 1, max = 3),
        # opts_sizing(rescale = FALSE),
        opts_selection(type = "single", css = ""),
        opts_hover_inv(css = "opacity:0.25"),
        opts_hover(css = "")
      )
    ) %>% return(.)
    
  })
  
  
  output$gene_expr = renderGirafe({
    shiny::validate(need(input$gene_expr_gene%in%rna_genes, "" ))
    plot_gene_expr()
  })
  
  #############
  ## Mapping ##
  #############
  
  plot_mapping <- reactive({
    
    ## START TEST ##
    # input$mapping_sample <- "E8.5_rep1"
    ## END TEST ##
  
    # Fetch data    
    cell_metadata_sample.dt <- cell_metadata.dt[sample==input$mapping_sample] %>% 
      .[,celltype:=stringr::str_replace_all(celltype,c("/"=" ","-"=" ","_"=" "))]
      
    # Define dot size  
    size.values <- c(0.30, 0.30)
    names(size.values) <- c(input$mapping_sample, "Atlas")
    
    # Define dot alpha  
    alpha.values <- c(0.70, 0.70)
    names(alpha.values) <- c(input$mapping_sample, "Atlas")
    
    # Define dot colours  
    colour.values <- c("red", "lightgrey")
    names(colour.values) <- c(input$mapping_sample, "Atlas")
    
    # Subset atlas for faster plotting
    if (input$mapping_subset_cells) {
      umap_reference_subset.dt <- rbind(
        umap_reference.dt[cell%in%unique(cell_metadata_sample.dt$closest.cell)],
        umap_reference.dt[!cell%in%unique(cell_metadata_sample.dt$closest.cell),.SD[sample.int(n=.N, size=1e4)]]
      )
    } else {
      umap_reference_subset.dt <- rbind(
        umap_reference.dt[cell%in%unique(cell_metadata_sample.dt$closest.cell)],
        umap_reference.dt[!cell%in%unique(cell_metadata_sample.dt$closest.cell),.SD[sample.int(n=.N, size=5e4)]]
      )
    }
    
    to.plot <- umap_reference_subset.dt %>% 
      .[,index:=match(cell, cell_metadata_sample.dt$closest.cell)] %>% 
      .[,mapped:=as.factor(!is.na(index))] %>% 
      .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",input$mapping_sample))] %>%
      setorder(mapped) 
    
    # p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas", rasterise = TRUE, subset = TRUE, legend = FALSE)
    
    p.mapping <- ggplot(to.plot, aes(x=V1, y=V2)) +
      geom_point_interactive(aes(size=mapped, alpha=mapped, colour=mapped, tooltip=celltype, data_id=celltype)) +
      scale_size_manual(values = size.values) +
      scale_alpha_manual(values = alpha.values) +
      scale_colour_manual(values = colour.values) +
      guides(colour = guide_legend(override.aes = list(size=6))) +
      theme_classic() +
      ggplot_theme_NoAxes() +
      theme(
        legend.position = "top",
        legend.title = element_blank()
      )
    
    # Plot celltype proportions
    celltype_proportions.dt <- cell_metadata_sample.dt %>%
      .[,N:=.N,by="sample"] %>%
      .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","stage","celltype")]
    
    # tmp <- celltype_colours[names(celltype_colours) %in% unique(celltype_proportions.dt$celltype)]
    tmp <- celltype_colours; names(tmp) <- stringr::str_replace_all(names(tmp),c("/"=" ","-"=" ","_"=" "))
    celltype_proportions.dt[,celltype:=factor(celltype, levels=rev(names(tmp)))]
    
    p.celltype_proportions <- ggplot(celltype_proportions.dt, aes(x=celltype, y=N)) +
      geom_bar_interactive(aes(tooltip=celltype, data_id=celltype, fill=celltype), stat="identity", color="black") +
      scale_fill_manual(values=tmp, drop=FALSE) +
      scale_x_discrete(drop=FALSE) +
      coord_flip() +
      labs(y="Number of cells") +
      theme_bw() +
      theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color="black", size=rel(0.9)),
        axis.title.x = element_text(color="black", size=rel(0.9)),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=rel(0.9), color="black"),
        axis.text.x = element_text(size=rel(1), color="black")
      )
    
    girafe(
      code = print(p.mapping+p.celltype_proportions + plot_layout(nrow=1)),
      width_svg = 10, height_svg = 6,
      options = list(
        # opts_zoom(min = 1, max = 3),
        # opts_sizing(rescale = FALSE),
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer;fill:magenta;stroke:magenta;color:magenta;fill-opacity:1")
      )
    ) %>% return(.)
    
    
  })
  
  output$mapping_plot = renderGirafe({
    shiny::validate(need(input$mapping_sample%in%samples,""))
    plot_mapping()
  })
  
  ################################################################
  ## chromatin accessibility (cells vs metacells vs pseudobulk) ##
  ################################################################
  
  observeEvent(input$chr_acc_feature_chr, {
    tmp <- grep(sprintf("%s:",input$chr_acc_feature_chr),atac_peaks, value=T)
    updateSelectizeInput(session = session, "chr_acc_feature", choices = tmp, server = TRUE, selected = tmp[1])
  })
  
  observeEvent(input$chr_acc_feature_type, {
    if (input$chr_acc_feature_type=="peaks") {
      tmp <- grep(sprintf("%s:",input$chr_acc_feature_chr),atac_peaks,value=T)
      updateSelectizeInput(session = session, inputId = 'chr_acc_feature', choices = tmp, server = TRUE, selected = tmp[1])
    } else if (input$chr_acc_feature_type=="tss") {
      updateSelectizeInput(session = session, inputId = 'chr_acc_feature', choices = atac_genes, server = TRUE, selected = "Foxa2")
    }
  })
  
  plot_chromatin_accessibility <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$chr_acc_celltypes <- celltypes
    # input$chr_acc_feature_type <- "tss"
    # input$chr_acc_feature_chr <- "chr18"
    # input$chr_acc_feature <- "chr18:32542271-32542871"
    # input$chr_acc_feature <- "Foxa2"
    ## END TEST ##
    
    if (input$chr_acc_feature_type=="peaks") {
      atac_cells.array <- atac_peaks_cells.array
      atac_metacells.array <- atac_peaks_metacells.array
      atac_pseudobulk.array <- atac_peaks_pseudobulk.mtx
      atac_pseudobulk_replicates.array <- atac_peaks_pseudobulk_replicates.mtx
    } else if (input$chr_acc_feature_type=="tss") {
      atac_metacells.array <- atac_genes_metacells.array
      atac_cells.array <- atac_genes_cells.array
      atac_pseudobulk.array <- atac_genes_pseudobulk.mtx
      atac_pseudobulk_replicates.array <- atac_genes_pseudobulk_replicates.mtx
    }
    
    # Plot cells
    cells_acc.dt <- data.table(
      acc = as.numeric(atac_cells.array[input$chr_acc_feature,]>0),
      celltype = cell_metadata.dt[colnames(atac_cells.array),celltype]
    ) %>% .[celltype%in%input$chr_acc_celltypes] %>% 
      .[,.(acc=mean(acc)), by="celltype"] %>%
      .[,celltype:=factor(celltype,levels=input$chr_acc_celltypes)]
    
    p.cells <- ggplot(cells_acc.dt, aes(x=celltype, y=acc, fill=celltype)) +
      geom_bar_interactive(aes(tooltip=celltype, data_id=celltype), stat="identity", color="black") +
      scale_fill_manual(values=celltype_colours[input$chr_acc_celltypes]) +
      labs(x="",y="", title=sprintf("%s accessibility (cells)",input$chr_acc_feature)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    # Plot metacells
    metacells_acc.dt <- data.table(
      acc = as.numeric(atac_metacells.array[input$chr_acc_feature,]),
      celltype = metacell_metadata.dt[colnames(atac_metacells.array),celltype],
      class = "metacells"
    ) %>% .[celltype%in%input$chr_acc_celltypes] %>% .[,celltype:=factor(celltype,levels=input$chr_acc_celltypes)]
    
    p.metacells <- ggplot(metacells_acc.dt, aes(x=celltype, y=acc, fill=celltype, tooltip=celltype, data_id=celltype)) +
      geom_violin_interactive(scale="width", alpha=0.75) +
      geom_boxplot_interactive(width=0.5, outlier.shape=NA, alpha=0.75) +
      stat_summary(fun.data = function(x) { return(c(y = max(metacells_acc.dt$acc)+1, label = length(x)))}, geom = "text", size=2.5) +
      scale_fill_manual(values=celltype_colours[input$chr_acc_celltypes]) +
      labs(x="",y="", title=sprintf("%s accessibility (metacells)",input$chr_acc_feature)) +
      guides(x = guide_axis(angle = 90)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.75)),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )

    # Plot pseudobulk
    pseudobulk_acc.dt <- data.table(
      acc = atac_pseudobulk_replicates.array[input$chr_acc_feature,],
      sample = colnames(atac_pseudobulk_replicates.array),
      class = "pseudobulk"
    ) %>% 
      .[,celltype:=strsplit(sample,"-") %>% map_chr(1)] %>%
      .[celltype%in%input$chr_acc_celltypes] %>% .[,celltype:=factor(celltype,levels=input$chr_acc_celltypes)]
    
    tmp <- pseudobulk_acc.dt[,.(acc=mean(acc), sd=sd(acc)), by="celltype"]
    
    p.pseudobulk <- ggplot(pseudobulk_acc.dt, aes(x=celltype, y=acc, fill=celltype)) +
      geom_bar_interactive(aes(tooltip=celltype, data_id=celltype), stat="identity", color="black", alpha=0.9, data=tmp) +
      geom_jitter(size=1.5, alpha=0.9, width=0.08, shape=21) +
      geom_errorbar(aes(ymin=acc-sd, ymax=acc+sd), width=0.25, alpha=1, size=0.6, data=tmp) +
      stat_summary(fun.data = function(x) { return(c(y = max(pseudobulk_acc.dt$acc)+1, label = length(x)))}, geom = "text", size=2.5) +
      scale_fill_manual(values=celltype_colours[input$chr_acc_celltypes]) +
      labs(x="",y="", title=sprintf("%s accessibility (pseudobulk with bootstrapped replicates)",input$chr_acc_feature)) +
      guides(x = guide_axis(angle = 90)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(colour="black",size=rel(0.9)),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    girafe(
      code = print(p.cells+p.metacells+p.pseudobulk + plot_layout(nrow=3)),
      width_svg = 8, height_svg = 9,
      options = list(
        # opts_zoom(min = 1, max = 3),
        # opts_sizing(rescale = FALSE),
        opts_selection(type = "single", css = ""),
        opts_hover_inv(css = "opacity:0.25"),
        opts_hover(css = "")
      )
    ) %>% return(.)
    
  })
  
  
  output$chr_acc = renderGirafe({
    shiny::validate(need(input$chr_acc_feature%in%c(atac_peaks,atac_genes), "" ))
    shiny::validate(need(input$chr_acc_feature_chr%in%chr_mm10, "" ))
    plot_chromatin_accessibility()
  })
  

  #####################
  ## RNA vs chromVAR ##
  #####################
  
  # updateSelectizeInput(session = session, inputId = 'rna_vs_chromvar_tf', choices = TFs_all, server = TRUE, selected = "FOXA2")
  observeEvent(input$rna_vs_chromvar_motif_annotation, {
    updateSelectizeInput(session = session, "rna_vs_chromvar_tf", choices = TFs[[input$rna_vs_chromvar_motif_annotation]], selected = "FOXA2")
    motifs <- motif2gene.list[[input$rna_vs_chromvar_motif_annotation]][gene==input$rna_vs_chromvar_tf,motif]
    updateSelectizeInput(session = session, "rna_vs_chromvar_motif", choices = motifs, selected = motifs[1])
  })
  
  observeEvent(input$rna_vs_chromvar_tf, {
    motifs <- motif2gene.list[[input$rna_vs_chromvar_motif_annotation]][gene==input$rna_vs_chromvar_tf,motif]
    updateSelectizeInput(session = session, "rna_vs_chromvar_motif", choices = motifs, selected = motifs[1])
  })
  
  # Update TF selection upon click
  observeEvent(input$rna_vs_chromvar_plot_selected, {
    # updateTextInput(session = session, "rna_vs_chromvar_tf", value = input$rna_vs_chromvar_plot_selected)
    updateSelectizeInput(session = session, "rna_vs_chromvar_tf", selected = input$rna_vs_chromvar_plot_selected)
    motifs <- motif2gene.list[[input$rna_vs_chromvar_motif_annotation]][gene==input$rna_vs_chromvar_plot_selected,motif]
    updateSelectizeInput(session = session, "rna_vs_chromvar_motif", choices = motifs, selected = motifs[1])
  })
  
  
  
  plot_rna_vs_chromvar = reactive({
    
    ## START TEST ##
    # input$rna_vs_chromvar_tf <- "FOXA2"
    # input$rna_vs_chromvar_motif <- "FOXA2_573" # "FOXA2_355"
    # input$rna_vs_chromvar_motif_annotation <- "JASPAR"
    # input$chromvar_algorithm <- "chromvar_multiome_positive"
    # input$rna_vs_chromvar_cor_range <- c(-1,1)
    # input$rna_vs_chromvar_scatterplot_add_text <- TRUE
    ## END TEST ##
    
    # Fetch data
    if (input$chromvar_algorithm=="chromvar") {
      tmp <- "rna_vs_chromvar"
    } else if (input$chromvar_algorithm=="chromvar_multiome_positive") {
      tmp <- "rna_vs_chromvar_multiome"
    } else if (input$chromvar_algorithm=="chromvar_multiome_all") {
      print("Not implemented yet")
    }
    cor_rna_chromvar.dt <- fread(paste0(data_folder,sprintf("/%s/%s/cor_rna_vs_chromvar_pseudobulk.txt.gz",tmp,input$rna_vs_chromvar_motif_annotation))) %>% 
      .[,log_pval:=-log10(padj_fdr)]
    rna_chromvar.dt <- fread(paste0(data_folder,sprintf("/%s/%s/rna_vs_chromvar_pseudobulk.txt.gz",tmp,input$rna_vs_chromvar_motif_annotation))) %>%
      .[gene==input$rna_vs_chromvar_tf & motif==input$rna_vs_chromvar_motif]
    
    ## Volcano plot ##
    
    xlim_min <- input$rna_vs_chromvar_cor_range[1]
    xlim_max <- input$rna_vs_chromvar_cor_range[2]
    
    to.plot <- cor_rna_chromvar.dt %>% .[r>=input$rna_vs_chromvar_cor_range[1] & r<=input$rna_vs_chromvar_cor_range[2]]
    ylim_min <- min(to.plot$log_pval,na.rm=T)
    ylim_max <- max(to.plot$log_pval,na.rm=T)
    
    # negative_hits <- to.plot[sig==TRUE & r<0 & r>input$rna_vs_chromvar_cor_range[1],gene]
    # positive_hits <- to.plot[sig==TRUE & r>0 & r<input$rna_vs_chromvar_cor_range[2],gene]
    
    # Select dot size
    nhits <- nrow(to.plot)
    if (nhits>=500) { dot_size <- 3 } else if (nhits>100 & nhits<500) { dot_size <- 5 } else if (nhits<=100) { dot_size <- 8 }
    
    p.volcano <- ggplot(to.plot, aes(x=r, y=log_pval)) +
      geom_segment(x=0, xend=0, y=0, yend=ylim_max, color="orange", size=0.25) +
      geom_jitter_interactive(aes(fill=log_pval, alpha=log_pval, tooltip=gene, data_id=gene, onclick=gene), width=0.03, height=0.03, shape=21, size=dot_size) + 
      scale_fill_gradient(low = "gray80", high = "red") +
      scale_alpha_continuous(range=c(0.25,1)) +
      scale_x_continuous(limits=c(xlim_min-0.15,xlim_max+0.15)) +
      scale_y_continuous(limits=c(ylim_min,ylim_max+1)) +
      annotate("text", x=0, y=ylim_max+1, size=5, label=sprintf("(%d)", nhits)) +
      # annotate("text", x=xlim_min-0.05, y=ylim_max+2, size=4, label=sprintf("%d (-)",length(negative_hits))) +
      # annotate("text", x=xlim_max+0.05, y=ylim_max+2, size=4, label=sprintf("%d (+)",length(positive_hits))) +
      labs(x="Pearson correlation\n(RNA expression vs chromVAR z-score)", y=expression(paste("-log"[10],"(p.value)"))) +
      theme_classic() +
      theme(
        plot.margin = margin(t = 0, r = 25, b = 0, l = 0, unit = "pt"),
        axis.text = element_text(size=rel(1.1), color='black'),
        axis.title = element_text(size=rel(1.25), color='black'),
        legend.position="none"
      )
    
    if (cor_rna_chromvar.dt[gene==input$rna_vs_chromvar_tf,r]>0) { nudge_x = -0.55 } else { nudge_x = 0.55  }
    p.volcano <- p.volcano + 
      ggrepel::geom_text_repel(data=cor_rna_chromvar.dt[gene==input$rna_vs_chromvar_tf & motif==input$rna_vs_chromvar_motif], aes(x=r, y=log_pval, label=gene), size=6, 
                               nudge_x = nudge_x,
                               box.padding = 0.5,
                               nudge_y = 1,
                               segment.curvature = -0.1,
                               segment.ncp = 3,
                               segment.angle = 20,
                               arrow = arrow(length = unit(0.02, "npc"))
      )
    
    ## Scatterplot ##
    
    rna_chromvar.dt %>% .[,tooltip:=sprintf("%s\nRNA expression = %.2f\nchromVAR = %.2f",celltype,expr,chromvar_zscore)]
    
    p.scatter <- ggplot(rna_chromvar.dt, aes(x=chromvar_zscore, y=expr, fill=celltype)) +
      geom_point_interactive(aes(tooltip=tooltip, data_id = celltype), size=6, shape=21) +
      # stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
      celltype_palette_fill +
      labs(x="chromVAR z-score", y="Gene expression", title=input$rna_vs_chromvar_tf) +
      theme_classic() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=rel(1.5), color="black"),
        axis.title = element_text(size=rel(1.10), color="black"),
        axis.text = element_text(size=rel(1), color="black")
      )
    
    if (input$rna_vs_chromvar_scatterplot_add_text) {
      p.scatter <- p.scatter + ggrepel::geom_text_repel(aes(label=celltype), size=3, max.overlaps=Inf)
    }
    
    ## PAGA ##
    celltypes <- sapply(net.paga$val,"[[","vertex.names")
    alphas <- rep(0.9,length(celltypes)); names(alphas) <- celltypes
    sizes <- rep(6,length(celltypes)); names(sizes) <- celltypes
    p.paga <- ggnet2_interactive(
      net = net.paga,
      mode = c("x", "y"),
      color = celltype_colours[celltypes],
      node.alpha = alphas,
      node.size = sizes,    
      edge.size = 0.15,
      edge.color = "grey",
      label = TRUE,
      label.size = 2.5
    )
    
    # p.paga <- ggraph(igraph.paga.tbl, x = x, y = y) +
    #   geom_edge_link(edge_colour = "grey66", edge_alpha=0.8, edge_width=0.25) +
    #   geom_node_point(aes(color=celltype), size=6, alpha=0.8) +
    #   geom_node_text(aes(label=celltype), size=2.5) +
    #   scale_color_manual(values=celltype_colours) +
    #   theme_classic(base_size=14) +
    #   theme(
    #     axis.line = element_blank(), 
    #     axis.text = element_blank(),
    #     axis.ticks = element_blank(), 
    #     axis.title = element_blank(),
    #     legend.position = "none"
    #   )
    
    # Define cell type order
    cellype.order <- rownames(connectivity.mtx)
    
    # Plot expression
    expr.values <- rna_chromvar.dt[gene==input$rna_vs_chromvar_tf,expr]; names(expr.values) <- rna_chromvar.dt[gene==input$rna_vs_chromvar_tf,celltype]
    igraph.paga.tbl <- igraph.paga.tbl %>% activate(nodes) %>% mutate(expr=expr.values[cellype.order])
    
    p.paga.rna <- ggraph(igraph.paga.tbl, x = x, y = y) +
      geom_edge_link(edge_colour = "grey66", edge_alpha=1, edge_width=0.5) +
      geom_node_point(aes(fill=expr), size=7, shape=21) +
      scale_fill_gradient(low = "white", high = "darkgreen") +
      theme_classic(base_size=14) +
      theme(
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.position = "none"
      )
    
    # Plot chromvar
    chromvar.values <- rna_chromvar.dt$chromvar_zscore; names(chromvar.values) <- rna_chromvar.dt$celltype
    igraph.paga.tbl <- igraph.paga.tbl %>% activate(nodes) %>% mutate(chromvar=chromvar.values[cellype.order])
    
    p.paga.acc <- ggraph(igraph.paga.tbl, x = x, y = y) +
      geom_edge_link(edge_colour = "grey66", edge_alpha=1, edge_width=0.5) +
      geom_node_point(aes(fill=chromvar), size=7, shape=21) +
      scale_fill_gradient(low = "white", high = "purple") +
      theme_classic(base_size=14) +
      theme(
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.position = "none"
      )
    
    # htmlwidget call
     girafe(
      # ggobj = plot_grid(p1, p2, nrow=2, scale=0.95, rel_heights = c(2/5,3/5)), 
      code = print(((p.volcano)|(p.scatter/(p.paga+p.paga.rna+p.paga.acc))) + plot_layout(width=c(1,2))),
      width_svg = 14, height_svg = 9,
      options = list(
        opts_zoom(min = 1, max = 5),
        opts_sizing(rescale = TRUE),
        # opts_selection(type = "single", css = "cursor:pointer;r:1.5%"),
        # opts_hover_inv(css = "opacity:0.65;cursor:pointer;r:.7%"),
        # opts_hover_inv(css = "opacity:0.65;cursor:pointer;r:.7%"),
        # opts_hover(css = "cursor:pointer;r:1.5%")
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
  })
  
  
  output$rna_vs_chromvar_plot = renderGirafe({
    shiny::validate(need(input$rna_vs_chromvar_tf%in%TFs_all,""))
    return(plot_rna_vs_chromvar())
  })
  
  # output$rna_vs_chromvar_print <- renderPrint({
  #   if (length(input$rna_vs_chromvar_tf)>0) {
  #     tmp <- cor_rna_chromvar.list[[input$rna_vs_chromvar_motif_annotation]] %>%
  #       .[r>=input$rna_vs_chromvar_cor_range[1] & r<=input$rna_vs_chromvar_cor_range[2]] %>%
  #       .[,abs_r:=abs(r)] %>% setorder(-abs_r) %>% head(n=15)
  #     cat("Top 15 TFs with largest absolute RNA vs chromVAR+ correlation (within the selected range):\n")
  #     cat(sprintf("%s: r=%.2f\t", tmp$gene, round(tmp$r,2)))
  #   }
  # })
  
  output$rna_vs_chromvar_plot_motif <- renderPlot({
    shiny::validate(need(input$rna_vs_chromvar_tf%in%TFs_all,""))
    shiny::validate(need(input$rna_vs_chromvar_motif%in%motifs_all,""))
    # motif <- motif2gene.list[[input$rna_vs_chromvar_motif_annotation]][gene==input$rna_vs_chromvar_tf & motif==input$rna_vs_chromvar_motif,motif]
    tmp <- tf_motifs_pwm.list[[input$rna_vs_chromvar_motif_annotation]][["motifs"]][[input$rna_vs_chromvar_motif]]
    if (input$rna_vs_chromvar_motif_annotation=="JASPAR") {
      m <- toPWM(tmp, type="prob") %>% as.matrix
    } else if (input$rna_vs_chromvar_motif_annotation=="CISBP") {
      m <- (2**as.matrix(tmp))*0.25  # this is not 100% accurate..
    }
    
    ggseqlogo(m) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
      )
  }, bg="transparent", execOnResize = TRUE)
  
  
  ###########################
  ## Differential analysis ##
  ###########################
  
  observeEvent(input$diff_modality, {
    if (input$diff_modality=="RNA") {
      updateSelectizeInput(session = session, inputId = 'diff_feature', choices = rna_genes, server = TRUE, selected = "Foxa2")
    } else if (input$diff_modality=="ATAC") {
      tmp <- grep(sprintf("%s:",input$diff_atac_chr),atac_peaks,value=T)
      updateSelectizeInput(session = session, inputId ='diff_feature', choices = tmp, server = TRUE, selected = tmp[1])
    }
  })
  observeEvent(input$diff_atac_chr, {
    if (input$diff_modality=="ATAC") {
      tmp <- grep(sprintf("%s:",input$diff_atac_chr),atac_peaks,value=T)
      updateSelectizeInput(session = session, inputId = 'diff_feature', choices = tmp, server = TRUE, selected = tmp[1])
    }
  })
  
  # Update feature upon click
  observeEvent(input$diff_plot_selected, {
    if (input$diff_modality=="RNA") {
      updateSelectizeInput(session = session, inputId = 'diff_feature', choices = rna_genes, server = TRUE, selected = input$diff_plot_selected)
    } else if (input$diff_modality=="ATAC") {
      updateSelectizeInput(session = session, inputId ='diff_atac_chr', choices = chr_mm10, server = TRUE, selected = strsplit(input$diff_plot_selected,":") %>% map_chr(1))
      updateSelectizeInput(session = session, inputId ='diff_feature', choices = tmp, server = TRUE, selected = input$diff_plot_selected)
    }
  })


  plot_differential <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$diff_celltypeA <- "Gut"
    # input$diff_celltypeB <- "Neural_crest"
    # input$diff_modality <- "ATAC"
    # input$diff_atac_chr <- "chr18"
    # input$diff_feature <- "chr18:32542271-32542871"
    # input$diff_resolution <- "Cells"
    # input$diff_range <- c(-8,8)
    # input$min_log_pval <- 25
    ## END TEST ##
    
    ## Volcano ##

    # Fetch data
    if (input$diff_modality=="RNA") {
      if (input$diff_resolution=="Cells") {
        diff.dt <- fread(file.path(data_folder,sprintf("differential/rna/cells/%s_vs_%s.txt.gz",input$diff_celltypeA,input$diff_celltypeB))) %>%
          .[,c("groupA_N","groupB_N"):=NULL] %>% setnames("gene","feature")
      } else if (input$diff_resolution=="Metacells") {
        diff.dt <- fread(file.path(data_folder,sprintf("differential/rna/metacells/%s_vs_%s.txt.gz",input$diff_celltypeA,input$diff_celltypeB))) %>%
          .[,c("groupA_N","groupB_N"):=NULL] %>% setnames("gene","feature")
      } else if (input$diff_resolution=="Pseudobulk") {
        diff.dt <- fread(file.path(data_folder,sprintf("differential/rna/pseudobulk/%s_vs_%s.txt.gz",input$diff_celltypeA,input$diff_celltypeB))) %>% 
          setnames("gene","feature")
      }
      diff.dt <- diff.dt[feature%in%rna_genes]
    } else if (input$diff_modality=="ATAC") {
      if (input$diff_resolution=="Cells") {
        diff.dt <- fread(file.path(data_folder,sprintf("differential/atac/cells/%s_vs_%s.txt.gz",input$diff_celltypeA,input$diff_celltypeB))) %>%
          .[,c("groupA_N","groupB_N","AUC"):=NULL] %>% setnames(c("feature","logFC","padj_fdr"))
      } else if (input$diff_resolution=="Metacells") {
        diff.dt <- fread(file.path(data_folder,sprintf("differential/atac/metacells/%s_vs_%s.txt.gz",input$diff_celltypeA,input$diff_celltypeB))) %>%
          .[,c("groupA_N","groupB_N","mean_groupA","mean_groupB"):=NULL]
      } else if (input$diff_resolution=="Pseudobulk") {
        diff.dt <- fread(file.path(data_folder,sprintf("differential/atac/pseudobulk/%s_vs_%s.txt.gz",input$diff_celltypeA,input$diff_celltypeB))) %>%
          .[,c("groupA_N","groupB_N","mean_groupA","mean_groupB"):=NULL]
      }
      diff.dt <- diff.dt[feature%in%atac_peaks]
    }
    
    if (!input$diff_feature%in%diff.dt$feature) {
      diff.dt <- rbind(diff.dt, data.table(feature=input$diff_feature, logFC=0, padj_fdr=1))
    }
    
    to.plot <- diff.dt %>% 
      .[logFC>=input$diff_range[1] & logFC<=input$diff_range[2]] %>% 
      .[is.na(logFC),logFC:=0] %>%
      .[is.na(padj_fdr),padj_fdr:=1] %>%
      .[,log_pval:=-log10(padj_fdr+1e-150)]
      
    
    # xlim_min <- min(to.plot$logFC); xlim_max <- max(to.plot$logFC)
    if (input$diff_modality=="ATAC" & input$diff_resolution=="Cells") {
      xlim_min <- -1; xlim_max <- 1; margin_text <- 0.10; margin_dots <- 0.05
    } else {
      xlim_min <- -8; xlim_max <- 8; margin_text <- 2; margin_dots <- 0.15
    }
    to.plot[logFC<=xlim_min,logFC:=xlim_min]; to.plot[logFC>=xlim_max,logFC:=xlim_max]
    
    # negative_hits <- to.plot[sig==TRUE & r<0 & r>input$rna_vs_chromvar_cor_range[1],gene]
    # positive_hits <- to.plot[sig==TRUE & r>0 & r<input$rna_vs_chromvar_cor_range[2],gene]
    
    # to.plot.subset <- rbind(
    #   to.plot[log_pval>=10],
    #   to.plot[log_pval<=10][sample.int(.N,5e3)]
    # )
    to.plot.subset <- to.plot[log_pval>=input$min_log_pval]
    if (!input$diff_feature%in%to.plot.subset$feature) {
      to.plot.subset <- rbind(to.plot.subset,to.plot[feature==input$diff_feature])
    }
    
    ylim_min <- min(to.plot.subset$log_pval,na.rm=T); ylim_max <- max(to.plot$log_pval,na.rm=T)
    
    p.volcano <- ggplot(to.plot.subset, aes(x=logFC, y=log_pval)) +
      geom_segment(x=0, xend=0, y=input$min_log_pval, yend=ylim_max, color="orange", size=0.25) +
      geom_jitter_interactive(aes(fill=log_pval, tooltip=feature, data_id=feature, alpha=log_pval, size=log_pval, onclick=feature), width=0.03, height=0.03, shape=21) + 
      geom_point_interactive(aes(tooltip=feature, data_id=feature, onclick=feature), size=7, fill="green", shape=21, data=to.plot.subset[feature==input$diff_feature]) + 
      geom_hline(yintercept=input$min_log_pval, linetype="dashed") +
      ggrepel::geom_text_repel(aes(label=feature), size=ifelse(input$diff_modality=="RNA",5,4), data=to.plot[feature==input$diff_feature]) +
      scale_fill_gradient(low = "gray80", high = "red") +
      scale_alpha_continuous(range=c(0.25,1)) +
      scale_size_continuous(range=c(0.15,3.5)) +
      scale_x_continuous(limits=c(xlim_min-margin_dots,xlim_max+margin_dots)) +
      scale_y_continuous(limits=c(ylim_min,ylim_max+3)) +
      annotate("text", x=0, y=ylim_max+3, size=5, label=sprintf("(%d)", nrow(to.plot.subset))) +
      # coord_cartesian(xlim=c(xlim_min,xlim_max)) +
      # annotate("text", x=xlim_min-0.05, y=ylim_max+2, size=4, label=sprintf("%d (-)",length(negative_hits))) +
      # annotate("text", x=xlim_max+0.05, y=ylim_max+2, size=4, label=sprintf("%d (+)",length(positive_hits))) +
      annotate("text", x=xlim_min+margin_text, y=1, size=5, label=sprintf("Higher in %s",input$diff_celltypeA)) +
      annotate("text", x=xlim_max-margin_text, y=1, size=5, label=sprintf("Higher in %s",input$diff_celltypeB)) +
      labs(x=ifelse(input$diff_modality=="RNA","Differential expression","Differential accessibility"), y=expression(paste("-log"[10],"(p.value)"))) +
      theme_classic() +
      theme(
        axis.text = element_text(size=rel(1.25), color='black'),
        axis.title = element_text(size=rel(1.25), color='black'),
        legend.position="none"
      )
    
    ## Plot expression/accessibility values ##
    
    # Fetch data
    if (input$diff_modality=="RNA") {
      if (input$diff_resolution=="Cells") {
        expr.dt <- data.table(
          expr = as.numeric(rna_expr_cells.array[input$diff_feature,]),
          celltype = cell_metadata.dt[colnames(rna_expr_cells.array),celltype]
        )
      } else if (input$diff_resolution=="Metacells") {
        expr.dt <- data.table(
          expr = as.numeric(rna_expr_metacells.array[input$diff_feature,]),
          celltype = metacell_metadata.dt[colnames(rna_expr_metacells.array),celltype]
        )
      } else if (input$diff_resolution=="Pseudobulk") {
        expr.dt <- data.table(
          expr = as.numeric(rna_expr_pseudobulk_replicates.mtx[input$diff_feature,]),
          sample = colnames(rna_expr_pseudobulk_replicates.mtx)
        ) %>% .[,celltype:=strsplit(sample,"-") %>% map_chr(1)] 
      }
    } else if (input$diff_modality=="ATAC") {
      if (input$diff_resolution=="Cells") {
        expr.dt <- data.table(
          expr = as.numeric(atac_peaks_cells.array[input$diff_feature,]),
          celltype = cell_metadata.dt[colnames(atac_peaks_cells.array),celltype]
        )
      } else if (input$diff_resolution=="Metacells") {
        expr.dt <- data.table(
          expr = as.numeric(atac_peaks_metacells.array[input$diff_feature,]),
          celltype = metacell_metadata.dt[colnames(atac_peaks_metacells.array),celltype]
        )
      } else if (input$diff_resolution=="Pseudobulk") {
        expr.dt <- data.table(
          expr = as.numeric(atac_peaks_pseudobulk_replicates.mtx[input$diff_feature,]),
          sample = colnames(atac_peaks_pseudobulk_replicates.mtx)
        ) %>% .[,celltype:=strsplit(sample,"-") %>% map_chr(1)]
      }
    }
    
    expr.dt <- expr.dt %>% 
      .[celltype%in%c(input$diff_celltypeA,input$diff_celltypeB)] %>% 
      .[,celltype:=factor(celltype,levels=c(input$diff_celltypeA,input$diff_celltypeB))]
    
      
    if (input$diff_resolution%in%c("Cells","Metacells")) {
      
      p.expr <- ggplot(expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
        geom_jitter(size=2, width=0.05, alpha=0.5, shape=21) +
        geom_violin(scale="width", alpha=0.40) +
        geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.70) +
        # stat_summary(fun.data = function(x) { return(c(y = max(expr.dt$expr)+0.5, label = sprintf("N=%s",length(x))))}, geom = "text", size=4) +
        stat_summary(fun.data = function(x) { return(c(y = max(expr.dt$expr)+0.5, label = length(x)))}, geom = "text", size=5) +
        scale_fill_manual(values=celltype_colours[c(input$diff_celltypeA,input$diff_celltypeB)]) +
        labs(x="", y=sprintf("%s %s",input$diff_feature,ifelse(input$diff_modality=="RNA","expression","accessibility")), title=input$diff_feature) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust=0.5, size=rel(1.25)),
          axis.text.x = element_text(colour="black",size=rel(1.25)),
          axis.text.y = element_text(colour="black",size=rel(1.25)),
          axis.title.y = element_text(colour="black",size=rel(1.25)),
          axis.ticks.x = element_blank(),
          legend.position = "none"
        )
    } else if (input$diff_resolution=="Pseudobulk") {
      
      tmp <- expr.dt[,.(expr=mean(expr), sd=sd(expr)), by="celltype"]
      
      p.expr <- ggplot(expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
        geom_bar(stat="identity", color="black", alpha=0.9, data=tmp) +
        geom_jitter(size=2.5, alpha=0.9, width=0.08, shape=21) +
        geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), width=0.25, alpha=1, size=0.6, data=tmp) +
        stat_summary(fun.data = function(x) { return(c(y = max(expr.dt$expr)+0.5, label = length(x)))}, geom = "text", size=5) +
        scale_fill_manual(values=celltype_colours[c(input$diff_celltypeA,input$diff_celltypeB)]) +
        labs(x="",y=sprintf("%s %s",input$diff_feature,ifelse(input$diff_modality=="RNA","expression","accessibility")), title=input$diff_feature) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust=0.5, size=rel(1.25)),
          axis.text.x = element_text(colour="black",size=rel(1.25)),
          axis.text.y = element_text(colour="black",size=rel(1.25)),
          axis.title.y = element_text(colour="black",size=rel(1.25)),
          axis.ticks.x = element_blank(),
          legend.position = "none"
        )
    }
    
    # girafe call
    girafe(
      code = print(p.volcano + p.expr + plot_layout(ncol=2, widths=c(2,1))),
      width_svg = 14, height_svg = 9,
      options = list(
        opts_zoom(min = 1, max = 5),
        opts_sizing(rescale = TRUE),
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer")
      )
    ) %>% return(.)
    
  })
  
  output$diff_plot <- renderGirafe({
    shiny::validate(need(input$diff_celltypeA%in%celltypes,"Please select celltype A"))
    shiny::validate(need(input$diff_celltypeB%in%celltypes,"Please select celltype B"))
    shiny::validate(need(input$diff_celltypeA!=input$diff_celltypeB,"Celltype A and B must be different"))
    shiny::validate(need(input$diff_modality%in%c("RNA","ATAC"),"Please select data modality"))
    shiny::validate(need(input$diff_resolution%in%c("Cells","Metacells","Pseudobulk"),"Please select data resolution (cell, metacell, pseudobulk)"))
    if (input$diff_modality=="ATAC") {
      shiny::validate(need(input$diff_atac_chr%in%chr_mm10,"Please select chromosome"))
      shiny::validate(need(input$diff_feature%in%atac_peaks),"Please select a feature from our annotation")
    } else if (input$diff_modality=="RNA") {
      shiny::validate(need(input$diff_feature%in%rna_genes,"Please select a feature from our annotation"))
    }
    plot_differential()
  })
    
  ################
  ## TF markers ##
  ################
  
  plot_celltype_tf_markers <- reactive({
    
    ## START TEST ##
    # input$celltype_tf_markers_celltype <- "Gut"
    # input$celltype_tf_markers_range <- c(0,1)
    # input$celltype_tf_markers_tf <- "FOXA2"
    # input$celltype_tf_markers_motif_annotation <- "CISBP"
    # input$celltype_tf_markers_modality <- "RNA"
    ## END TEST ##
    
    # Fetch data
    if (input$celltype_tf_markers_modality=="rna") {
      tf_markers_celltype.dt <- tf_markers_rna.dt[celltype==input$celltype_tf_markers_celltype]
      tf_markers_gene.dt <- tf_markers_rna.dt[gene==input$celltype_tf_markers_tf]
    } else if (input$celltype_tf_markers_modality=="chromvar") {
      tf_markers_celltype.dt <- tf_markers_chromvar.dt[celltype==input$celltype_tf_markers_celltype]
      tf_markers_gene.dt <- tf_markers_chromvar.dt[gene==input$celltype_tf_markers_tf]
    }
    
    to.plot.text <- tf_markers_celltype.dt[tf_marker_score>=0.75]
    
    p.celltype <- ggplot(tf_markers_celltype.dt[tf_marker_score>0], aes(x=gene, y=tf_marker_score)) +
      geom_point_interactive(aes(size=tf_marker_score, tooltip=gene, data_id=gene), shape=21, color="black", fill="gray95", alpha=0.9) +
      ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=3, max.overlaps=Inf) +
      # geom_text(data=to.plot.text, aes(label=gene)) +
      scale_size_continuous(range = c(0.05,3.5)) +
      labs(x="Transcription Factor", y="TF marker score", title=input$celltype_tf_markers_celltype) +
      guides(size="none") +
      theme_classic() +
      theme(
        plot.title = element_text(size=rel(1.25), hjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color="black", size=rel(1)),
        axis.title = element_text(color="black", size=rel(1.25))
      )
    
    p.gene <- ggplot(tf_markers_gene.dt, aes(x=celltype, y=tf_marker_score)) +
      geom_bar_interactive(aes(fill=celltype, tooltip=celltype, data_id=celltype), stat = 'identity', color="black") +
      ggrepel::geom_text_repel(aes(label=gsub("_"," ",celltype)), size=5, max.overlaps=Inf, data=tf_markers_gene.dt[tf_marker_score>0.5]) +
      scale_fill_manual(values=celltype_colours) +
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by=0.25)) +
      coord_polar() + 
      theme_bw() +
      labs(x="", y="", title=input$celltype_tf_markers_tf) +
      theme(
        plot.title = element_text(size=rel(1.25), hjust=0.5),
        legend.position = "none",
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(size=rel(1.5)),
        panel.border = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line=element_blank()
      )
    
    # htmlwidget call
    girafe(
      code = print(p.celltype|p.gene) + plot_layout(width=c(1,1)),
      width_svg = 14, height_svg = 9,
      options = list(
        opts_zoom(min = 1, max = 5),
        opts_sizing(rescale = TRUE),
        # opts_selection(type = "single", css = "cursor:pointer;r:1.5%"),
        # opts_hover_inv(css = "opacity:0.65;cursor:pointer;r:.7%"),
        # opts_hover_inv(css = "opacity:0.65;cursor:pointer;r:.7%"),
        # opts_hover(css = "cursor:pointer;r:1.5%")
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
    
  })
  
  output$celltype_tf_markers_plot <- renderGirafe({
    shiny::validate(need(input$celltype_tf_markers_celltype%in%celltypes,"Please select a celltype"))
    shiny::validate(need(input$celltype_tf_markers_tf%in%TFs_all,"Please select a celltype"))
    # shiny::validate(need(input$celltype_tf_markers_motif_annotation%in%c("CISBP","JASPAR"),"Please select a motif annotation"))
    shiny::validate(need(input$celltype_tf_markers_modality%in%c("chromvar","rna"),"Please select a data resolution (cell, metacell, pseudobulk)"))
    plot_celltype_tf_markers()
  })
  
  # Update TF upon click
  observeEvent(input$celltype_tf_markers_plot_selected, {
    updateTextInput(session = session, "celltype_tf_markers_tf", value = input$celltype_tf_markers_plot_selected)
  })
  
  # Update celltype upon click
  observeEvent(input$celltype_tf_markers_plot_selected, {
    updateTextInput(session = session, "celltype_tf_markers_celltype", value = input$celltype_tf_markers_plot_selected)
  })

  
  #########################
  ## NMP differentiation ##
  #########################
  
  ## GRN ##
  
  plot_nmp_grn <- reactive({
  
    ## START TEST ##
    # input <- list()
    # input$nmp_diff_colourby <- "gene_expression"
    # input$nmp_diff_feature <- "T" #"chr18:32542271-32542871"
    # input$nmp_grn_colourby <- "celltype_expression"
    # input$nmp_grn_colourby_celltype <- "NMP"
    ## END TEST ##
      
    # to.plot <- insilico_chip_stats[tf%in%input$tf_chip]
    net <- nmp_grn.igraph
    
    color_palette <- grDevices::colorRamp(c("gray62", "purple"))( (1:100)/100 )
    color_palette <- cbind(color_palette, seq(100, 255, length.out = 100))
      
    if (input$nmp_grn_colourby=="eigenvalue_centrality") {
      eigenvalue_centrality_scores <- igraph::eigen_centrality(net)$vector[names(V(net))]
      color_nodes <- colourvalues::colour_values(eigenvalue_centrality_scores, palette = color_palette)
      names(color_nodes) <- names(V(net))
      V(net)$color <- color_nodes[names(V(net))]
    } else if (input$nmp_grn_colourby=="degree_centrality") {
      degree_centrality_scores <- igraph::degree(net)[names(V(net))]
      color_nodes <- colourvalues::colour_values(degree_centrality_scores, palette = color_palette)
      names(color_nodes) <- names(V(net))
      V(net)$color <- color_nodes[names(V(net))]
    } else if (input$nmp_grn_colourby=="celltype") {
      expr.mtx <- rna_expr_pseudobulk.mtx[stringr::str_to_title(names(V(net))),nmp_celltypes]
      rownames(expr.mtx) <- toupper(rownames(expr.mtx))
      color_nodes <- celltype_colours[colnames(expr.mtx)[apply(expr.mtx,1,which.max)]]
      names(color_nodes) <- rownames(expr.mtx)
      V(net)$color <- color_nodes[names(V(net))]
    } else if (input$nmp_grn_colourby=="celltype_expression") {
      expr <- rna_expr_pseudobulk.mtx[stringr::str_to_title(names(V(net))),input$nmp_grn_colourby_celltype] %>% minmax.normalisation
      # color_nodes <- colourvalues::colour_values(expr, palette = "viridis")
      color_nodes <- colourvalues::colour_values(expr, palette = color_palette)
      V(net)$color <- color_nodes
    }
    
    # Modify edge attributes
    E(net)[which(E(net)$weight<0)]$color <- "darkblue"  # Colour negative correlation edges as blue
    E(net)[which(E(net)$weight>0)]$color <- "darkred"   # Colour positive correlation edges as re
    
    # Modify vertex attributes
    # V(net)$group <- V(net)$class
    # V(net)$shape <- stringr::str_replace_all(V(net)$class,c("gene"="triangle","TF"="circle"))
    
    visIgraph(net, randomSeed=42) %>%
      visIgraphLayout(randomSeed=42, physics = TRUE) %>%
      # visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visNodes(shadow = TRUE, opacity=1, font = list(color="black", size=30)) %>%
      visEdges(width = 0.1, color = "black", smooth = FALSE) %>%
      # visGroups(groupname = "TF", shape = "circle", size=35, font = list(color="black", size=30), color = list(background="#63B8FF", hover="#4876FF", border="#63B8FF")) %>%
      # visGroups(groupname = "gene", shape = "triangle", size=20, font = list(color="black", size=35), color = list(background="#EE6363", hover="red", border="#EE6363")) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        minVelocity = 2,
        # maxVelocity = 2,
        forceAtlas2Based = list(gravitationalConstant = -300),
        stabilization = TRUE # By default, vis.js computes coordinates dynamically and waits for stabilization before rendering
      ) %>%
      visInteraction(
        hover = TRUE,
        dragNodes = TRUE,
        dragView = TRUE,
        zoomView = TRUE
      ) %>%
      visLegend(width = 0.1, position = "right", main = "Legend")
  })
  
  output$nmp_grn_plot <- renderVisNetwork({
    # shiny::validate(need(input$tf_chip%in%TFs, "Please select one or multiple TFs from our annotation." ))
    return(plot_nmp_grn())
  })
  
    
  ## Trajectory and box plots ##
  
  observeEvent(input$nmp_diff_colourby, {
    if (input$nmp_diff_colourby=="gene_expression") {
      updateSelectizeInput(session = session, inputId = 'nmp_diff_feature', choices = rna_genes, server = TRUE, selected = "Mesp1")
    } else if (input$nmp_diff_colourby=="chromatin_accessibility") {
      tmp <- grep(sprintf("%s:",input$nmp_diff_atac_chr),atac_peaks,value=T)
      updateSelectizeInput(session = session, inputId = 'nmp_diff_feature', choices = c(tmp,"chr7:79789147-79789747"), server = TRUE, selected = "chr7:79789147-79789747")
    }
  })
  observeEvent(input$nmp_diff_atac_chr, {
    if (input$nmp_diff_colourby=="chromatin_accessibility") {
      tmp <- grep(sprintf("%s:",input$nmp_diff_atac_chr),atac_peaks,value=T)
      updateSelectizeInput(session = session, inputId = 'nmp_diff_feature', choices = c(tmp,"chr7:79789147-79789747"), server = TRUE, selected = "chr7:79789147-79789747")
    }
  })
  
  plot_nmp_differentiation <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$nmp_diff_colourby <- "chromatin_accessibility"
    # input$nmp_diff_feature <-  "chr18:32542271-32542871"
    ## END TEST ##
    
    ## Fetch data ##
    
    to.plot <- nmp_metadata.dt %>% .[celltype%in%nmp_celltypes]
    
    if (input$nmp_diff_colourby == "gene_expression") {
      expr.dt <- data.table(
        metacell = colnames(nmp_rna_expr.array),
        color = as.numeric(nmp_rna_expr.array[input$nmp_diff_feature,])
      )
      to.plot <- to.plot %>% merge(expr.dt,by="metacell") %>% setorder(color)
      plot_title <- sprintf("%s expression",input$nmp_diff_feature)
    } else if (input$nmp_diff_colourby == "chromatin_accessibility") {
      expr.dt <- data.table(
        metacell = colnames(nmp_atac_peaks.array),
        color = as.numeric(nmp_atac_peaks.array[input$nmp_diff_feature,])
      )
      to.plot <- to.plot %>% merge(expr.dt,by="metacell") %>% setorder(color)
      plot_title <- sprintf("%s accessibility",input$nmp_diff_feature)
    } else {
      to.plot$color <- to.plot[[input$nmp_diff_colourby]]
    }
    
    ## Plot trajectory ##
    
    p.trajectory <- ggplot(to.plot, aes(x=V1, y=V2, color=color)) +
      geom_point_interactive(aes(tooltip = celltype, data_id = celltype), size=2.25, alpha=0.9) +
      theme_classic() +
      ggplot_theme_NoAxes() +
      theme(
        legend.position = "none"
      )
    
    # Modify legends    
    if (input$nmp_diff_colourby%in%c("celltype","stage","sample")) {
      p.trajectory <- p.trajectory + guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))
    } else {
      p.trajectory <- p.trajectory + theme(legend.title = element_blank(), legend.position = "top")
      
    }
    
    # Define palette
    palette <- switch(input$nmp_diff_colourby, 
                      "celltype" = celltype_palette, 
                      "stage" = stage_palette, 
                      "sample" = sample_palette, 
                      "gene_expression" = rna_palette, 
                      "chromatin_accessibility" = atac_palette
    )
    p.trajectory <- p.trajectory + palette
    
    ## Plot RNA expression or chromatin accessibility per celltype ##

    if (input$nmp_diff_colourby%in%c("gene_expression","chromatin_accessibility")) {
      to.plot <- nmp_metadata.dt[,c("metacell","celltype")] %>% 
        .[celltype%in%nmp_celltypes] %>% 
        merge(expr.dt,by="metacell") %>%
        .[,celltype:=factor(celltype,levels=nmp_celltypes)]
      
      p.expr <- ggplot(to.plot, aes(x=celltype, y=color, fill=celltype, tooltip=celltype, data_id = celltype)) +
        geom_jitter_interactive(size=2, width=0.05, alpha=0.5, shape=21) +
        geom_violin_interactive(scale="width", alpha=0.40) +
        geom_boxplot_interactive(width=0.5, outlier.shape=NA, alpha=0.70) +
        stat_summary(fun.data = function(x) { return(c(y = max(to.plot$color)+0.5, label = length(x)))}, geom = "text", size=2) +
        scale_fill_manual(values=celltype_colours) +
        labs(x="", y="", title=plot_title) +
        # guides(x = guide_axis(angle = 90)) +
        theme_classic() +
        theme(
          # plot.title = element_text(hjust=0.5, size=rel(0.9)),
          plot.title = element_blank(),
          axis.text.x = element_text(colour="black",size=rel(1)),
          # axis.text.y = element_text(colour="black",size=rel(1.25)),
          axis.title.y = element_text(colour="black",size=rel(1)),
          axis.ticks.x = element_blank(),
          legend.position = "none"
        )
      
      p <- p.trajectory + p.expr + plot_layout(ncol=2, widths=c(1.25,1))
    } else {
      p <- p.trajectory
    }
      
    girafe(
      code = print(p),
      width_svg = 7, height_svg = 5,
      options = list( 
        opts_sizing(rescale = FALSE),
        # opts_selection(type = "single", css = "cursor:pointer;fill:magenta;color:magenta"),
        opts_selection(type = "single", css = ""),
        # opts_hover_inv(css = "opacity:0.45;"),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
    
  })
  
  output$nmp_diff_plot = renderGirafe({
    if (input$nmp_diff_colourby=="chromatin_accessibility") shiny::validate(need(input$nmp_diff_atac_chr%in%chr_mm10,"Please select chromosome"))
    if (input$nmp_diff_colourby%in%c("gene_expression","chromatin_accessibility")) shiny::validate(need(input$nmp_diff_feature%in%c(atac_peaks,atac_genes,rna_genes),"Please select a feature from our annotation"))
    plot_nmp_differentiation()
  })
  
  ########################
  ## In silico ChIP-seq ##
  ########################
  
  plot_insilico_chip = reactive({
    
    ## START TEST ##
    # input$insilico_chip_tf <- c("FOXA2","FOXC1","FOXB1")
    # input$insilico_chip_motif_annotation <- "JASPAR"
    # input$insilico_chip_min_score <- 0.25
    ## END TEST ##
    
    tmp <- file.path(data_folder,sprintf("insilico_chip/%s/%s.bed.gz",tolower(input$insilico_chip_motif_annotation),input$insilico_chip_tf))
    tfs.to.plot <- input$insilico_chip_tf[file.exists(tmp)]  # if (any(!file.exists(tmp))) {
    
    virtual_chip.dt <- tfs.to.plot %>% map(function(i) {
      fread(file.path(data_folder,sprintf("insilico_chip/%s/%s.bed.gz",tolower(input$insilico_chip_motif_annotation),i))) %>%
        # setnames(c("chr","start","end","score","correlation_score","max_accessibility_score","motif_score","motif_counts")) %>%
        setnames(c("chr","start","end","score")) %>%
        .[!is.na(score)] %>%
        .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
        .[,c("chr","start","end"):=NULL] %>%
        .[,tf:=i] %>%
        return
    }) %>% rbindlist %>% .[,tf:=factor(tf,levels=tfs.to.plot)]
    
    ## Plot minimum score vs number of peaks ##
    
    to.plot <- insilico_chip_stats.dt[tf%in%tfs.to.plot]
    
    p1 <- ggplot(to.plot, aes_string(x="min_score", y="log2_N", color="tf")) +
      geom_line_interactive(size=1.5, aes(data_id=tf)) +
      labs(y="Number of predicted binding sites (log2)", x="Minimum score") +
      geom_vline(xintercept=input$insilico_chip_min_score, linetype="dashed") +
      scale_x_continuous(breaks=seq(0,0.75,by=0.10)) +
      theme_classic() +
      theme(
        axis.text = element_text(color="black"),
        legend.position = "top",
        legend.title = element_blank()
      )
    
    ## Plot fraction of positively correlated peaks ##
    
    to.plot <- virtual_chip.dt[abs(score)>=input$insilico_chip_min_score] %>% 
      .[,sign:=factor(c("-","+"))[(sign(score)>0)+1]] %>%
      .[,.N,by=c("tf","sign")]
    
    p2 <- ggplot(to.plot, aes(x=tf, y=N, fill=sign)) +
      geom_bar_interactive(stat="identity", aes(data_id=tf), color="black") +
      # scale_fill_discrete(drop=FALSE) + 
      scale_x_discrete(drop=FALSE) +
      labs(x="", y="Number of putative binding sites") +
      scale_fill_manual(values=c("+"="#d8b365","-"="#5ab4ac")) +
      # scale_fill_brewer(palette="Dark2") +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        axis.text.y = element_text(size=rel(0.75), color="black"),
        axis.text.x = element_text(size=rel(1), color="black")
      )
    
    # plot_grid(plotlist=list(p1,p2), scale=0.95, ncol=2, rel_widths = c(1/2,1/2)) %>% return(.)
    
    girafe(
      code = print(p1+p2),
      width_svg = 12, height_svg = 8,
      options = list( 
        opts_sizing(rescale = FALSE),
        opts_hover_inv(css = "opacity:0.2;"),
        opts_hover(css = "")
      )
    ) %>% return(.)
  })
  
  output$insilico_chip_plot = renderGirafe({
    shiny::validate(need(input$insilico_chip_motif_annotation%in%c("CISBP","JASPAR"), "Please select motif annotation" ))
    shiny::validate(need(input$insilico_chip_tf%in%TFs_all, "Please select one or multiple TFs from our annotation." ))
    return(plot_insilico_chip())
  })
  
  ## Download ##
  
  roots1 <- c(wd = file.path(data_folder,sprintf('insilico_chip')))
  shinyFileChoose(input, 'insilico_chip_files', roots = roots1, filetypes=c('gz'))
  
  output$insilico_chip_download_print <- renderPrint({
    if (length(input$insilico_chip_files)>1) {
      tmp <- sapply(input$insilico_chip_files[["files"]], function(x) x[[3]]) %>% unname
      cat(paste(tmp,collapse="\n"))
    }
  })
  
  output$insilico_chip_download <- downloadHandler(
    filename = function() {
      files.to.download <- as.character(parseFilePaths(roots1, input$insilico_chip_files)$datapath)
      if (length(files.to.download)>1) {
        "insilico_chip.zip"
      } else {
        input$insilico_chip_files[[1]][[1]][[2]]
      }
    },
    content = function(file) {
      files.to.download <- as.character(parseFilePaths(roots1, input$insilico_chip_files)$datapath)
      if (length(files.to.download)>1) {
        utils::zip(file, files=files.to.download, flags = "-jr9X")
        # tar(file, files.to.download)
      } else {
        file.copy(from=files.to.download, to=file)
      }
    },
    contentType = "application/zip"
  )
  
}
