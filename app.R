## Install and load required packages ----
required_packages <- c(
  "shiny", "Seurat", "ggplot2", "EnhancedVolcano", "dplyr", "epitools",
  "pwr", "DT", "MAST", "DESeq2", "tibble", "shinyjs", "clusterProfiler",
  "slingshot", "scales", "viridis", "RColorBrewer", "metap", "RANN",
  "reshape2"
)

cran_pkgs <- setdiff(required_packages, c("MAST", "DESeq2", "clusterProfiler", "slingshot"))
bioc_pkgs <- intersect(required_packages, c("MAST", "DESeq2", "clusterProfiler", "slingshot"))

options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

if (length(bioc_pkgs)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}

# Increase maximum upload size to 6000 MB
options(shiny.maxRequestSize = 6000 * 1024^2)  # 6000 MB
options(shiny.launch.browser = TRUE)
################### Define UI for the app ########
ui <- fluidPage(
  titlePanel("Seurat Data Visualization and Differential Gene Expression Analysis by Robert Vanagas: Version 5.4"),
  
  sidebarLayout(
    sidebarPanel(
      # ==== Data Upload and Loading 
      fileInput("seurat_file", "Upload already clustered Seurat Object (.rds):",
                accept = c(".rds")),
      h4("Object Requirements"),
      h6("1. Normalization + Transformation already performed (Ex by SCTransform)
     Check this by running head(seurat_object$Idents)."),
      h6("2. Clustering and UMAP already performed (Ex by FindClusters and RunUMAP functions)."),
      h6("3. If multiple tissues/samples in object, the metadata needs a column specifying what tissue each cell belongs to."),
      actionButton("load_data", "Load/Reset Object"),
      hr(),
      
      # ==== NEW: switch active identities 
      selectInput(
        "identity_column",
        "Use this metadata column as current identities:",
        choices = NULL            # filled after object is loaded
      ),
      actionButton("apply_identity", "Apply"),
      hr(),
      
      # ------------------- Tissue Grouping 
      h4("Tissue Grouping"),
      selectInput("tissue_column", "Select Metadata Column for Tissue Information:",
                  choices = NULL),
      
      numericInput(
        "num_groups",
        "Number of Tissue Groups:",
        value = 2, min = 2, step = 1
      ),
      uiOutput("groups_ui"),
      actionButton("assign_tissue_groups", "Assign Tissue Groups"),
      hr(),
      
      # ------------------- Sub-clustering 
      h4("Subclustering"),
      h6("If you want to re-cluster one or multiple clusters, please use the subset option instead and download the subsetted object.
     Perform re-clustering outside of the app and then upload it for analysis."),
      uiOutput("cluster_selection_ui"),
      numericInput("resolution", "Resolution:", value = 0.1, step = 0.1),
      uiOutput("graph_selection_ui"),
      actionButton("subcluster_button", "Subcluster"),
      hr(),
      
      # ------------------- Rename Clusters 
      h4("Rename Clusters"),
      selectInput("cluster_to_rename", "Select Cluster to Rename:",
                  choices = NULL),
      textInput("new_cluster_name", "New Cluster Name:", value = ""),
      actionButton("rename_cluster_button", "Rename Cluster"),
      hr(),
      
      # ------------------- Merge Clusters 
      h4("Merge Clusters"),
      selectInput("merge_from", "Select Cluster to Merge From:",
                  choices = NULL),
      selectInput("merge_into", "Select Cluster to Merge Into:",
                  choices = NULL),
      actionButton("merge_clusters", "Merge Clusters"),
      hr(),
      
      # ---- Universal downloader ----------------------------------------------------
      hr(),
      selectInput(
        "download_what",
        "Download plots and tables:",
        choices = character(0)          # start empty – will be filled from the server
      ),
      downloadButton("download_selected", "Save selected"),
      
      # ------------------- Subset / Save Object 
      h4("Save Seurat Object"),
      selectizeInput("clusters_to_subset",
                     "Select Clusters to Include in Subset. Select all clusters to save the current object entirely:",
                     choices = NULL, multiple = TRUE),
      textInput("subset_filename",
                "Filename for Subsetted Seurat Object (without extension):",
                value = "subset_seurat_object"),
      downloadButton("download_subset_seurat", "Download Subsetted Seurat Object")
    ),
    
    
    mainPanel(
      tabsetPanel(
        # **UMAP Plot Tab**
        tabPanel("UMAP",
                 h4("UMAP Settings"),
                 # Dynamic UI for Split By
                 uiOutput("split_by_option_ui"),
                 h6("BUG: If niches do not appear in Group by: Rename one of the clusters and it should update"),
                 selectInput("group_by", "Group By:",
                             choices = NULL),  # Will be updated dynamically
                 selectizeInput("cells_highlight", "Cells to Highlight (Cluster IDs):",
                                choices = NULL, multiple = TRUE),
                 numericInput("sizes_highlight", "Highlight Size:", value = 0.5),
                 numericInput("xlim1", "X-axis Min:", value = NA),
                 numericInput("xlim2", "X-axis Max:", value = NA),
                 numericInput("ylim1", "Y-axis Min:", value = NA),
                 numericInput("ylim2", "Y-axis Max:", value = NA),
                 actionButton("plot_umap", "Plot UMAP"),
                 hr(),
                 plotOutput("umap_plot", height = "800px")  # Increased height for better visibility
        ),
        
        # **Feature Plot Tab**
        tabPanel("Feature Plot",
                 h4("Feature Plot Settings"),
                 selectizeInput("feature_genes", "Select Genes:",
                                choices = NULL, multiple = TRUE),
                 # Dynamic UI for Split By
                 uiOutput("split_by_feature_ui"),
                 actionButton("plot_feature", "Plot Feature"),
                 hr(),
                 div(style = 'height: 800px; overflow-y: scroll;',
                     uiOutput("feature_plots_ui")
                 )
        ),
        
        # **Volcano Plot Tab**
        tabPanel("Volcano Plot",
                 h4("Volcano Plot Settings 
     Don't forget to change view to show all genes in Volcano plot!"),
                 
                 # New Group By Setting
                 selectInput("volcano_group_by", "Group By:",
                             choices = NULL, selected = "ident"),  # Will be updated dynamically
                 
                 # Dynamic UI for Ident 1 based on Group By selection
                 uiOutput("volcano_ident1_ui"),
                 
                 # Dynamic UI for Ident 2 based on Group By selection
                 uiOutput("volcano_ident2_ui"),
                 
                 numericInput("xlim_left", "X-axis Min:", value = -10),
                 numericInput("xlim_right", "X-axis Max:", value = 10),
                 numericInput("ylim_top", "Y-axis Max:", value = 25),
                 numericInput("p_cutoff", "p-value Cutoff:", value = 0.05),
                 numericInput("fc_cutoff", "Log2 FC Cutoff:", value = 1.5),
                 
                 # **Added Inputs for P-value Type and Test Use**
                 selectInput("volcano_pval_type", "P-value Type:",
                             choices = c("Adjusted p-value" = "p_val_adj"),
                             selected = "p_val_adj"),
                 
                 selectInput("volcano_test_use", "Test Used:",
                             choices = c("wilcox", "wilcox_limma", "bimod", "roc", "t",
                                         "negbinom", "poisson", "LR", "MAST", "DESeq2"),
                             selected = "wilcox"),
                 
                 # **Added Input for Number of Entries to Show**
                 numericInput("volcano_table_entries", "Number of Entries to Show in Table:",
                              value = 10, min = 1),
                 
                 actionButton("plot_volcano", "Plot Volcano"),
                 hr(),
                 plotOutput("volcano_plot", height = "800px"),
                 DT::dataTableOutput("volcano_table")  # **Added Output for Volcano Table**
        ),
        
        
        # **Heatmap Tab**
        tabPanel("Heatmap",
                 h4("Heatmap Settings"),
                 
                 # New Group By Setting
                 selectInput("heatmap_group_by", "Group By:",
                             choices = NULL, selected = "ident"),  # Will be updated dynamically
                 
                 numericInput("top_n_markers", "Number of Top Markers per Cluster:",
                              value = 5, min = 1),
                 numericInput("logfc_threshold", "Log2 FC Threshold:", value = 1.5),
                 
                 selectizeInput("additional_genes", "Select Additional Genes to add at bottom of Heatmap", choices = NULL, multiple = TRUE),
                 
                 # **Added Inputs for Test Use and P-value Type**
                 selectInput("heatmap_test_use", "Test Used:",
                             choices = c("wilcox", "wilcox_limma", "bimod", "roc", "t",
                                         "negbinom", "poisson", "LR", "MAST", "DESeq2"),
                             selected = "wilcox"),
                 
                 selectInput("heatmap_pval_type", "P-value Type:",
                             choices = c("Adjusted p-value" = "p_val_adj"),
                             selected = "p_val_adj"),
                 
                 # **Added Input for Number of Entries to Show**
                 numericInput("heatmap_table_entries", "Number of Entries to Show in Table:",
                              value = 10, min = 1),
                 
                 actionButton("plot_heatmap", "Plot Heatmap"),
                 hr(),
                 plotOutput("heatmap_plot", height = "800px"),
                 DT::dataTableOutput("heatmap_table")  # **Added Output for Heatmap Table**
        ),
        
        
        
        
        # **Spatial Plot Tab**
        tabPanel("Spatial Plot",
                 h4("Spatial Plot Settings"),
                 h6("The purpose of this tab is to quickly inspect that your Seurat object has Spatial data in it"),
                 # Dynamic UI to select tissue/image if multiple are available
                 uiOutput("spatial_tissue_ui"),
                 # Select grouping variable
                 selectInput("spatial_group_by", "Group By:",
                             choices = NULL, selected = "seurat_clusters"),
                 # Point size
                 numericInput("spatial_size", "Point Size:", value = 1, min = 0.1, step = 0.1),
                 # Alpha transparency
                 numericInput("spatial_alpha", "Alpha Transparency:", value = 1, min = 0, max = 1, step = 0.1),
                 # Border color
                 textInput("spatial_border_color", "Border Color:", value = "black"),
                 # Toggle axes
                 checkboxInput("spatial_axes", "Show Axes", value = FALSE),
                 actionButton("plot_spatial", "Plot Spatial"),
                 hr(),
                 plotOutput("spatial_plot", height = "800px")  # Increased height for better visibility
        ),
        
        # **Spatial Niche Analysis Tab**
        tabPanel("Spatial Niche Analysis",
                 h4("Build Niche Assay"),
                 h6("This builds niche assay based on spatial clustering. To identify the clusters, go to UMAP, Volcano or Heatmap tabs and select Group by: niches"),
                 sidebarLayout(
                   sidebarPanel(
                     uiOutput("niche_fov_ui"),
                     numericInput("niche_k", "Number of Niches (niches.k):", value = 2, min = 1),
                     numericInput("niche_neighbors_k", "Number of Neighbors (neighbors.k):", value = 20, min = 1),
                     actionButton("build_niche_assay", "Build Niche Assay"),
                     hr(),
                     h5("Rename Niches"),
                     uiOutput("niche_selection_ui"),
                     textInput("new_niche_name", "New Niche Name:", value = ""),
                     actionButton("rename_niche_button", "Rename Niche"),
                     hr(),
                   ),
                   mainPanel(
                     plotOutput("niche_imagedimplot", height = "800px"),
                     DT::dataTableOutput("niche_table")
                   )
                 )
        ),
        
        
        
        
        # **Neighborhood Analysis Tab**
        tabPanel("Neighborhood Analysis",
                 h4("Neighborhood Analysis Settings"),
                 h6("Be patient! This can take a few minutes for large datasets"),
                 # Select clusters to analyze
                 selectizeInput("neighbor_clusters", "Select Cluster(s) to Analyze:",
                                choices = NULL, multiple = TRUE),
                 
                 # Option to select niches
                 checkboxInput("neighbor_use_niches", "Restrict Cells to Selected Niches", value = FALSE),
                 
                 # Conditional UI to select niches
                 conditionalPanel(
                   condition = "input.neighbor_use_niches == true",
                   selectizeInput("neighbor_niches", "Select Niche(s):",
                                  choices = NULL, multiple = TRUE)
                 ),
                 
                 # Number of cells to use
                 numericInput("neighbor_num_cells", "Number of Cells per Niche to Use (0 for all):",
                              value = 0, min = 0, step = 1),
                 
                 # Specify k values
                 numericInput("neighbor_k_min", "Minimum k (Number of Neighbors):",
                              value = 1, min = 1, step = 1),
                 numericInput("neighbor_k_max", "Maximum k (Number of Neighbors):",
                              value = 40, min = 1, step = 1),
                 
                 # Plotting options
                 textInput("neighbor_plot_title", "Plot Title:",
                           value = "Neighborhood Composition of Selected Clusters"),
                 textInput("neighbor_x_label", "X-axis Label:",
                           value = "Number of Nearest Neighbors (k)"),
                 textInput("neighbor_y_label", "Y-axis Label:",
                           value = "Proportion of Neighboring Cells"),
                 selectInput("neighbor_theme", "Select Plot Theme:",
                             choices = c("Minimal" = "theme_minimal",
                                         "Classic" = "theme_classic",
                                         "Grey" = "theme_grey",
                                         "Dark" = "theme_dark",
                                         "Light" = "theme_light"),
                             selected = "theme_minimal"),
                 checkboxInput("neighbor_legend", "Show Legend", value = TRUE),
                 
                 # Additional plotting options
                 numericInput("neighbor_line_size", "Line Thickness:",
                              value = 1, min = 0.1, step = 0.1),
                 numericInput("neighbor_axis_text_size", "Axis Text Size:",
                              value = 12, min = 1, step = 1),
                 numericInput("neighbor_axis_label_size", "Axis Label Size:",
                              value = 14, min = 1, step = 1),
                 
                 # Select image to use
                 uiOutput("neighbor_image_ui"),  # Dynamic UI to display image choices
                 
                 # Action buttons
                 actionButton("run_neighbor_analysis", "Run Analysis"),
                 hr(),
                 
                 # Output plot
                 plotOutput("neighbor_plot", height = "800px")
        ),
        
        
        
        # **Violin Plot Tab**
        tabPanel("Violin Plot",
                 h4("Violin Plot Settings"),
                 selectInput("vln_group_by", "Group By:",
                             choices = NULL, selected = "ident"),
                 selectizeInput("vln_features", "Select Genes:",
                                choices = NULL, multiple = TRUE),
                 numericInput("vln_point_size", "Point Size:", value = 0.1, min = 0.01, step = 0.01),
                 # Removed width input
                 actionButton("plot_vln", "Plot Violin"),
                 hr(),
                 plotOutput("vln_plot", height = "800px")  # Increased height for better visibility
        ),
        
        # **Dot Plot Tab**
        tabPanel("Dot Plot",
                 h4("Dot Plot Settings"),
                 selectInput("dot_group_by", "Group By:",
                             choices = NULL, selected = "ident"),
                 selectizeInput("dot_features", "Select Genes:",
                                choices = NULL, multiple = TRUE),
                 numericInput("dot_scale_size", "Scale Size:", value = 8, min = 1, step = 1),
                 numericInput("dot_scale_color", "Scale Color Intensity:", value = 6, min = 1, step = 1),
                 actionButton("plot_dot", "Plot Dot"),
                 hr(),
                 plotOutput("dot_plot", height = "800px")  # Increased height for better visibility
        ),
        
        
        
        
        
        # **Cluster Composition Tab UI**
        tabPanel("Cluster Composition",
                 sidebarLayout(
                   sidebarPanel(
                     h4("Cluster Composition Settings"),
                     
                     # **Tissue or Tissue Group Selection**
                     radioButtons("composition_group_type", "Group By:",
                                  choices = c("Individual Tissue" = "individual",
                                              "Tissue Group" = "group"),
                                  selected = "individual"),
                     
                     # Conditional UI for selecting individual tissues
                     conditionalPanel(
                       condition = "input.composition_group_type == 'individual'",
                       selectizeInput("composition_tissues", "Select Tissue(s):",
                                      choices = NULL,  # To be updated dynamically
                                      multiple = TRUE)
                     ),
                     
                     # Conditional UI for selecting multiple tissue groups
                     conditionalPanel(
                       condition = "input.composition_group_type == 'group'",
                       selectizeInput("composition_tissue_group", "Select Tissue Group(s):",
                                      choices = NULL,  # To be updated dynamically
                                      multiple = TRUE,
                                      options = list(maxItems = NULL))  # Allows unlimited selections
                     ),
                     
                     # **Cluster Selection (Optional)**
                     selectizeInput("composition_clusters", "Select Cluster(s) to Include (Leave empty for all):",
                                    choices = NULL, multiple = TRUE),
                     
                     # **Plot Options**
                     checkboxInput("show_percentage", "Display Percentage Values", value = TRUE),
                     
                     # **Color Palette Selection**
                     selectInput("composition_color_palette", "Select Color Palette:",
                                 choices = c("Set1", "Set2", "Set3", "Paired", "Pastel1", "Pastel2",
                                             "Dark2", "Accent",
                                             "viridis", "magma", "plasma", "inferno", "cividis"),  # Added viridis palettes
                                 selected = "Set1"),
                     

                   ),
                   
                   mainPanel(
                     h4("Cluster Composition Stacked Bar Plot"),
                     plotOutput("composition_plot", height = "600px")
                   )
                 )
        ),
        
        
        # **Cluster Comparison Tab**
        tabPanel("Cluster Comparison",
                 h4("Compare Differential Gene Expression Within a Cluster"),
                 h6("Compare gene expression within a selected cell cluster across different samples/tissues/conditions to uncover sample/tissue/condition-specific changes"),
                 sidebarLayout(
                   sidebarPanel(
                     # **Volcano Plot Section**
                     h5("Volcano Plot Settings"),
                     
                     # Select Cluster for Volcano Plot
                     selectInput("comp_cluster_volcano", "Select Cluster for Volcano Plot:",
                                 choices = NULL, selected = NULL),
                     
                     # Comparison Mode for Volcano Plot
                     radioButtons("comp_mode_volcano", "Comparison Mode:",
                                  choices = c("Tissue Groups" = "groups",
                                              "Specific Tissues" = "specific"),
                                  selected = "groups"),
                     
                     # Conditional UI for Specific Tissues in Volcano Plot
                     conditionalPanel(
                       condition = "input.comp_mode_volcano == 'specific'",
                       selectInput("comp_volcano_tissue1", "Select Tissue 1: (Right side)",
                                   choices = NULL, selected = NULL),
                       selectInput("comp_volcano_tissue2", "Select Tissue 2: (Left side)",
                                   choices = NULL, selected = NULL)
                     ),
                     
                     # **Conditional UI for Tissue Groups in Volcano Plot**
                     conditionalPanel(
                       condition = "input.comp_mode_volcano == 'groups'",
                       selectInput("comp_volcano_group1", "Select Tissue Group 1: (Right side)",
                                   choices = NULL, selected = NULL),
                       selectInput("comp_volcano_group2", "Select Tissue Group 2: (Left side)",
                                   choices = NULL, selected = NULL)
                     ),
                     
                     # DE Analysis Parameters for Volcano Plot
                     numericInput("comp_pval_volcano", "Adjusted P-value Cutoff:",
                                  value = 0.05, min = 0.0001, max = 1, step = 0.001),
                     numericInput("comp_logfc_volcano", "Log2 Fold Change Cutoff:",
                                  value = 1, min = 0, max = 5, step = 0.1),
                     
                     # **Added Inputs for P-value Type and Test Used in Volcano Plot**
                     selectInput("comp_volcano_pval_type", "P-value Type:",
                                 choices = c("Adjusted p-value" = "p_val_adj",
                                             "Raw p-value" = "p_val"),
                                 selected = "p_val_adj"),
                     
                     selectInput("comp_volcano_test_use", "Test Used:",
                                 choices = c("wilcox", "wilcox_limma", "bimod", "roc", "t",
                                             "negbinom", "poisson", "LR", "MAST", "DESeq2"),
                                 selected = "wilcox"),
                     
                     # **Added Input for Number of Entries to Show in Volcano Table**
                     numericInput("comp_volcano_table_entries", "Number of Entries to Show in Table:",
                                  value = 10, min = 1),
                     
                     # Button to Generate Volcano Plot
                     actionButton("generate_volcano", "Generate Volcano Plot"),
                     
                     hr(),
                     
                     # **Heatmap Section**
                     h5("Heatmap Settings"),
                     
                     # Select Cluster for Heatmap
                     selectInput("comp_cluster_heatmap", "Select Cluster for Heatmap:",
                                 choices = NULL, selected = NULL),
                     
                     # Comparison Mode for Heatmap
                     radioButtons("comp_mode_heatmap", "Comparison Mode:",
                                  choices = c("Tissue Groups" = "groups",
                                              "All Tissues" = "all"),
                                  selected = "groups"),
                     
                     # DE Analysis Parameters for Heatmap
                     numericInput("comp_pval_heatmap", "Adjusted P-value Cutoff:",
                                  value = 0.05, min = 0.0001, max = 1, step = 0.001),
                     numericInput("comp_logfc_heatmap", "Log2 Fold Change Cutoff when calculating DE",
                                  value = 1, min = -10, max = 10, step = 0.1),
                     
                     # **Added Inputs for Test Use and P-value Type in Heatmap**
                     selectInput("comp_heatmap_test_use", "Test Used:",
                                 choices = c("wilcox", "wilcox_limma", "bimod", "roc", "t",
                                             "negbinom", "poisson", "LR", "MAST", "DESeq2"),
                                 selected = "wilcox"),
                     
                     selectInput("comp_heatmap_pval_type", "P-value Type:",
                                 choices = c("Adjusted p-value" = "p_val_adj",
                                             "Raw p-value" = "p_val"),
                                 selected = "p_val_adj"),
                     
                     # **Added Input for Number of Entries to Show in Heatmap Table**
                     numericInput("comp_heatmap_table_entries", "Number of Entries to Show in Table:",
                                  value = 10, min = 1),
                     
                     # **Added Inputs for Heatmap Customization**
                     numericInput("comp_top_n_markers", "Number of Top Markers per Cluster:",
                                  value = 5, min = 1),
                     numericInput("comp_logfc_threshold_heatmap", "Log2 FC Threshold to be included in heatmap:", value = 1.5, min = -10, max = 10, step = 0.1),
                     
                     selectizeInput("comp_additional_genes_heatmap", "Select Additional Genes to add at bottom of Heatmap",
                                    choices = NULL, multiple = TRUE),
                     
                     # Button to Generate Heatmap
                     actionButton("generate_heatmap", "Generate Heatmap"),
                     
                     hr(),
                     
                   ),
                   
                   mainPanel(
                     tabsetPanel(
                       # **Volcano Plot Tab**
                       tabPanel("Volcano Plot",
                                plotOutput("comp_volcano_plot", height = "600px"),
                                br(),
                                DT::dataTableOutput("comp_de_table_volcano")
                       ),
                       
                       # **Heatmap Tab**
                       tabPanel("Heatmap",
                                plotOutput("comp_heatmap_plot", height = "600px"),
                                br(),
                                DT::dataTableOutput("comp_de_table_heatmap")
                       )
                     )
                   )
                 )
        ),
      )
    )
  )
)

# Define server logic for the app
server <- function(input, output, session) {
  
  ########################## SECTION 1.  General server functionality #######################
  
  # Reactive value that stores the Seurat object 
  seurat_obj <- reactiveVal(NULL)
  
  # Reactive values that store plots / tables for download 
  reactive_values <- reactiveValues(
    umap_plot                         = NULL,
    volcano_plot                      = NULL,
    heatmap_plot                      = NULL,
    FindAllMarkers_table              = NULL,
    FindMarkers_table                 = NULL,
    conserved_volcano_plot            = NULL,
    conserved_heatmap_plot            = NULL,
    FindConservedMarkers_table        = NULL,
    FindConservedMarkers_heatmap_table= NULL,
    neighbor_plot                     = NULL
  )
  
  # --------------------------------------------------------------------------- #
  #  helper: refresh every UI control that shows cluster levels                 #
  # --------------------------------------------------------------------------- #
  refreshClusterInputs <- function(obj) {
    updated_levels <- levels(obj)
    
    single_selects <- c(
      "cluster_select", "ident1",             # various pickers
      "cluster_to_rename", "stats_cluster",
      "merge_from", "merge_into",
      "comp_cluster_volcano", "comp_cluster_heatmap"
    )
    multi_selects  <- c(
      "clusters_to_subset", "cells_highlight",
      "neighbor_clusters",  "composition_clusters",
      "stats_combined_clusters", "comparison_clusters"
    )
    
    for (id in single_selects)
      updateSelectInput(session, id, choices = updated_levels)
    
    for (id in multi_selects)
      updateSelectizeInput(session, id, choices = updated_levels, server = TRUE)
    
    # controls that need extra items in the choices ----------------------------
    updateSelectInput(session, "ident2",
                      choices = c("None", updated_levels))
    updateSelectizeInput(session, "stats_denominator_clusters",
                         choices  = c("All available cells" = "all", updated_levels),
                         selected = "all", server = TRUE)
  }
  
  # --------------------------------------------------------------------------- #
  #  observer: user clicks **Apply** to switch the active identities            #
  # --------------------------------------------------------------------------- #
  observeEvent(input$apply_identity, {
    req(seurat_obj())
    
    chosen_col <- input$identity_column
    obj        <- seurat_obj()
    
    if (!chosen_col %in% colnames(obj@meta.data)) {
      showNotification("Selected column is not in the metadata.", type = "error")
      return(NULL)
    }
    
    # 1) make that column the Idents
    Idents(obj) <- obj@meta.data[[chosen_col]]
    
    # 2) keep the synthetic 'ident' column in sync so plots using `group.by="ident"` still work
    obj$ident <- as.character(Idents(obj))
    
    # 3) save back
    seurat_obj(obj)
    
    # 4) refresh every cluster-picker UI element
    refreshClusterInputs(obj)
    
    showNotification(paste("Current identities set to", chosen_col), type = "message")
  })
  
  # --------------------------------------------------------------------------- #
  #  reactive values for tissue grouping (used in SECTION 2) 
  tissue_groups <- reactiveValues(groups = NULL)
  
  
  #  load an .rds Seurat file and initialise the app
  observeEvent(input$load_data, {
    req(input$seurat_file)
    showModal(modalDialog("Loading data...", footer = NULL))
    
    # try reading the file 
    loaded_seurat <- tryCatch({
      readRDS(input$seurat_file$datapath)
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error loading Seurat object:", e$message), type = "error")
      return(NULL)
    })
    req(loaded_seurat)
    
    # guarantee required metadata columns 
    if (!"seurat_clusters" %in% colnames(loaded_seurat@meta.data)) {
      if (!is.null(Idents(loaded_seurat))) {
        loaded_seurat$seurat_clusters <- as.character(Idents(loaded_seurat))
      } else {
        removeModal()
        showNotification("No clustering found – run FindClusters first.", type = "error")
        return(NULL)
      }
    }
    if (!"stim" %in% colnames(loaded_seurat@meta.data)) {
      if ("orig.ident" %in% colnames(loaded_seurat@meta.data)) {
        loaded_seurat$stim <- loaded_seurat$orig.ident
      } else {
        removeModal()
        showNotification("'orig.ident' metadata column missing.", type = "error")
        return(NULL)
      }
    }
    
    # keep synthetic 'ident' column in sync
    loaded_seurat$ident <- as.character(Idents(loaded_seurat))
    
    # store the object 
    seurat_obj(loaded_seurat)
    removeModal()
    
    # ---------- NEW: populate the identity-switch drop-down -------------------
    updateSelectInput(
      session, "identity_column",
      choices  = colnames(seurat_obj()@meta.data),
      selected = "ident"
    )
    
    # every picker that needs cluster levels
    refreshClusterInputs(seurat_obj())
    
    # -------------- graph list for sub-clustering -----------------------------
    observeEvent(seurat_obj(), {
      req(seurat_obj())
      available_graphs <- names(seurat_obj()@graphs)
      output$graph_selection_ui <- renderUI({
        selectInput("graph_to_use", "Select Graph for Subclustering:",
                    choices = available_graphs, selected = "RNA_snn")
      })
    })
    
    # -------------- assorted other pickers ------------------------------------
    metadata_columns <- c("ident", colnames(seurat_obj()@meta.data))
    updateSelectInput(session, "group_by",
                      choices = metadata_columns, selected = "ident")
    
    all_genes <- rownames(seurat_obj())
    updateSelectizeInput(session, "feature_genes",
                         choices = all_genes, server = TRUE)
    
    updateSelectInput(session, "tissue_column",
                      choices  = colnames(seurat_obj()@meta.data),
                      selected = ifelse("stim" %in% colnames(seurat_obj()@meta.data), "stim", NULL))
    
    # neighbourhood analysis – image choices
    output$neighbor_image_ui <- renderUI({
      image_choices <- names(seurat_obj()@images)
      if (length(image_choices) == 0) {
        helpText("No spatial images available in the Seurat object.")
      } else {
        selectInput("neighbor_image", "Select Image:", choices = image_choices)
      }
    })
  })
  
  # --------------------------------------------------------------------------- #
  #  misc outputs used in UI 
  output$current_clusters <- renderPrint({
    req(seurat_obj())
    levels(seurat_obj())
  })
  
  output$cluster_selection_ui <- renderUI({
    req(seurat_obj())
    selectInput("cluster_select", "Select Cluster to Subcluster:",
                choices = levels(seurat_obj()))
  })
  
  ########################## SECTION 2.  Tissue group functionality #######################
  
  
  # Dynamic UI for Tissue Groups
  output$groups_ui <- renderUI({
    req(input$num_groups)          # Ensure the number of groups is specified
    req(seurat_obj())              # Ensure the Seurat object is loaded
    req(input$tissue_column)       # Ensure the tissue column is selected
    
    num <- input$num_groups
    group_ui_list <- lapply(1:num, function(i) {
      tagList(
        textInput(
          inputId = paste0("group", i, "_name"), 
          label = paste("Name for Group", i, ":"), 
          value = paste("Group", i)
        ),
        selectizeInput(
          inputId = paste0("group", i, "_tissues"), 
          label = paste("Select Tissues for Group", i, ":"), 
          choices = unique(seurat_obj()@meta.data[[input$tissue_column]]),
          multiple = TRUE
        ),
        hr()
      )
    })
    do.call(tagList, group_ui_list)
  })
  
  # Assign Tissue Groups based on user input
  observeEvent(input$assign_tissue_groups, {
    req(seurat_obj())
    req(input$num_groups)
    req(input$tissue_column)
    
    num <- input$num_groups
    group_names <- sapply(1:num, function(i) input[[paste0("group", i, "_name")]])
    group_tissues <- lapply(1:num, function(i) input[[paste0("group", i, "_tissues")]])
    
    # Validate that all group names are provided
    if (any(group_names == "")) {
      showNotification("Please provide names for all groups.", type = "error")
      return(NULL)
    }
    
    # Validate that tissues are assigned to only one group
    all_assigned_tissues <- unlist(group_tissues)
    unique_assigned_tissues <- unique(all_assigned_tissues)
    if (length(all_assigned_tissues) != length(unique_assigned_tissues)) {
      overlapping_tissues <- all_assigned_tissues[duplicated(all_assigned_tissues)]
      showNotification(
        paste("Tissues cannot be assigned to multiple groups:", 
              paste(overlapping_tissues, collapse = ", ")), 
        type = "error"
      )
      return(NULL)
    }
    
    # Optional: Validate that all tissues are assigned to a group
    # Uncomment the following lines if you want to ensure that all tissues are grouped
    # all_tissues <- unique(seurat_obj()@meta.data[[input$tissue_column]])
    # unassigned_tissues <- setdiff(all_tissues, unique_assigned_tissues)
    # if (length(unassigned_tissues) > 0) {
    #   showNotification(
    #     paste("The following tissues are not assigned to any group:", 
    #           paste(unassigned_tissues, collapse = ", ")), 
    #     type = "warning"
    #   )
    # }
    
    # Assign tissue_group based on selections
    updated_seurat <- seurat_obj()
    updated_seurat$tissue_group <- "Other"  # Default category
    
    for (i in 1:num) {
      group_name <- group_names[i]
      tissues <- group_tissues[[i]]
      updated_seurat$tissue_group[updated_seurat@meta.data[[input$tissue_column]] %in% tissues] <- group_name
    }
    
    # Update the Seurat object
    seurat_obj(updated_seurat)
    
    # Store group names in reactive values
    tissue_groups$groups <- group_names
    
    # Update group_by options to include the new tissue_group
    metadata_columns <- c("ident", colnames(seurat_obj()@meta.data))
    updateSelectInput(session, "group_by",
                      choices = metadata_columns, selected = "ident")
    
    showNotification("Tissue groups assigned successfully.", type = "message")
  })
  
  ################ SECTION 3. Subclustering functionality ########################
  
  # Sub-clustering 
  observeEvent(input$subcluster_button, {
    req(seurat_obj())
    
    # --- user settings ----------------------------------------------------------
    cluster_to_sub <- input$cluster_select
    resolution     <- input$resolution
    graph_name     <- input$graph_to_use
    
    # --- run Seurat::FindSubCluster --------------------------------------------
    obj <- seurat_obj()
    obj <- FindSubCluster(
      object         = obj,
      cluster        = cluster_to_sub,
      graph.name     = graph_name,
      subcluster.name = "sub.cluster",   # ⬅ creates ‘sub.cluster’ in meta.data
      resolution     = resolution
    )
    
    # make the new sub-clusters the active identities
    obj <- SetIdent(obj, value = "sub.cluster")
    
    # keep synthetic ‘ident’ column in sync so existing plots keep working
    obj$ident <- as.character(Idents(obj))
    
    # --- save back & refresh every UI element that lists clusters --------------
    seurat_obj(obj)
    
    # (A)  picker that switches active identities  – **NEW**
    updateSelectInput(session, "identity_column",
                      choices  = colnames(seurat_obj()@meta.data),
                      selected = "sub.cluster")
    
    # (B)  all inputs whose choices are the current cluster levels
    refreshClusterInputs(obj)
    
    # (C)  picker for group-by options in plots
    metadata_columns <- c("ident", colnames(obj@meta.data))
    updateSelectInput(session, "group_by",
                      choices  = metadata_columns,
                      selected = "ident")
    
    # (optional) show a quick message
    showNotification("Sub-clustering finished – ‘sub.cluster’ set as Idents",
                     type = "message")
  })
  
  
  
  
  
  # Rename Clusters
  observeEvent(input$rename_cluster_button, {
    req(seurat_obj())
    old_name <- input$cluster_to_rename
    new_name <- input$new_cluster_name
    if (new_name != "") {
      immune.subsets.13.cl <- seurat_obj()
      # Get current identities
      current_idents <- Idents(immune.subsets.13.cl)
      # Rename cluster
      levels(current_idents)[levels(current_idents) == old_name] <- new_name
      # Set the new identities
      Idents(immune.subsets.13.cl) <- current_idents
      # **Add this line to synchronize 'ident' with active identities**
      immune.subsets.13.cl$ident <- as.character(Idents(immune.subsets.13.cl))
      # Update the reactive seurat object
      seurat_obj(immune.subsets.13.cl)
      # Update cluster choices in all relevant dropdowns
      updateSelectInput(session, "cluster_to_rename",
                        choices = levels(seurat_obj()))
      updateSelectInput(session, "cluster_select",
                        choices = levels(seurat_obj()))
      updateSelectInput(session, "ident1",
                        choices = levels(seurat_obj()))
      updateSelectInput(session, "ident2",
                        choices = c("None", levels(seurat_obj())))
      updateSelectInput(session, "stats_cluster",
                        choices = levels(seurat_obj()))
      updateSelectizeInput(session, "cells_highlight",
                           choices = levels(seurat_obj()), server = TRUE)
      
      # Update clusters for subsetting
      updateSelectizeInput(session, "clusters_to_subset",
                           choices = levels(seurat_obj()), server = TRUE)
      
      # Update comparison clusters for Statistics
      updateSelectizeInput(session, "comparison_clusters",
                           choices = levels(seurat_obj()), server = TRUE)
      
      # Update group_by options
      metadata_columns <- c("ident", colnames(seurat_obj()@meta.data))
      updateSelectInput(session, "group_by",
                        choices = metadata_columns, selected = "ident")
      # Update current clusters display
      output$current_clusters <- renderPrint({
        levels(seurat_obj())
      })
    }
  })
  
  
  # **Dynamic UI for Split By in UMAP**
  output$split_by_option_ui <- renderUI({
    req(seurat_obj())
    req(input$tissue_column)
    choices <- c("None", input$tissue_column, "tissue_group")
    radioButtons("split_by_option", "Split By:",
                 choices = choices,
                 selected = "None")
  })
  
  # **Dynamic UI for Split By in Feature Plot**
  output$split_by_feature_ui <- renderUI({
    req(seurat_obj())
    req(input$tissue_column)
    choices <- c("None", input$tissue_column, "tissue_group")
    radioButtons("split_by_feature", "Split By:",
                 choices = choices,
                 selected = "None")
  })
  ########################## SECTION 4.  UMAP plot functionality ##################
  
  # UMAP Plot
  observeEvent(input$plot_umap, {
    req(seurat_obj())
    output$umap_plot <- renderPlot({
      split_by <- if (input$split_by_option == "None") NULL else input$split_by_option
      group_by <- input$group_by
      xlim_vals <- c(input$xlim1, input$xlim2)
      ylim_vals <- c(input$ylim1, input$ylim2)
      cells_to_highlight <- input$cells_highlight
      sizes_highlight <- input$sizes_highlight
      
      # Handle split_by correctly
      p <- DimPlot(seurat_obj(), label = TRUE, group.by = group_by, split.by = split_by)
      if (length(cells_to_highlight) > 0) {
        cells.highlight <- WhichCells(seurat_obj(), idents = cells_to_highlight)
        p <- DimPlot(seurat_obj(), label = TRUE, group.by = group_by, split.by = split_by,
                     cells.highlight = cells.highlight, sizes.highlight = sizes_highlight)
      }
      if (!any(is.na(xlim_vals))) {
        p <- p + xlim(xlim_vals)
      }
      if (!any(is.na(ylim_vals))) {
        p <- p + ylim(ylim_vals)
      }
      # Store the plot object
      reactive_values$umap_plot <- p
      print(p)
    })
  })
  
  
  ########################## SECTION 5.  Feature plot functionality #######################
  
  # Feature Plot
  observeEvent(input$plot_feature, {
    req(seurat_obj())
    features <- input$feature_genes
    split_by_input <- input$split_by_feature  # Capture the input value
    output$feature_plots_ui <- renderUI({
      plot_output_list <- lapply(features, function(feature) {
        plotname <- paste("feature_plot_", feature, sep = "")
        plotOutput(plotname)
      })
      do.call(tagList, plot_output_list)
    })
    # For each feature, render the plot
    lapply(features, function(feature) {
      output[[paste("feature_plot_", feature, sep = "")]] <- renderPlot({
        split_by <- if (split_by_input == "None") NULL else split_by_input
        FeaturePlot(seurat_obj(), features = feature, split.by = split_by)
      })
    })
  })
  

  
  ########################## SECTION 6.  Volcano Plot functionality #######################
  
  # Dynamic UI for Group By in Volcano Plot
  observe({
    req(seurat_obj())
    
    # Update Group By options to include metadata columns, including 'niches'
    metadata_columns <- c("ident", colnames(seurat_obj()@meta.data))
    updateSelectInput(session, "volcano_group_by",
                      choices = metadata_columns,
                      selected = "ident")
    
    # Update Group By options for Heatmap as well
    updateSelectInput(session, "heatmap_group_by",
                      choices = metadata_columns,
                      selected = "ident")
  })
  
  # Dynamic UI for Ident 1 based on Group By selection in Volcano Plot
  output$volcano_ident1_ui <- renderUI({
    req(seurat_obj())
    req(input$volcano_group_by)
    
    group_by_col <- input$volcano_group_by
    unique_groups <- unique(seurat_obj()@meta.data[[group_by_col]])
    
    selectInput("ident1", "Ident 1: (Right side of volcano plot)",
                choices = unique_groups, multiple = TRUE)
  })
  
  # Dynamic UI for Ident 2 based on Group By selection in Volcano Plot
  output$volcano_ident2_ui <- renderUI({
    req(seurat_obj())
    req(input$volcano_group_by)
    
    group_by_col <- input$volcano_group_by
    unique_groups <- unique(seurat_obj()@meta.data[[group_by_col]])
    
    selectInput("ident2", "Ident 2: (Left side of volcano plot)",
                choices = c("Test against all populations", unique_groups), selected = "Test against all populations", multiple = TRUE)
  })
  
  # Enhanced Volcano Plot
  observeEvent(input$plot_volcano, {
    req(seurat_obj())
    req(input$ident1)
    
    output$volcano_plot <- renderPlot({
      ident1 <- input$ident1
      ident2 <- input$ident2
      
      group_by_col <- input$volcano_group_by
      
      # Determine if ident2 is set to test against all
      if ("Test against all populations" %in% ident2 || is.null(ident2)) {
        ident2 <- NULL
      }
      
      # Perform differential expression with selected test.use
      markers <- FindMarkers(seurat_obj(), ident.1 = ident1, ident.2 = ident2,
                             group.by = group_by_col,
                             test.use = input$volcano_test_use)
      
      # Determine the p-value column based on user selection
      p_val_column <- input$volcano_pval_type
      
      # Ensure that the selected p-value column exists
      if (!(p_val_column %in% colnames(markers))) {
        showNotification(paste("Selected p-value type", p_val_column, "does not exist in the markers table."), type = "error")
        return(NULL)
      }
      
      # Use the selected p-value column
      y_val <- p_val_column
      
      # Enhanced Volcano requires specific column names, so rename accordingly
      markers_renamed <- markers
      colnames(markers_renamed)[which(colnames(markers_renamed) == y_val)] <- "p_val_plot"
      if ("avg_log2FC" %in% colnames(markers_renamed)) {
        colnames(markers_renamed)[which(colnames(markers_renamed) == "avg_log2FC")] <- "log2FC_plot"
      } else if ("avg_logFC" %in% colnames(markers_renamed)) {
        # Handle different naming conventions
        colnames(markers_renamed)[which(colnames(markers_renamed) == "avg_logFC")] <- "log2FC_plot"
      } else {
        showNotification("Cannot find log2 fold change column in markers.", type = "error")
        return(NULL)
      }
      
      # Define plot limits
      xlim_range <- c(input$xlim_left, input$xlim_right)
      ylim_top <- input$ylim_top
      pCutoff <- input$p_cutoff
      FCcutoff <- input$fc_cutoff
      
      # Create Enhanced Volcano Plot
      p <- EnhancedVolcano(markers_renamed, lab = rownames(markers_renamed),
                           x = 'log2FC_plot', y = 'p_val_plot',
                           xlim = xlim_range,
                           ylim = c(0, ylim_top),
                           pCutoff = pCutoff,
                           FCcutoff = FCcutoff,
                           title = paste('Enhanced Volcano Plot:', paste(ident1, collapse=", "), 
                                         'vs', ifelse(is.null(ident2), 'All', paste(ident2, collapse=", "))),
                           subtitle = 'Differential Gene Expression',
                           caption = 'Log2 fold change vs Selected p-value',
                           pointSize = 2.0, labSize = 3.0,
                           legendLabels = c('NS', 'Log2 FC', 'p-value', 'p-value & Log2 FC'),
                           legendLabSize = 12, legendIconSize = 4.0)
      # Store the plot object
      reactive_values$volcano_plot <- p
      
      # Store the FindMarkers table for interactive display
      reactive_values$FindMarkers_table <- markers
      
      print(p)
    })
  })
  
  # Render Interactive Volcano Table
  output$volcano_table <- DT::renderDataTable({
    req(reactive_values$FindMarkers_table)
    
    # Select the p-value type based on user input
    pval_type <- input$volcano_pval_type
    
    table_to_show <- reactive_values$FindMarkers_table
    
    # Round numeric columns for better readability
    table_to_show <- table_to_show %>%
      mutate(
        avg_log2FC = if("avg_log2FC" %in% colnames(.)) round(avg_log2FC, 3) else NA,
        p_val = signif(p_val, 3),
        p_val_adj = signif(p_val_adj, 3)
      )
    
    # Add gene names as a separate column
    table_to_show <- cbind(Gene = rownames(table_to_show), table_to_show)
    
    # Select the appropriate p-value column
    if (pval_type %in% colnames(table_to_show)) {
      table_to_show <- table_to_show %>%
        dplyr::select(Gene, everything(), all_of(pval_type))
    } else {
      showNotification(paste("Selected p-value type", pval_type, "does not exist in the markers table."), type = "error")
      return(NULL)
    }
    
    # Rename the p-value column for clarity
    colnames(table_to_show)[ncol(table_to_show)] <- "P-value"
    
    DT::datatable(table_to_show, 
                  options = list(
                    pageLength = input$volcano_table_entries,  # User-specified number of entries
                    scrollX = TRUE,
                    dom = 'Bfrtip',  # Include buttons in the table
                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                  ),
                  rownames = FALSE,
                  extensions = 'Buttons',
                  filter = 'top',
                  caption = 'Differential Expression Results')
  })
  

  
  ########################## SECTION 7.  Heatmap functionality ##############################
  
  # ---- (Dynamic UI for Group-By is still created earlier) 
  
  # ---------- Generate Heatmap 
  observeEvent(input$plot_heatmap, {
    req(seurat_obj())
    
    ## -- 1. Differential expression 
    top_n           <- input$top_n_markers
    logfc_threshold <- input$logfc_threshold
    test_use        <- input$heatmap_test_use
    group_by_col    <- input$heatmap_group_by
    
    Immunemarkers <- FindAllMarkers(
      object   = seurat_obj(),
      test.use = test_use,
      group.by = group_by_col
    )
    
    reactive_values$FindAllMarkers_table <- Immunemarkers  # for the data table
    
    top_markers <- Immunemarkers %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC >  logfc_threshold |
                      avg_log2FC < -logfc_threshold) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = top_n) %>%
      ungroup()
    
    selected_features <- unique(c(top_markers$gene,
                                  input$additional_genes))
    
    ## -- 2. Ensure all requested genes are in scale.data
    da            <- DefaultAssay(seurat_obj())
    scaled_genes  <- rownames(GetAssayData(seurat_obj(), assay = da,
                                           slot  = "scale.data"))
    missing_feats <- setdiff(selected_features, scaled_genes)
    
    obj_for_heat <- if (length(missing_feats) > 0) {
      message("Scaling ", length(missing_feats), " missing gene(s) for heatmap …")
      ScaleData(seurat_obj(),
                features = missing_feats,
                assay    = da,
                verbose  = FALSE)
    } else {
      seurat_obj()
    }
    
    ## -- 3. Render heatmap 
    output$heatmap_plot <- renderPlot({
      p <- DoHeatmap(
        object   = subset(obj_for_heat, downsample = 1000),
        features = selected_features,
        group.by = group_by_col
      ) +
        NoLegend() +
        theme(text = element_text(size = 8))
      
      reactive_values$heatmap_plot <- p
      print(p)
    })
  })
  
  # ---------- Interactive Heatmap Table 
  output$heatmap_table <- DT::renderDataTable({
    req(reactive_values$FindAllMarkers_table)
    
    pval_type <- input$heatmap_pval_type
    tbl       <- reactive_values$FindAllMarkers_table %>%
      mutate(
        avg_log2FC = round(avg_log2FC, 3),
        p_val      = signif(p_val, 3),
        p_val_adj  = signif(p_val_adj, 3)
      ) %>%
      tibble::rownames_to_column("Gene")
    
    if (!(pval_type %in% colnames(tbl))) {
      showNotification(paste("Selected p-value type", pval_type,
                             "does not exist in the markers table."),
                       type = "error")
      return(NULL)
    }
    
    colnames(tbl)[which(colnames(tbl) == pval_type)] <- "P-value"
    
    DT::datatable(tbl,
                  options = list(
                    pageLength = input$heatmap_table_entries,
                    scrollX    = TRUE,
                    dom        = 'Bfrtip',
                    buttons    = c('copy','csv','excel','pdf','print')
                  ),
                  rownames   = FALSE,
                  extensions = 'Buttons',
                  filter     = 'top',
                  caption    = 'Marker Genes Results')
  })
  

  
  

  
  ########################## SECTION 9.  DOWNLOAD SUBSET functionality ####################
  # Download Subsetted Seurat Object
  output$download_subset_seurat <- downloadHandler(
    filename = function() {
      paste0(input$subset_filename, ".rds")
    },
    content = function(file) {
      req(seurat_obj())
      selected_clusters <- input$clusters_to_subset
      if (length(selected_clusters) == 0) {
        showNotification("Please select at least one cluster to subset.", type = "error")
        return(NULL)
      }
      immune.subsets.13.cl <- seurat_obj()
      # Subset the Seurat object based on selected clusters
      subset_seurat_obj <- subset(immune.subsets.13.cl, idents = selected_clusters)
      # Save the subsetted Seurat object to the specified file
      saveRDS(subset_seurat_obj, file = file)
    }
  )
  
  ########################## SECTION 10. Visualisation functionality #######################
  
  # Reactive value to store the spatial images names
  spatial_images <- reactiveVal(NULL)
  
  # Update Spatial Plot UI and related inputs when data is loaded
  observeEvent(seurat_obj(), {
    req(seurat_obj())
    
    # Check if the Seurat object has spatial data
    if (length(seurat_obj()@images) > 0) {
      spatial_imgs <- names(seurat_obj()@images)
      spatial_images(spatial_imgs)
      
      # Update the spatial group_by selection input
      updateSelectInput(session, "spatial_group_by",
                        choices = c("ident", colnames(seurat_obj()@meta.data)),
                        selected = "ident")
      
      # Update the spatial tissue selection input
      updateSelectInput(session, "spatial_tissue",
                        choices = spatial_imgs,
                        selected = spatial_imgs[1])
    } else {
      spatial_images(NULL)
      
      # Notify the user that no spatial data is available
      showNotification("The loaded Seurat object does not contain spatial data.", type = "warning")
      
      # Disable Spatial Plot tab inputs
      output$spatial_tissue_ui <- renderUI({
        helpText("No spatial data available for this Seurat object.")
      })
    }
  })
  
  ########################## SECTION 11. ImageDimPlot functionality #######################
  
  # **Dynamic UI for Spatial Tissue Selection**
  output$spatial_tissue_ui <- renderUI({
    if (!is.null(spatial_images())) {
      selectInput("spatial_tissue", "Select Tissue/Image:", choices = spatial_images())
    } else {
      helpText("No spatial data available for this Seurat object.")
    }
  })
  
  # **Render Spatial Plot Using ImageDimPlot**
  observeEvent(input$plot_spatial, {
    req(seurat_obj())
    
    if (is.null(spatial_images())) {
      showNotification("Spatial data is not available in the loaded Seurat object.", type = "error")
      return(NULL)
    }
    
    # Ensure that the selected tissue/image exists
    if (!(input$spatial_tissue %in% spatial_images())) {
      showNotification("Selected tissue/image does not exist.", type = "error")
      return(NULL)
    }
    
    output$spatial_plot <- renderPlot({
      req(input$spatial_tissue)
      req(input$spatial_group_by)
      
      # Generate the ImageDimPlot without manual color customization
      p <- tryCatch({
        ImageDimPlot(
          object = seurat_obj(),
          fov = input$spatial_tissue,
          group.by = input$spatial_group_by,
          size = input$spatial_size,
          alpha = input$spatial_alpha,
          border.color = input$spatial_border_color,
          axes = input$spatial_axes,
          coord.fixed = TRUE,
          flip_xy = FALSE
        ) + 
          ggtitle(paste("Spatial DimPlot for", input$spatial_group_by, "in", input$spatial_tissue)) +
          theme(plot.title = element_text(hjust = 0.5))
      }, error = function(e) {
        showNotification(paste("Error generating spatial plot:", e$message), type = "error")
        return(NULL)
      })
      
      # Store the plot object for download
      if (!is.null(p)) {
        reactive_values$spatial_plot <- p
        print(p)
      }
    })
  })
  

  
  ########################## SECTION 12. SPATIAL NICHE ASSAY ###########################
  # Output UI for fov selection
  output$niche_fov_ui <- renderUI({
    req(seurat_obj())
    if (length(seurat_obj()@images) > 0) {
      selectInput("niche_fov", "Select Field of View (fov):", choices = names(seurat_obj()@images))
    } else {
      helpText("No spatial data available.")
    }
  })
  
  # Observe Build Niche Assay button
  observeEvent(input$build_niche_assay, {
    req(seurat_obj())
    
    # Check if spatial data is available
    if (length(seurat_obj()@images) == 0) {
      showNotification("No spatial data available in the Seurat object.", type = "error")
      return(NULL)
    }
    
    # Retrieve input values
    fov <- input$niche_fov
    niches_k <- input$niche_k
    neighbors_k <- input$niche_neighbors_k
    
    # Save the current active assay
    current_assay <- DefaultAssay(seurat_obj())
    
    # Run BuildNicheAssay
    seurat_obj_updated <- tryCatch({
      BuildNicheAssay(
        object = seurat_obj(),
        fov = fov,
        assay = "niche",
        group.by = "seurat_clusters",
        niches.k = niches_k,
        neighbors.k = neighbors_k
      )
    }, error = function(e) {
      showNotification(paste("Error building niche assay:", e$message), type = "error")
      return(NULL)
    })
    
    # If BuildNicheAssay failed, exit the observer
    if (is.null(seurat_obj_updated)) {
      return(NULL)
    }
    
    # Reset the active assay to the previous one
    DefaultAssay(seurat_obj_updated) <- current_assay
    
    # Update the reactive Seurat object
    seurat_obj(seurat_obj_updated)
    
    # Generate the ImageDimPlot
    output$niche_imagedimplot <- renderPlot({
      ImageDimPlot(
        object = seurat_obj(),
        group.by = "niches",
        fov = fov,
        size = 1.5,
        dark.background = FALSE
      )
    })
    
    # Generate the interactive table
    output$niche_table <- DT::renderDataTable({
      table_data <- as.data.frame(table(Idents(seurat_obj()), seurat_obj()$niches))
      colnames(table_data) <- c("Cluster", "Niche", "Count")
      DT::datatable(
        table_data,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
        ),
        rownames = FALSE,
        extensions = 'Buttons',
        filter = 'top',
        caption = 'Cluster vs. Niche Composition'
      )
    })
    
    # Update the niche selection input for renaming
    output$niche_selection_ui <- renderUI({
      niches <- unique(seurat_obj()$niches)
      selectInput("niche_to_rename", "Select Niche to Rename:", choices = niches)
    })
  })
  
  # Observe Rename Niche button
  observeEvent(input$rename_niche_button, {
    req(seurat_obj())
    
    old_niche <- input$niche_to_rename
    new_niche <- input$new_niche_name
    
    if (new_niche != "") {
      # Retrieve the Seurat object
      seurat_object <- seurat_obj()
      
      # Modify the niches in the Seurat object
      niche_data <- seurat_object$niches
      niche_data[niche_data == old_niche] <- new_niche
      seurat_object$niches <- niche_data
      
      # Update the reactive Seurat object
      seurat_obj(seurat_object)
      
      # Regenerate the ImageDimPlot
      output$niche_imagedimplot <- renderPlot({
        fov <- input$niche_fov
        ImageDimPlot(
          object = seurat_obj(),
          group.by = "niches",
          fov = fov,
          size = 1.5,
          dark.background = FALSE
        )
      })
      
      # Update the niche table
      output$niche_table <- DT::renderDataTable({
        table_data <- as.data.frame(table(Idents(seurat_obj()), seurat_obj()$niches))
        colnames(table_data) <- c("Cluster", "Niche", "Count")
        DT::datatable(
          table_data,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
          ),
          rownames = FALSE,
          extensions = 'Buttons',
          filter = 'top',
          caption = 'Cluster vs. Niche Composition'
        )
      })
      
      # Update the niche selection input
      output$niche_selection_ui <- renderUI({
        niches <- unique(seurat_obj()$niches)
        selectInput("niche_to_rename", "Select Niche to Rename:", choices = niches)
      })
    }
  })
  

  
  
  ########################## SECTION 13. Neighborhood Analysis #########################
  # Update cluster choices when Seurat object is loaded or updated
  observe({
    req(seurat_obj())
    
    # Update clusters
    updateSelectizeInput(session, "neighbor_clusters",
                         choices = levels(Idents(seurat_obj())),
                         server = TRUE)
  })
  
  # Update niche choices when niches are created or modified
  observe({
    req(seurat_obj())
    
    # Check if 'niches' column exists in meta.data
    if ("niches" %in% colnames(seurat_obj()@meta.data)) {
      niches_available <- unique(seurat_obj()$niches)
      niches_available <- niches_available[!is.na(niches_available)]
      
      updateSelectizeInput(session, "neighbor_niches",
                           choices = niches_available,
                           server = TRUE)
    } else {
      updateSelectizeInput(session, "neighbor_niches",
                           choices = NULL,
                           server = TRUE)
    }
  })
  
  # Update image choices for Neighborhood Analysis
  output$neighbor_image_ui <- renderUI({
    req(seurat_obj())
    image_choices <- names(seurat_obj()@images)
    if (length(image_choices) == 0) {
      helpText("No spatial images available in the Seurat object.")
    } else {
      selectInput("neighbor_image", "Select Image:",
                  choices = image_choices)
    }
  })
  
  # Observe the Run Analysis button
  observeEvent(input$run_neighbor_analysis, {
    req(seurat_obj())
    req(input$neighbor_clusters)
    
    # Get the selected clusters
    selected_clusters <- input$neighbor_clusters
    
    # Determine if niches are to be used
    if (input$neighbor_use_niches) {
      req(input$neighbor_niches)
      if (!"niches" %in% colnames(seurat_obj()@meta.data)) {
        showNotification("Niche information is not available. Please build the niche assay first.", type = "error")
        return(NULL)
      }
      selected_niches <- input$neighbor_niches
    } else {
      selected_niches <- NULL  # No niche restriction
    }
    
    num_cells <- input$neighbor_num_cells
    k_min <- input$neighbor_k_min
    k_max <- input$neighbor_k_max
    image_name <- input$neighbor_image
    
    # Extract spatial coordinates
    if (length(seurat_obj()@images) == 0) {
      showNotification("Seurat object does not contain spatial data.", type = "error")
      return(NULL)
    }
    
    # Ensure that the selected image exists
    if (!(image_name %in% names(seurat_obj()@images))) {
      showNotification("Selected image does not exist.", type = "error")
      return(NULL)
    }
    
    # Extract coordinates from the selected image
    coords_all <- GetTissueCoordinates(seurat_obj(), image = image_name)
    coords_all <- as.data.frame(coords_all)
    
    # Ensure 'cell' column exists and set as row names
    if(!"cell" %in% colnames(coords_all)) {
      coords_all$cell <- rownames(coords_all)
    }
    rownames(coords_all) <- coords_all$cell
    
    # Remove any duplicated row names
    coords_all <- coords_all[!duplicated(rownames(coords_all)), ]
    
    # Ensure 'x' and 'y' are numeric
    coords_all$x <- as.numeric(coords_all$x)
    coords_all$y <- as.numeric(coords_all$y)
    
    # Remove rows with NA values (if any)
    coords_all <- coords_all[complete.cases(coords_all), ]
    
    # Initialize a list to store results
    results_list <- list()
    
    if (is.null(selected_niches)) {
      # No niche restriction
      # Get cells in selected clusters
      selected_cells <- WhichCells(seurat_obj(), idents = selected_clusters)
      
      # Ensure selected_cells are present in coords_all
      selected_cells <- intersect(selected_cells, rownames(coords_all))
      
      # Handle case where there are no cells in the selected clusters
      if (length(selected_cells) == 0) {
        showNotification(paste("No cells found in cluster(s):", paste(selected_clusters, collapse = ", ")), type = "error")
        return(NULL)
      }
      
      # Determine number of cells to use
      if (num_cells <= 0 || num_cells >= length(selected_cells)) {
        cells_subset <- selected_cells
      } else {
        cells_subset <- sample(selected_cells, num_cells)
      }
      
      # Subset coordinates for the selected cells
      coords_subset <- coords_all[cells_subset, ]
      
      # Define the range of k values
      k_values <- seq(k_min, k_max)
      
      # Get the cluster assignments for all cells
      all_clusters <- Idents(seurat_obj())
      cluster_levels <- levels(all_clusters)  # Ensure consistent ordering
      
      # Loop over each k value
      for(k in k_values) {
        # Calculate nearest neighbors
        nn_result <- nn2(
          data = coords_all[, c("x", "y")],
          query = coords_subset[, c("x", "y")],
          k = k + 1  # +1 to exclude the cell itself
        )
        
        # Extract neighbor indices excluding the first column (the cell itself)
        neighbor_indices <- nn_result$nn.idx[, -1, drop = FALSE]  # Ensures matrix format even if k=1
        
        # Convert neighbor indices to cell IDs
        neighbor_cell_ids <- matrix(rownames(coords_all)[neighbor_indices], nrow = nrow(neighbor_indices))
        
        # Retrieve cluster assignments for neighbors
        neighbor_clusters <- matrix(all_clusters[neighbor_cell_ids], nrow = nrow(neighbor_cell_ids))
        
        # Calculate counts of each cluster for each cell
        cluster_counts <- apply(neighbor_clusters, 1, function(x) table(factor(x, levels = cluster_levels)))
        
        # Convert to matrix (clusters x cells)
        if (is.vector(cluster_counts)) {
          cluster_counts <- matrix(cluster_counts, nrow = length(cluster_counts), dimnames = list(names(cluster_counts), NULL))
        } else {
          cluster_counts <- as.matrix(cluster_counts)
        }
        
        # Sum counts across all cells in the subset
        total_counts <- rowSums(cluster_counts)
        
        # Calculate proportions
        cluster_proportions <- total_counts / sum(total_counts)
        
        # Store the proportions with cluster names
        cluster_proportions <- as.data.frame(t(cluster_proportions))
        cluster_proportions$k <- k
        
        # Add to results
        results_list[[as.character(k)]] <- cluster_proportions
      }
      
      # Combine results
      combined_results <- do.call(rbind, results_list)
      
      # Melt the data frame for plotting
      results_melted <- reshape2::melt(combined_results, id.vars = "k", variable.name = "cluster", value.name = "proportion")
      
      # Ensure 'cluster' is treated as a factor with levels in 'cluster_levels'
      results_melted$cluster <- factor(results_melted$cluster, levels = cluster_levels)
      
      # Plot the results
      plot_theme <- switch(input$neighbor_theme,
                           "theme_minimal" = theme_minimal(),
                           "theme_classic" = theme_classic(),
                           "theme_grey" = theme_grey(),
                           "theme_dark" = theme_dark(),
                           "theme_light" = theme_light(),
                           theme_minimal())  # Default theme
      
      p <- ggplot(results_melted, aes(x = k, y = proportion, color = cluster)) +
        geom_line(size = input$neighbor_line_size) +
        labs(
          title = input$neighbor_plot_title,
          x = input$neighbor_x_label,
          y = input$neighbor_y_label,
          color = "Cluster"
        ) +
        plot_theme +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = input$neighbor_axis_label_size),
          axis.text = element_text(size = input$neighbor_axis_text_size),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)
        )
      
      # Optionally hide legend
      if (!input$neighbor_legend) {
        p <- p + theme(legend.position = "none")
      }
      
    } else {
      # Niche restriction
      # Initialize a list to store results for each niche
      niche_results_list <- list()
      
      # Loop over selected niches
      for (niche in selected_niches) {
        # Filter cells based on selected clusters and niche
        cells_in_niche <- WhichCells(seurat_obj(), expression = niches == niche)
        
        # Cells in selected clusters within the niche
        cells_in_clusters <- WhichCells(seurat_obj(), idents = selected_clusters)
        selected_cells <- intersect(cells_in_clusters, cells_in_niche)
        
        # Ensure selected_cells are present in coords_all
        selected_cells <- intersect(selected_cells, rownames(coords_all))
        
        # Handle case where there are no cells in the selected clusters and niche
        if (length(selected_cells) == 0) {
          showNotification(paste("No cells found in cluster(s):", paste(selected_clusters, collapse = ", "),
                                 "within niche:", niche), type = "warning")
          next  # Skip to the next niche
        }
        
        # Determine number of cells to use
        if (num_cells <= 0 || num_cells >= length(selected_cells)) {
          cells_subset <- selected_cells
        } else {
          cells_subset <- sample(selected_cells, num_cells)
        }
        
        # Subset coordinates for the selected cells
        coords_subset <- coords_all[cells_subset, ]
        
        # Define the range of k values
        k_values <- seq(k_min, k_max)
        
        # Get the cluster assignments for all cells
        all_clusters <- Idents(seurat_obj())
        cluster_levels <- levels(all_clusters)  # Ensure consistent ordering
        
        # Loop over each k value
        for(k in k_values) {
          # Calculate nearest neighbors
          nn_result <- nn2(
            data = coords_all[, c("x", "y")],
            query = coords_subset[, c("x", "y")],
            k = k + 1  # +1 to exclude the cell itself
          )
          
          # Extract neighbor indices excluding the first column (the cell itself)
          neighbor_indices <- nn_result$nn.idx[, -1, drop = FALSE]  # Ensures matrix format even if k=1
          
          # Convert neighbor indices to cell IDs
          neighbor_cell_ids <- matrix(rownames(coords_all)[neighbor_indices], nrow = nrow(neighbor_indices))
          
          # Retrieve cluster assignments for neighbors
          neighbor_clusters <- matrix(all_clusters[neighbor_cell_ids], nrow = nrow(neighbor_cell_ids))
          
          # Calculate counts of each cluster for each cell
          cluster_counts <- apply(neighbor_clusters, 1, function(x) table(factor(x, levels = cluster_levels)))
          
          # Convert to matrix (clusters x cells)
          if (is.vector(cluster_counts)) {
            cluster_counts <- matrix(cluster_counts, nrow = length(cluster_counts), dimnames = list(names(cluster_counts), NULL))
          } else {
            cluster_counts <- as.matrix(cluster_counts)
          }
          
          # Sum counts across all cells in the subset
          total_counts <- rowSums(cluster_counts)
          
          # Calculate proportions
          cluster_proportions <- total_counts / sum(total_counts)
          
          # Store the proportions with cluster names
          cluster_proportions <- as.data.frame(t(cluster_proportions))
          cluster_proportions$k <- k
          cluster_proportions$niche <- niche
          
          # Add to niche results
          niche_results_list[[paste0(niche, "_", k)]] <- cluster_proportions
        }
      }
      
      # Combine results from all niches
      if (length(niche_results_list) == 0) {
        showNotification("No data available for the selected clusters and niches.", type = "error")
        return(NULL)
      }
      
      combined_results <- do.call(rbind, niche_results_list)
      
      # Melt the data frame for plotting
      results_melted <- reshape2::melt(combined_results, id.vars = c("k", "niche"), variable.name = "cluster", value.name = "proportion")
      
      # Ensure 'cluster' is treated as a factor with levels in 'cluster_levels'
      results_melted$cluster <- factor(results_melted$cluster, levels = cluster_levels)
      
      # Plot the results using ggplot2 with customizations
      plot_theme <- switch(input$neighbor_theme,
                           "theme_minimal" = theme_minimal(),
                           "theme_classic" = theme_classic(),
                           "theme_grey" = theme_grey(),
                           "theme_dark" = theme_dark(),
                           "theme_light" = theme_light(),
                           theme_minimal())  # Default theme
      
      p <- ggplot(results_melted, aes(x = k, y = proportion, color = cluster)) +
        geom_line(size = input$neighbor_line_size) +
        facet_wrap(~ niche) +
        labs(
          title = input$neighbor_plot_title,
          x = input$neighbor_x_label,
          y = input$neighbor_y_label,
          color = "Cluster"
        ) +
        plot_theme +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = input$neighbor_axis_label_size),
          axis.text = element_text(size = input$neighbor_axis_text_size),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)
        )
      
      # Optionally hide legend
      if (!input$neighbor_legend) {
        p <- p + theme(legend.position = "none")
      }
    }
    
    # Store the plot in reactive_values for download
    reactive_values$neighbor_plot <- p
    
    # Render the plot
    output$neighbor_plot <- renderPlot({
      print(p)
    })
    
  })

  ########################## SECTION 14. Violin plot functionality #######################
  
  # Update UI elements for Violin Plot when data is loaded
  observeEvent(seurat_obj(), {
    req(seurat_obj())
    
    # Update group_by options for Violin Plot
    updateSelectInput(session, "vln_group_by",
                      choices = c("ident", colnames(seurat_obj()@meta.data)),
                      selected = "seurat_clusters")
    
    # Update feature gene selection for Violin Plot
    all_genes <- rownames(seurat_obj())
    updateSelectizeInput(session, "vln_features",
                         choices = all_genes, server = TRUE)
  })
  
  # Generate Violin Plot
  observeEvent(input$plot_vln, {
    req(seurat_obj())
    req(input$vln_features)
    req(input$vln_group_by)
    
    output$vln_plot <- renderPlot({
      p <- VlnPlot(object = seurat_obj(),
                   features = input$vln_features,
                   group.by = input$vln_group_by,
                   pt.size = input$vln_point_size,
                   # width = input$vln_violin_width,  # Removed to avoid unused argument
                   combine = TRUE) + 
        theme(plot.title = element_text(hjust = 0.5))
      
      # Store the plot object for download
      reactive_values$vln_plot <- p
      
      print(p)
    })
  })
  

  ########################## SECTION 15. DOT PLOT functionality #######################
  
  # Update UI elements for Dot Plot when data is loaded
  observeEvent(seurat_obj(), {
    req(seurat_obj())
    
    # Update group_by options for Dot Plot
    updateSelectInput(session, "dot_group_by",
                      choices = c("ident", colnames(seurat_obj()@meta.data)),
                      selected = "seurat_clusters")
    
    # Update feature gene selection for Dot Plot
    all_genes <- rownames(seurat_obj())
    updateSelectizeInput(session, "dot_features",
                         choices = all_genes, server = TRUE)
  })
  
  # Generate Dot Plot
  observeEvent(input$plot_dot, {
    req(seurat_obj())
    req(input$dot_features)
    req(input$dot_group_by)
    
    output$dot_plot <- renderPlot({
      p <- DotPlot(object = seurat_obj(),
                   features = input$dot_features,
                   group.by = input$dot_group_by) + 
        scale_size(range = c(1, input$dot_scale_size)) +
        scale_color_gradient(low = "lightgray", high = "blue", 
                             limits = c(0, input$dot_scale_color)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5))
      
      # Store the plot object for download
      reactive_values$dot_plot <- p
      
      print(p)
    })
  })
  

  ########################## SECTION 16. CLUSTER COMPARISON SERVER LOGIC ###########
  
  # Initialize reactive values to store DE results and plots
  reactive_values$comp_de_results_volcano <- NULL
  reactive_values$comp_de_results_heatmap <- NULL
  reactive_values$comp_volcano_plot <- NULL
  reactive_values$comp_heatmap_plot <- NULL
  
  # Function to Retrieve All Unique Tissue Groups
  get_all_tissue_groups <- reactive({
    req(seurat_obj())
    if ("tissue_group" %in% colnames(seurat_obj()@meta.data)) {
      unique(seurat_obj()@meta.data$tissue_group)
    } else {
      character(0)
    }
  })
  
  # Update Cluster and Tissue Group Selection Inputs Upon Loading Seurat Object or Assigning Tissue Groups
  observe({
    req(seurat_obj())
    
    # Update Cluster Choices for Volcano Plot
    updateSelectInput(session, "comp_cluster_volcano",
                      choices = levels(Idents(seurat_obj())),
                      selected = levels(Idents(seurat_obj()))[1])
    
    # Update Cluster Choices for Heatmap
    updateSelectInput(session, "comp_cluster_heatmap",
                      choices = levels(Idents(seurat_obj())),
                      selected = levels(Idents(seurat_obj()))[1])
    
    # Retrieve all tissue groups
    all_tissue_groups <- get_all_tissue_groups()
    
    # Update Tissue Group Selection Choices for Volcano Plot
    updateSelectInput(session, "comp_volcano_group1",
                      choices = all_tissue_groups,
                      selected = all_tissue_groups[1])
    updateSelectInput(session, "comp_volcano_group2",
                      choices = all_tissue_groups,
                      selected = all_tissue_groups[2])
    
    # Update Specific Tissues Selection Inputs for Volcano Plot
    updateSelectInput(session, "comp_volcano_tissue1",
                      choices = unique(seurat_obj()@meta.data[[input$tissue_column]]),
                      selected = unique(seurat_obj()@meta.data[[input$tissue_column]])[1])
    updateSelectInput(session, "comp_volcano_tissue2",
                      choices = unique(seurat_obj()@meta.data[[input$tissue_column]]),
                      selected = unique(seurat_obj()@meta.data[[input$tissue_column]])[2])
    
    # Update Additional Genes Selection for Heatmap
    all_genes <- rownames(seurat_obj())
    updateSelectizeInput(session, "comp_additional_genes_heatmap",
                         choices = all_genes, server = TRUE)
  })
  
  # **Volcano Plot Generation**
  observeEvent(input$generate_volcano, {
    req(seurat_obj())
    req(input$comp_cluster_volcano)
    
    # Disable download buttons until plot is generated
    disable("download_comp_volcano")
    disable("download_comp_table_volcano")
    
    # Determine Comparison Mode
    comparison_mode <- input$comp_mode_volcano
    
    # Subset to Selected Cluster for Volcano Plot
    seurat_volcano <- subset(seurat_obj(), idents = input$comp_cluster_volcano)
    
    # Set Idents Based on Comparison Mode
    if (comparison_mode == "groups") {
      # Use Tissue Groups as Idents
      validate(
        need(!is.null(seurat_volcano@meta.data$tissue_group),
             "Tissue groups are not defined in the metadata.")
      )
      Idents(seurat_volcano) <- seurat_volcano@meta.data$tissue_group
      ident1 <- input$comp_volcano_group1
      ident2 <- input$comp_volcano_group2
      
      # Validate that selected groups exist
      validate(
        need(ident1 %in% levels(Idents(seurat_volcano)), "Selected Tissue Group 1 does not exist."),
        need(ident2 %in% levels(Idents(seurat_volcano)), "Selected Tissue Group 2 does not exist.")
      )
      
    } else if (comparison_mode == "specific") {
      # Use Specific Tissues as Idents
      ident1 <- input$comp_volcano_tissue1
      ident2 <- input$comp_volcano_tissue2
      
      # Validate that selected tissues exist
      validate(
        need(all(c(ident1, ident2) %in% seurat_volcano@meta.data[[input$tissue_column]]),
             "One or both selected tissues do not exist in the selected cluster.")
      )
      
      Idents(seurat_volcano) <- seurat_volcano@meta.data[[input$tissue_column]]
    }
    
    # Perform Differential Expression Analysis with User-Selected Test
    de_markers_volcano <- tryCatch({
      FindMarkers(seurat_volcano,
                  ident.1 = ident1,
                  ident.2 = ident2,
                  test.use = input$comp_volcano_test_use,
                  min.pct = 0.25,
                  logfc.threshold = input$comp_logfc_volcano)
    }, error = function(e) {
      showNotification(paste("Error during differential expression analysis:", e$message), type = "error")
      return(NULL)
    })
    
    req(de_markers_volcano)
    
    # Determine the P-value column based on user selection
    p_val_column <- input$comp_volcano_pval_type
    
    # Ensure that the selected p-value column exists
    if (!(p_val_column %in% colnames(de_markers_volcano))) {
      showNotification(paste("Selected p-value type", p_val_column, "does not exist in the DE results."), type = "error")
      return(NULL)
    }
    
    # Filter DE Results Based on Cutoffs
    de_filtered_volcano <- de_markers_volcano %>%
      filter(get(p_val_column) < input$comp_pval_volcano, abs(avg_log2FC) > input$comp_logfc_volcano)
    
    # Store DE Results
    reactive_values$comp_de_results_volcano <- de_filtered_volcano
    
    # Prepare Data for EnhancedVolcano
    volcano_data <- de_markers_volcano
    volcano_data$Gene <- rownames(volcano_data)
    
    # Generate Enhanced Volcano Plot
    p_volcano <- EnhancedVolcano(volcano_data,
                                 lab = volcano_data$Gene,
                                 x = 'avg_log2FC',
                                 y = p_val_column,
                                 pCutoff = input$comp_pval_volcano,
                                 FCcutoff = input$comp_logfc_volcano,
                                 title = paste("Volcano Plot: Cluster", input$comp_cluster_volcano),
                                 subtitle = paste(ident1, "vs", ident2),
                                 legendPosition = 'right',
                                 legendLabSize = 12,
                                 legendIconSize = 4.0)
    
    # Store the plot object
    reactive_values$comp_volcano_plot <- p_volcano
    
    # Render Volcano Plot
    output$comp_volcano_plot <- renderPlot({
      req(reactive_values$comp_volcano_plot)
      reactive_values$comp_volcano_plot
    })
    
    # Render DE Genes Table for Volcano Plot
    output$comp_de_table_volcano <- DT::renderDataTable({
      req(reactive_values$comp_de_results_volcano)
      de_table_volcano <- reactive_values$comp_de_results_volcano
      
      # Convert rownames to a column named "Gene"
      de_table_volcano <- tibble::rownames_to_column(de_table_volcano, var = "Gene") %>%
        arrange(get(input$comp_volcano_pval_type))
      
      DT::datatable(de_table_volcano,
                    options = list(
                      pageLength = input$comp_volcano_table_entries,  # User-specified number of entries
                      scrollX = TRUE,
                      dom = 'Bfrtip',  # Include buttons in the table
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                    ),
                    rownames = FALSE,
                    extensions = 'Buttons',
                    filter = 'top',
                    caption = 'Differentially Expressed Genes (Volcano Plot)')
    })
    
    # Enable download buttons after generating the Volcano Plot and Table
    enable("download_comp_volcano")
    enable("download_comp_table_volcano")
  })
  
  # ------------- Cluster-comparison Heatmap
  observeEvent(input$generate_heatmap, {
    req(seurat_obj())
    req(input$comp_cluster_heatmap)
    
    disable("download_comp_heatmap")
    disable("download_comp_table_heatmap")
    
    # ---- subset to the chosen cluster 
    seurat_heatmap <- subset(seurat_obj(), idents = input$comp_cluster_heatmap)
    
    # ---- set Idents according to comparison mode 
    if (input$comp_mode_heatmap == "groups") {
      validate(need(!is.null(seurat_heatmap$tissue_group),
                    "Tissue groups are not defined in the metadata."))
      Idents(seurat_heatmap) <- seurat_heatmap$tissue_group
    } else {
      Idents(seurat_heatmap) <- seurat_heatmap[[input$tissue_column]][, 1]
    }
    
    # ---- differential expression 
    de_markers_heatmap <- FindAllMarkers(
      seurat_heatmap,
      only.pos        = FALSE,
      test.use        = input$comp_heatmap_test_use,
      min.pct         = 0.25,
      logfc.threshold = input$comp_logfc_heatmap
    )
    
    pval_col <- input$comp_heatmap_pval_type
    validate(need(pval_col %in% colnames(de_markers_heatmap),
                  "Chosen p-value column is missing."))
    
    de_filtered <- de_markers_heatmap %>%
      filter(get(pval_col) < input$comp_pval_heatmap,
             abs(avg_log2FC) > input$comp_logfc_heatmap)
    
    reactive_values$comp_de_results_heatmap <- de_filtered
    
    # ---- assemble gene list 
    top_markers <- de_filtered %>%
      group_by(cluster) %>%
      top_n(n = input$comp_top_n_markers, wt = avg_log2FC) %>%
      ungroup()
    
    gene_list <- unique(c(top_markers$gene,
                          input$comp_additional_genes_heatmap))
    validate(need(length(gene_list) > 0,
                  "No genes meet the specified thresholds."))
    
    # ---- make sure they are scaled
    da            <- DefaultAssay(seurat_heatmap)
    scaled_genes  <- rownames(GetAssayData(seurat_heatmap, assay = da,
                                           slot  = "scale.data"))
    missing_feats <- setdiff(gene_list, scaled_genes)
    
    if (length(missing_feats) > 0) {
      message("Scaling ", length(missing_feats), " gene(s)…")
      seurat_heatmap <- ScaleData(seurat_heatmap,
                                  features = missing_feats,
                                  assay    = da,
                                  verbose  = FALSE)
    }
    
    # ---- plot 
    p_heatmap <- DoHeatmap(seurat_heatmap, features = gene_list) +
      ggtitle(paste("Heatmap – Cluster", input$comp_cluster_heatmap)) +
      NoLegend()
    
    reactive_values$comp_heatmap_plot <- p_heatmap
    
    output$comp_heatmap_plot <- renderPlot({ p_heatmap })
    
    # ---- table 
    output$comp_de_table_heatmap <- DT::renderDataTable({
      de_filtered %>%
        tibble::rownames_to_column("Gene") %>%
        arrange(get(pval_col)) %>%
        DT::datatable(options = list(
          pageLength = input$comp_heatmap_table_entries,
          scrollX    = TRUE,
          dom        = 'Bfrtip',
          buttons    = c('copy','csv','excel','pdf','print')
        ),
        rownames  = FALSE,
        extensions = 'Buttons',
        filter     = 'top',
        caption    = 'Differentially Expressed Genes (Heatmap)')
    })
    
    enable("download_comp_heatmap")
    enable("download_comp_table_heatmap")
  })
  

  ########################## SECTION 17. Merge Clusters Functionality ###########
  
  # Update Merge Cluster Selection Inputs Whenever Clusters Change
  observe({
    req(seurat_obj())
    clusters <- levels(seurat_obj())
    
    updateSelectInput(session, "merge_from",
                      choices = clusters,
                      selected = NULL)
    
    updateSelectInput(session, "merge_into",
                      choices = clusters,
                      selected = NULL)
  })
  
  # Helper Function to Synchronize 'ident' with Active Identities
  synchronize_ident <- function(obj) {
    obj$ident <- as.character(Idents(obj))
    return(obj)
  }
  
  # Observe Event to Perform Cluster Merge
  observeEvent(input$merge_clusters, {
    req(seurat_obj())
    merge_from <- input$merge_from
    merge_into <- input$merge_into
    
    # Validation Checks
    if (is.null(merge_from) || merge_from == "") {
      showNotification("Please select a cluster to merge from.", type = "error")
      return(NULL)
    }
    
    if (is.null(merge_into) || merge_into == "") {
      showNotification("Please select a cluster to merge into.", type = "error")
      return(NULL)
    }
    
    if (merge_from == merge_into) {
      showNotification("Source and target clusters must be different.", type = "error")
      return(NULL)
    }
    
    # Confirm Merge Action
    showModal(modalDialog(
      title = "Confirm Cluster Merge",
      paste("Are you sure you want to merge cluster", merge_from, "into cluster", merge_into, "?"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_merge", "Confirm")
      )
    ))
  })
  
  # Perform the Merge After Confirmation
  observeEvent(input$confirm_merge, {
    removeModal()
    
    req(seurat_obj())
    merge_from <- input$merge_from
    merge_into <- input$merge_into
    
    # Retrieve the current Seurat object
    obj <- seurat_obj()
    
    # Get current identities
    current_idents <- Idents(obj)
    
    # Perform the merge by setting the identities of 'merge_from' to 'merge_into'
    new_idents <- as.character(current_idents)
    new_idents[new_idents == merge_from] <- merge_into
    Idents(obj) <- new_idents
    
    # **Synchronize the 'ident' column with the new active identities**
    obj <- synchronize_ident(obj)
    
    # Update the Seurat object
    seurat_obj(obj)
    
    # Notify the user
    showNotification(paste("Cluster", merge_from, "has been merged into cluster", merge_into, "."), type = "message")
    
    # Update all cluster-related UI elements to reflect the new cluster levels
    updated_clusters <- levels(seurat_obj())
    
    # List of all cluster-related input IDs to update
    cluster_inputs <- c("cluster_select", "ident1", "ident2", "cluster_to_rename",
                        "stats_cluster", "merge_from", "merge_into")
    
    for (input_id in cluster_inputs) {
      updateSelectInput(session, input_id,
                        choices = updated_clusters)
    }
    
    # Update selectizeInputs that allow multiple selections
    updateSelectizeInput(session, "clusters_to_subset",
                         choices = updated_clusters, server = TRUE)
    
    updateSelectizeInput(session, "comparison_clusters",
                         choices = updated_clusters, server = TRUE)
    
    updateSelectizeInput(session, "cells_highlight",
                         choices = updated_clusters, server = TRUE)
    
    # Update any other UI elements that depend on cluster levels as needed
    # For example:
    metadata_columns <- c("ident", colnames(seurat_obj()@meta.data))
    updateSelectInput(session, "group_by",
                      choices = metadata_columns, selected = "ident")
    
    # Refresh current clusters display if applicable
    output$current_clusters <- renderPrint({
      levels(seurat_obj())
    })
  })
  
  ########################## SECTION 18. CLUSTER COMPOSITION #########################
  
  # 18-A  Refresh the three pickers whenever a Seurat object is loaded 
  observeEvent(seurat_obj(), {
    req(seurat_obj())
    
    ## -- individual tissues 
    updateSelectizeInput(
      session, "composition_tissues",
      choices = unique(seurat_obj()@meta.data[[input$tissue_column]]),
      server  = TRUE
    )
    
    ## -- tissue groups
    if ("tissue_group" %in% colnames(seurat_obj()@meta.data)) {
      tg <- unique(seurat_obj()@meta.data$tissue_group)
      updateSelectizeInput(session, "composition_tissue_group",
                           choices = tg, selected = tg[1], server = TRUE)
    } else {
      updateSelectizeInput(session, "composition_tissue_group",
                           choices = "No tissue groups defined",
                           selected = "No tissue groups defined", server = TRUE)
    }
    
    ## -- clusters (incl. sub-clusters)
    updateSelectizeInput(
      session, "composition_clusters",
      choices = levels(Idents(seurat_obj())), server = TRUE
    )
  })
  
  ## 18-B  Render the stacked-bar composition plot
  output$composition_plot <- renderPlot({
    req(seurat_obj())
    
    ## ---- 1. subset rows according to user choice 
    if (input$composition_group_type == "individual") {
      req(input$composition_tissues)
      plot_data <- seurat_obj()@meta.data %>%
        filter(.data[[input$tissue_column]] %in% input$composition_tissues) %>%
        mutate(Group = .data[[input$tissue_column]])
    } else {                                           # tissue groups
      req(input$composition_tissue_group)
      plot_data <- seurat_obj()@meta.data %>%
        filter(tissue_group %in% input$composition_tissue_group) %>%
        mutate(Group = tissue_group)
    }
    
    ## ---- 2. add cluster identities (length-matched) 
    plot_data <- plot_data %>%
      mutate(Current_Ident = Idents(seurat_obj())[rownames(.)])
    
    ## ---- 3. optional cluster filter 
    if (length(input$composition_clusters)) {
      plot_data <- plot_data %>%
        filter(Current_Ident %in% input$composition_clusters)
    }
    
    validate(need(nrow(plot_data) > 0,
                  "No data available for the selected options."))
    
    ## ---- 4. summarise counts & proportions 
    summary_data <- plot_data %>%
      group_by(Group, Current_Ident) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Group) %>%
      mutate(Proportion = Count / sum(Count)) %>%
      ungroup()
    
    ## ---- 5. build a colour palette 
    palette_name <- input$composition_color_palette
    n_clu        <- length(unique(summary_data$Current_Ident))
    
    if (palette_name %in% rownames(brewer.pal.info)) {
      maxc  <- brewer.pal.info[palette_name, "maxcolors"]
      colours <- if (n_clu > maxc)
        colorRampPalette(brewer.pal(maxc, palette_name))(n_clu)
      else
        brewer.pal(max(n_clu, 3), palette_name)
    } else if (palette_name %in% c("viridis","magma","plasma","inferno","cividis")) {
      colours <- viridis(n_clu, option = palette_name)
    } else {
      colours <- hue_pal()(n_clu)
    }
    
    ## ---- 6. draw the plot 
    p <- ggplot(summary_data,
                aes(x = Group, y = Proportion, fill = Current_Ident)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = colours, name = "Cluster") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(title = "Cluster Composition Across Groups",
           x = ifelse(input$composition_group_type == "individual",
                      "Tissue", "Tissue Group"),
           y = "Proportion of Cells") +
      theme_prism() +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1),
            plot.title   = element_text(hjust = 0.5))
    
    if (input$show_percentage) {
      p <- p + geom_text(aes(label = paste0(round(Proportion*100, 1), "%")),
                         position = position_stack(vjust = 0.5),
                         size = 5, colour = "white")
    }
    reactive_values$composition_plot <- p
    print(p)
  })

  ########################## SECTION 19. Plot/table DOWNLOAD HANDLER ##############
  
  # 19-A  A master list that maps internal object names → nice labels
  download_choices_all <- list(
    "UMAP plot"                 = "umap_plot",
    "Feature plots (PDF)"       = "feature_plots",
    "Volcano plot"              = "volcano_plot",
    "Volcano table (CSV)"       = "FindMarkers_table",
    "Heatmap plot"              = "heatmap_plot",
    "Heatmap table (CSV)"       = "FindAllMarkers_table",
    "Spatial plot"              = "spatial_plot",
    "Neighbour plot"            = "neighbor_plot",
    "Violin plot"               = "vln_plot",
    "Dot plot"                  = "dot_plot",
    "Cluster-cmp Volcano plot"  = "comp_volcano_plot",
    "Cluster-cmp Heatmap plot"  = "comp_heatmap_plot",
    "Cluster-cmp Volcano table" = "comp_de_results_volcano",
    "Cluster-cmp Heatmap table" = "comp_de_results_heatmap",
    "Cluster Composition plot"  = "composition_plot"
  )
  
  ## 19-B  When *any* item in reactive_values changes, refresh the drop-down
  observe({
    rv_list <- reactiveValuesToList(reactive_values)           # converts to normal list
    available_keys <- names(rv_list)[!vapply(rv_list, is.null, logical(1))]
    
    ## ‘feature_plots’ is a flag we set when the Feature-Plot button is pressed
    if (isTRUE(reactive_values$feature_plots)) available_keys <- c(available_keys, "feature_plots")
    
    current_choices <- download_choices_all[download_choices_all %in% available_keys]
    
    updateSelectInput(
      session, "download_what",
      choices = current_choices,
      selected = if (length(current_choices)) current_choices[[1]] else character(0)
    )
  })
  
  ## 19-C  Set the ‘feature_plots’ flag when the user generates Feature Plots
  observeEvent(input$plot_feature, {
    reactive_values$feature_plots <- TRUE
  })
  
  ## 19-D  Clear everything when a new object is loaded / reset
  observeEvent(input$load_data, {
    isolate({
      for (nm in names(reactiveValuesToList(reactive_values)))
        reactive_values[[nm]] <- NULL
      reactive_values$feature_plots <- NULL
    })
  })
  
  ## 19-E  Single download handler
  output$download_selected <- downloadHandler(
    
    ## --- file name
    filename = function() {
      choice <- input$download_what
      ext <- switch(choice,
                    "feature_plots"            = ".pdf",
                    "FindMarkers_table"        = ".csv",
                    "FindAllMarkers_table"     = ".csv",
                    "comp_de_results_volcano"  = ".csv",
                    "comp_de_results_heatmap"  = ".csv",
                    ".png"                                 # default for plots
      )
      paste0(choice, "_", Sys.Date(), ext)
    },
    
    ## --- file content
    content = function(file) {
      
      choice <- input$download_what
      obj    <- reactive_values[[choice]]
      
      ## special: multipage Feature-Plot PDF
      if (choice == "feature_plots") {
        validate(need(length(input$feature_genes) > 0,
                      "Select genes and click 'Plot Feature' first."))
        pdf(file)
        split_by <- if (input$split_by_feature == "None") NULL else input$split_by_feature
        for (feat in input$feature_genes)
          print(FeaturePlot(seurat_obj(), features = feat, split.by = split_by))
        dev.off()
        return(invisible())
      }
      
      ## tables
      if (is.data.frame(obj))
        return(write.csv(obj, file, row.names = FALSE))
      
      ## ggplots 
      if (inherits(obj, c("gg","gtable","patchwork"))) {
        width  <- ifelse(grepl("heatmap", choice), 12, 16)
        height <- ifelse(grepl("heatmap", choice), 16, 12)
        ggsave(file, plot = obj, width = width, height = height, dpi = 300)
        return(invisible())
      }
      
      ## fallback 
      showNotification("I don't know how to save this object.", type = "error")
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
