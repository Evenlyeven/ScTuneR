suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(patchwork)))

## ===== define a function of the analysis ===== ##
scTuneR <- function(seurat_obj, output_dir, sub_analysis = TRUE, splitby, integration = TRUE, pc_range, res_range) {
  ## ----- handle input file ----- ##
  if (sub_analysis) {
    # Expecting a single Seurat object file (.rds or .RData)
    if (grepl("\\.rds$", seurat_obj)) {
      seurat_obj <- readRDS(seurat_obj)
    } else if (grepl("\\.RData$", seurat_obj)) {
      temp_env <- new.env()
      load(seurat_obj, envir = temp_env)
      
      # Look for Seurat object in the environment
      seurat_vars <- ls(temp_env)
      seurat_objs <- Filter(function(x) inherits(temp_env[[x]], "Seurat"), seurat_vars)
      
      if (length(seurat_objs) == 0) {
        stop("No Seurat object found in the RData file.")
      } else if (length(seurat_objs) > 1) {
        warning("Multiple Seurat objects found. Using the first one: ", seurat_objs[1])
      }
      
      seurat_obj <- temp_env[[seurat_objs[1]]]
      message("Using Seurat object from .RData: ", seurat_objs[1])
    } else {
      stop("For sub-analysis, the input must be a .rds or .RData file containing a Seurat object.")
    }
    
    # Validate
    if (!inherits(seurat_obj, "Seurat")) {
      stop("The loaded object is not a Seurat object.")
    }
    
  } else {
    # Not a sub-analysis: expecting a list of Seurat objects
    if (grepl("\\.rds$", seurat_obj)) {
      seurat_obj <- readRDS(seurat_obj)
    } else if (grepl("\\.RData$", seurat_obj)) {
      temp_env <- new.env()
      load(seurat_obj, envir = temp_env)
      
      # Look for list of Seurat objects
      seurat_vars <- ls(temp_env)
      seurat_list_objs <- Filter(function(x)
        is.list(temp_env[[x]]) &&
          all(sapply(temp_env[[x]], function(obj) inherits(obj, "Seurat"))),
        seurat_vars
      )
      
      if (length(seurat_list_objs) == 0) {
        stop("No list of Seurat objects found in the RData file.")
      } else if (length(seurat_list_objs) > 1) {
        warning("Multiple Seurat object lists found. Using the first one: ", seurat_list_objs[1])
      }
      
      seurat_obj <- temp_env[[seurat_list_objs[1]]]
      message("Using list of Seurat objects from .RData: ", seurat_list_objs[1])
    } else {
      stop("For standard analysis, the input must be a .rds or .RData file containing a list of Seurat objects.")
    }
    
    # Validate list
    if (!is.list(seurat_obj) || !all(sapply(seurat_obj, function(obj) inherits(obj, "Seurat")))) {
      stop("The loaded object is not a list of Seurat objects.")
    }
    
    # ensure it's named
    if (is.null(names(seurat_obj)) || any(names(seurat_obj) == "")) {
      stop("The list of Seurat objects must be named (e.g., with sample IDs).")
    }
  }
  
  ## ----- handle output directory ----- ##
  # Wrap output in a timestamped folder, so repeated runs do not overwrite results unless desired
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- file.path(output_dir, paste0("ScTuneR_", timestamp))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## ----- perform sub-analysis if requested ----- ##
  if (sub_analysis) {
    # Perform sub-analysis on the Seurat object
    message("Performing sub-analysis on the Seurat object...")
    
    seurat_obj <- SplitObject(seurat_obj, split.by = splitby)
    message("Seurat object split by '", splitby, "'.")
  } else {
    message("Skipping sub-analysis as per user request.")
  }
  
  ## ----- normalization ----- ##
  message("Normalizing the Seurat object with SCTransformation...")
  #at this point, seurat_obj is a names list of seurat objects
  for (i in seq_along(names(seurat_obj))){
    DefaultAssay(seurat_obj[[i]]) <- "RNA"
    seurat_obj[[i]] %<>% SCTransform(
      vst.flavor = "v2",
      variable.features.n = 3000,
      return.only.var.genes = FALSE
    )
  }
  
  saveRDS(seurat_obj, file = file.path(output_dir, "SCTed_seu_obj_ls.rds"))
  message("Normalization completed and saved to ", file.path(output_dir, "SCTed_seu_obj_ls.rds"))
  
  ## ----- prep for integration ----- ##
  #we will need to merge regardless of whether we will integrate or not
  merged <- merge(seurat_obj[[1]], y = seurat_obj[2:length(seurat_obj)],
                  add.cell.ids = names(seurat_obj))
  
  #prep for DE analysis and visualization
  merged %<>% PrepSCTFindMarkers(assay = "SCT")
  
  #variable features
  VariableFeatures(merged) <- unique(unlist(lapply(seurat_obj[1:length(seurat_obj)], VariableFeatures)))
  
  #run PCA
  merged %<>% RunPCA()
  saveRDS(merged, file = file.path(output_dir, "merged_seu_obj.rds"))
  message("Merged Seurat object saved to ", file.path(output_dir, "merged_seu_obj.rds"))
  
  #elbow plot of PCs
  p <- ElbowPlot(merged, ndims = 50) +
    labs(title = "Elow plot of PCs")
  
  png(filename = file.path(output_dir, "merged_PCA.png"), width = 800, height = 300)
  plot(p)
  dev.off()
  message("PCA completed and elbow plot saved to ", file.path(output_dir, "merged_PCA.png"))
  
  ## ----- merge only or integration ----- ##
  if (integration) {
    message("Integrating the Seurat object...")
    merged %<>% IntegrateLayers(method = HarmonyIntegration,
                                normalization.method = "SCT")
    #the default pc used in harmony is 50, at least in this seurat wrapper
    saveRDS(merged, file = file.path(output_dir, "integrated_seu_obj.rds"))
    message("Integration completed with Harmony.")
    message("Integrated Seurat object saved to ", file.path(output_dir, "integrated_seu_obj.rds"))
  } else {
    message("Skipping integration as per user request.")
  }
  
  ## ----- exploration of PCs and resolutions ----- ##
  message("Exploring PCs and resolutions...")
  i_range <- as.numeric(unlist(strsplit(pc_range, ",")))#"30,40,50" → c("30", "40", "50") → c(30, 40, 50)
  j_range <- as.numeric(unlist(strsplit(res_range, ",")))
  
  #for dimensional reduction, if integration TRUE -> reduction = "harmony", FALSE -> "pca" 
  if (integration) {
    reduction_method <- "harmony"
  } else {
    reduction_method <- "pca"
  }
  seu_m <- merged
  plots_list <- list()
  for (i in i_range){
    for (j in j_range){
      set.seed(10086)
      
      seu_m.c <- seu_m %>%
        FindNeighbors(dims = 1:i, reduction = reduction_method) %>%
        FindClusters(resolution = j) %>%
        RunUMAP(dims = 1:i, reduction = reduction_method)
      
      plot_obj <- DimPlot(seu_m.c, reduction = "umap", raster = FALSE, label = TRUE, repel = TRUE) +
        labs(caption = paste0("PC=", as.character(i), " res=", as.character(j)))
      plots_list[[paste0("PC", i, "res", j)]] <- plot_obj
    }
  }
  p <- wrap_plots(plots_list, ncol = length(j_range))
  png(filename = file.path(output_dir, "UMAP_PCres.png"), width = 600*length(j_range), height = 500*length(i_range))
  print(p)
  dev.off()
  message("UMAP plots for different PCs and resolutions saved to ", file.path(output_dir, "UMAP_PCres.png"))
  message("GOODBYE! Analysis completed successfully. Results saved in: ", output_dir)
}

## ===== define options for the script ===== ##
description_text <- 'This script performs flexible scRNA-seq analysis with optional subcluster analysis and Harmony-based integration.

REQUIRED:
--seurat_obj      Path to a Seurat object (.rds or .RData) or a list of Seurat objects

OPTIONAL:
--output_dir      Output directory [default = current working directory]
--sub_analysis    Whether to split merged object by sample [default = TRUE]
--splitby         Metadata column to split on [default = "orig.ident"]
--integration     Whether to run Harmony integration [default = TRUE]
--pc_range        PCs to explore, comma-separated [default = "30,40,50"]
--res_range       Resolutions to explore, comma-separated [default = "0.2,0.4,0.6"]

Example usage:
Rscript subTuner.R --seurat_obj mydata.rds --pc_range "20,30" --res_range "0.2,0.6"'

option_list <- list(
  make_option(c("--seurat_obj"), type = "character", default = NULL,#default=NULL + check input chunk below made this option required
              help = "Path to the Seurat object file (.rds or .RData) or a named list of Seurat objects."),
  make_option(c("--output_dir"), type = "character", default = ".", 
              help = "Directory to save the output results. Default is current working directory."),
  make_option(c("--sub_analysis"), type = "logical", default = TRUE,
              help = "Perform sub-analysis on the Seurat object. Default is TRUE."),
  make_option(c("--splitby"), type = "character", default = "orig.ident",
              help = "Variable to split the Seurat object by (e.g., 'orig.ident'). Default is 'orig.ident'."),
  make_option(c("--integration"), type = "logical", default = TRUE,
              help = "Perform Harmony integration on the merged Seurat object. Default is TRUE."),
  make_option(c("--pc_range"), type = "character", default = "30,40,50",
              help = "Range of PCs to explore, e.g., '30,40,50'. Default is '30,40,50'."),
  make_option(c("--res_range"), type = "character", default = "0.2,0.4,0.6",
              help = "Range of resolutions to explore, e.g., '0.2,0.4,0.6'. Default is '0.2,0.4,0.6'.")
)

opt_parser <- OptionParser(option_list = option_list, description = description_text)
opt <- parse_args(opt_parser)

## ----- check the input parameters ----- ##
if (is.null(opt$seurat_obj)) {
  stop("Please provide a Seurat object file or a named list of Seurat objects using --seurat_obj.")
}

## ----- call the function and run the analysis ----- ##
scTuneR(
  seurat_obj = opt$seurat_obj,
  output_dir = opt$output_dir,
  sub_analysis = opt$sub_analysis,
  splitby = opt$splitby,
  integration = opt$integration,
  pc_range = opt$pc_range,
  res_range = opt$res_range
)