suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(patchwork)))

## ===== define a function of the analysis ===== ##
scTuneR <- function(seurat_obj,
                    output_dir,
                    sub_analysis = TRUE,
                    splitby,
                    integration = TRUE,
                    integration_method = "harmony",
                    pc_range,
                    res_range,
                    regress_cc = FALSE,
                    npcs = 50,
                    resume_from = NULL) {
  ## ----- validate integration_method ----- ##
  # Both methods take their batch grouping from the SCT models carried in the merged
  # object -- i.e. from the samples that were SCTransformed separately below, which is
  # the split defined by `splitby` (or by the names of the input list). Neither method
  # reads a metadata column: IntegrateLayers derives the groups itself via
  # CreateIntegrationGroups(), and HarmonyIntegration hard-codes vars_use to those
  # groups (it has no group.by.vars argument to override).
  #   "harmony" corrects the PCA embedding in place, writing reduction "harmony".
  #   "rpca"    finds reciprocal-PCA anchors between the SCT models and corrects the
  #             same PCA embedding, writing reduction "integrated.rpca".
  valid_integration <- c("harmony", "rpca")
  integration_method <- tolower(integration_method)
  if (!integration_method %in% valid_integration) {
    stop("integration_method must be one of: ",
         paste(valid_integration, collapse = ", "))
  }
  # reduction each method writes into the object, used again in Stage 4
  integration_reduction <- c(harmony = "harmony", rpca = "integrated.rpca")[[integration_method]]

  ## ----- validate resume_from ----- ##
  # resume_from lets the user re-enter the pipeline at a checkpoint by pointing
  # --seurat_obj at a previously saved .rds, so a failed step does not require
  # restarting from the raw input. Valid checkpoints:
  #   "sct"        -> --seurat_obj is SCTed_seu_obj_ls.rds (named list of SCT'd objects)
  #   "merged"     -> --seurat_obj is merged_seu_obj.rds  (single object, PCA done)
  #   "integrated" -> --seurat_obj is integrated_seu_obj.rds (single object, integration done;
  #                   pass the same --integration_method that produced it)
  valid_resume <- c("sct", "merged", "integrated")
  if (!is.null(resume_from)) {
    resume_from <- tolower(resume_from)
    if (!resume_from %in% valid_resume) {
      stop("resume_from must be one of: ",
           paste(valid_resume, collapse = ", "))
    }
    message("Resuming pipeline from checkpoint: '", resume_from, "'.")
  }

  ## ----- parse and validate the exploration grid up front ----- ##
  # Parsed here rather than in Stage 4 so a bad --pc_range / --res_range / --npcs
  # combination fails in seconds, instead of after SCTransform, merge, integration
  # and two multi-GB .rds writes have already run.
  #non-numeric entries become NA and are reported by the checks below, so the coercion
  #warning is suppressed to keep the error message clean
  i_range <- suppressWarnings(as.numeric(unlist(strsplit(pc_range, ","))))  #"30,40,50" -> c(30, 40, 50)
  j_range <- suppressWarnings(as.numeric(unlist(strsplit(res_range, ","))))
  if (length(i_range) == 0 || anyNA(i_range) || any(i_range < 2)) {
    stop("pc_range must be a comma-separated list of numbers >= 2, e.g. '30,40,50'. Got: ",
         pc_range)
  }
  if (length(j_range) == 0 || anyNA(j_range) || any(j_range <= 0)) {
    stop("res_range must be a comma-separated list of positive numbers, e.g. '0.2,0.4'. Got: ",
         res_range)
  }
  #the PCs requested for exploration cannot exceed the dims computed by RunPCA/integration
  if (max(i_range) > npcs) {
    stop("pc_range requests ", max(i_range),
         " PCs but only ", npcs,
         " were computed. Increase --npcs or lower --pc_range.")
  }

  ## ----- helper to load a resume checkpoint (.rds or .RData) ----- ##
  # expect_list = TRUE  -> a named list of Seurat objects (the "sct" checkpoint)
  # expect_list = FALSE -> a single Seurat object ("merged"/"integrated" checkpoints)
  load_resume_obj <- function(path, expect_list) {
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
      obj <- readRDS(path)
    } else if (grepl("\\.RData$", path, ignore.case = TRUE)) {
      temp_env <- new.env()
      load(path, envir = temp_env)
      vars <- ls(temp_env)
      if (length(vars) == 0) {
        stop("No object found in resume file: ", path)
      }
      obj <- temp_env[[vars[1]]]
    } else {
      stop("Resume file must be a .rds or .RData file: ", path)
    }

    if (expect_list) {
      if (!is.list(obj) ||
          !all(sapply(obj, function(x) inherits(x, "Seurat")))) {
        stop("Expected a named list of Seurat objects in: ", path)
      }
      if (is.null(names(obj)) || any(names(obj) == "")) {
        stop("The list of Seurat objects must be named (e.g., with sample IDs): ",
             path)
      }
    } else {
      if (!inherits(obj, "Seurat")) {
        stop("Expected a single Seurat object in: ", path)
      }
    }
    obj
  }

  # Anchor features are chosen in Stage 2. Declared here because Stage 2 is skipped when
  # resuming from the "merged"/"integrated" checkpoints, and Stage 3 falls back to the
  # VariableFeatures already stored on the checkpoint object in that case.
  integration_features <- NULL

  ## ----- handle output directory ----- ##
  # Wrap output in a timestamped folder, so repeated runs do not overwrite results unless desired
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- file.path(output_dir, paste0("ScTuneR_", timestamp))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ## ===== Stage 1: load raw input, split, and normalize ===== ##
  # Skipped entirely when resuming from any checkpoint.
  if (is.null(resume_from)) {
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
        seurat_objs <- Filter(function(x)
          inherits(temp_env[[x]], "Seurat"), seurat_vars)

        if (length(seurat_objs) == 0) {
          stop("No Seurat object found in the RData file.")
        } else if (length(seurat_objs) > 1) {
          warning("Multiple Seurat objects found. Using the first one: ",
                  seurat_objs[1])
        }

        seurat_obj <- temp_env[[seurat_objs[1]]]
        message("Using Seurat object from .RData: ", seurat_objs[1])
      } else {
        stop(
          "For sub-analysis, the input must be a .rds or .RData file containing a Seurat object."
        )
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
            all(sapply(temp_env[[x]], function(obj)
              inherits(obj, "Seurat"))),
          seurat_vars)

        if (length(seurat_list_objs) == 0) {
          stop("No list of Seurat objects found in the RData file.")
        } else if (length(seurat_list_objs) > 1) {
          warning("Multiple Seurat object lists found. Using the first one: ",
                  seurat_list_objs[1])
        }

        seurat_obj <- temp_env[[seurat_list_objs[1]]]
        message("Using list of Seurat objects from .RData: ",
                seurat_list_objs[1])
      } else {
        stop(
          "For standard analysis, the input must be a .rds or .RData file containing a list of Seurat objects."
        )
      }

      # Validate list
      if (!is.list(seurat_obj) ||
          !all(sapply(seurat_obj, function(obj)
            inherits(obj, "Seurat")))) {
        stop("The loaded object is not a list of Seurat objects.")
      }

      # ensure it's named
      if (is.null(names(seurat_obj)) ||
          any(names(seurat_obj) == "")) {
        stop("The list of Seurat objects must be named (e.g., with sample IDs).")
      }
    }

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
    for (i in seq_along(names(seurat_obj))) {
      DefaultAssay(seurat_obj[[i]]) <- "RNA"

      if (regress_cc) {
        message("Scoring cell cycle genes (human) for sample: ",
                names(seurat_obj)[i])
        seurat_obj[[i]] <- NormalizeData(seurat_obj[[i]])  # needed for CellCycleScoring
        seurat_obj[[i]] <- CellCycleScoring(
          seurat_obj[[i]],
          s.features   = cc.genes$s.genes,
          g2m.features = cc.genes$g2m.genes
        )
        seurat_obj[[i]] %<>% SCTransform(
          vst.flavor = "v2",
          variable.features.n = 3000,
          return.only.var.genes = FALSE,
          vars.to.regress = c("S.Score", "G2M.Score")
        )
      } else {
        seurat_obj[[i]] %<>% SCTransform(
          vst.flavor = "v2",
          variable.features.n = 3000,
          return.only.var.genes = FALSE
        )
      }
    }

    saveRDS(seurat_obj, file = file.path(output_dir, "SCTed_seu_obj_ls.rds"))
    message(
      "Normalization completed and saved to ",
      file.path(output_dir, "SCTed_seu_obj_ls.rds")
    )
  }

  ## ===== Stage 2: merge + PCA ===== ##
  # Runs on a full pipeline or when resuming from the "sct" checkpoint.
  if (is.null(resume_from) || resume_from == "sct") {
    if (identical(resume_from, "sct")) {
      # --seurat_obj points at SCTed_seu_obj_ls.rds (named list of SCT'd objects)
      seurat_obj <- load_resume_obj(seurat_obj, expect_list = TRUE)
      message("Loaded SCTed Seurat object list from checkpoint.")
    }

    ## ----- prep for integration ----- ##
    #we will need to merge regardless of whether we will integrate or not
    merged <- merge(seurat_obj[[1]],
                    y = seurat_obj[2:length(seurat_obj)],
                    add.cell.ids = names(seurat_obj))

    #prep for DE analysis and visualization
    merged %<>% PrepSCTFindMarkers(assay = "SCT")

    #variable features: rank the per-sample SCT variable features by how many
    #samples share them and keep the top 3000, rather than taking the union. This
    #keeps the PCA below driven by genes that are variable across samples instead
    #of sample-specific ones, and it is the feature set passed explicitly to
    #IntegrateLayers in Stage 3 as the genes RPCA finds anchors on. Note that
    #setting VariableFeatures() alone is NOT enough: for an SCTAssay,
    #IntegrateLayers overwrites `features` with SelectSCTIntegrationFeatures()
    #unless `features` is supplied, so it has to be handed over there directly.
    integration_features <- SelectIntegrationFeatures(
      object.list = seurat_obj,
      nfeatures = 3000,
      assay = rep("SCT", length(seurat_obj))
    )

    VariableFeatures(merged) <- integration_features

    #run PCA on those features (npcs controls how many PCs are computed, so a wider
    #pc_range can be explored downstream; it is also the number of dims Harmony
    #corrects below, so raise it only when a larger pc_range actually needs it)
    merged %<>% RunPCA(
      assay = "SCT",
      features = integration_features,
      npcs = npcs
    )
    saveRDS(merged, file = file.path(output_dir, "merged_seu_obj.rds"))
    message("Merged Seurat object saved to ",
            file.path(output_dir, "merged_seu_obj.rds"))

    #elbow plot of PCs (cap ndims at npcs so ElbowPlot never exceeds computed PCs)
    p <- ElbowPlot(merged, ndims = min(50, npcs)) +
      labs(title = "Elbow plot of PCs")

    png(
      filename = file.path(output_dir, "merged_PCA.png"),
      width = 800,
      height = 300
    )
    plot(p)
    dev.off()
    message("PCA completed and elbow plot saved to ",
            file.path(output_dir, "merged_PCA.png"))
  }

  ## ===== Stage 3: merge only or integration ===== ##
  # Runs on a full pipeline or when resuming from "sct"/"merged".
  if (is.null(resume_from) || resume_from %in% c("sct", "merged")) {
    if (identical(resume_from, "merged")) {
      # --seurat_obj points at merged_seu_obj.rds (single object with PCA done)
      merged <- load_resume_obj(seurat_obj, expect_list = FALSE)
      message("Loaded merged Seurat object from checkpoint.")
    }

    if (integration) {
      message("Integrating the Seurat object with '",
              integration_method,
              "'...")

      #Stage 2 is skipped when resuming from "merged", so fall back to the feature set
      #already stored on the checkpoint object.
      if (is.null(integration_features)) {
        integration_features <- VariableFeatures(merged)
      }
      if (length(integration_features) == 0) {
        stop("No variable features found on the merged object; cannot pick anchor features.")
      }
      #report the groups the correction will actually use, since neither method takes a
      #batch argument (levels() of an SCTAssay lists its per-sample SCT models)
      sct_models <- levels(merged[[DefaultAssay(merged)]])
      message("Batch grouping comes from the ", length(sct_models),
              " SCT models in the merged object: ",
              paste(sct_models, collapse = ", "))

      if (integration_method == "harmony") {
        #HarmonyIntegration ignores both `npcs` and `features`, and always corrects every
        #dimension of `orig.reduction` ("pca"), so how many dims Harmony sees is set
        #by RunPCA(npcs = npcs) above, not by anything passed here. Its batch variable is
        #not settable: it hard-codes vars_use to the SCT-model groups reported above.
        merged %<>% IntegrateLayers(
          method = HarmonyIntegration,
          normalization.method = "SCT",
          new.reduction = integration_reduction
        )
      } else {
        #RPCA anchors are found between the separately-SCTransformed samples that the
        #merge brought in. `features` is passed explicitly so the anchors use the same
        #3000 genes the PCA above was built on; left to itself IntegrateLayers would
        #call SelectSCTIntegrationFeatures() and ignore VariableFeatures(merged).
        #dims and k.filter are left at RPCAIntegration's defaults (anchors in dims 1:30,
        #k.filter NA, i.e. no anchor filtering -- note FindIntegrationAnchors called
        #directly defaults to k.filter 200 instead). This does not limit pc_range: the
        #corrected embedding keeps every dimension of the PCA it corrects
        #(dims.to.integrate defaults to all `npcs` of them), and the per-cell anchor
        #weights are computed in that full space, not in the 30 anchor dims. Only the
        #choice of which cell pairs become anchors is made in dims 1:30.
        merged %<>% IntegrateLayers(
          method = RPCAIntegration,
          normalization.method = "SCT",
          features = integration_features,
          new.reduction = integration_reduction
        )
      }

      saveRDS(merged, file = file.path(output_dir, "integrated_seu_obj.rds"))
      message("Integration completed with '",
              integration_method,
              "'; embedding stored in reduction '",
              integration_reduction,
              "'.")
      message(
        "Integrated Seurat object saved to ",
        file.path(output_dir, "integrated_seu_obj.rds")
      )
    } else {
      message("Skipping integration as per user request.")
    }
  }

  ## ===== Stage 4: exploration of PCs and resolutions ===== ##
  if (identical(resume_from, "integrated")) {
    # --seurat_obj points at integrated_seu_obj.rds (single object, integration done)
    merged <- load_resume_obj(seurat_obj, expect_list = FALSE)
    message("Loaded integrated Seurat object from checkpoint.")
  }

  #i_range / j_range were parsed and checked against npcs at the top of the function
  message("Exploring PCs and resolutions...")

  #for dimensional reduction, if integration TRUE -> the reduction written by the
  #chosen method ("harmony" or "integrated.rpca"), FALSE -> "pca"
  if (integration) {
    reduction_method <- integration_reduction
  } else {
    reduction_method <- "pca"
  }

  #when resuming from the "integrated" checkpoint the object was built by whichever
  #method was used then, so fail early (and informatively) on a mismatch
  if (!reduction_method %in% Reductions(merged)) {
    stop("Reduction '", reduction_method,
         "' not found in the object. Available: ",
         paste(Reductions(merged), collapse = ", "),
         ". Check --integration / --integration_method.")
  }

  #the corrected embedding may hold fewer dims than npcs (e.g. an object carried
  #over from an earlier run), so check against what is actually there
  n_dims_avail <- ncol(Embeddings(merged, reduction = reduction_method))
  if (max(i_range) > n_dims_avail) {
    stop("pc_range requests ", max(i_range),
         " dims but reduction '", reduction_method,
         "' only has ", n_dims_avail, ".")
  }

  seu_m <- merged
  plots_list <- list()
  for (i in i_range) {
    for (j in j_range) {
      set.seed(10086)

      seu_m.c <- seu_m %>%
        FindNeighbors(dims = 1:i, reduction = reduction_method) %>%
        FindClusters(resolution = j) %>%
        RunUMAP(dims = 1:i, reduction = reduction_method)

      plot_obj <- DimPlot(
        seu_m.c,
        reduction = "umap",
        raster = FALSE,
        label = TRUE,
        repel = TRUE
      ) +
        labs(caption = paste0(
          "PC=", as.character(i),
          " res=", as.character(j),
          " reduction=", reduction_method
        ))
      plots_list[[paste0("PC", i, "res", j)]] <- plot_obj
    }
  }
  p <- wrap_plots(plots_list, ncol = length(j_range))
  png(
    filename = file.path(output_dir, "UMAP_PCres.png"),
    width = 600 * length(j_range),
    height = 500 * length(i_range)
  )
  print(p)
  dev.off()
  message(
    "UMAP plots for different PCs and resolutions saved to ",
    file.path(output_dir, "UMAP_PCres.png")
  )
  #remind the user which reduction the follow-up ScHerdeR run has to be pointed at
  message("Downstream: run ScHerdeR with --reduction '", reduction_method, "'.")
  message("GOODBYE! Analysis completed successfully. Results saved in: ",
          output_dir)
}

## ===== define options for the script ===== ##
description_text <- 'This script performs flexible scRNA-seq analysis with optional subcluster analysis and batch integration (Harmony or RPCA).

REQUIRED:
--seurat_obj      Path to a Seurat object (.rds or .RData) or a list of Seurat objects.
                  When --resume_from is set, this must point to the corresponding
                  checkpoint .rds file (see --resume_from).

OPTIONAL:
--output_dir      Output directory [default = current working directory]
--sub_analysis    Whether to split merged object by sample [default = TRUE]
--splitby         Metadata column to split on [default = "orig.ident"]
--integration     Whether to run integration [default = TRUE]
--integration_method  Integration method to use [default = "harmony"]. One of:
                    "harmony"  Harmony correction of the PCA embedding. Writes
                               reduction "harmony".
                    "rpca"     Seurat reciprocal-PCA anchor integration, correcting the
                               same PCA embedding. Writes reduction "integrated.rpca".
                  Both methods take their batch grouping from the separately-
                  SCTransformed samples (the split set by --splitby, or the names of
                  the input list); neither reads a metadata column, and the Harmony
                  batch variable is not settable. The run logs the groups it used.
                  The resulting reduction name is what a follow-up ScHerdeR run
                  must be given via its --reduction option.
--npcs            Number of PCs computed by RunPCA [default = 50]. Must be >= the
                  largest value in --pc_range. Harmony corrects all of them, so raise
                  this only when a larger --pc_range needs it (e.g. --npcs 100).
--pc_range        PCs to explore, comma-separated [default = "30,40,50"]
--res_range       Resolutions to explore, comma-separated [default = "0.2,0.4,0.6"]
--regress_cc      Regress out cell cycle scores during SCTransform [default = FALSE]
                  NOTE: Human only. Uses Seurat built-in cc.genes (S and G2M gene sets).
                  Recommended when cell cycle is a confounding factor rather than
                  a biological question of interest.
--resume_from     Resume from a previously saved checkpoint instead of starting from
                  the raw input [default = NULL, i.e. full pipeline]. One of:
                    "sct"        --seurat_obj = SCTed_seu_obj_ls.rds (skip load/split/SCTransform)
                    "merged"     --seurat_obj = merged_seu_obj.rds   (skip through PCA)
                    "integrated" --seurat_obj = integrated_seu_obj.rds (skip through integration)
                                 Pass the same --integration_method used to build it.

Example usage:
Rscript scTuner.R --seurat_obj mydata.rds --pc_range "20,30" --res_range "0.2,0.6"
Rscript scTuner.R --seurat_obj mydata.rds --integration_method rpca --pc_range "20,30"
Rscript scTuner.R --seurat_obj SCTed_seu_obj_ls.rds --resume_from sct --pc_range "20,30"'

option_list <- list(
  make_option(
    c("--seurat_obj"),
    type = "character",
    default = NULL,
    #default=NULL + check input chunk below made this option required
    help = "Path to the Seurat object file (.rds or .RData) or a named list of Seurat objects."
  ),
  make_option(
    c("--output_dir"),
    type = "character",
    default = ".",
    help = "Directory to save the output results. Default is current working directory."
  ),
  make_option(
    c("--sub_analysis"),
    type = "logical",
    default = TRUE,
    help = "Perform sub-analysis on the Seurat object. Default is TRUE."
  ),
  make_option(
    c("--splitby"),
    type = "character",
    default = "orig.ident",
    help = "Variable to split the Seurat object by (e.g., 'orig.ident'). Default is 'orig.ident'."
  ),
  make_option(
    c("--integration"),
    type = "logical",
    default = TRUE,
    help = "Perform integration on the merged Seurat object. Default is TRUE."
  ),
  make_option(
    c("--integration_method"),
    type = "character",
    default = "harmony",
    help = "Integration method: 'harmony' (reduction 'harmony') or 'rpca' (Seurat
                    reciprocal PCA, reduction 'integrated.rpca'). Both group cells by
                    the separately-SCTransformed samples set by --splitby.
                    Default is 'harmony'."
  ),
  make_option(
    c("--npcs"),
    type = "integer",
    default = 50,
    help = "Number of PCs computed by RunPCA. Must be >= the largest value in
                    --pc_range. Harmony corrects all of them, so raise this only when
                    a larger --pc_range needs it. Default is 50."
  ),
  make_option(
    c("--pc_range"),
    type = "character",
    default = "30,40,50",
    help = "Range of PCs to explore, e.g., '30,40,50'. Default is '30,40,50'."
  ),
  make_option(
    c("--res_range"),
    type = "character",
    default = "0.2,0.4,0.6",
    help = "Range of resolutions to explore, e.g., '0.2,0.4,0.6'. Default is '0.2,0.4,0.6'."
  ),
  make_option(
    c("--regress_cc"),
    type = "logical",
    default = FALSE,
    help = "Regress out cell cycle scores (S.Score, G2M.Score) during SCTransform.
                    Human only — uses Seurat built-in cc.genes. Default is FALSE."
  ),
  make_option(
    c("--resume_from"),
    type = "character",
    default = NULL,
    help = "Resume from a saved checkpoint: 'sct' (SCTed_seu_obj_ls.rds),
                    'merged' (merged_seu_obj.rds), or 'integrated' (integrated_seu_obj.rds).
                    --seurat_obj must point to the matching file. Default is NULL (full pipeline)."
  )
)

opt_parser <- OptionParser(option_list = option_list, description = description_text)
opt <- parse_args(opt_parser)

## ----- check the input parameters ----- ##
if (is.null(opt$seurat_obj)) {
  stop(
    "Please provide a Seurat object file or a named list of Seurat objects using --seurat_obj."
  )
}

## ----- call the function and run the analysis ----- ##
scTuneR(
  seurat_obj = opt$seurat_obj,
  output_dir = opt$output_dir,
  sub_analysis = opt$sub_analysis,
  splitby = opt$splitby,
  integration = opt$integration,
  integration_method = opt$integration_method,
  pc_range = opt$pc_range,
  res_range = opt$res_range,
  regress_cc  = opt$regress_cc,
  npcs = opt$npcs,
  resume_from = opt$resume_from
)
