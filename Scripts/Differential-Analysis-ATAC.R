library(ArchR)

set.seed(1)

#########################
#########################

#' This function computes a pairwise differential test between two groups
#' @param se summarized experiment
#' @param compareCol a columnname in colData(se) where the comparisons will be made from
#' @param topGroup vector of names to compare top (ie. top / bottom) from colData(se)[,compareCol]
#' @param bottomGroup vector of names to compare bottom (ie. top / bottom) from colData(se)[,compareCol]
#' @param method method for differential testing : edgeR
#' @param comparisonName name of comparison
#' @param lfcThreshold log2FoldChange Threshold to test default is 0
#' @param cooksCutoff see DESeq2 ??results but cooksCutoff is for filtering default is Inf
#' @param fitType for DESeq2 params parametric, mean, local for dispersion fit
#' @param threads number of threads to use
#' @param BPPARAM this function will attempt to use sock for windows and multicoreparam if it does not work just set this to Windows and sock will automatically be ran
#' @param housekeeping a vector of integers specifying which rows to be used as normalization available for DESeq2 and Voom
#'
pairwiseTest <- function(se, compareCol = "Group", topGroup, bottomGroup, method = "edgeR", compareTop = NULL,
  comparisonName = "Differential Test", fitType = "parametric", threads = 1, BPPARAM = NULL, round = FALSE,
  housekeeping = NULL, lfcThreshold = 0, cooksCutoff = Inf){

  stopifnot("counts" %in% assayNames(se))
  assays(se) <- assays(se)[c(which(assayNames(se) %in% "counts"), which(assayNames(se) %ni% "counts"))]

  #Identify Top and Bottom Idx
  stopifnot(!any(topGroup %in% bottomGroup)) #topGroup and bottomGroup can not intersect
  topGroupIdx <- which(paste0(colData(se)[,compareCol]) %in% topGroup)
  bottomGroupIdx <- which(paste0(colData(se)[,compareCol]) %in% bottomGroup)
  stopifnot(length(topGroupIdx) > 1) #topGroup needs 2 samples
  stopifnot(length(bottomGroupIdx) > 1) #bottomGroup needs 2 samples
  message(sprintf("Comparing\nTop Samples: %s\nvs\nBottom Samples: %s", paste(colnames(se)[topGroupIdx],collapse=" "), paste(colnames(se)[bottomGroupIdx],collapse=" ")))

  if(!is.null(compareTop)){
    message("Subsetting the top features from Top and Bottom!")
    topFeaturesIdx <- head(order(rowMeans(edgeR::cpm(assay(se[,c(topGroupIdx)]))),decreasing=TRUE),compareTop)
    bottomFeaturesIdx <- head(order(rowMeans(edgeR::cpm(assay(se[,c(bottomGroupIdx)]))),decreasing=TRUE),compareTop)
    keepIdx <- unique(c(topFeaturesIdx,bottomFeaturesIdx))
    se <- se[keepIdx,]
  }

  #Keep only Columns important for comparison
  se <- se[,c(topGroupIdx, bottomGroupIdx)]
  colData(se)[,compareCol,drop=FALSE] <- paste0(colData(se)[,compareCol])

  #Round to Integer
  if(round){
    message("Rounding to integer for pairwise Testing!")
    assay(se) <- round(assay(se))
  }

  #Design
  design <- as.formula(paste0("~",compareCol))

  #Run Differential Testing
  message("Finding Initial Differential Parameters based on method of choosing...")
  diffObj <- pairwiseTest_helper(se = se, method = method, design = design, 
  	fitType = fitType, threads = threads, housekeeping = housekeeping, BPPARAM = BPPARAM)

  #Pairwise Test
  diffTest <- diffObj$diffTestResults
  method <- diffObj$diffTestMethod
  rm(diffObj)
  gc()

  if(toupper(method) == "EDGER"){

    suppressPackageStartupMessages(require(edgeR))
    message(sprintf("Running Pairwise test using edgeR on comparison %s...", comparisonName))

    #Check input
    stopifnot(all(topGroup %in% colnames(diffTest$design)))
    stopifnot(all(bottomGroup %in% colnames(diffTest$design)))

    #Make Contrasts
    contrast <- rep(0, ncol(diffTest$design))
    contrast[which(colnames(diffTest$design) %in% topGroup)] <- 1/length(which(colnames(diffTest$design) %in% topGroup))
    contrast[which(colnames(diffTest$design) %in% bottomGroup)] <- -1/length(which(colnames(diffTest$design) %in% bottomGroup))

    #Test
    if(lfcThreshold!=0){ #Test a lfc
      e <- glmQLFTest(diffTest, contrast = contrast) %>% glmTreat(., contrast = contrast, lfc = lfcThreshold)
          topTags(n = nrow(diffTest$counts), sort.by = "none") %>% data.frame()
    }else{
      e <- glmQLFTest(diffTest, contrast = contrast) %>% 
        topTags(n = nrow(diffTest$counts), sort.by = "none") %>% data.frame()
    }

    message("Organizing Output...")
    out <- data.frame(row.names = row.names(e),
                    log2Mean = e$logCPM,
                    log2FoldChange = e$logFC,
                    lfcSE = NA, #edgeR does not give this intentionally
                    pval = e$PValue,
                    padj = e$FDR)

  }else{

    stop("Unrecognized Differential Method")

  }

  message("Returning Pairwise Summarized Experiment...")
  #if a padj value is NA set to 1 ie non significant and if 0 convert to the least non-zero because -log10(fdr) plots become unhappy
  out[is.na(out$padj),"padj"] <- 1
  out[which(out$padj==0),"padj"] <- min(out[which(out$padj!=0),"padj"])
  seDiff <- SummarizedExperiment(assays = SimpleList(differential = as.matrix(out)), rowRanges = rowRanges(se))
  seDiff@metadata$compName <- comparisonName
  seDiff@metadata$compGroups <- list(topGroup = topGroup, bottomGroup = bottomGroup)
       
  return(seDiff)

}

pairwiseTest_helper <- function(se, method = "edgeR", design = ~Group, 
	fitType = "parametric", threads = NULL, BPPARAM = NULL, housekeeping = NULL){

  #Make Model Matrix
  design0 <- as.formula(paste0(paste0(as.character(design),collapse="")," + 0"))
  modelMatrix <- model.matrix(design0, data = colData(se))
    
  if(length(all.vars(design)) == 1){
    colnames(modelMatrix) <- gsub(all.vars(design),"",colnames(modelMatrix))
  }
  
  #Check Housekeeping
  if(!is.null(housekeeping)){
    if(!all(housekeeping %in% seq_len(nrow(se)))){
      stop("Housekeeping features provided are not within seq_len(nrow(se))! Please check input!")
    }
  }

 if(toupper(method) == "EDGER"){

    mainMethod <- "edgeR"

    message("Running EdgeR TMM Pipeline...")
    suppressPackageStartupMessages(require(edgeR))

    out <- DGEList(assay(se)) %>% calcNormFactors(., method = "TMM") %>% 
      estimateDisp(., design = modelMatrix, robust = TRUE) %>%
      glmQLFit(., design = modelMatrix)


  }else{

    stop("Unrecognized Differential Method")

  }

  return(list(design = design, diffTestMethod = mainMethod, diffTestResults = out))

}

#########################
#########################

#Input Summarized Experiment
se <- readRDS("output/ATAC_Analysis/Main/ATAC-Summarized-Experiment.rds")

#Differential Test
diffTest <- pairwiseTest(se, compareCol = "Group",
 topGroup = unique(colData(se)$Group)[1], bottomGroup = unique(colData(se)$Group)[2], method = "edger")

#Save
saveRDS(diffTest, "Differential-Test.rds")


