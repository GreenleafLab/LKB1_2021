library(ArchR)

set.seed(1)

groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowMeans(mat[, which(groups == x), drop = F], 
                na.rm = na.rm)
        }
        else {
            rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

#Input Summarized Experiment
se <- readRDS("output/ATAC_Analysis/Main/ATAC-Summarized-Experiment.rds")

#Get Variable Peaks
gM <- groupMeans(assay(se), colData(se)$Group)
rV <- matrixStats::rowVars(gM)
varPeaks <- head(order(rV,decreasing=TRUE), 10000)

#Heatmap Averaged Groups
mat <- assays(se[varPeaks,])[["log2norm"]]
mat <- groupMeans(mat, colData(se)[,"Group"])
mat <- rowZscores(mat) 
mat[mat > ceiling] <- ceiling
mat[mat < -ceiling] <- -ceiling

#KMeans
nk <- 6 #Number of KMeans
km <- kmeans(mat, nk, 100, 10)

ht1 <- ArchR:::.ArchRHeatmap(
    mat = mat,
    scale = FALSE,
    limits = c(-1,1),
    split = km$cluster,
    clusterRows = TRUE,
    clusterCols = TRUE,
    labelCols = TRUE,
    showColDendrogram = TRUE,
    draw = FALSE,
    useRaster = TRUE
)

#Plot
pdf("~/Desktop/Variable-Peaks-Heatmap.pdf", width = 8, height = 8)
ht1
dev.off()